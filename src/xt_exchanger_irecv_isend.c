/**
 * @file xt_exchanger_irecv_isend.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/ppm_xfuncs.h"
#include "xt_config_internal.h"
#include "xt/xt_mpi.h"
#include "xt/xt_request_msgs.h"
#include "xt_request_msgs_internal.h"
#include "xt_mpi_internal.h"
#include "xt_exchanger_internal.h"
#include "xt_exchanger_irecv_isend.h"
#include "xt_exchanger_irecv_isend_internal.h"
#include "xt_exchanger_simple_base.h"
#include "xt_exchanger_simple_base_internal.h"

enum {
#ifdef _OPENMP
  nvtab = 2,
#else
  nvtab = 1,
#endif
};

struct xt_exchanger_vtable
xt_exchanger_irecv_isend_vtable[nvtab];

static void
xt_exchanger_irecv_isend_s_exchange(
  Xt_exchanger exchanger,
  const void *src_data, void *dst_data);

static void
xt_exchanger_irecv_isend_a_exchange(
  Xt_exchanger exchanger,
  const void *src_data, void *dst_data,
  Xt_request *request);

#ifdef _OPENMP
static void
xt_exchanger_irecv_isend_s_exchange_omp(
  Xt_exchanger exchanger, const void *src_data, void *dst_data);

static void
xt_exchanger_irecv_isend_a_exchange_omp(Xt_exchanger exchanger,
                                        const void *src_data, void *dst_data,
                                        Xt_request *request);
#endif

static Xt_exchanger_thread_share
xt_exchanger_irecv_isend_create_thread_share(
  Xt_exchanger exchanger, Xt_config config);

static void
xt_exchanger_irecv_isend_destroy_thread_share(
  Xt_exchanger exchanger, Xt_exchanger_thread_share thread_share);

void
xt_exchanger_irecv_isend_init(void)
{
  xt_exchanger_irecv_isend_vtable[0]
    = xt_exchanger_simple_base_vtable;
#ifdef _OPENMP
  xt_exchanger_irecv_isend_vtable[1]
    = xt_exchanger_simple_base_vtable;
#endif
  xt_exchanger_irecv_isend_vtable[0].s_exchange
    = xt_exchanger_irecv_isend_s_exchange;
  xt_exchanger_irecv_isend_vtable[0].a_exchange
    = xt_exchanger_irecv_isend_a_exchange;
  xt_exchanger_irecv_isend_vtable[0].create_thread_share
    = xt_exchanger_irecv_isend_create_thread_share;
  xt_exchanger_irecv_isend_vtable[0].destroy_thread_share
    = xt_exchanger_irecv_isend_destroy_thread_share;
#ifdef _OPENMP
  xt_exchanger_irecv_isend_vtable[1].s_exchange
    = xt_exchanger_irecv_isend_s_exchange_omp;
  xt_exchanger_irecv_isend_vtable[1].a_exchange
    = xt_exchanger_irecv_isend_a_exchange_omp;
  xt_exchanger_irecv_isend_vtable[1].create_thread_share
    = xt_exchanger_irecv_isend_create_thread_share;
  xt_exchanger_irecv_isend_vtable[1].destroy_thread_share
    = xt_exchanger_irecv_isend_destroy_thread_share;
#endif
}



/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

static void
msgs_to_reqs(const void *src_data, void *dst_data,
             int nsend, int nrecv,
             const struct Xt_msg_param *send_msgs,
             const struct Xt_msg_param *recv_msgs,
             Xt_exchanger_simple_base base,
             MPI_Request *requests)
{
  int tag_offset = base->tag_offset;
  MPI_Comm comm = base->comm;
  for (int i = 0; i < nrecv; ++i)
    xt_mpi_call(MPI_Irecv(
                  (unsigned char *)dst_data + recv_msgs[i].displacement,
                  recv_msgs[i].count, recv_msgs[i].datatype, recv_msgs[i].rank,
                  tag_offset + xt_mpi_tag_exchange_msg, comm,
                  requests+i), comm);
#if XT_ENABLE_ALTERNATIVE_EXCHANGER_SEND_MODES
  typedef int (*ifp)(const void *buf, int count, MPI_Datatype datatype,
                     int dest, int tag, MPI_Comm comm, MPI_Request *request);
  enum Xt_exch_send_comm_mode send_mode
    = XT_CONFIG_FLAGS_GET_EXCH_SEND_COMM_MODE(base->config_flags);
  ifp a_send = (ifp)xt_send_func_tab[0][send_mode];
#else
#  define a_send MPI_Isend
#endif
  for (int i = 0; i < nsend; ++i)
    xt_mpi_call(a_send(
                  CAST_MPI_SEND_BUF((const unsigned char *)src_data
                                    + send_msgs[i].displacement),
                  send_msgs[i].count, send_msgs[i].datatype, send_msgs[i].rank,
                  tag_offset + xt_mpi_tag_exchange_msg, comm,
                  requests+nrecv+i), comm);
#if XT_ENABLE_ALTERNATIVE_EXCHANGER_SEND_MODES
#  undef a_send
#endif
}

static void
xt_exchanger_irecv_isend_s_exchange(
  Xt_exchanger exchanger, const void *src_data, void *dst_data)
{
  Xt_exchanger_simple_base base = (Xt_exchanger_simple_base)exchanger;
  int nsend = base->nmsg[xt_direction_send],
    nrecv = base->nmsg[xt_direction_recv];
  MPI_Request *requests
    = xmalloc((size_t)(nrecv + nsend) * sizeof (*requests));

  msgs_to_reqs(src_data, dst_data, nsend, nrecv,
               xt_exchanger_simple_base_get_send_msg(base),
               xt_exchanger_simple_base_get_recv_msg(base),
               base, requests);
  xt_mpi_call(MPI_Waitall(nrecv + nsend, requests, MPI_STATUSES_IGNORE),
              base->comm);
  free(requests);
}

#ifdef _OPENMP
static void
xt_exchanger_irecv_isend_a_exchange_mt(const void *src_data, void *dst_data,
                                       int nsend, int nrecv,
                                       const struct Xt_msg_param *send_msgs,
                                       const struct Xt_msg_param *recv_msgs,
                                       Xt_exchanger_simple_base base,
                                       Xt_exchanger_thread_share shared_req) {
  MPI_Request *requests
    = xt_request_msgs_get_req_ptr((Xt_request)shared_req);
  int num_threads = omp_get_num_threads(),
    tid = omp_get_thread_num();
  int start_send = (nsend * tid) / num_threads,
    nsend_ = (nsend * (tid+1)) / num_threads - start_send,
    start_recv = (nrecv * tid) / num_threads,
    nrecv_ = (nrecv * (tid+1)) / num_threads - start_recv;

  msgs_to_reqs(src_data, dst_data, nsend_, nrecv_,
               send_msgs+start_send, recv_msgs+start_recv,
               base, requests+start_recv+start_send);
}

static void
xt_exchanger_irecv_isend_s_exchange_mt(const void *src_data, void *dst_data,
                                       int nsend, int nrecv,
                                       const struct Xt_msg_param *send_msgs,
                                       const struct Xt_msg_param *recv_msgs,
                                       Xt_exchanger_simple_base base,
                                       Xt_exchanger_thread_share shared_req)
{
  int num_threads = omp_get_num_threads(),
    tid = omp_get_thread_num();
  xt_exchanger_irecv_isend_a_exchange_mt(
    src_data, dst_data, nsend, nrecv,
    send_msgs, recv_msgs, base, shared_req);
  MPI_Request *requests
    = xt_request_msgs_get_req_ptr((Xt_request)shared_req);
  int start_send = (nsend * tid) / num_threads,
    nsend_ = (nsend * (tid+1)) / num_threads - start_send,
    start_recv = (nrecv * tid) / num_threads,
    nrecv_ = (nrecv * (tid+1)) / num_threads - start_recv,
    start_req = start_send+start_recv,
    nreq_ = nrecv_+nsend_;
  xt_mpi_call(MPI_Waitall(nreq_, requests+start_req,
                          MPI_STATUSES_IGNORE), base->comm);
}

static void
xt_exchanger_irecv_isend_s_exchange_omp(
  Xt_exchanger exchanger, const void *src_data, void *dst_data)
{
  size_t num_threads_limit
    = xt_exchanger_simple_base_get_num_omp_threads(exchanger);
  if (num_threads_limit) {
    Xt_exchanger_simple_base base = (Xt_exchanger_simple_base)exchanger;
    int nsend = base->nmsg[xt_direction_send],
      nrecv = base->nmsg[xt_direction_recv];
    Xt_exchanger_thread_share shared_req
      = (Xt_exchanger_thread_share)xt_request_msgs_alloc(nsend+nrecv, base->comm,
                                                         &xt_default_config);
#pragma omp parallel num_threads(num_threads_limit)
    xt_exchanger_irecv_isend_s_exchange_mt(
      src_data, dst_data, nsend, nrecv,
      xt_exchanger_simple_base_get_send_msg(base),
      xt_exchanger_simple_base_get_recv_msg(base),
      base, shared_req);
    free(shared_req);
  }
}

#endif

static void
xt_exchanger_irecv_isend_a_exchange(
  Xt_exchanger exchanger,
  const void *src_data, void *dst_data,
  Xt_request *request)
{
  struct Xt_config_ conf = xt_default_config;
  xt_config_set_redist_mthread_mode(&conf, XT_MT_NONE);
  Xt_exchanger_simple_base base = (Xt_exchanger_simple_base)exchanger;
  int nsend = base->nmsg[xt_direction_send],
    nrecv = base->nmsg[xt_direction_recv];
  Xt_request requests_ = *request
    = xt_request_msgs_alloc(nrecv + nsend, base->comm, &conf);

  MPI_Request *requests
    = xt_request_msgs_get_req_ptr(requests_);
  msgs_to_reqs(src_data, dst_data, nsend, nrecv,
               xt_exchanger_simple_base_get_send_msg(base),
               xt_exchanger_simple_base_get_recv_msg(base),
               base, requests);
}

#ifdef _OPENMP
static void
xt_exchanger_irecv_isend_a_exchange_omp(Xt_exchanger exchanger,
                                        const void *src_data, void *dst_data,
                                        Xt_request *request) {

  size_t num_threads_limit
    = xt_exchanger_simple_base_get_num_omp_threads(exchanger);
  if (num_threads_limit) {
    struct Xt_config_ conf = xt_default_config;
    xt_config_set_redist_mthread_mode(&conf, XT_MT_OPENMP);
    Xt_exchanger_simple_base base = (Xt_exchanger_simple_base)exchanger;
    int nsend = base->nmsg[xt_direction_send],
      nrecv = base->nmsg[xt_direction_recv];
    Xt_exchanger_thread_share shared_req
      = (Xt_exchanger_thread_share)(
        *request = xt_request_msgs_alloc(nrecv + nsend, base->comm, &conf));
#pragma omp parallel num_threads(num_threads_limit)
    xt_exchanger_irecv_isend_a_exchange_mt(
      src_data, dst_data, nsend, nrecv,
      xt_exchanger_simple_base_get_send_msg(base),
      xt_exchanger_simple_base_get_recv_msg(base),
      base, shared_req);
  } else
    *request = XT_REQUEST_NULL;
}
#endif


static Xt_exchanger_thread_share
xt_exchanger_irecv_isend_create_thread_share(
  Xt_exchanger exchanger, Xt_config config)
{
  Xt_exchanger_simple_base base = (Xt_exchanger_simple_base)exchanger;
  int nsend = base->nmsg[xt_direction_send],
    nrecv = base->nmsg[xt_direction_recv];
  Xt_request request = xt_request_msgs_alloc(nsend + nrecv, base->comm, config);
  xt_request_msgs_set_auto_free(request, 0);
  return (Xt_exchanger_thread_share)request;
}

static void
xt_exchanger_irecv_isend_destroy_thread_share(
  Xt_exchanger exchanger, Xt_exchanger_thread_share thread_share)
{
  (void)exchanger;
  /* assert that thread_share contains no in-flight requests */
  free(thread_share);
}

Xt_exchanger
xt_exchanger_irecv_isend_new(int nsend, int nrecv,
                             const struct Xt_msg_param *send_msgs,
                             const struct Xt_msg_param *recv_msgs,
                             MPI_Comm comm, int tag_offset,
                             Xt_config config) {

  assert(nsend >= 0 && nrecv >= 0);
  /** note: tag_offset + xt_mpi_tag_exchange_msg must not
   *        be used on @a comm by any other part of the program during the
   *        lifetime of the created exchanger object
   */
  Xt_exchanger_simple_base exchanger =
    xt_exchanger_simple_base_alloc((size_t)nsend+(size_t)nrecv, 0);
  int mthread_mode = xt_config_get_redist_mthread_mode(config);
  xt_exchanger_simple_base_init(
    exchanger, nsend, nrecv, send_msgs, recv_msgs,
    comm, tag_offset,
    config);
  exchanger->vtable = xt_exchanger_irecv_isend_vtable+mthread_mode;
  return (Xt_exchanger)exchanger;
}


/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
