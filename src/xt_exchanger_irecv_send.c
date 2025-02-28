/**
 * @file xt_exchanger_irecv_send.c
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

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#if defined OMPI_MAJOR_VERSION &&                                       \
  OMPI_MAJOR_VERSION == 4 && OMPI_MINOR_VERSION == 0 && OMPI_RELEASE_VERSION <= 1
#  include <pthread.h>
#  include <stdio.h>
#endif
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_exchanger.h"
#include "xt_exchanger_internal.h"
#include "xt_config_internal.h"
#include "xt_mpi_internal.h"
#include "xt_request_msgs_internal.h"
#include "xt_exchanger_irecv_send.h"
#include "xt_exchanger_irecv_send_internal.h"
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
xt_exchanger_irecv_send_vtable[nvtab];

static void
xt_exchanger_irecv_send_s_exchange(
  Xt_exchanger exchanger,
  const void *src_data, void *dst_data);

#ifdef _OPENMP
static void
xt_exchanger_irecv_send_s_exchange_omp(
  Xt_exchanger exchanger, const void *src_data, void *dst_data);
#endif

static Xt_exchanger_thread_share
xt_exchanger_irecv_send_create_thread_share(
  Xt_exchanger exchanger, Xt_config config);

static void
xt_exchanger_irecv_send_destroy_thread_share(
  Xt_exchanger exchanger, Xt_exchanger_thread_share thread_share);

void
xt_exchanger_irecv_send_init(void)
{
#ifdef _OPENMP
  xt_exchanger_irecv_send_vtable[1] =
#endif
  xt_exchanger_irecv_send_vtable[0]
    = xt_exchanger_simple_base_vtable;
  xt_exchanger_irecv_send_vtable[0].s_exchange
    = xt_exchanger_irecv_send_s_exchange;
  xt_exchanger_irecv_send_vtable[0].a_exchange
    = (xt_exchanger_a_exchange_func)0;
  xt_exchanger_irecv_send_vtable[0].create_thread_share
    = xt_exchanger_irecv_send_create_thread_share;
  xt_exchanger_irecv_send_vtable[0].destroy_thread_share
    = xt_exchanger_irecv_send_destroy_thread_share;
#ifdef _OPENMP
  xt_exchanger_irecv_send_vtable[1].s_exchange
    = xt_exchanger_irecv_send_s_exchange_omp;
  xt_exchanger_irecv_send_vtable[1].a_exchange
    = (xt_exchanger_a_exchange_func)0;
  xt_exchanger_irecv_send_vtable[1].create_thread_share
    = xt_exchanger_irecv_send_create_thread_share;
  xt_exchanger_irecv_send_vtable[1].destroy_thread_share
    = xt_exchanger_irecv_send_destroy_thread_share;
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
xt_exchanger_irecv_send_s_exchange__(const void *src_data, void *dst_data,
                                     int nsend, int nrecv,
                                     const struct Xt_msg_param *send_msgs,
                                     const struct Xt_msg_param *recv_msgs,
                                     Xt_exchanger_simple_base base,
                                     MPI_Request *recv_request)
{
  int tag_offset = base->tag_offset;
  MPI_Comm comm = base->comm;
  for (int i = 0; i < nrecv; ++i) {
    void *dst_data_ = (unsigned char *)dst_data + recv_msgs[i].displacement;
    xt_mpi_call(MPI_Irecv(dst_data_, recv_msgs[i].count, recv_msgs[i].datatype,
                          recv_msgs[i].rank,
                          tag_offset + xt_mpi_tag_exchange_msg, comm,
                          recv_request+i), comm);
  }
#if XT_ENABLE_ALTERNATIVE_EXCHANGER_SEND_MODES
  typedef int (*ifp)(const void *buf, int count, MPI_Datatype datatype,
                     int dest, int tag, MPI_Comm comm);
  int send_mode
    = XT_CONFIG_FLAGS_GET_EXCH_SEND_COMM_MODE(base->config_flags);
  ifp send_op = (ifp)xt_send_func_tab[1][send_mode];
#else
#  define send_op MPI_Send
#endif
  for (int i = 0; i < nsend; ++i) {
    const void *src_data_
      = (const unsigned char *)src_data + send_msgs[i].displacement;
    xt_mpi_call(send_op(CAST_MPI_SEND_BUF(src_data_), send_msgs[i].count,
                        send_msgs[i].datatype,
                        send_msgs[i].rank,
                        tag_offset + xt_mpi_tag_exchange_msg, comm),
                comm);
  }
#if XT_ENABLE_ALTERNATIVE_EXCHANGER_SEND_MODES
#  undef send_op
#endif
  xt_mpi_call(MPI_Waitall(nrecv, recv_request, MPI_STATUSES_IGNORE), comm);
}

static void
xt_exchanger_irecv_send_s_exchange(
  Xt_exchanger exchanger, const void *src_data, void *dst_data)
{
  Xt_exchanger_simple_base base = (Xt_exchanger_simple_base)exchanger;
  int nsend = base->nmsg[xt_direction_send],
    nrecv = base->nmsg[xt_direction_recv];
  MPI_Request *requests
    = xmalloc((size_t)nrecv * sizeof (*requests));

  xt_exchanger_irecv_send_s_exchange__(
    src_data, dst_data, nsend, nrecv,
    xt_exchanger_simple_base_get_send_msg(base),
    xt_exchanger_simple_base_get_recv_msg(base),
    base, requests);

  free(requests);
}

#ifdef _OPENMP
static void
xt_exchanger_irecv_send_s_exchange_mt(const void *src_data, void *dst_data,
                                      int nsend, int nrecv,
                                      const struct Xt_msg_param *send_msgs,
                                      const struct Xt_msg_param *recv_msgs,
                                      Xt_exchanger_simple_base base,
                                      Xt_exchanger_thread_share shared_req)
{
  MPI_Request *recv_request
    = xt_request_msgs_get_req_ptr((Xt_request)shared_req);
  int num_threads = omp_get_num_threads(),
    tid = omp_get_thread_num();
  int start_send = (nsend * tid) / num_threads,
    nsend_ = (nsend * (tid+1)) / num_threads - start_send,
    start_recv = (nrecv * tid) / num_threads,
    nrecv_ = (nrecv * (tid+1)) / num_threads - start_recv;
  xt_exchanger_irecv_send_s_exchange__(src_data,
                                       dst_data,
                                       nsend_, nrecv_,
                                       send_msgs+start_send,
                                       recv_msgs+start_recv,
                                       base,
                                       recv_request+start_recv);
}

static void
xt_exchanger_irecv_send_s_exchange_omp(
  Xt_exchanger exchanger, const void *src_data, void *dst_data)
{

  struct Xt_config_ conf = xt_default_config;
  xt_config_set_redist_mthread_mode(&conf, XT_MT_OPENMP);
  Xt_exchanger_simple_base base = (Xt_exchanger_simple_base)exchanger;
  int nsend = base->nmsg[xt_direction_send],
    nrecv = base->nmsg[xt_direction_recv];
  Xt_exchanger_thread_share shared_req
    = (Xt_exchanger_thread_share)xt_request_msgs_alloc(
      nrecv, base->comm, &conf);
#pragma omp parallel
  xt_exchanger_irecv_send_s_exchange_mt(
    src_data, dst_data, nsend, nrecv,
    xt_exchanger_simple_base_get_send_msg(base),
    xt_exchanger_simple_base_get_recv_msg(base),
    base, shared_req);
  free(shared_req);
}

#endif

static Xt_exchanger_thread_share
xt_exchanger_irecv_send_create_thread_share(
  Xt_exchanger exchanger, Xt_config config)
{
  Xt_exchanger_simple_base base = (Xt_exchanger_simple_base)exchanger;
  Xt_request req = xt_request_msgs_alloc(
    base->nmsg[xt_direction_recv], base->comm, config);
  xt_request_msgs_set_auto_free(req, 0);
  return (Xt_exchanger_thread_share)req;
}

static void
xt_exchanger_irecv_send_destroy_thread_share(
  Xt_exchanger exchanger, Xt_exchanger_thread_share thread_share)
{
  (void)exchanger;
  /* assert that exchanger matches and thread_share contains no in-flight requests */
  free(thread_share);
}

#if defined OMPI_MAJOR_VERSION &&                                       \
  OMPI_MAJOR_VERSION == 4 && OMPI_MINOR_VERSION == 0 && OMPI_RELEASE_VERSION <= 1
static void
xt_ompi_ssend_diag(void)
{
  fputs("warning: with OpenMPI 4.0 and 4.0.1 a deadlock may occur\n"
        "when MPI_Ssend is called and the destination rank also does after\n"
        "posting a matching MPI_Irecv.\n"
        "For this reason, the requested combination of a synchronous mode,\n"
        "blocking exchanger is not fulfilled and standard mode will be used.\n"
        "*********************************\n"
        "PLEASE CONSIDER UPDATING YOUR MPI INSTALLATION, IF POSSIBLE\n"
        "*********************************\n",
        stderr);
}
#endif

Xt_exchanger
xt_exchanger_irecv_send_new(int nsend, int nrecv,
                            const struct Xt_msg_param * send_msgs,
                            const struct Xt_msg_param * recv_msgs,
                            MPI_Comm comm, int tag_offset, Xt_config config)
{
  /** note: tag_offset + xt_mpi_tag_exchange_msg must not
   *        be used on @a comm by any other part of the program during the
   *        lifetime of the created exchanger object
   */
  Xt_exchanger_simple_base exchanger =
    xt_exchanger_simple_base_alloc((size_t)nsend+(size_t)nrecv, 0);
  int mthread_mode = xt_config_get_redist_mthread_mode(config);
  xt_exchanger_simple_base_init(
    exchanger, nsend, nrecv, send_msgs, recv_msgs,
    comm, tag_offset, config);
#if defined OMPI_MAJOR_VERSION &&                                       \
  OMPI_MAJOR_VERSION == 4 && OMPI_MINOR_VERSION == 0 && OMPI_RELEASE_VERSION <= 1
  /*
   * OpenMPI 4.0 and 4.0.1 seem to have a bug in not making progress
   * when MPI_Ssend is called and the destination rank also does after
   * posting a matching MPI_Irecv. Hence prohibit creation of
   * synchronous exchanger here.
   */
  int config_flags = config->flags, send_mode
    = XT_CONFIG_FLAGS_GET_EXCH_SEND_COMM_MODE(config_flags);
  if (send_mode == xt_exch_send_synchronous) {
    exchanger->config_flags
      = XT_CONFIG_FLAGS_SET_EXCH_SEND_COMM_MODE(config_flags,
                                                xt_exch_send_standard);
    int rank;
    xt_mpi_call(MPI_Comm_rank(comm, &rank), comm);
    if (rank == 0) {
      static pthread_once_t ssend_diag_control = PTHREAD_ONCE_INIT;
      pthread_once(&ssend_diag_control, xt_ompi_ssend_diag);
    }
  }
#endif
  exchanger->vtable = xt_exchanger_irecv_send_vtable+mthread_mode;
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
