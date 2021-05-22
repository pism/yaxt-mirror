/**
 * @file xt_exchanger_mix_isend_irecv.c
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
 * URL: https://doc.redmine.dkrz.de/yaxt/html/
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

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt/xt_request_msgs.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt_exchanger.h"
#include "xt_exchanger_mix_isend_irecv.h"

static Xt_exchanger
xt_exchanger_mix_isend_irecv_copy(Xt_exchanger exchanger,
                                  MPI_Comm newComm, int new_tag_offset);
static void xt_exchanger_mix_isend_irecv_delete(Xt_exchanger exchanger);
static void xt_exchanger_mix_isend_irecv_s_exchange(Xt_exchanger exchanger,
                                                    const void * src_data,
                                                    void * dst_data);
static void xt_exchanger_mix_isend_irecv_a_exchange(
  Xt_exchanger exchanger, const void * src_data, void * dst_data,
  Xt_request *request);
static int
xt_exchanger_mix_isend_irecv_get_msg_ranks(Xt_exchanger exchanger,
                                           enum xt_msg_direction direction,
                                           int *restrict *ranks);

static MPI_Datatype
xt_exchanger_mix_isend_irecv_get_MPI_Datatype(Xt_exchanger exchanger,
                                              int rank,
                                              enum xt_msg_direction direction,
                                              bool do_dup);

static const struct xt_exchanger_vtable exchanger_mix_isend_irecv_vtable = {
  .copy = xt_exchanger_mix_isend_irecv_copy,
  .delete = xt_exchanger_mix_isend_irecv_delete,
  .s_exchange = xt_exchanger_mix_isend_irecv_s_exchange,
  .a_exchange = xt_exchanger_mix_isend_irecv_a_exchange,
  .get_msg_ranks = xt_exchanger_mix_isend_irecv_get_msg_ranks,
  .get_MPI_Datatype = xt_exchanger_mix_isend_irecv_get_MPI_Datatype,
};

typedef struct Xt_exchanger_mix_isend_irecv_ * Xt_exchanger_mix_isend_irecv;

struct mix_msg {
  struct Xt_redist_msg data;
#if SIZEOF_MPI_DATATYPE == 2 * SIZEOF_INT
#  define MSG_DIR(msg) ((enum xt_msg_direction)((msg).data.padding))
#  define type data.padding
#else
  enum xt_msg_direction type;
#  define MSG_DIR(msg) ((msg).type)
#endif
};

struct Xt_exchanger_mix_isend_irecv_ {

  const struct xt_exchanger_vtable * vtable;

  int n, tag_offset;
  MPI_Comm comm;
  struct mix_msg msgs[];
};

static Xt_exchanger_mix_isend_irecv
xt_exchanger_mix_isend_irecv_alloc(size_t nmsg)
{
  Xt_exchanger_mix_isend_irecv exchanger;
  size_t header_size = sizeof (*exchanger),
    body_size = sizeof (struct mix_msg) * nmsg;
  exchanger = xmalloc(header_size + body_size);
  exchanger->n = (int)nmsg;
  exchanger->vtable = &exchanger_mix_isend_irecv_vtable;
  return exchanger;
}

Xt_exchanger
xt_exchanger_mix_isend_irecv_new(int nsend, int nrecv,
                                 const struct Xt_redist_msg *send_msgs,
                                 const struct Xt_redist_msg *recv_msgs,
                                 MPI_Comm comm, int tag_offset) {

  /** note: tag_offset + xt_mpi_tag_exchange_msg must not
   *        be used on @a comm by any other part of the program during the
   *        lifetime of the created exchanger object
   */

  assert((nsend >= 0) & (nrecv >= 0));
  size_t nmsg = (size_t)nsend + (size_t)nrecv;
  Xt_exchanger_mix_isend_irecv exchanger
    = xt_exchanger_mix_isend_irecv_alloc(nmsg);
  exchanger->comm = comm;
  exchanger->tag_offset = tag_offset;
  struct mix_msg *restrict msgs = exchanger->msgs;
  xt_redist_msgs_strided_copy((size_t)nsend, send_msgs, sizeof (send_msgs[0]),
                              &(msgs[0].data), sizeof (msgs[0]), comm);
  for (size_t i = 0; i < (size_t)nsend; ++i)
    msgs[i].type = SEND;
  xt_redist_msgs_strided_copy((size_t)nrecv, recv_msgs, sizeof (recv_msgs[0]),
                              &(msgs[nsend].data), sizeof (msgs[0]), comm);
  for (size_t i = 0; i < (size_t)nrecv; ++i)
    msgs[i + (size_t)nsend].type = RECV;

  xt_exchanger_internal_optimize(nmsg, msgs, sizeof(*msgs), comm);

  for (size_t i = 1; i < nmsg; ++i) {

    if (msgs[i-1].data.rank == msgs[i].data.rank && MSG_DIR(msgs[i]) == SEND) {

      struct mix_msg temp = msgs[i-1];
      msgs[i-1] = msgs[i];
      msgs[i] = temp;
      i++;
    }
  }

  return (Xt_exchanger)exchanger;
}

static Xt_exchanger
xt_exchanger_mix_isend_irecv_copy(Xt_exchanger exchanger,
                                  MPI_Comm new_comm, int new_tag_offset)
{
  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;
  size_t nmsg = (size_t)exchanger_msr->n;
  Xt_exchanger_mix_isend_irecv exchanger_copy
    = xt_exchanger_mix_isend_irecv_alloc(nmsg);
  exchanger_copy->comm = new_comm;
  exchanger_copy->tag_offset = new_tag_offset;
  struct mix_msg *restrict new_msgs = exchanger_copy->msgs,
    *restrict orig_msgs = exchanger_msr->msgs;
  xt_redist_msgs_strided_copy(nmsg, &orig_msgs->data, sizeof (*orig_msgs),
                              &new_msgs->data, sizeof (*new_msgs),
                              new_comm);
  for (size_t i = 0; i < nmsg; ++i)
    new_msgs[i].type = orig_msgs[i].type;
  return (Xt_exchanger)exchanger_copy;
}

enum { max_on_stack_req = 16 };

static void xt_exchanger_mix_isend_irecv_delete(Xt_exchanger exchanger) {

  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;

  size_t nmsg = (size_t)exchanger_msr->n;
  struct mix_msg *restrict msgs = exchanger_msr->msgs;

  xt_redist_msgs_strided_destruct(nmsg, &msgs[0].data, exchanger_msr->comm,
                                  sizeof (*msgs));
  free(exchanger_msr);
}

static void xt_exchanger_mix_isend_irecv_s_exchange(Xt_exchanger exchanger,
                                                    const void * src_data,
                                                    void * dst_data) {

  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;

  if (exchanger_msr->n > 0) {
    size_t nmsg = (size_t)exchanger_msr->n;
    MPI_Comm comm = exchanger_msr->comm;
    struct mix_msg *restrict msgs = exchanger_msr->msgs;
    int tag_offset = exchanger_msr->tag_offset;
    MPI_Request req_buf[max_on_stack_req];
    MPI_Request *requests
      = nmsg <= max_on_stack_req
      ? req_buf : xmalloc(nmsg * sizeof (*requests));
    for (size_t i = 0; i < nmsg; ++i) {
      typedef int (*ifp)(void *buf, int count, MPI_Datatype datatype, int dest,
                         int tag, MPI_Comm comm, MPI_Request *request);
      ifp op = MSG_DIR(msgs[i]) == SEND ? (ifp)MPI_Isend : (ifp)MPI_Irecv;
      void *data = MSG_DIR(msgs[i]) == SEND ? (void *)src_data : dst_data;
      xt_mpi_call(op(data, 1, msgs[i].data.datatype,
                     msgs[i].data.rank,
                     tag_offset + xt_mpi_tag_exchange_msg,
                     comm, requests+i), comm);
    }
    xt_mpi_call(MPI_Waitall((int)nmsg, requests, MPI_STATUSES_IGNORE), comm);
    if (requests != req_buf)
      free(requests);
  }
}

static void xt_exchanger_mix_isend_irecv_a_exchange(
  Xt_exchanger exchanger, const void * src_data, void * dst_data,
  Xt_request *request) {

  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;

  Xt_request requests = XT_REQUEST_NULL;

  if (exchanger_msr->n > 0) {
    size_t nmsg = (size_t)exchanger_msr->n;
    MPI_Comm comm = exchanger_msr->comm;
    struct mix_msg *restrict msgs = exchanger_msr->msgs;
    int tag_offset = exchanger_msr->tag_offset;
    MPI_Request req_buf[max_on_stack_req];
    MPI_Request *tmp_requests
      = nmsg <= max_on_stack_req
      ? req_buf : xmalloc(nmsg * sizeof (*tmp_requests));
    for (size_t i = 0; i < nmsg; ++i) {
      typedef int (*ifp)(void *buf, int count, MPI_Datatype datatype, int dest,
                         int tag, MPI_Comm comm, MPI_Request *request);
      ifp op = MSG_DIR(msgs[i]) == SEND ? (ifp)MPI_Isend : (ifp)MPI_Irecv;
      void *data = MSG_DIR(msgs[i]) == SEND ? (void *)src_data : dst_data;
      xt_mpi_call(op(data, 1, msgs[i].data.datatype,
                     msgs[i].data.rank,
                     tag_offset + xt_mpi_tag_exchange_msg,
                     comm, tmp_requests+i), comm);
    }
    requests = xt_request_msgs_new((int)nmsg, tmp_requests, comm);
    if (tmp_requests != req_buf)
      free(tmp_requests);
  }

  *request = requests;
}

static int
xt_exchanger_mix_isend_irecv_get_msg_ranks(Xt_exchanger exchanger,
                                           enum xt_msg_direction direction,
                                           int *restrict *ranks)
{
  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;
  size_t nmsg = 0, nmsg_all = (size_t)exchanger_msr->n;
  const struct mix_msg *restrict msgs = exchanger_msr->msgs;
  for (size_t i = 0; i < nmsg_all; ++i)
    nmsg += MSG_DIR(msgs[i]) == direction;
  int *restrict ranks_ = *ranks;
  if (!ranks_)
    ranks_ = *ranks = xmalloc(nmsg * sizeof (*ranks_));
  for (size_t i = 0, j = (size_t)-1; i < nmsg_all; ++i)
    if (MSG_DIR(msgs[i]) == direction)
      ranks_[++j] = msgs[i].data.rank;
  return (int)nmsg;
}

static MPI_Datatype
xt_exchanger_mix_isend_irecv_get_MPI_Datatype(Xt_exchanger exchanger,
                                              int rank,
                                              enum xt_msg_direction direction,
                                              bool do_dup)
{
  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;
  size_t nmsg = (size_t)exchanger_msr->n;
  struct mix_msg *restrict msgs = exchanger_msr->msgs;
  MPI_Datatype datatype_copy = MPI_DATATYPE_NULL;
  for (size_t i = 0; i < nmsg; ++i)
    if (MSG_DIR(msgs[i]) == direction && msgs[i].data.rank == rank) {
      if (do_dup)
        xt_mpi_call(MPI_Type_dup(msgs[i].data.datatype, &datatype_copy),
                    exchanger_msr->comm);
      else
        datatype_copy = msgs[i].data.datatype;
      break;
    }
  return datatype_copy;
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
