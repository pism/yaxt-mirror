/**
 * @file xt_exchanger_mix_isend_irecv.c
 *
 * @copyright Copyright  (C)  2014 Jörg Behrens <behrens@dkrz.de>
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
 * URL: https://redmine.dkrz.de/doc/yaxt/html/index.html
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

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt_exchanger.h"
#include "xt_exchanger_mix_isend_irecv.h"

static void xt_exchanger_mix_isend_irecv_delete(Xt_exchanger exchanger);
static void xt_exchanger_mix_isend_irecv_s_exchange(Xt_exchanger exchanger,
                                                    const void * src_data,
                                                    void * dst_data);

static const struct xt_exchanger_vtable exchanger_mix_isend_irecv_vtable = {
  .delete = xt_exchanger_mix_isend_irecv_delete,
  .s_exchange = xt_exchanger_mix_isend_irecv_s_exchange,
};

typedef struct Xt_exchanger_mix_isend_irecv_ * Xt_exchanger_mix_isend_irecv;

struct mix_msg {
  struct Xt_redist_msg data;
  enum {SEND, RECV} type;
};

struct Xt_exchanger_mix_isend_irecv_ {

  const struct xt_exchanger_vtable * vtable;

  int n, tag_offset;
  struct mix_msg * msgs;
  MPI_Comm comm;
};

Xt_exchanger
xt_exchanger_mix_isend_irecv_new(int nsend, int nrecv,
                                 struct Xt_redist_msg * send_msgs,
                                 struct Xt_redist_msg * recv_msgs,
                                 MPI_Comm comm) {

  Xt_exchanger_mix_isend_irecv exchanger = xmalloc(1 * sizeof(*exchanger));

  exchanger->vtable = &exchanger_mix_isend_irecv_vtable;

  struct Xt_redist_msg * temp_msgs;

  exchanger->comm = xt_mpi_comm_smart_dup(comm, &exchanger->tag_offset);
  exchanger->n = nsend + nrecv;
  exchanger->msgs = xmalloc((size_t)exchanger->n * sizeof(*(exchanger->msgs)));
  temp_msgs = xt_redist_msgs_copy((size_t)nsend, send_msgs, exchanger->comm);
  for (int i = 0; i < nsend; ++i) {
    exchanger->msgs[i].data = temp_msgs[i];
    exchanger->msgs[i].type = SEND;
  }
  free(temp_msgs);
  temp_msgs = xt_redist_msgs_copy((size_t)nrecv, recv_msgs, exchanger->comm);
  for (int i = 0; i < nrecv; ++i) {
    exchanger->msgs[i+nsend].data = temp_msgs[i];
    exchanger->msgs[i+nsend].type = RECV;
  }
  free(temp_msgs);

  xt_exchanger_internal_optimize((size_t)exchanger->n, exchanger->msgs,
                                 sizeof(*(exchanger->msgs)), exchanger->comm);

  for (int i = 1; i < exchanger->n; ++i) {

    if ((exchanger->msgs[i-1].data.rank == exchanger->msgs[i].data.rank) &&
        (exchanger->msgs[i].type == SEND)) {

      struct mix_msg temp = exchanger->msgs[i-1];
      exchanger->msgs[i-1] = exchanger->msgs[i];
      exchanger->msgs[i] = temp;
      i++;
    }
  }

  return (Xt_exchanger)exchanger;
}

static void xt_exchanger_mix_isend_irecv_delete(Xt_exchanger exchanger) {

  Xt_exchanger_mix_isend_irecv exchanger_msr =
    (Xt_exchanger_mix_isend_irecv)exchanger;

  size_t nmsg = (size_t)exchanger_msr->n;

  xt_redist_msgs_strided_destruct(nmsg, &exchanger_msr->msgs->data,
                                  exchanger_msr->comm,
                                  sizeof (*exchanger_msr->msgs));
  free(exchanger_msr->msgs);

  xt_mpi_comm_smart_dedup(&exchanger_msr->comm, exchanger_msr->tag_offset);
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
    MPI_Request *requests = xmalloc(nmsg * sizeof (*requests));
    for (size_t i = 0; i < nmsg; ++i)
      if (msgs[i].type == SEND)
        xt_mpi_call(MPI_Isend((void *)src_data, 1, msgs[i].data.datatype,
                              msgs[i].data.rank,
                              tag_offset + xt_mpi_tag_exchange_msg,
                              comm, requests+i), comm);
      else
        xt_mpi_call(MPI_Irecv(dst_data, 1, msgs[i].data.datatype,
                              msgs[i].data.rank,
                              tag_offset + xt_mpi_tag_exchange_msg,
                              comm, requests+i), comm);

    xt_mpi_call(MPI_Waitall((int)nmsg, requests, MPI_STATUSES_IGNORE), comm);
    free(requests);
  }
}
