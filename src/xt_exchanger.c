/**
 * @file xt_exchanger.c
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

#include <string.h>
#include <stdlib.h>

#include "xt/xt_mpi.h"
#include "xt_exchanger.h"
#include "xt_exchanger_irecv_isend.h"
#include "xt_exchanger_mix_isend_irecv.h"
#include "xt_exchanger_irecv_isend_packed.h"
#include "xt_exchanger_neigh_alltoall.h"

Xt_exchanger (*xt_exchanger_default_constructor)
  (int nsend, int nrecv, const struct Xt_redist_msg *send_msgs,
   const struct Xt_redist_msg *recv_msgs, MPI_Comm comm, int tag_offset) =
    xt_exchanger_mix_isend_irecv_new;
    // xt_exchanger_irecv_isend_new;
    // xt_exchanger_neigh_alltoall_new;

Xt_exchanger
xt_exchanger_copy(Xt_exchanger exchanger, MPI_Comm new_comm, int new_tag_offset)
{
  return exchanger->vtable->copy(exchanger, new_comm, new_tag_offset);
}

void xt_exchanger_delete(Xt_exchanger exchanger) {

  exchanger->vtable->delete(exchanger);
}

void xt_exchanger_s_exchange(Xt_exchanger exchanger, const void *src_data, void *dst_data) {

  exchanger->vtable->s_exchange(exchanger, src_data, dst_data);
}

void xt_exchanger_a_exchange(Xt_exchanger exchanger,
                             const void *src_data, void *dst_data,
                             Xt_request *request) {

  exchanger->vtable->a_exchange(exchanger, src_data, dst_data, request);
}

static int compare_exchange_messages (const void * msg_a, const void * msg_b)
{
  int rank_a = *(const int *)msg_a;
  int rank_b = *(const int *)msg_b;
  return ( rank_a - rank_b );
}

void xt_exchanger_internal_optimize(size_t n, void * msgs, size_t msg_type_size,
                                    MPI_Comm comm) {

  int comm_size, comm_rank, is_inter;
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  xt_mpi_call(MPI_Comm_test_inter(comm, &is_inter), comm);
  int (*get_comm_size)(MPI_Comm, int *)
    = is_inter ? MPI_Comm_remote_size : MPI_Comm_size;
  xt_mpi_call(get_comm_size(comm, &comm_size), comm);
  /* In order to avoid congestion of messages, the order of send and receive
   * messages is changed. This is done by sorting the messages according to the
   * rank of the respective message partner. Before the sorting to ranks that
   * are smaller or equal to the local rank the size of the communicator is
   * added.
   * example: process 5 is supposed to communicate with processes: 9, 5, 2, 6, 1
   * 1. add comm_size(10): 9, 15, 12, 6, 11
   * 2. sort: 6, 9, 11, 12, 15 -> final order: 6, 9, 1, 2, 5
   */
  unsigned char *msgs_ = (unsigned char *)msgs;
  /* apply offset for ranks <= comm_rank */
  for (size_t i = 0; i < n; ++i) {
    int *p = (int *)(void *)(msgs_ + i * msg_type_size),
      r = *p;
    if (r <= comm_rank) *p = r + comm_size;
  }
  qsort(msgs, n, msg_type_size, compare_exchange_messages);
  /* undo offset */
  for (size_t i = 0; i < n; ++i) {
    int *p = (int *)(void *)(msgs_ + i * msg_type_size),
      r = *p;
    if (r >= comm_size) *p = r - comm_size;
  }
}


int
xt_exchanger_get_msg_ranks(Xt_exchanger exchanger,
                           enum xt_msg_direction direction,
                           int *restrict *ranks)
{
  return exchanger->vtable->get_msg_ranks(exchanger, direction, ranks);
}

MPI_Datatype
xt_exchanger_get_MPI_Datatype(Xt_exchanger exchanger,
                              int rank, enum xt_msg_direction direction) {
  return exchanger->vtable->get_MPI_Datatype(exchanger, rank, direction);
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
