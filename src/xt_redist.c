/**
 * @file xt_redist.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
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

#include <stdlib.h>
#include "xt/xt_core.h"
#include "xt/xt_redist.h"
#include "xt/xt_mpi.h"
#include "core/ppm_xfuncs.h"
#include "xt_redist_internal.h"

static void
redist_s_exchange_internal_irecv_send(void *src_data, void *dst_data,
                                      int nsend, int nrecv,
                                      struct Xt_redist_msg * send_msgs,
                                      struct Xt_redist_msg * recv_msgs,
                                      MPI_Comm comm);
static void
redist_s_exchange_internal_irecv_isend(void *src_data, void *dst_data,
                                       int nsend, int nrecv,
                                       struct Xt_redist_msg * send_msgs,
                                       struct Xt_redist_msg * recv_msgs,
                                       MPI_Comm comm);

void (*xt_redist_s_exchange_internal)(
  void *src_data, void *dst_data, int nsend, int nrecv,
  struct Xt_redist_msg * send_msgs, struct Xt_redist_msg * recv_msgs,
  MPI_Comm comm) =
    // redist_s_exchange_internal_irecv_send;
    redist_s_exchange_internal_irecv_isend;

void xt_redist_delete(Xt_redist redist) {

  redist->vtable->delete(redist);
}

void xt_redist_s_exchange(Xt_redist redist, int num_arrays,
                          void **src_data, void **dst_data) {

  redist->vtable->s_exchange(redist, num_arrays, src_data, dst_data);
}

void xt_redist_s_exchange1(Xt_redist redist, void *src_data, void *dst_data) {

  redist->vtable->s_exchange1(redist, src_data, dst_data);
}

MPI_Datatype xt_redist_get_send_MPI_Datatype(Xt_redist redist, int rank) {

  return redist->vtable->get_send_MPI_Datatype(redist, rank);
}

MPI_Datatype xt_redist_get_recv_MPI_Datatype(Xt_redist redist, int rank) {

  return redist->vtable->get_recv_MPI_Datatype(redist, rank);
}

MPI_Comm xt_redist_get_MPI_Comm(Xt_redist redist) {

  return redist->vtable->get_MPI_Comm(redist);
}

static void
redist_s_exchange_internal_irecv_send(void *src_data, void *dst_data,
                                      int nsend, int nrecv,
                                      struct Xt_redist_msg * send_msgs,
                                      struct Xt_redist_msg * recv_msgs,
                                      MPI_Comm comm) {

  MPI_Request *recv_request
    = xmalloc((size_t)nrecv * sizeof (*recv_request));

  for (int i = 0; i < nrecv; ++i)
    xt_mpi_call(MPI_Irecv(dst_data, 1, recv_msgs[i].datatype,
                          recv_msgs[i].rank, 0, comm,
                          recv_request+i), comm);

  for (int i = 0; i < nsend; ++i)
    xt_mpi_call(MPI_Send(src_data, 1, send_msgs[i].datatype,
                         send_msgs[i].rank, 0, comm), comm);

  xt_mpi_call(MPI_Waitall(nrecv, recv_request, MPI_STATUSES_IGNORE), comm);

  free(recv_request);
}

static void
redist_s_exchange_internal_irecv_isend(void *src_data, void *dst_data,
                                       int nsend, int nrecv,
                                       struct Xt_redist_msg * send_msgs,
                                       struct Xt_redist_msg * recv_msgs,
                                       MPI_Comm comm) {

  /* prevent warning about unused static function */
  (void)redist_s_exchange_internal_irecv_send;
  MPI_Request *requests
    = xmalloc((size_t)(nrecv + nsend) * sizeof (*requests));

  for (int i = 0; i < nrecv; ++i)
    xt_mpi_call(MPI_Irecv(dst_data, 1, recv_msgs[i].datatype,
                          recv_msgs[i].rank, 0, comm,
                          requests+i), comm);

  for (int i = 0; i < nsend; ++i)
    xt_mpi_call(MPI_Isend(src_data, 1, send_msgs[i].datatype,
                          send_msgs[i].rank, 0, comm, requests+nrecv+i),
                          comm);

  xt_mpi_call(MPI_Waitall(nrecv + nsend, requests, MPI_STATUSES_IGNORE), comm);

  free(requests);
}
