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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#if ! defined TLS && defined HAVE_PTHREAD
#include <pthread.h>
#endif

#include <stdlib.h>
#include "xt/xt_core.h"
#include "xt/xt_redist.h"
#include "xt/xt_mpi.h"
#include "core/ppm_xfuncs.h"
#include "xt_redist_internal.h"

static void
redist_s_exchange_internal_irecv_send(const void *src_data, void *dst_data,
                                      int nsend, int nrecv,
                                      struct Xt_redist_msg * send_msgs,
                                      struct Xt_redist_msg * recv_msgs,
                                      MPI_Comm comm);
static void
redist_s_exchange_internal_irecv_isend(const void *src_data, void *dst_data,
                                       int nsend, int nrecv,
                                       struct Xt_redist_msg * send_msgs,
                                       struct Xt_redist_msg * recv_msgs,
                                       MPI_Comm comm);

void (*xt_redist_s_exchange_internal)(
  const void *src_data, void *dst_data, int nsend, int nrecv,
  struct Xt_redist_msg * send_msgs, struct Xt_redist_msg * recv_msgs,
  MPI_Comm comm) =
    // redist_s_exchange_internal_irecv_send;
    redist_s_exchange_internal_irecv_isend;

void xt_redist_delete(Xt_redist redist) {

  redist->vtable->delete(redist);
}

void xt_redist_s_exchange(Xt_redist redist, int num_arrays,
                          const void **src_data, void **dst_data) {

  redist->vtable->s_exchange(redist, num_arrays, src_data, dst_data);
}

void xt_redist_s_exchange1(Xt_redist redist, const void *src_data, void *dst_data) {

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
redist_s_exchange_internal_irecv_send(const void *src_data, void *dst_data,
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
    xt_mpi_call(MPI_Send((void *)src_data, 1, send_msgs[i].datatype,
                         send_msgs[i].rank, 0, comm), comm);

  xt_mpi_call(MPI_Waitall(nrecv, recv_request, MPI_STATUSES_IGNORE), comm);

  free(recv_request);
}

static void
redist_s_exchange_internal_irecv_isend(const void *src_data, void *dst_data,
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
    xt_mpi_call(MPI_Isend((void *)src_data, 1, send_msgs[i].datatype,
                          send_msgs[i].rank, 0, comm, requests+nrecv+i),
                          comm);

  xt_mpi_call(MPI_Waitall(nrecv + nsend, requests, MPI_STATUSES_IGNORE), comm);

  free(requests);
}

#if ! defined TLS && defined HAVE_PTHREAD

static pthread_key_t tls_qsort_data_key;

void xt_redist_init(void) {
  int ierror = pthread_key_create(&tls_qsort_data_key, NULL);
  if (ierror) {
    fprintf(stderr, "%s: error creating pthread key: %s\n",
            __func__, strerror(ierror));
    die("Failed pthread_key_create in initialization!");
  }
}

#else

static TLS int comm_rank, comm_size;

#endif

static int compare_exchange_messages (const void * msg_a, const void * msg_b)
{
  int rank_a = *(int*)msg_a;
  int rank_b = *(int*)msg_b;
#if ! defined TLS && defined HAVE_PTHREAD
  int comm_rank, comm_size;
  {
    int *qsort_data = pthread_getspecific(tls_qsort_data_key);
    comm_size = qsort_data[0];
    comm_rank = qsort_data[1];
  }
#endif
  if (rank_a <= comm_rank) rank_a += comm_size;
  if (rank_b <= comm_rank) rank_b += comm_size;

  return ( rank_a - rank_b );
}

void xt_redist_internal_optimize(int nsend, int nrecv,
                                 void * send_msgs,
                                 void * recv_msgs,
                                 size_t msg_type_size,
                                 MPI_Comm comm) {

#if ! defined TLS && defined HAVE_PTHREAD
  int qsort_data[2];
  MPI_Comm_rank(comm, qsort_data + 1);
  MPI_Comm_size(comm, qsort_data + 0);
  pthread_setspecific(tls_qsort_data_key, qsort_data);
#else
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);
#endif
  /* In order to avoid congestion of messages, the order of send and receive
   * messages is changed. This is done by sorting the messages according to the
   * rank of the respective message partner. Before the sorting to ranks that
   * are smaller or equal to the local rank the size of the communicator is
   * added.
   * example: process 5 is supposed to communicate with processes: 9, 5, 2, 6, 1
   * 1. add comm_size(10): 9, 15, 12, 6, 11
   * 2. sort: 6, 9, 11, 12, 15 -> final order: 6, 9, 1, 2, 5
   */
  qsort(send_msgs, (size_t)nsend, msg_type_size, compare_exchange_messages);
  qsort(recv_msgs, (size_t)nrecv, msg_type_size, compare_exchange_messages);
}
