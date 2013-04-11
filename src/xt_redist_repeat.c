/**
 * @file xt_redist_repeat.c
 *
 * @copyright Copyright  (C)  2012 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#include <assert.h>
#include <stdlib.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt/xt_redist_repeat.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt_redist_internal.h"

static void
redist_repeat_delete(Xt_redist redist);

static void
redist_repeat_s_exchange(Xt_redist redist, unsigned num_arrays,
                         void **src_data, void **dst_data);

static void
redist_repeat_s_exchange1(Xt_redist redist,
                          void *src_data, void *dst_data);

static MPI_Datatype
redist_repeat_get_send_MPI_Datatype(Xt_redist redist, int rank);

static MPI_Datatype
redist_repeat_get_recv_MPI_Datatype(Xt_redist redist, int rank);

static MPI_Comm
redist_repeat_get_MPI_Comm(Xt_redist redist);

static const struct xt_redist_vtable redist_repeat_vtable = {
  .delete                = redist_repeat_delete,
  .s_exchange            = redist_repeat_s_exchange,
  .s_exchange1           = redist_repeat_s_exchange1,
  .get_send_MPI_Datatype = redist_repeat_get_send_MPI_Datatype,
  .get_recv_MPI_Datatype = redist_repeat_get_recv_MPI_Datatype,
  .get_MPI_Comm          = redist_repeat_get_MPI_Comm
};

struct redist_repeat_msg {

  int rank;
  MPI_Datatype datatype;
};

typedef struct Xt_redist_repeat_ *Xt_redist_repeat;

struct Xt_redist_repeat_ {

  const struct xt_redist_vtable *vtable;

  int ndst, nsrc;
  struct redist_repeat_msg * send_msgs;
  struct redist_repeat_msg * recv_msgs;

  MPI_Comm comm;
};

static void
generate_msg_infos(struct redist_repeat_msg ** msgs, int * nmsgs,
                   MPI_Aint extent, int * displacements, Xt_redist redist,
                   unsigned num_repetitions, MPI_Comm comm,
                   MPI_Datatype (*get_MPI_datatype)(Xt_redist,int)) {

  size_t msgs_array_size = 0;

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  MPI_Datatype datatype;

  assert(*nmsgs >= 0);
  size_t num_messages = (size_t)*nmsgs;
  struct redist_repeat_msg *p = *msgs;
  for (int i = 0; i < comm_size; ++i) {

    datatype = get_MPI_datatype(redist, i);

    if (datatype != MPI_DATATYPE_NULL) {

      MPI_Aint curr_lb, curr_extent, curr_true_lb, curr_true_extent;
      MPI_Datatype datatype_with_extent;

      // adjust extent of datatype to match the displacements
      xt_mpi_call(MPI_Type_get_extent(datatype, &curr_lb, &curr_extent), comm);
      xt_mpi_call(MPI_Type_get_true_extent(datatype, &curr_true_lb,
                                           &curr_true_extent), comm);
      if (curr_true_extent > extent)
        Xt_abort(comm, "ERROR: new datatype extent is too small (Xt_redist_repeat)\n",
                 __FILE__, __LINE__);
      xt_mpi_call(MPI_Type_create_resized(datatype, curr_lb, extent,
                                          &datatype_with_extent), comm);

      ENSURE_ARRAY_SIZE(p, msgs_array_size, num_messages+1);

      p[num_messages].rank = i;
      p[num_messages].datatype =
        xt_mpi_generate_datatype(displacements, num_repetitions,
                                 datatype_with_extent, comm);
      ++num_messages;

      MPI_Type_free(&datatype_with_extent);
      MPI_Type_free(&datatype);
    }
  }

  if (num_messages > 0)
    p = xrealloc(p, num_messages * sizeof(*p));
  *msgs = p;
  *nmsgs = (int)num_messages;
}

Xt_redist xt_redist_repeat_new(Xt_redist redist, MPI_Aint src_extent,
                               MPI_Aint dst_extent, unsigned num_repetitions,
                               int displacements[num_repetitions]) {

  Xt_redist_repeat redist_rep = xmalloc(1 * sizeof(*redist_rep));

  redist_rep->vtable = &redist_repeat_vtable;
  redist_rep->ndst = 0;
  redist_rep->nsrc = 0;
  redist_rep->send_msgs = NULL;
  redist_rep->recv_msgs = NULL;

  MPI_Comm comm = xt_redist_get_MPI_Comm(redist);
  xt_mpi_call(MPI_Comm_dup(comm, &(redist_rep->comm)), comm);

  if (num_repetitions < 1)
    Xt_abort(comm, "ERROR: invalid number of repetitions (Xt_redist_repeat)\n",
             __FILE__, __LINE__);


  generate_msg_infos(&(redist_rep->send_msgs), &(redist_rep->nsrc), src_extent,
                     displacements, redist, num_repetitions, redist_rep->comm,
                     xt_redist_get_send_MPI_Datatype);

  generate_msg_infos(&(redist_rep->recv_msgs), &(redist_rep->ndst), dst_extent, 
                     displacements, redist, num_repetitions, redist_rep->comm,
                     xt_redist_get_recv_MPI_Datatype);

  return (Xt_redist)redist_rep;
}

static inline Xt_redist_repeat
xrcs(void *redist)
{
  return (Xt_redist_repeat)redist;
}

static void
redist_repeat_s_exchange(Xt_redist redist,
                         unsigned num_arrays,
                         void **src_data,
                         void **dst_data)
{
  Xt_redist_repeat redist_rep = xrcs(redist);
  if (num_arrays == 1U)
    redist_repeat_s_exchange1(redist, src_data[0], dst_data[0]);
  else
    Xt_abort(redist_rep->comm, "ERROR: multi-array s_exchange is not"
             " implemented for this xt_redist type "
             "(Xt_redist_repeat)", __FILE__, __LINE__);
}

static void
redist_repeat_delete(Xt_redist redist) {

  Xt_redist_repeat redist_rep = xrcs(redist);

  for (int i = 0; i < redist_rep->nsrc; ++i)
    MPI_Type_free(&(redist_rep->send_msgs[i].datatype));
  free(redist_rep->send_msgs);

  for (int i = 0; i < redist_rep->ndst; ++i)
    MPI_Type_free(&(redist_rep->recv_msgs[i].datatype));
  free(redist_rep->recv_msgs);

  xt_mpi_call(MPI_Comm_free(&(redist_rep->comm)), MPI_COMM_WORLD);

  free(redist_rep);
}

static MPI_Datatype
copy_msg_dt(int nmsg, const struct redist_repeat_msg msg[nmsg],
            MPI_Comm comm, int rank)
{
  MPI_Datatype datatype_copy = MPI_DATATYPE_NULL;

  for (int i = 0; i < nmsg; ++i)
    if (msg[i].rank == rank) {
      xt_mpi_call(MPI_Type_dup(msg[i].datatype, &datatype_copy), comm);
      break;
    }

  return datatype_copy;
}


static MPI_Datatype
redist_repeat_get_send_MPI_Datatype(Xt_redist redist, int rank) {

  Xt_redist_repeat redist_rep = xrcs(redist);
  return copy_msg_dt(redist_rep->nsrc, redist_rep->send_msgs,
                     redist_rep->comm, rank);
}

static MPI_Datatype
redist_repeat_get_recv_MPI_Datatype(Xt_redist redist, int rank) {

  Xt_redist_repeat redist_rep = xrcs(redist);
  return copy_msg_dt(redist_rep->ndst, redist_rep->recv_msgs,
                     redist_rep->comm, rank);
}

static void
redist_repeat_s_exchange1(Xt_redist redist,
                                     void *src_data, void *dst_data) {

  Xt_redist_repeat redist_rep = xrcs(redist);

  MPI_Request * recv_requests;

  recv_requests = xmalloc((size_t)redist_rep->ndst * sizeof(*recv_requests));

  for (int i = 0; i < redist_rep->ndst; ++i)
    xt_mpi_call(MPI_Irecv(dst_data, 1, redist_rep->recv_msgs[i].datatype,
                          redist_rep->recv_msgs[i].rank, 0, redist_rep->comm,
                          recv_requests+i), redist_rep->comm);

  for (int i = 0; i < redist_rep->nsrc; ++i)
    xt_mpi_call(MPI_Send(src_data, 1, redist_rep->send_msgs[i].datatype,
                         redist_rep->send_msgs[i].rank, 0, redist_rep->comm),
                    redist_rep->comm);

  xt_mpi_call(MPI_Waitall(redist_rep->ndst, recv_requests,
                          MPI_STATUSES_IGNORE), redist_rep->comm);

  free(recv_requests);
}

static MPI_Comm
redist_repeat_get_MPI_Comm(Xt_redist redist) {

  Xt_redist_repeat redist_rep = xrcs(redist);

  return redist_rep->comm;
}
