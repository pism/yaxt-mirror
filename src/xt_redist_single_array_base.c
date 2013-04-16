/**
 * @file xt_redist_single_array_base.c
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
#include <stdio.h>
#include <assert.h>

#include <mpi.h>

#include "xt/xt_mpi.h"
#include "xt_redist_single_array_base.h"
#include "xt_redist_internal.h"
#include "xt/xt_xmap.h"
#include "xt/xt_idxlist.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

static void
redist_sab_delete(Xt_redist redist);

static void
redist_sab_s_exchange(Xt_redist redist, unsigned num_arrays,
                      void **src_data, void **dst_data);

static void
redist_sab_s_exchange1(Xt_redist redist, void *src_data, void *dst_data);

static MPI_Datatype
redist_sab_get_send_MPI_Datatype(Xt_redist redist, int rank);

static MPI_Datatype
redist_sab_get_recv_MPI_Datatype(Xt_redist redist, int rank);

static MPI_Comm
redist_sab_get_MPI_Comm(Xt_redist redist);

static const struct xt_redist_vtable redist_sab_vtable = {
  .delete                = redist_sab_delete,
  .s_exchange            = redist_sab_s_exchange,
  .s_exchange1           = redist_sab_s_exchange1,
  .get_send_MPI_Datatype = redist_sab_get_send_MPI_Datatype,
  .get_recv_MPI_Datatype = redist_sab_get_recv_MPI_Datatype,
  .get_MPI_Comm          = redist_sab_get_MPI_Comm
};

typedef struct Xt_redist_sab_ *Xt_redist_sab;

struct Xt_redist_sab_ {

  const struct xt_redist_vtable *vtable;

  int nsend, nrecv;

  struct Xt_redist_msg * send_msgs;
  struct Xt_redist_msg * recv_msgs;

  MPI_Comm comm;
};

Xt_redist xt_redist_single_array_base_new(int nsend, int nrecv,
                                          struct Xt_redist_msg * send_msgs,
                                          struct Xt_redist_msg * recv_msgs,
                                          MPI_Comm comm) {

  Xt_redist_sab redist = xmalloc(sizeof (*redist));

  redist->vtable = &redist_sab_vtable;
  redist->nsend = nsend;
  redist->nrecv = nrecv;
  redist->send_msgs = send_msgs;
  redist->recv_msgs = recv_msgs;
  redist->comm = comm;

  return (Xt_redist)redist;
}

static inline Xt_redist_sab
xrsab(void *redist)
{
  return (Xt_redist_sab)redist;
}

static void
redist_sab_delete(Xt_redist redist) {

  Xt_redist_sab redist_sab = xrsab(redist);

  int i;

  for (i = 0; i < redist_sab->nsend; ++i)
    if (redist_sab->send_msgs[i].datatype != MPI_DATATYPE_NULL)
      MPI_Type_free(&(redist_sab->send_msgs[i].datatype));
  free(redist_sab->send_msgs);

  for (i = 0; i < redist_sab->nrecv; ++i)
    if (redist_sab->recv_msgs[i].datatype != MPI_DATATYPE_NULL)
      MPI_Type_free(&(redist_sab->recv_msgs[i].datatype));
  free(redist_sab->recv_msgs);

  if ((redist_sab->comm != MPI_COMM_WORLD) &&
      (redist_sab->comm != MPI_COMM_SELF))
    MPI_Comm_free(&(redist_sab->comm));

  free(redist_sab);
}

static void
redist_sab_s_exchange(Xt_redist redist, unsigned num_arrays,
                      void **src_data, void **dst_data)
{
  Xt_redist_sab redist_rep = xrsab(redist);
  if (num_arrays == 1U)
    redist_sab_s_exchange1(redist, src_data[0], dst_data[0]);
  else
    Xt_abort(redist_rep->comm, "ERROR: multi-array s_exchange is not"
             " implemented for this xt_redist type "
             "(Xt_redist_single_array_base)", __FILE__, __LINE__);
}

static void
redist_sab_s_exchange1(Xt_redist redist, void *src_data, void *dst_data) {

  Xt_redist_sab redist_sab = xrsab(redist);

  MPI_Request *recv_request
    = xmalloc((size_t)redist_sab->nrecv * sizeof (*recv_request));

  int i;

  for (i = 0; i < redist_sab->nrecv; ++i)
    xt_mpi_call(MPI_Irecv(dst_data, 1, redist_sab->recv_msgs[i].datatype,
                          redist_sab->recv_msgs[i].rank, 0, redist_sab->comm,
                          recv_request+i), redist_sab->comm);

  for (i = 0; i < redist_sab->nsend; ++i)
    xt_mpi_call(MPI_Send(src_data, 1, redist_sab->send_msgs[i].datatype,
                         redist_sab->send_msgs[i].rank, 0, redist_sab->comm),
                redist_sab->comm);

  xt_mpi_call(MPI_Waitall(redist_sab->nrecv, recv_request,
                              MPI_STATUSES_IGNORE), redist_sab->comm);

  free(recv_request);
}

static MPI_Datatype
copy_msg_dt(int nmsg, const struct Xt_redist_msg msg[nmsg],
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
redist_sab_get_send_MPI_Datatype(Xt_redist redist, int rank) {

  Xt_redist_sab redist_sab = xrsab(redist);
  return copy_msg_dt(redist_sab->nsend, redist_sab->send_msgs,
                     redist_sab->comm, rank);
}

static MPI_Datatype
redist_sab_get_recv_MPI_Datatype(Xt_redist redist, int rank) {

  Xt_redist_sab redist_sab = xrsab(redist);
  return copy_msg_dt(redist_sab->nrecv, redist_sab->recv_msgs,
                     redist_sab->comm, rank);
}

static MPI_Comm
redist_sab_get_MPI_Comm(Xt_redist redist) {

  Xt_redist_sab redist_sab = xrsab(redist);

  return redist_sab->comm;
}
