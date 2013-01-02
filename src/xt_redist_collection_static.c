/**
 * @file xt_redist_collection_static.c
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
#include "xt/xt_redist_collection_static.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt_redist_internal.h"

static void
redist_collection_static_delete(Xt_redist redist);

static void
redist_collection_static_s_exchange(Xt_redist redist, void **src_data,
                                    unsigned num_src_arrays, void **dst_data,
                                    unsigned num_dst_arrays);

static void
redist_collection_static_s_exchange1(Xt_redist redist, void *src_data, void *dst_data);

static MPI_Datatype
redist_collection_static_get_send_MPI_Datatype(Xt_redist redist, int rank);

static MPI_Datatype
redist_collection_static_get_recv_MPI_Datatype(Xt_redist redist, int rank);

static const struct xt_redist_vtable redist_collection_static_vtable = {
  .delete                = redist_collection_static_delete,
  .s_exchange            = redist_collection_static_s_exchange,
  .s_exchange1           = redist_collection_static_s_exchange1,
  .get_send_MPI_Datatype = redist_collection_static_get_send_MPI_Datatype,
  .get_recv_MPI_Datatype = redist_collection_static_get_recv_MPI_Datatype
};

struct redist_collection_static_msg {

  int rank;
  MPI_Datatype datatype;
};

struct Xt_redist_collection_static {

  const struct xt_redist_vtable *vtable;

  int ndst, nsrc;
  struct redist_collection_static_msg * send_msgs;
  struct redist_collection_static_msg * recv_msgs;

  MPI_Comm comm;
};

static MPI_Datatype
generate_datatype(MPI_Aint * displacements, int * block_lengths,
                  unsigned num_redists, MPI_Datatype * datatypes,
                  MPI_Comm comm) {

  MPI_Datatype datatype;

  unsigned num_datatypes = 0;

  for (unsigned i = 0; i < num_redists; ++i)
    if (datatypes[i] != MPI_DATATYPE_NULL)
      ++num_datatypes;

  MPI_Datatype * datatypes_;
  MPI_Aint * displacements_;

  if (num_datatypes != num_redists) {

    datatypes_ = xmalloc(num_datatypes * sizeof(*datatypes_));
    displacements_ = xmalloc(num_datatypes * sizeof(*displacements));

    num_datatypes = 0;

    for (unsigned i = 0; i < num_redists; ++i) {
      if (datatypes[i] != MPI_DATATYPE_NULL) {

        datatypes_[num_datatypes] = datatypes[i];
        displacements_[num_datatypes] = displacements[i];
        ++num_datatypes;
      }
    }
  } else {

    datatypes_ = datatypes;
    displacements_ = displacements;
  }

  assert(num_datatypes <= INT_MAX);
  xt_mpi_call(MPI_Type_create_struct((int)num_datatypes, block_lengths,
                                     displacements_, datatypes_, &datatype),
              comm);

  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  if (num_datatypes != num_redists) {
    free(datatypes_);
    free(displacements_);
  }

  return datatype;
}

static void
generate_msg_infos(struct redist_collection_static_msg ** msgs, int * nmsgs,
                   MPI_Aint * displacements, Xt_redist * redists,
                   unsigned num_redists, MPI_Comm comm,
                   MPI_Datatype (*get_MPI_datatype)(Xt_redist,int)) {

  size_t msgs_array_size = 0;

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int block_lengths[num_redists];
  MPI_Datatype datatypes[num_redists];

  for (unsigned i = 0; i < num_redists; ++i)
    block_lengths[i] = 1;

  assert(*nmsgs >= 0);
  size_t num_messages = (size_t)*nmsgs;
  struct redist_collection_static_msg *p = *msgs;
  for (int i = 0; i < comm_size; ++i) {

    for (unsigned j = 0; j < num_redists; ++j)
      datatypes[j] = get_MPI_datatype(redists[j], i);

    for (unsigned j = 0; j < num_redists; ++j) {

      if (datatypes[j] != MPI_DATATYPE_NULL) {

        ENSURE_ARRAY_SIZE(p, msgs_array_size, num_messages+1);

        p[num_messages].rank = i;
        p[num_messages].datatype
          = generate_datatype(displacements, block_lengths,
                              num_redists, datatypes, comm);
        ++num_messages;

        break;
      }
    }

    for (unsigned j = 0; j < num_redists; ++j)
      if (datatypes[j] != MPI_DATATYPE_NULL)
        xt_mpi_call(MPI_Type_free(datatypes+j), comm);
  }

  if (num_messages > 0)
    p = xrealloc(p, num_messages * sizeof(*p));
  *msgs = p;
  assert(num_messages <= INT_MAX);
  *nmsgs = (int)num_messages;
}

Xt_redist xt_redist_collection_static_new(Xt_redist * redists, unsigned num_redists,
                                          MPI_Aint * src_displacements,
                                          MPI_Aint * dst_displacements, MPI_Comm comm) {

  struct Xt_redist_collection_static * redist_coll;

  redist_coll = xmalloc(1 * sizeof(*redist_coll));

  redist_coll->vtable = &redist_collection_static_vtable;
  redist_coll->ndst = 0;
  redist_coll->nsrc = 0;
  redist_coll->send_msgs = NULL;
  redist_coll->recv_msgs = NULL;
  xt_mpi_call(MPI_Comm_dup(comm, &(redist_coll->comm)), comm);

  generate_msg_infos(&(redist_coll->send_msgs), &(redist_coll->nsrc), 
                     src_displacements, redists, num_redists, redist_coll->comm,
                     xt_redist_get_send_MPI_Datatype);

  generate_msg_infos(&(redist_coll->recv_msgs), &(redist_coll->ndst),
                     dst_displacements, redists, num_redists, redist_coll->comm,
                     xt_redist_get_recv_MPI_Datatype);

  return (Xt_redist)redist_coll;
}

static void
redist_collection_static_s_exchange(Xt_redist redist,
                                    void **XT_UNUSED(src_data),
                                    unsigned XT_UNUSED(num_src_arrays),
                                    void **XT_UNUSED(dst_data),
                                    unsigned XT_UNUSED(num_dst_arrays)) {

  struct Xt_redist_collection_static *redist_coll
    = (struct Xt_redist_collection_static *)redist;

  Xt_abort(redist_coll->comm, "ERROR: s_exchange is not implemented for"
           " this xt_redist type (xt_redist_collection_static)", __FILE__, __LINE__);
}

static void
redist_collection_static_delete(Xt_redist redist) {

  struct Xt_redist_collection_static *redist_coll
    = (struct Xt_redist_collection_static *)redist;

  int i;

  for (i = 0; i < redist_coll->nsrc; ++i)
    MPI_Type_free(&(redist_coll->send_msgs[i].datatype));
  free(redist_coll->send_msgs);

  for (i = 0; i < redist_coll->ndst; ++i)
    MPI_Type_free(&(redist_coll->recv_msgs[i].datatype));
  free(redist_coll->recv_msgs);

  xt_mpi_call(MPI_Comm_free(&(redist_coll->comm)), MPI_COMM_WORLD);

  free(redist_coll);
}

static MPI_Datatype
redist_collection_static_get_send_MPI_Datatype(Xt_redist redist, int rank) {

  struct Xt_redist_collection_static * redist_coll;

  redist_coll = (struct Xt_redist_collection_static *)redist;

  MPI_Datatype datatype_copy;

  datatype_copy = MPI_DATATYPE_NULL;

  for (int i = 0; i < redist_coll->nsrc; ++i)
    if (redist_coll->send_msgs[i].rank == rank) {
      xt_mpi_call(MPI_Type_dup(redist_coll->send_msgs[i].datatype, &datatype_copy),
                      redist_coll->comm);
      break;
    }

  return datatype_copy;
}

static MPI_Datatype
redist_collection_static_get_recv_MPI_Datatype(Xt_redist redist, int rank) {

  struct Xt_redist_collection_static * redist_coll;

  redist_coll = (struct Xt_redist_collection_static *)redist;

  MPI_Datatype datatype_copy;

  datatype_copy = MPI_DATATYPE_NULL;

  for (int i = 0; i < redist_coll->ndst; ++i)
    if (redist_coll->recv_msgs[i].rank == rank) {
      xt_mpi_call(MPI_Type_dup(redist_coll->recv_msgs[i].datatype, &datatype_copy),
                      redist_coll->comm);
      break;
    }

  return datatype_copy;
}

static void
redist_collection_static_s_exchange1(Xt_redist redist, void *src_data, void *dst_data) {

  struct Xt_redist_collection_static *redist_coll
    = (struct Xt_redist_collection_static *) redist;

  MPI_Request * recv_requests;

  recv_requests = xmalloc((size_t)redist_coll->ndst * sizeof(*recv_requests));

  for (int i = 0; i < redist_coll->ndst; ++i)
    xt_mpi_call(MPI_Irecv(dst_data, 1, redist_coll->recv_msgs[i].datatype,
                          redist_coll->recv_msgs[i].rank, 0, redist_coll->comm,
                          recv_requests+i), redist_coll->comm);

  for (int i = 0; i < redist_coll->nsrc; ++i)
    xt_mpi_call(MPI_Send(src_data, 1, redist_coll->send_msgs[i].datatype,
                             redist_coll->send_msgs[i].rank, 0, redist_coll->comm),
                    redist_coll->comm);

  xt_mpi_call(MPI_Waitall(redist_coll->ndst, recv_requests,
                              MPI_STATUSES_IGNORE), redist_coll->comm);

  free(recv_requests);
}
