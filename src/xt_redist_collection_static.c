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
#include "xt_redist_single_array_base.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt_redist_internal.h"


static MPI_Datatype
generate_datatype(unsigned num_redists,
                  MPI_Aint displacements[num_redists],
                  MPI_Datatype datatypes[num_redists],
                  int block_lengths[num_redists],
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
generate_msg_infos(struct Xt_redist_msg ** msgs, int * nmsgs,
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
  struct Xt_redist_msg *p = *msgs;
  for (int i = 0; i < comm_size; ++i) {

    int non_empty_xfer = 0;
    for (unsigned j = 0; j < num_redists; ++j)
      non_empty_xfer |= (datatypes[j] = get_MPI_datatype(redists[j], i))
        != MPI_DATATYPE_NULL;

    if (non_empty_xfer)
    {
      ENSURE_ARRAY_SIZE(p, msgs_array_size, num_messages+1);

      p[num_messages].rank = i;
      p[num_messages].datatype
        = generate_datatype(num_redists, displacements, datatypes,
                            block_lengths, comm);
      ++num_messages;

      for (unsigned j = 0; j < num_redists; ++j)
        if (datatypes[j] != MPI_DATATYPE_NULL)
          xt_mpi_call(MPI_Type_free(datatypes+j), comm);
    }
  }

  if (num_messages > 0)
    p = xrealloc(p, num_messages * sizeof(*p));
  *msgs = p;
  *nmsgs = (int)num_messages;
}

static void check_comms(Xt_redist * redists, unsigned num_redists,
                        MPI_Comm comm) {

  int result;

  for (unsigned i = 0; i < num_redists; ++i) {

    xt_mpi_call(MPI_Comm_compare(xt_redist_get_MPI_Comm(redists[i]),
                                 comm, &result), comm);

    if ((result != MPI_IDENT) && (result != MPI_CONGRUENT))
      Xt_abort(comm, "ERROR: MPI communicators do not match; cannot build "
               "redist_collection_static\n", __FILE__, __LINE__);
  }
}

Xt_redist
xt_redist_collection_static_new(Xt_redist * redists, unsigned num_redists,
                                MPI_Aint src_displacements[num_redists],
                                MPI_Aint dst_displacements[num_redists],
                                MPI_Comm comm) {

  int nsend, nrecv;
  struct Xt_redist_msg * send_msgs = NULL;
  struct Xt_redist_msg * recv_msgs = NULL;
  MPI_Comm new_comm;

  nsend = 0;
  nrecv = 0;
  send_msgs = NULL;
  recv_msgs = NULL;
  xt_mpi_call(MPI_Comm_dup(comm, &new_comm), comm);

  check_comms(redists, num_redists, comm);

  generate_msg_infos(&send_msgs, &nsend, src_displacements, redists,
                     num_redists, new_comm, xt_redist_get_send_MPI_Datatype);

  generate_msg_infos(&recv_msgs, &nrecv, dst_displacements, redists,
                     num_redists, new_comm, xt_redist_get_recv_MPI_Datatype);

  return xt_redist_single_array_base_new(nsend, nrecv, send_msgs, recv_msgs,
                                         new_comm);
}
