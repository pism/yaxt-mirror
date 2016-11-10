/**
 * @file xt_redist_repeat.c
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
#include <limits.h>
#include <stdlib.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/xt_redist_repeat.h"
#include "xt_redist_single_array_base.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt_redist_internal.h"

static void
generate_msg_infos(struct Xt_redist_msg **msgs, int *nmsgs,
                   MPI_Aint extent, const int *displacements, Xt_redist redist,
                   int num_repetitions, MPI_Comm comm,
                   enum xt_msg_direction direction,
                   MPI_Datatype (*get_MPI_datatype)(Xt_redist,int)) {

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  assert(*nmsgs >= 0);
  size_t num_messages = (size_t)*nmsgs;
  int *restrict ranks = NULL;
  size_t num_ranks
    = (size_t)xt_redist_get_msg_ranks(redist, direction, (int **)&ranks);
  struct Xt_redist_msg *restrict p
    = xrealloc(*msgs, sizeof (*p) * (num_messages + num_ranks));
  for (size_t i = 0; i < num_ranks; ++i) {
    MPI_Datatype datatype = get_MPI_datatype(redist, ranks[i]);
    MPI_Aint curr_lb, curr_extent;
    MPI_Datatype datatype_with_extent;

    // adjust extent of datatype to match the displacements
    xt_mpi_call(MPI_Type_get_extent(datatype, &curr_lb, &curr_extent), comm);
    xt_mpi_call(MPI_Type_create_resized(datatype, curr_lb, extent,
                                        &datatype_with_extent), comm);

    p[num_messages + i].rank = ranks[i];
    p[num_messages + i].datatype
      = xt_mpi_generate_datatype(displacements, num_repetitions,
                                 datatype_with_extent, comm);
    MPI_Type_free(&datatype_with_extent);
    MPI_Type_free(&datatype);
  }
  free(ranks);
  *msgs = p;
  *nmsgs = (int)(num_messages + num_ranks);
}

Xt_redist xt_redist_repeat_asym_new(Xt_redist redist, MPI_Aint src_extent,
                                    MPI_Aint dst_extent, int num_repetitions,
                                    const int src_displacements[num_repetitions],
                                    const int dst_displacements[num_repetitions]) {

  int nsend = 0, nrecv = 0;
  struct Xt_redist_msg *send_msgs = NULL, *recv_msgs = NULL;
  int tag_offset;
  MPI_Comm comm
    = xt_mpi_comm_smart_dup(xt_redist_get_MPI_Comm(redist), &tag_offset);

  if (num_repetitions < 1)
    Xt_abort(comm, "ERROR: invalid number of repetitions (Xt_redist_repeat)\n",
             __FILE__, __LINE__);


  generate_msg_infos(&send_msgs, &nsend, src_extent,
                     src_displacements, redist, num_repetitions, comm,
                     SEND, xt_redist_get_send_MPI_Datatype);

  generate_msg_infos(&recv_msgs, &nrecv, dst_extent,
                     dst_displacements, redist, num_repetitions, comm,
                     RECV, xt_redist_get_recv_MPI_Datatype);

  Xt_redist result
    = xt_redist_single_array_base_new(nsend, nrecv, send_msgs, recv_msgs, comm);
  xt_mpi_comm_smart_dedup(&comm, tag_offset);
  return result;
}

Xt_redist xt_redist_repeat_new(Xt_redist redist, MPI_Aint src_extent,
                               MPI_Aint dst_extent, int num_repetitions,
                               const int displacements[num_repetitions]) {

  return xt_redist_repeat_asym_new(redist, src_extent,
                                   dst_extent, num_repetitions,
                                   displacements,
                                   displacements);
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
