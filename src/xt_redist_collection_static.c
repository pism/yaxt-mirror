/**
 * @file xt_redist_collection_static.c
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
#include <stdbool.h>
#include <stdlib.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/xt_redist_collection_static.h"
#include "xt/xt_redist_single_array_base.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt_redist_internal.h"
#include "xt_config_internal.h"


static size_t
generate_msg_infos(struct Xt_redist_msg ** msgs,
                   const MPI_Aint *displacements,
                   const Xt_redist *redists,
                   size_t num_redists, MPI_Comm comm,
                   enum xt_msg_direction direction)
{
  size_t num_ranks[num_redists];
  int *restrict ranks[num_redists];

  size_t num_messages = xt_redist_agg_msg_count(num_redists, direction,
                                                redists, num_ranks, ranks);
  /* get lists of ranks to send/receive message to/from */
  size_t rank_pos[num_redists];
  for (size_t j = 0; j < num_redists; ++j)
    rank_pos[j] = 0;
  struct Xt_redist_msg *p = NULL;
  if (num_messages) {
    MPI_Datatype datatypes[num_redists];
    int block_lengths[num_redists];
    for (size_t i = 0; i < num_redists; ++i)
      block_lengths[i] = 1;
    p = xmalloc(sizeof (*p) * num_messages);
    for (size_t i = 0; i < num_messages; ++i) {
      int min_rank = INT_MAX;
      for (size_t j = 0; j < num_redists; ++j)
        if (rank_pos[j] < num_ranks[j] && ranks[j][rank_pos[j]] < min_rank)
          min_rank = ranks[j][rank_pos[j]];

      for (size_t j = 0; j < num_redists; ++j)
        datatypes[j] =
          (rank_pos[j] < num_ranks[j] && ranks[j][rank_pos[j]] == min_rank)
          ? xt_redist_get_MPI_Datatype(redists[j], min_rank, direction, false)
          : MPI_DATATYPE_NULL;

      p[i].rank = min_rank;
      p[i].datatype
        = xt_create_compound_datatype(num_redists, displacements, datatypes,
                                      block_lengths, comm);
      for (size_t j = 0; j < num_redists; ++j) {
        rank_pos[j]
          += (rank_pos[j] < num_ranks[j] && ranks[j][rank_pos[j]] == min_rank);
      }
    }
  }
  free(ranks[0]);
  *msgs = p;
  return num_messages;
}

Xt_redist
xt_redist_collection_static_new(Xt_redist *redists, int num_redists,
                                const MPI_Aint src_displacements[num_redists],
                                const MPI_Aint dst_displacements[num_redists],
                                MPI_Comm comm)
{
  return xt_redist_collection_static_custom_new(
    redists, num_redists, src_displacements, dst_displacements, comm,
    (Xt_config)&xt_default_config);
}

Xt_redist
xt_redist_collection_static_custom_new(
  Xt_redist * redists, int num_redists,
  const MPI_Aint src_displacements[num_redists],
  const MPI_Aint dst_displacements[num_redists],
  MPI_Comm comm, Xt_config config) {
  // ensure that yaxt is initialized
  assert(xt_initialized());

  struct Xt_redist_msg *send_msgs, *recv_msgs;

  int tag_offset;
  MPI_Comm new_comm = xt_mpi_comm_smart_dup(comm, &tag_offset);

  xt_redist_check_comms(redists, num_redists, comm);

  size_t num_redists_ = num_redists >= 0 ? (size_t)num_redists : 0;
  size_t nsend
    = generate_msg_infos(&send_msgs, src_displacements, redists, num_redists_,
                         new_comm, SEND);

  size_t nrecv
    = generate_msg_infos(&recv_msgs, dst_displacements, redists, num_redists_,
                         new_comm, RECV);

  Xt_redist redist_collection =
    xt_redist_single_array_base_custom_new(
      (int)nsend, (int)nrecv, send_msgs, recv_msgs, new_comm, config);

  xt_redist_msgs_free(nsend, send_msgs, new_comm);
  xt_redist_msgs_free(nrecv, recv_msgs, new_comm);
  xt_mpi_comm_smart_dedup(&new_comm, tag_offset);
  return redist_collection;
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
