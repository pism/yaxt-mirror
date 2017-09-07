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
#include "xt_redist_single_array_base.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt_redist_internal.h"
#include "xt/xt_sort.h"


static size_t
generate_msg_infos(struct Xt_redist_msg ** msgs,
                   const MPI_Aint *displacements, Xt_redist *redists,
                   size_t num_redists, MPI_Comm comm,
                   enum xt_msg_direction direction,
                   MPI_Datatype (*get_MPI_datatype)(Xt_redist,int)) {

  int block_lengths[num_redists];
  MPI_Datatype datatypes[num_redists];

  for (size_t i = 0; i < num_redists; ++i)
    block_lengths[i] = 1;

  size_t num_ranks[num_redists], rank_pos[num_redists];
  int *restrict ranks[num_redists];
  bool ranks_left = false;
  /* get lists of ranks to send/receive message to/from */
  for (size_t j = 0; j < num_redists; ++j) {
    num_ranks[j]
      = (size_t)xt_redist_get_msg_ranks(redists[j], direction,
                                        (int **)(ranks + j));
    /* sort list */
    xt_sort_int(ranks[j], num_ranks[j]);
    ranks_left |= (num_ranks[j] > 0);
    rank_pos[j] = 0;
  }

  size_t num_messages = ranks_left
    ? xt_ranks_uniq_count(num_redists, num_ranks, (const int **)ranks) : 0;
  struct Xt_redist_msg *p = NULL;
  if (num_messages) {
    p = xmalloc(sizeof (*p) * num_messages);
    for (size_t i = 0; i < num_messages; ++i) {
      int min_rank = INT_MAX;
      for (size_t j = 0; j < num_redists; ++j)
        if (rank_pos[j] < num_ranks[j] && ranks[j][rank_pos[j]] < min_rank)
          min_rank = ranks[j][rank_pos[j]];

      for (size_t j = 0; j < num_redists; ++j)
        datatypes[j] =
          (rank_pos[j] < num_ranks[j] && ranks[j][rank_pos[j]] == min_rank)
          ? get_MPI_datatype(redists[j], min_rank) : MPI_DATATYPE_NULL;

      p[i].rank = min_rank;
      p[i].datatype
        = xt_create_compound_datatype(num_redists, displacements, datatypes,
                                      block_lengths, comm);
      for (size_t j = 0; j < num_redists; ++j) {
        if (datatypes[j] != MPI_DATATYPE_NULL)
          xt_mpi_call(MPI_Type_free(datatypes+j), comm);
        rank_pos[j]
          += (rank_pos[j] < num_ranks[j] && ranks[j][rank_pos[j]] == min_rank);
      }
    }
  }
  for (size_t j = 0; j < num_redists; ++j)
    free(ranks[j]);
  *msgs = p;
  return num_messages;
}

Xt_redist
xt_redist_collection_static_new(Xt_redist * redists, int num_redists,
                                const MPI_Aint src_displacements[num_redists],
                                const MPI_Aint dst_displacements[num_redists],
                                MPI_Comm comm) {

  struct Xt_redist_msg *send_msgs, *recv_msgs;

  int tag_offset;
  MPI_Comm new_comm = xt_mpi_comm_smart_dup(comm, &tag_offset);

  xt_redist_check_comms(redists, num_redists, comm);

  size_t num_redists_ = num_redists >= 0 ? (size_t)num_redists : 0;
  size_t nsend
    = generate_msg_infos(&send_msgs, src_displacements, redists, num_redists_,
                         new_comm, SEND, xt_redist_get_send_MPI_Datatype);

  size_t nrecv
    = generate_msg_infos(&recv_msgs, dst_displacements, redists, num_redists_,
                         new_comm, RECV, xt_redist_get_recv_MPI_Datatype);

  Xt_redist redist_collection =
    xt_redist_single_array_base_new((int)nsend, (int)nrecv,
                                    send_msgs, recv_msgs, new_comm);

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
