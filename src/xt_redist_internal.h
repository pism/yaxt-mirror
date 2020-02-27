/**
 * @file xt_redist_internal.h
 * @brief redistribution of data, non-public declarations
 *
 * contains declaration the redistribution data structure, which
 * is derived from one or more xt_xmaps
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
#ifndef XT_REDIST_INTERNAL_H
#define XT_REDIST_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <mpi.h>

#include "xt/xt_redist.h"
#include "xt/xt_request.h"

enum xt_msg_direction {SEND, RECV};

struct xt_redist_vtable {

  Xt_redist (*copy)(Xt_redist);
  void (*delete)(Xt_redist);
  void (*s_exchange)(Xt_redist, int, const void **, void **);
  void (*a_exchange)(Xt_redist, int, const void **, void **, Xt_request *);
  void (*s_exchange1)(Xt_redist, const void *, void *);
  void (*a_exchange1)(Xt_redist, const void *, void *, Xt_request *);
  MPI_Datatype (*get_msg_MPI_Datatype)(Xt_redist, int, enum xt_msg_direction);
  int (*get_num_msg)(Xt_redist, enum xt_msg_direction);
  int (*get_msg_ranks)(Xt_redist, enum xt_msg_direction, int *restrict *);
  MPI_Comm (*get_MPI_Comm)(Xt_redist);
};

struct Xt_redist_ {
  const struct xt_redist_vtable *vtable;
};

void
xt_redist_msgs_strided_copy(size_t n,
                            const struct Xt_redist_msg *restrict src,
                            size_t src_stride,
                            struct Xt_redist_msg *restrict dst,
                            size_t dst_stride,
                            MPI_Comm comm);

void xt_redist_msgs_strided_destruct(size_t n, struct Xt_redist_msg *msgs,
                                     MPI_Comm comm, size_t ofs_stride);

static inline void
xt_redist_msgs_free(size_t n, struct Xt_redist_msg *msgs, MPI_Comm comm)
{
  xt_redist_msgs_strided_destruct(n, msgs, comm, sizeof (*msgs));
  free(msgs);
}

void xt_redist_check_comms(Xt_redist *redists, int num_redists, MPI_Comm comm);

int xt_redist_get_msg_ranks(Xt_redist redist, enum xt_msg_direction direction,
                            int *restrict *ranks);

MPI_Datatype xt_redist_get_MPI_Datatype(Xt_redist redist, int rank,
                                        enum xt_msg_direction direction);

size_t
xt_ranks_uniq_count(size_t num_rank_sets,
                    const size_t num_ranks[num_rank_sets],
                    const int *const ranks[num_rank_sets]);

MPI_Datatype
xt_create_compound_datatype(size_t num_redists,
                            const MPI_Aint displacements[num_redists],
                            const MPI_Datatype datatypes[num_redists],
                            const int block_lengths[num_redists],
                            MPI_Comm comm);

unsigned
xt_redist_agg_msg_count(size_t num_redists, enum xt_msg_direction direction,
                        const Xt_redist redists[num_redists],
                        size_t num_ranks[num_redists],
                        int *restrict ranks[num_redists]);

#endif
/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
