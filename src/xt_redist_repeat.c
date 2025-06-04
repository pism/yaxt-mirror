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
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#include "xt_mpi_ddt_wrap.h"
#include "xt/xt_redist_repeat.h"
#include <xt_redist_single_array_base_internal.h>
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt_redist_internal.h"
#include "xt_config_internal.h"
#include "xt_exchanger_internal.h"



static const char filename[] = "xt_redist_repeat.c";

static void
generate_msg_infos(struct Xt_msg_param *restrict send_msgs,
                   struct Xt_msg_param *restrict recv_msgs,
                   Xt_exchanger_iter exch_iter,
                   int num_repetitions,
                   const MPI_Aint (*hdisp)[num_repetitions],
                   MPI_Comm comm)
{
  size_t send_msg_idx = 0, recv_msg_idx = 0;
  do {
    MPI_Datatype datatype
      = xt_exchanger_iter_get_MPI_Datatype(exch_iter);
    int direction
      = xt_exchanger_iter_get_direction(exch_iter);
    struct Xt_msg_param *msg = direction == xt_direction_send
      ? send_msgs + send_msg_idx++ : recv_msgs + recv_msg_idx++;
    msg->displacement = xt_exchanger_iter_get_displacement(exch_iter);
    int count = xt_exchanger_iter_get_count(exch_iter);
    msg->rank = xt_exchanger_iter_get_rank(exch_iter);
    msg->count = 1;
    Xt_Type_create_hindexed_block(num_repetitions, count,
                                  hdisp[direction],
                                  datatype, &msg->datatype, comm);
    xt_mpi_call(MPI_Type_commit(&msg->datatype), comm);
  } while (xt_exchanger_iter_next(exch_iter));
}

Xt_redist
xt_redist_repeat_asym_new(Xt_redist redist, MPI_Aint src_extent,
                          MPI_Aint dst_extent, int num_repetitions,
                          const int src_displacements[num_repetitions],
                          const int dst_displacements[num_repetitions])
{
  return xt_redist_repeat_asym_custom_new(
    redist, src_extent, dst_extent, num_repetitions,
    src_displacements, dst_displacements, (Xt_config)&xt_default_config);
}

Xt_redist
xt_redist_repeat_asym_custom_new(Xt_redist redist, MPI_Aint src_extent,
                                 MPI_Aint dst_extent, int num_repetitions,
                                 const int src_displacements[num_repetitions],
                                 const int dst_displacements[num_repetitions],
                                 Xt_config config)
{
  // ensure that yaxt is initialized
  assert(xt_initialized());

  int tag_offset;
  MPI_Comm comm
    = xt_mpi_comm_smart_dup(xt_redist_get_MPI_Comm(redist), &tag_offset);

  if (num_repetitions < 1)
    Xt_abort(comm, "ERROR: invalid number of repetitions (Xt_redist_repeat)\n",
             filename, __LINE__);

  int nmsg[2];
  for (int i = 0; i < 2; ++i)
    nmsg[i] = redist->vtable->get_num_msg(redist, (enum xt_msg_direction)i);
  size_t num_messages = (size_t)nmsg[SEND] + (size_t)nmsg[RECV];
  enum { auto_limit = 32 };
  MPI_Aint hdisp_auto[2][auto_limit],
    (*restrict hdisp)[num_repetitions] = hdisp_auto;
  struct Xt_msg_param msgs_auto[auto_limit],
    *msgs = msgs_auto;
  void *buf = NULL;
  if (num_repetitions > auto_limit
      || num_messages > auto_limit) {
    size_t msgs_alloc_size
      = (num_messages > auto_limit ? num_messages * sizeof (*msgs) : (size_t)0),
      hdisp_alloc_size
      = (num_repetitions > auto_limit
         ? (size_t)num_repetitions * sizeof (MPI_Aint) * 2
         : (size_t)0),
      alloc_size = msgs_alloc_size + hdisp_alloc_size;
    buf = xmalloc(alloc_size);
    if (num_repetitions > auto_limit)
      hdisp = buf;
    if (num_messages > auto_limit)
      msgs = (void *)((unsigned char *)buf + hdisp_alloc_size);
  }

  Xt_exchanger exchanger = xt_redist_get_exchanger(redist);
  for (size_t j = 0; j < (size_t)num_repetitions; ++j)
    hdisp[xt_direction_send][j] = (MPI_Aint)src_displacements[j] * src_extent;
  for (size_t j = 0; j < (size_t)num_repetitions; ++j)
    hdisp[xt_direction_recv][j] = (MPI_Aint)dst_displacements[j] * dst_extent;
  struct Xt_config_ config_ = *config;
  config_.flags |= exch_no_dt_dup;
  Xt_exchanger_iter exch_iter
    = xt_exchanger_get_iterator(exchanger, &config_);
  generate_msg_infos(msgs, msgs+nmsg[0], exch_iter,
                     num_repetitions,
                     (const MPI_Aint (*)[(size_t)num_repetitions])hdisp, comm);
  xt_exchanger_iter_delete(exch_iter);

  Xt_redist result = xt_redist_single_array_base_ofs_custom_new(
    nmsg[SEND], nmsg[RECV], msgs, msgs+nmsg[SEND], comm, &config_);
  if (num_repetitions > auto_limit
      || num_messages > auto_limit)
    free(buf);
  xt_mpi_comm_smart_dedup(&comm, tag_offset);
  return result;
}

Xt_redist xt_redist_repeat_new(Xt_redist redist, MPI_Aint src_extent,
                               MPI_Aint dst_extent, int num_repetitions,
                               const int displacements[num_repetitions]) {
  return xt_redist_repeat_asym_custom_new(
    redist, src_extent, dst_extent, num_repetitions, displacements,
    displacements, (Xt_config)&xt_default_config);
}

Xt_redist xt_redist_repeat_custom_new(Xt_redist redist, MPI_Aint src_extent,
                                      MPI_Aint dst_extent, int num_repetitions,
                                      const int displacements[num_repetitions],
                                      Xt_config config)
{
  return xt_redist_repeat_asym_custom_new(
    redist, src_extent, dst_extent, num_repetitions, displacements,
    displacements, config);
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
