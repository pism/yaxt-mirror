/**
 * @file xt_redist_p2p.c
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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <mpi.h>

#include "xt/xt_mpi.h"
#include "xt/xt_redist_p2p.h"
#include "xt_redist_internal.h"
#include "xt_redist_single_array_base.h"
#include "xt/xt_xmap.h"
#include "xt/xt_idxlist.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

static MPI_Datatype
generate_datatype(int const * transfer_pos, int num_transfer_pos,
                  int *offsets, MPI_Datatype base_datatype, MPI_Comm comm) {

  MPI_Datatype type;

  int const * displ;
  int * tmp_displ = NULL;

  if (offsets != NULL) {

    tmp_displ = xmalloc((size_t)num_transfer_pos * sizeof(int));

    for (int i = 0; i < num_transfer_pos; ++i)
      tmp_displ[i] = offsets[transfer_pos[i]];

      displ = tmp_displ;

  } else
    displ = transfer_pos;

  type = xt_mpi_generate_datatype(displ, num_transfer_pos, base_datatype, comm);

  free(tmp_displ);

  return type;
}

static void
generate_msg_infos(int num_msgs, Xt_xmap_iter iter, int *offsets,
                   MPI_Datatype base_datatype, struct Xt_redist_msg ** msgs,
                   MPI_Comm comm) {

  if (num_msgs <= 0) {
    *msgs = NULL;
    return;
  }

  *msgs = xmalloc((size_t)num_msgs * sizeof(**msgs));

  struct Xt_redist_msg * curr_msg;

  curr_msg = *msgs;

  do {

    int const * curr_transfer_pos;
    int curr_num_transfer_pos;

    curr_transfer_pos = xt_xmap_iterator_get_transfer_pos(iter);
    curr_num_transfer_pos = xt_xmap_iterator_get_num_transfer_pos(iter);

    curr_msg->datatype
      = generate_datatype(curr_transfer_pos, curr_num_transfer_pos,
                          offsets, base_datatype, comm);
    curr_msg->rank = xt_xmap_iterator_get_rank(iter);

    curr_msg++;

  } while (xt_xmap_iterator_next(iter));
}

Xt_redist xt_redist_p2p_off_new(Xt_xmap xmap, int *src_offsets,
                                int *dst_offsets, MPI_Datatype datatype) {

  int nsend, nrecv;
  struct Xt_redist_msg * send_msgs = NULL;
  struct Xt_redist_msg * recv_msgs = NULL;
  MPI_Comm comm;

  {
    MPI_Comm xmap_comm = xt_xmap_get_communicator(xmap);
    xt_mpi_call(MPI_Comm_dup(xmap_comm, &comm), xmap_comm);
  }

  nrecv = xt_xmap_get_num_sources(xmap);
  Xt_xmap_iter dst_iter = xt_xmap_get_in_iterator(xmap);
  generate_msg_infos(nrecv, dst_iter, dst_offsets, datatype, &recv_msgs,
                     comm);
  if (dst_iter) xt_xmap_iterator_delete(dst_iter);

  nsend = xt_xmap_get_num_destinations(xmap);
  Xt_xmap_iter src_iter = xt_xmap_get_out_iterator(xmap);
  generate_msg_infos(nsend, src_iter, src_offsets, datatype, &send_msgs,
                     comm);
  if (src_iter) xt_xmap_iterator_delete(src_iter);

  return xt_redist_single_array_base_new(nsend, nrecv, send_msgs, recv_msgs,
                                         comm);
}

/* ====================================================================== */

struct ext_disp
{
  int disp, ext_idx;
};

static inline struct ext_disp
pos2disp(int pos, int num_ext, struct Xt_offset_ext extents[],
         int psum_ext_size[])
{
  int j = 0;
  /* FIXME: use bsearch if linear search is too slow, i.e. num_ext >> 1000 */
  /* what extent covers the pos'th position? */
  while (j < num_ext && pos >= psum_ext_size[j + 1])
    ++j;
  int disp = extents[j].start + (pos - psum_ext_size[j]) * extents[j].stride;
  return (struct ext_disp){ .disp = disp, .ext_idx = j };
}

static inline struct ext_disp
pos2disp2(int pos, int num_ext, struct Xt_offset_ext extents[],
          int psum_ext_size[], int start_ext)
{
  int j = start_ext;
  if (pos < psum_ext_size[j + 1] && pos >= psum_ext_size[j])
    ;
  else if (pos < psum_ext_size[j + 1])
  {
    j = 0;
    while (j < start_ext && pos >= psum_ext_size[j + 1])
      ++j;
  }
  else
    while (j < num_ext && pos >= psum_ext_size[j + 1])
      ++j;
  int disp = extents[j].start + (pos - psum_ext_size[j]) * extents[j].stride;
  return (struct ext_disp){ .disp = disp, .ext_idx = j };
}

static MPI_Datatype
generate_ext_datatype(int num_transfer_pos, const int transfer_pos[],
                      int num_ext, struct Xt_offset_ext extents[],
                      int psum_ext_size[],
                      MPI_Datatype base_datatype, MPI_Comm comm)
{
  if (num_transfer_pos > 0)
  {
    size_t blocks_size = 8, block_ofs = 0;
    int *displ = NULL, *blocklen = NULL;
    int i = 0;
    do
    {
      blocks_size *= 2;
      displ = xrealloc(displ, blocks_size * sizeof (*displ));
      blocklen = xrealloc(blocklen, blocks_size * sizeof (*blocklen));
      do
      {
        int j = 1;
        struct ext_disp pos = pos2disp(transfer_pos[i],
                                       num_ext, extents, psum_ext_size);
        int block_disp = pos.disp;
        while (i + j < num_transfer_pos
               && block_disp + j
               == (pos = pos2disp2(transfer_pos[i + j],
                                   num_ext, extents, psum_ext_size,
                                   pos.ext_idx)).disp)
          ++j;
        displ[block_ofs] = block_disp;
        blocklen[block_ofs] = j;
        i = i + j;
      } while (++block_ofs < blocks_size && i < num_transfer_pos);
    } while (i < num_transfer_pos);
    MPI_Datatype type
      = xt_mpi_generate_datatype_block(displ, blocklen, (int)block_ofs,
                                       base_datatype, comm);
    free(displ);
    free(blocklen);
    return type;
  }
  else
    return MPI_DATATYPE_NULL;
}

static void
generate_ext_msg_infos(int num_msgs, Xt_xmap_iter iter,
                       int num_ext,
                       struct Xt_offset_ext extents[],
                       MPI_Datatype base_datatype,
                       struct Xt_redist_msg **msgs,
                       MPI_Comm comm)
{
  if (num_msgs <= 0) {
    *msgs = NULL;
    return;
  }

  struct Xt_redist_msg *curr_msg =
    *msgs = xmalloc((size_t)num_msgs * sizeof(**msgs));

  /* partial sums of ext sizes */
  int *psum_ext_size
    = xmalloc(((size_t)num_ext + 1) * sizeof (psum_ext_size[0]));
  psum_ext_size[0] = 0;
  for (size_t i = 0; i < (size_t)num_ext; ++i)
    psum_ext_size[i + 1] = psum_ext_size[i] + extents[i].size;

  do {

    const int *curr_transfer_pos = xt_xmap_iterator_get_transfer_pos(iter);
    int curr_num_transfer_pos = xt_xmap_iterator_get_num_transfer_pos(iter);

    curr_msg->datatype
      = generate_ext_datatype(curr_num_transfer_pos, curr_transfer_pos,
                              num_ext, extents, psum_ext_size,
                              base_datatype, comm);
    curr_msg->rank = xt_xmap_iterator_get_rank(iter);

    curr_msg++;

  } while (xt_xmap_iterator_next(iter));
  free(psum_ext_size);
}

Xt_redist xt_redist_p2p_ext_new(Xt_xmap xmap,
                                int num_src_ext,
                                struct Xt_offset_ext src_extents[],
                                int num_dst_ext,
                                struct Xt_offset_ext dst_extents[],
                                MPI_Datatype datatype)
{
  int nsend, nrecv;
  struct Xt_redist_msg * send_msgs = NULL;
  struct Xt_redist_msg * recv_msgs = NULL;
  MPI_Comm comm;

  {
    MPI_Comm xmap_comm = xt_xmap_get_communicator(xmap);
    xt_mpi_call(MPI_Comm_dup(xmap_comm, &comm), xmap_comm);
  }

  nrecv = xt_xmap_get_num_sources(xmap);
  Xt_xmap_iter dst_iter = xt_xmap_get_in_iterator(xmap);
  generate_ext_msg_infos(nrecv, dst_iter, num_dst_ext, dst_extents,
                         datatype, &recv_msgs, comm);
  if (dst_iter) xt_xmap_iterator_delete(dst_iter);

  nsend = xt_xmap_get_num_destinations(xmap);
  Xt_xmap_iter src_iter = xt_xmap_get_out_iterator(xmap);
  generate_ext_msg_infos(nsend, src_iter, num_src_ext, src_extents,
                         datatype, &send_msgs, comm);
  if (src_iter) xt_xmap_iterator_delete(src_iter);

  return xt_redist_single_array_base_new(nsend, nrecv, send_msgs, recv_msgs,
                                         comm);
}

/* ====================================================================== */

static void
aux_gen_simple_block_offsets(int *block_offsets, int *block_sizes,
                             int block_num) {

  if (block_num<1) return;
  block_offsets[0] = 0;
  for (int i = 1; i < block_num; ++i) {
    block_offsets[i] = block_offsets[i-1] + block_sizes[i-1];
  }
}

static MPI_Datatype
generate_block_datatype(int const * transfer_pos, int num_transfer_pos,
                        int *block_offsets, int *block_sizes,
                        MPI_Datatype base_datatype, MPI_Comm comm) {

  MPI_Datatype type;

  int *bdispl_vec;
  int *blen_vec;

  assert(block_sizes != NULL);

  bdispl_vec = xmalloc(2 * (size_t)num_transfer_pos * sizeof(*bdispl_vec));
  blen_vec = bdispl_vec + num_transfer_pos;
  assert(block_offsets);

  for (int i = 0; i < num_transfer_pos; ++i) {
    int j = transfer_pos[i];
    bdispl_vec[i] = block_offsets[j];
    blen_vec[i] = block_sizes[j];
  }

  type = xt_mpi_generate_datatype_block(bdispl_vec, blen_vec,
                                        num_transfer_pos,
                                        base_datatype, comm);

  free(bdispl_vec);

  return type;
}

static void
generate_block_msg_infos(int num_msgs, Xt_xmap_iter iter, int *block_offsets,
                         int *block_sizes, MPI_Datatype base_datatype,
                         struct Xt_redist_msg ** msgs, MPI_Comm comm) {

  if (num_msgs <= 0) {
    *msgs = NULL;
    return;
  }

  *msgs = xmalloc((size_t)num_msgs * sizeof(**msgs));

  struct Xt_redist_msg * curr_msg;

  curr_msg = *msgs;

  do {

    int const * curr_transfer_pos;
    int curr_num_transfer_pos;

    curr_transfer_pos = xt_xmap_iterator_get_transfer_pos(iter);
    curr_num_transfer_pos = xt_xmap_iterator_get_num_transfer_pos(iter);

    curr_msg->datatype
      = generate_block_datatype(curr_transfer_pos, curr_num_transfer_pos,
                                block_offsets, block_sizes, base_datatype,
                                comm);
    curr_msg->rank = xt_xmap_iterator_get_rank(iter);

    curr_msg++;

  } while (xt_xmap_iterator_next(iter));
}

Xt_redist
xt_redist_p2p_blocks_off_new(Xt_xmap xmap,
                             int *src_block_offsets, int *src_block_sizes,
                             int src_block_num,
                             int *dst_block_offsets, int *dst_block_sizes,
                             int dst_block_num,
                             MPI_Datatype datatype) {

  int nsend, nrecv;
  struct Xt_redist_msg * send_msgs = NULL;
  struct Xt_redist_msg * recv_msgs = NULL;
  MPI_Comm comm;

  if (!src_block_sizes)
    die("xt_redist_p2p_blocks_off_new: undefined src_block_sizes");
  if (!dst_block_sizes)
    die("xt_redist_p2p_blocks_off_new: undefined dst_block_sizes");

  MPI_Comm xmap_comm;

  xmap_comm = xt_xmap_get_communicator(xmap);

  xt_mpi_call(MPI_Comm_dup(xmap_comm, &comm), xmap_comm);

  nsend = xt_xmap_get_num_destinations(xmap);
  nrecv = xt_xmap_get_num_sources(xmap);

  int *aux_offsets = NULL;

  Xt_xmap_iter dst_iter, src_iter;

  dst_iter = xt_xmap_get_in_iterator(xmap);
  src_iter = xt_xmap_get_out_iterator(xmap);

  // dst part:
  int max_dst_pos = xt_xmap_get_max_dst_pos(xmap);
  if (dst_block_num < max_dst_pos)
    die("xt_redist_p2p_blocks_off_new: dst_block_num too small");

  if (dst_block_offsets)
    aux_offsets = dst_block_offsets;
  else {
    aux_offsets = xmalloc((size_t)dst_block_num * sizeof(*aux_offsets));
    aux_gen_simple_block_offsets(aux_offsets, dst_block_sizes, dst_block_num);
  }

  generate_block_msg_infos(nrecv, dst_iter, aux_offsets, dst_block_sizes,
                           datatype, &recv_msgs, comm);

  if (!dst_block_offsets) free(aux_offsets);

  // src part:
  int max_src_pos = xt_xmap_get_max_src_pos(xmap);
  if (src_block_num < max_src_pos)
    die("xt_redist_p2p_blocks_off_new: src_block_num too small");

  if (src_block_offsets)
    aux_offsets = src_block_offsets;
  else {
    aux_offsets = xmalloc((size_t)src_block_num * sizeof(*aux_offsets));
    aux_gen_simple_block_offsets(aux_offsets, src_block_sizes, src_block_num);
  }

  generate_block_msg_infos(nsend, src_iter, aux_offsets, src_block_sizes,
                           datatype, &send_msgs, comm);

  if (!src_block_offsets) free(aux_offsets);

  if (dst_iter) xt_xmap_iterator_delete(dst_iter);
  if (src_iter) xt_xmap_iterator_delete(src_iter);

  return xt_redist_single_array_base_new(nsend, nrecv, send_msgs, recv_msgs,
                                         comm);
}

Xt_redist xt_redist_p2p_blocks_new(Xt_xmap xmap,
                                   int *src_block_sizes, int src_block_num,
                                   int *dst_block_sizes, int dst_block_num,
                                   MPI_Datatype datatype) {

  return xt_redist_p2p_blocks_off_new(xmap,
                                      NULL, src_block_sizes, src_block_num,
                                      NULL, dst_block_sizes, dst_block_num,
                                      datatype);

}


Xt_redist xt_redist_p2p_new(Xt_xmap xmap, MPI_Datatype datatype) {

  return xt_redist_p2p_off_new(xmap, NULL, NULL, datatype);
}
