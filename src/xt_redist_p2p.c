/**
 * @file xt_redist_p2p.c
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>

#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/xt_redist_p2p.h"
#include "xt_redist_internal.h"
#include "xt/xt_redist_single_array_base.h"
#include "xt/xt_xmap.h"
#include "xt/xt_idxlist.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"
#include "xt_config_internal.h"

#include "xt_arithmetic_util.h"

#include "xt_mpi_ddt_cache.h"

/* the following two functions fullfil the same purpose as
 * xt_disp2ext and xt_disp2ext_count but work with an indirection  */
static size_t
xt_mdisp2ext_count(size_t disp_len, const int *disp, const int *pos)
{
  if (!disp_len) return 0;
  size_t i = 0;
  int cur_stride = 1, cur_size = 1;
  int last_disp = disp[pos[0]];
  for (size_t p = 1; p < disp_len; ++p) {
    int new_disp = disp[pos[p]];
    int new_stride = new_disp - last_disp;
    if (cur_size == 1) {
      cur_stride = new_stride;
      cur_size = 2;
    } else if (new_stride == cur_stride) {
      // cur_size >= 2:
      cur_size++;
    } else if (cur_size > 2 || (cur_size == 2 && cur_stride == 1) ) {
      // we accept small contiguous vectors (nstrides==2, stride==1)
      i++;
      cur_stride = 1;
      cur_size = 1;
    } else { // cur_size == 2, next offset doesn't match current stride
      // break up trivial vec:
      i++;
      cur_size = 2;
      cur_stride = new_stride;
    }
    last_disp = new_disp;
  }
  // tail cases:
  if (cur_size > 2 || (cur_size == 2 && cur_stride == 1)) {
    i++;
  } else if (cur_size == 2) {
    i+=2;
  } else { // cur_size == 1
    i++;
  }

  return i;
}

static size_t
xt_mdisp2ext(size_t disp_len, const int *disp, const int *pos,
             int *restrict v_disp, int *restrict v_size, int *restrict v_stride)
{
  if (disp_len<1) return 0;

  int cur_start = disp[pos[0]], cur_stride = 1, cur_size = 1;
  int last_disp = cur_start;
  size_t i = 0;
  for (size_t p = 1; p < disp_len; ++p) {
    int new_disp = disp[pos[p]];
    int new_stride = new_disp - last_disp;
    if (cur_size == 1) {
      cur_stride = new_stride;
      cur_size = 2;
    } else if (new_stride == cur_stride) {
      // cur_size >= 2:
      cur_size++;
    } else if (cur_size > 2 || (cur_size == 2 && cur_stride == 1) ) {
      // we accept small contiguous vectors (nstrides==2, stride==1)
      v_disp[i] = cur_start;
      v_stride[i] = cur_stride;
      v_size[i] = cur_size;
      i++;
      cur_start = new_disp;
      cur_stride = 1;
      cur_size = 1;
    } else { // cur_size == 2, next offset doesn't match current stride
      // break up trivial vec:
      v_disp[i] = cur_start;
      v_size[i] = 1;
      v_stride[i] = 1;
      i++;
      cur_start += cur_stride;
      cur_size = 2;
      cur_stride = new_stride;
    }
    last_disp = new_disp;
  }
  // tail cases:
  v_disp[i] = cur_start;
  if (cur_size != 2 || cur_stride == 1) {
    v_size[i] = cur_size;
    v_stride[i] = cur_stride;
  } else {
    v_size[i] = 1;
    v_stride[i] = 1;
    i++;
    v_disp[i] = cur_start + cur_stride;
    v_size[i] = 1;
    v_stride[i] = 1;
  }
  i++;
  return i;
}

static size_t
xt_mdisp2blk_count(size_t disp_len, const int *disp, const int *pos)
{
  size_t count = 0, i = 0;
  while (i < disp_len) {
    size_t j = 1;
    int blk_strt_pos = pos[disp[i]];
    while (i + j < disp_len && blk_strt_pos + (int)j == pos[disp[i + j]])
      ++j;
    ++count;
    i += j;
  } while (i < disp_len);
  return count;
}

static size_t
xt_mdisp2blk(size_t disp_len, const int *disp, const int *pos,
             int *restrict blk_disp, int *restrict blk_len)
{
  size_t new_count = 0;
  size_t i = 0;
  while (i < disp_len) {
    size_t j = 1;
    int blk_strt_pos = pos[disp[i]];
    while (i + j < disp_len && blk_strt_pos + (int)j == pos[disp[i + j]])
      ++j;
    blk_disp[new_count] = blk_strt_pos;
    blk_len[new_count++] = (int)j;
    i += j;
  }
  return new_count;
}


static void
xt_dt_stripe_bufs_ext_grow(
  size_t *work_buf_size,
  void **work_buf,
  int **dt_stripes_disp,
  int **restrict dt_stripes_size,
  int **restrict dt_stripes_stride,
  size_t num_current_elems,
  size_t needed_elems);

static inline size_t
xt_dt_stripe_alloc_size_ext(size_t needed_elems);


static MPI_Datatype
generate_datatype(const int *transfer_pos, int num_transfer_pos,
                  const int *offsets,
                  size_t *work_buf_size, void **work_buf,
                  struct Xt_mpi_strp_prs_params *params,
                  Xt_config config)
{
  MPI_Datatype dt;
  size_t num_transfer_pos_ = num_transfer_pos >= 0 ? (size_t)num_transfer_pos : 0;
  if (config->ddt_cache_limit >= -1) {
    int *v_disp = *work_buf, *v_size = *work_buf, *v_stride = *work_buf;
    size_t vlen = offsets != NULL
      ? xt_mdisp2ext_count(num_transfer_pos_, offsets, transfer_pos)
      : xt_disp2ext_count(num_transfer_pos_, transfer_pos);
    if (xt_dt_stripe_alloc_size_ext(vlen) > *work_buf_size) {
      xt_dt_stripe_bufs_ext_grow(
        work_buf_size, work_buf, &v_disp, &v_size, &v_stride,
        0, vlen);
    } else {
      v_disp = *work_buf;
      v_stride = v_disp + vlen;
      v_size = v_stride + vlen;
    }
    if (offsets != NULL)
      xt_mdisp2ext(num_transfer_pos_, offsets, transfer_pos,
                   v_disp, v_size, v_stride);
    else
      xt_disp2ext(num_transfer_pos_, transfer_pos,
                  v_disp, v_size, v_stride);
    dt = xt_mpi_parse_stripe(vlen, v_disp, v_size, v_stride, params);
  } else {
    size_t num_blk = offsets != NULL
      ? xt_mdisp2blk_count(num_transfer_pos_, transfer_pos, offsets)
      : xt_disp2blk_count(num_transfer_pos_, transfer_pos);
    int *restrict blocklengths = *work_buf;
    size_t needed_buf_size = num_blk * sizeof(*blocklengths)
      * (1U + ((num_blk != num_transfer_pos_) || offsets != NULL));
    if (needed_buf_size > *work_buf_size) {
      *work_buf = blocklengths = xrealloc(blocklengths, needed_buf_size);
      *work_buf_size = needed_buf_size;
    }
    const int *displ;
    if (num_blk != num_transfer_pos_ || offsets != NULL) {
      if (offsets != NULL)
        xt_mdisp2blk(num_transfer_pos_, transfer_pos, offsets,
                     blocklengths + num_blk, blocklengths);
      else
        xt_disp2blk(num_transfer_pos_, transfer_pos,
                    blocklengths + num_blk, blocklengths);
      displ = blocklengths + num_blk;
    } else {
      for (size_t i = 0; i < num_blk; ++i)
        blocklengths[i] = 1;
      displ = transfer_pos;
    }
    dt = xt_mpi_generate_datatype_block_custom(
      displ, blocklengths, (int)num_blk, params->old_type,
      params->comm, config);
  }
  return dt;
}

static void
generate_msg_infos(int num_msgs, Xt_xmap_iter iter, const int *offsets,
                   struct Xt_redist_msg *msgs,
                   struct Xt_mpi_strp_prs_params *params,
                   Xt_config config) {

  if (num_msgs > 0) {
    size_t work_buf_size = 0;
    void *work_buf = NULL;
    struct Xt_redist_msg *restrict curr_msg = msgs;
    do {

      const int *curr_transfer_pos = xt_xmap_iterator_get_transfer_pos(iter);
      int curr_num_transfer_pos = xt_xmap_iterator_get_num_transfer_pos(iter);

      curr_msg->datatype
        = generate_datatype(curr_transfer_pos, curr_num_transfer_pos,
                            offsets, &work_buf_size, &work_buf, params, config);
      curr_msg->rank = xt_xmap_iterator_get_rank(iter);

      curr_msg++;
    } while (xt_xmap_iterator_next(iter));
    free(work_buf);
  }
}

Xt_redist
xt_redist_p2p_off_new(Xt_xmap xmap, const int *src_offsets,
                      const int *dst_offsets, MPI_Datatype datatype)
{
  return xt_redist_p2p_off_custom_new(xmap, src_offsets, dst_offsets, datatype,
                                      (Xt_config)&xt_default_config);
}

Xt_redist
xt_redist_p2p_off_custom_new(Xt_xmap xmap, const int *src_offsets,
                             const int *dst_offsets, MPI_Datatype datatype,
                             Xt_config config) {
  // ensure that yaxt is initialized
  assert(xt_initialized());

  int nsend = xt_xmap_get_num_destinations(xmap),
    nrecv = xt_xmap_get_num_sources(xmap);
  size_t nmsg = (size_t)nsend + (size_t)nrecv;
  struct Xt_redist_msg *msgs = xmalloc(nmsg * sizeof (*msgs)),
    *send_msgs = msgs, *recv_msgs = msgs + nsend;
  int tag_offset;
  MPI_Comm comm
    = xt_mpi_comm_smart_dup(xt_xmap_get_communicator(xmap), &tag_offset);

  struct Xt_mpi_strp_prs_params params;
  xt_init_mpi_strp_prs_params(&params, datatype, config, comm);
  Xt_xmap_iter dst_iter = xt_xmap_get_in_iterator(xmap);
  generate_msg_infos(nrecv, dst_iter, dst_offsets, recv_msgs, &params, config);
  if (dst_iter) xt_xmap_iterator_delete(dst_iter);

  Xt_xmap_iter src_iter = xt_xmap_get_out_iterator(xmap);
  generate_msg_infos(nsend, src_iter, src_offsets, send_msgs, &params, config);
  if (src_iter) xt_xmap_iterator_delete(src_iter);
  Xt_mpi_ddt_cache_check_retention(&params.ddt_list, nmsg, msgs);
  xt_destroy_mpi_strp_prs_params(&params);

  struct Xt_config_ config_ = *config;
  config_.flags |= exch_no_dt_dup;

  Xt_redist result = xt_redist_single_array_base_custom_new(
    nsend, nrecv, send_msgs, recv_msgs, comm, &config_);

  free(msgs);
  xt_mpi_comm_smart_dedup(&comm, tag_offset);
  return result;
}

/* ====================================================================== */

static inline int
pos2disp(int pos, int num_ext, const int psum_ext_size[])
{
  int j = 0;
  /* FIXME: use bsearch if linear search is too slow, i.e. num_ext >> 1000 */
  /* what extent covers the pos'th position? */
  while (j < num_ext && pos >= psum_ext_size[j + 1])
    ++j;
  return j;
}

static inline int
pos2disp2(int pos, int num_ext,
          const int psum_ext_size[], int start_ext)
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
  return j;
}

#define XT_OFS_TYPE int
#define XT_EXT_TYPE struct Xt_offset_ext
#define XT_EXT_TAG ext
#define XT_MPI_PARSE_STRIPE xt_mpi_parse_stripe
#define XT_EXT_STRIDE_MASK isign_mask_current_pos_ext_size
#define XT_EXT_STRIDE_MASK_PREP
#include "xt_redist_p2p_ext.h"
#undef XT_OFS_TYPE
#undef XT_EXT_TYPE
#undef XT_EXT_TAG
#undef XT_MPI_PARSE_STRIPE
#undef XT_EXT_STRIDE_MASK
#undef XT_EXT_STRIDE_MASK_PREP

#define XT_OFS_TYPE MPI_Aint
#define XT_EXT_TYPE struct Xt_aoffset_ext
#define XT_EXT_TAG aext
#define XT_MPI_PARSE_STRIPE xt_mpi_parse_astripe
#define XT_EXT_STRIDE_MASK asign_mask_current_pos_ext_size
#define XT_EXT_STRIDE_MASK_PREP MPI_Aint asign_mask_current_pos_ext_size \
  = asign_mask(current_pos_ext.size)
#include "xt_redist_p2p_ext.h"
#undef XT_OFS_TYPE
#undef XT_EXT_TYPE
#undef XT_EXT_TAG
#undef XT_MPI_PARSE_STRIPE
#undef XT_EXT_STRIDE_MASK
#undef XT_EXT_STRIDE_MASK_PREP


/* ====================================================================== */

static inline void
aux_gen_simple_block_offsets(int block_offsets[], const int block_sizes[],
                             size_t num_blocks) {

  if (num_blocks > 0) {
    int accum = 0;
    for (size_t i = 0; i < num_blocks; ++i) {
      block_offsets[i] = accum;
      accum += block_sizes[i];
    }
  }
}

static MPI_Datatype
generate_block_datatype(const int *transfer_pos, int num_transfer_pos,
                        const int *block_offsets, const int *block_sizes,
                        struct Xt_mpi_strp_prs_params *params,
                        Xt_config config) {

  assert(block_sizes && block_offsets);

  bool use_stripe_parse = config->ddt_cache_limit >= -1;
  int *bdispl_vec = xmalloc((size_t)(2 + use_stripe_parse)
                            * (size_t)num_transfer_pos * sizeof(*bdispl_vec)),
    *blen_vec = bdispl_vec + num_transfer_pos,
    *stride_vec = blen_vec + num_transfer_pos;

  for (int i = 0; i < num_transfer_pos; ++i) {
    int j = transfer_pos[i];
    bdispl_vec[i] = block_offsets[j];
    blen_vec[i] = block_sizes[j];
  }

  MPI_Datatype dt;
  if (use_stripe_parse) {
    for (int i = 0; i < num_transfer_pos; ++i)
      stride_vec[i] = 1;
    dt = xt_mpi_parse_stripe(
      (size_t)num_transfer_pos, bdispl_vec, blen_vec, stride_vec, params);
  } else
    dt = xt_mpi_generate_datatype_block_custom(
      bdispl_vec, blen_vec, num_transfer_pos,
      params->old_type, params->comm, config);

  free(bdispl_vec);

  return dt;
}

static void
generate_block_msg_infos(int num_msgs, Xt_xmap_iter iter,
                         const int *block_offsets,
                         const int *block_sizes, int **aux_offsets,
                         size_t num_blocks,
                         struct Xt_mpi_strp_prs_params *params,
                         Xt_config config,
                         struct Xt_redist_msg *msgs) {

  if (num_msgs > 0) {

    const int *block_offsets_;
    if (block_offsets)
      block_offsets_ = block_offsets;
    else {
      block_offsets_ = *aux_offsets
        = xrealloc(*aux_offsets, num_blocks * sizeof(*block_offsets_));
      aux_gen_simple_block_offsets(*aux_offsets, block_sizes, num_blocks);
    }

    size_t ofs = 0;
    do {
      const int *curr_transfer_pos = xt_xmap_iterator_get_transfer_pos(iter);
      int curr_num_transfer_pos = xt_xmap_iterator_get_num_transfer_pos(iter);
      msgs[ofs].datatype
        = generate_block_datatype(curr_transfer_pos, curr_num_transfer_pos,
                                  block_offsets_, block_sizes, params, config);
      msgs[ofs].rank = xt_xmap_iterator_get_rank(iter);

      ofs++;
    } while (xt_xmap_iterator_next(iter));

  }
}

Xt_redist
xt_redist_p2p_blocks_off_new(Xt_xmap xmap,
                             const int *src_block_offsets,
                             const int *src_block_sizes,
                             int src_block_num,
                             const int *dst_block_offsets,
                             const int *dst_block_sizes,
                             int dst_block_num,
                             MPI_Datatype datatype)
{
  return xt_redist_p2p_blocks_off_custom_new(
    xmap, src_block_offsets, src_block_sizes, src_block_num, dst_block_offsets,
    dst_block_sizes, dst_block_num, datatype, (Xt_config)&xt_default_config);
}


Xt_redist
xt_redist_p2p_blocks_off_custom_new(Xt_xmap xmap,
                                    const int *src_block_offsets,
                                    const int *src_block_sizes,
                                    int src_block_num,
                                    const int *dst_block_offsets,
                                    const int *dst_block_sizes,
                                    int dst_block_num,
                                    MPI_Datatype datatype,
                                    Xt_config config)
{
  // ensure that yaxt is initialized
  assert(xt_initialized() && src_block_sizes && dst_block_sizes);

  int tag_offset;
  MPI_Comm comm
    = xt_mpi_comm_smart_dup(xt_xmap_get_communicator(xmap), &tag_offset);


  int nsend = xt_xmap_get_num_destinations(xmap),
    nrecv = xt_xmap_get_num_sources(xmap);

  size_t nmsg = ((size_t)nsend + (size_t)nrecv);
  struct Xt_redist_msg *msgs = xmalloc(nmsg * sizeof (*msgs));

  int *aux_offsets = NULL;

  Xt_xmap_iter dst_iter = xt_xmap_get_in_iterator(xmap),
    src_iter = xt_xmap_get_out_iterator(xmap);

  // dst part:
#ifndef NDEBUG
  int max_dst_pos = xt_xmap_get_max_dst_pos(xmap);
  if (dst_block_num < max_dst_pos)
    die("xt_redist_p2p_blocks_off_new: dst_block_num too small");
#endif
  struct Xt_mpi_strp_prs_params params;
  xt_init_mpi_strp_prs_params(&params, datatype, config, comm);
  generate_block_msg_infos(nrecv, dst_iter, dst_block_offsets, dst_block_sizes,
                           &aux_offsets, (size_t)dst_block_num,
                           &params, config, msgs);
  if (dst_iter) xt_xmap_iterator_delete(dst_iter);

  // src part:
#ifndef NDEBUG
  int max_src_pos = xt_xmap_get_max_src_pos(xmap);
  if (src_block_num < max_src_pos)
    die("xt_redist_p2p_blocks_off_new: src_block_num too small");
#endif
  generate_block_msg_infos(nsend, src_iter, src_block_offsets, src_block_sizes,
                           &aux_offsets, (size_t)src_block_num,
                           &params, config, msgs+nrecv);
  xt_destroy_mpi_strp_prs_params(&params);
  free(aux_offsets);

  if (src_iter) xt_xmap_iterator_delete(src_iter);

  struct Xt_config_ config_ = *config;
  config_.flags |= exch_no_dt_dup;

  Xt_redist result
    = xt_redist_single_array_base_custom_new(
      nsend, nrecv, msgs+nrecv, msgs, comm, &config_);

  free(msgs);
  xt_mpi_comm_smart_dedup(&comm, tag_offset);
  return result;
}

Xt_redist xt_redist_p2p_blocks_new(Xt_xmap xmap,
                                   const int *src_block_sizes,
                                   int src_block_num,
                                   const int *dst_block_sizes,
                                   int dst_block_num,
                                   MPI_Datatype datatype)
{
  return xt_redist_p2p_blocks_custom_new(
    xmap, src_block_sizes, src_block_num, dst_block_sizes, dst_block_num,
    datatype, (Xt_config)&xt_default_config);
}

Xt_redist
xt_redist_p2p_blocks_custom_new(Xt_xmap xmap,
                                const int *src_block_sizes, int src_block_num,
                                const int *dst_block_sizes, int dst_block_num,
                                MPI_Datatype datatype,
                                Xt_config config)
{
  return xt_redist_p2p_blocks_off_custom_new(
    xmap, NULL, src_block_sizes, src_block_num,
    NULL, dst_block_sizes, dst_block_num, datatype, config);
}


Xt_redist xt_redist_p2p_new(Xt_xmap xmap, MPI_Datatype datatype)
{
  return xt_redist_p2p_custom_new(xmap, datatype,
                                  (Xt_config)&xt_default_config);
}

Xt_redist
xt_redist_p2p_custom_new(Xt_xmap xmap, MPI_Datatype datatype, Xt_config config)
{
  return xt_redist_p2p_off_custom_new(xmap, NULL, NULL, datatype, config);
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
