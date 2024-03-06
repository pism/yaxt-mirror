/**
 * @file xt_idxsection_stripe_intersection.h
 * Contains code parts independent of whether section
 * is destination or source side of intersection computation.
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
// intersection between an idxsection and a set of stripes
size_t num_stripes = xt_idxstripes_get_num_index_stripes(idxstripes_list);

const struct Xt_stripe *stripes
  = xt_idxstripes_get_index_stripes_const(idxstripes_list);
size_t num_dimensions = (size_t)idxsection->ndim;
#ifndef __PGI
struct {
  struct Xt_idxsection_ h;
  struct dim_desc d[num_dimensions];
} sss;                        /* sorted section stack copy object */
idxsection_init_sorted_copy(idxsection, &sss.h);
struct dim_desc *dims = sss.d;
#else
struct Xt_idxsection_ *sssp
= xmalloc(sizeof (struct Xt_idxsection_)
          + sizeof (struct dim_desc) * num_dimensions);
idxsection_init_sorted_copy(idxsection, sssp);
struct dim_desc *dims = sssp->dims;
#endif
enum {
  size_t_bits = sizeof (size_t) * CHAR_BIT,
};
#ifndef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
size_t *restrict isect_stripes_bm
  = xcalloc((num_stripes+size_t_bits-1)/size_t_bits, sizeof (size_t));
#else
size_t num_indices_in_section
  = (size_t)idxsection->parent.num_indices,
  section_bm_size = (num_indices_in_section+size_t_bits-1)/size_t_bits,
  *restrict pos_used_bm = xcalloc(section_bm_size
                         + (num_stripes+size_t_bits-1)/size_t_bits,
                         sizeof (size_t)),
  *restrict isect_stripes_bm = pos_used_bm + section_bm_size;
#endif
size_t num_result_stripes_total = 0;
{
  /* start with non-appendable terminator */
  Xt_int accum_start = XT_INT_MIN;
  int buf_nstrides = 0;
  Xt_int buf_stride = -1;
  Xt_int index_continuation = XT_INT_MIN-1;
  for (size_t k = 0; k < num_stripes; ++k) {
    struct Xt_stripe query = stripes[k];
    size_t num_result_stripes = 0;
    for (int j = 0; j < query.nstrides; ++j) {
      Xt_int insert_index = (Xt_int)(query.start + query.stride * j),
        running_index = insert_index;
      bool out_of_bounds = false;
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
      int local_pos = 0;
#endif
      for (size_t i = 0; i < num_dimensions; ++i) {
        Xt_int abs_global_stride = dims[i].global_stride;
        Xt_int curr_global_position
          = (Xt_int)(running_index / abs_global_stride);
        running_index = (Xt_int)(running_index % abs_global_stride);
        out_of_bounds |= (curr_global_position < dims[i].local_start)
          | (curr_global_position >= dims[i].local_start + dims[i].local_size);
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
        int curr_local_pos = (int)(curr_global_position - dims[i].local_start);
        local_pos += curr_local_pos * (int)dims[i].local_stride;
#endif
      }
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
      out_of_bounds = out_of_bounds ||
        (pos_used_bm[local_pos/size_t_bits] >> (local_pos % size_t_bits))&1;
#endif
      if (!out_of_bounds) {
        any_match = true;
        if (insert_index == index_continuation || buf_nstrides == 1) {
          /* already existing stripe can be continued */
          ++buf_nstrides;
          if (insert_index != index_continuation)
            buf_stride = (Xt_int)(insert_index - accum_start);
        } else {
          buf_stride = 1;
          buf_nstrides = 1;
          accum_start = insert_index;
          ++num_result_stripes;
        }
        index_continuation = (Xt_int)(insert_index + buf_stride);
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
        pos_used_bm[local_pos/size_t_bits]
          |= ((size_t)1 << (local_pos % size_t_bits));
#endif
      }
    }
    num_result_stripes_total += num_result_stripes;
    isect_mask |= (size_t)any_match << (k % size_t_bits);
    if ((k+1) % size_t_bits == 0) {
      isect_stripes_bm[k/size_t_bits] = isect_mask;
      isect_mask = 0;
    }
  }
  if (num_stripes % size_t_bits)
    isect_stripes_bm[num_stripes/size_t_bits] = isect_mask;
}
Xt_idxlist result;
if (num_result_stripes_total) {
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
  memset(pos_used_bm, 0, (num_indices_in_section+size_t_bits-1)/size_t_bits *
         sizeof (size_t));
#endif
  struct Xt_stripes_alloc result_alloc
    = xt_idxstripes_alloc(num_result_stripes_total);
  struct Xt_stripe *restrict result_stripes = result_alloc.stripes;
  size_t num_result_stripes = (size_t)-1;
  Xt_int index_continuation = XT_INT_MIN-1;
  int buf_nstrides = 0;
  Xt_int buf_stride = -1;
  for (size_t k = 0; k < num_stripes; ++k) {
    if (isect_stripes_bm[k/size_t_bits] & ((size_t)1 << (k % size_t_bits)))
      ; else continue;
    struct Xt_stripe query = stripes[k];
    for (int j = 0; j < query.nstrides; ++j) {
      Xt_int insert_index = (Xt_int)(query.start + query.stride * j),
        running_index = insert_index;
      bool out_of_bounds = false;
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
      int local_pos = 0;
#endif
      for (size_t i = 0; i < num_dimensions; ++i) {
        Xt_int abs_global_stride = dims[i].global_stride;
        Xt_int curr_global_position
          = (Xt_int)(running_index / abs_global_stride);
        running_index = (Xt_int)(running_index % abs_global_stride);
        out_of_bounds |= (curr_global_position < dims[i].local_start)
          | (curr_global_position >= dims[i].local_start + dims[i].local_size);
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
        int curr_local_pos = (int)(curr_global_position - dims[i].local_start);
        local_pos += curr_local_pos * (int)dims[i].local_stride;
#endif
      }
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
      out_of_bounds = out_of_bounds ||
        (pos_used_bm[local_pos/size_t_bits] >> (local_pos % size_t_bits))&1;
#endif
      if (!out_of_bounds) {
        if (insert_index == index_continuation || buf_nstrides == 1) {
          /* already existing stripe can be expanded */
          ++buf_nstrides;
          if (insert_index != index_continuation)
            buf_stride = (Xt_int)(insert_index
                                  - result_stripes[num_result_stripes].start);
        } else {
          result_stripes[num_result_stripes].nstrides = buf_nstrides;
          result_stripes[num_result_stripes].stride = buf_stride;
          buf_stride = 1;
          buf_nstrides = 1;
          result_stripes[++num_result_stripes].start = insert_index;
        }
        index_continuation = (Xt_int)(insert_index + buf_stride);
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
        pos_used_bm[local_pos/size_t_bits]
          |= ((size_t)1 << (local_pos % size_t_bits));
#endif
      }
    }
  }
  result_stripes[num_result_stripes].stride = buf_stride;
  result_stripes[num_result_stripes].nstrides = buf_nstrides;
  assert(num_result_stripes_total == num_result_stripes+1);
  Xt_idxlist initial_result = xt_idxstripes_congeal(result_alloc);
  if (XT_CONFIG_GET_FORCE_NOSORT(config)
      || xt_idxlist_get_sorting(initial_result) == 1)
    result = initial_result;
  else {
    result = xt_idxlist_sorted_copy(initial_result);
    xt_idxlist_delete(initial_result);
  }
} else
  result = xt_idxempty_new();
#ifdef __PGI
free(sssp);
#endif
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
free(pos_used_bm);
#else
free(isect_stripes_bm);
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
