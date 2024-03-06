/**
 * @file xt_idxsection_stripe_iterate.h
 *
 * @brief Loop over stripes to compute intersection. This source is
 * meant to be included for multiple dimensionalities to facilitate
 * better unrolling and vectorization of the innermost loop.
 *
 * @copyright Copyright  (C)  2024 Jörg Behrens <behrens@dkrz.de>
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

  for (size_t k = 0; k < num_stripes; ++k) {
    struct Xt_stripe query = stripes[k];
    size_t num_result_stripes = 0;
    bool any_match = false;
    for (int j = 0; j < query.nstrides; ++j) {
      Xt_int insert_index = (Xt_int)(query.start + query.stride * j),
        running_index = insert_index;
      bool out_of_bounds = false;
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
      int local_pos = 0;
#endif
      for (size_t i = 0; i < NUM_DIMENSIONS; ++i) {
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

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
