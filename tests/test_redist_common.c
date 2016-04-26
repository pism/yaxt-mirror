/**
 * @file test_redist_common.c
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

#include <string.h>

#include "yaxt.h"
#include "tests.h"
#include "test_redist_common.h"
#include "core/ppm_xfuncs.h"


/*
 * build xmap for destination list containing all odd elements of
 * source list dimensioned 1 to src_num_indices
 */
Xt_xmap
build_odd_selection_xmap(int src_num_indices)
{
  enum {
    selection_stride = 2,
  };
  if (src_num_indices < 0)
    PUT_ERR("error: src_num_indices < 0");
  Xt_int *index_list = xmalloc((size_t)src_num_indices
                               * sizeof (index_list[0]));
  for (int i = 0; i < src_num_indices; ++i)
    index_list[i] = (Xt_int)(i + 1);
  Xt_idxlist src_idxlist = xt_idxvec_new(index_list, src_num_indices);
  int dst_num_indices
    = (int)((src_num_indices + selection_stride - 1) / selection_stride);
  for (int i = 0; i < dst_num_indices; ++i)
    index_list[i] = (Xt_int)(i * selection_stride + 1);
  Xt_idxlist dst_idxlist = xt_idxvec_new(index_list, dst_num_indices);
  free(index_list);
  Xt_xmap xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);
  return xmap;
}

int communicators_are_congruent(MPI_Comm comm1, MPI_Comm comm2) {

  int result;

  xt_mpi_call(MPI_Comm_compare(comm1, comm2, &result), MPI_COMM_WORLD);

  return ((result == MPI_IDENT) || (result == MPI_CONGRUENT));
}

void
check_redist(Xt_redist redist, const void *src,
             size_t dst_size, void *dst, const void *ref_dst_data)
{
  void * tmp_dst = malloc(dst_size);
  for (int i = 0; i < 2; ++i) {
    memcpy(tmp_dst, dst, dst_size);
    if (i == 0) {
      xt_redist_s_exchange1(redist, src, tmp_dst);
    } else {
      Xt_request request = xt_redist_a_exchange1(redist, src, tmp_dst);
      xt_request_wait(&request);
    }
    if (memcmp(tmp_dst, ref_dst_data, dst_size))
      PUT_ERR("error in xt_redist_s_exchange\n");
  }
  free(tmp_dst);
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
