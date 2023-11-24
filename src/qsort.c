/**
 * @file qsort.c
 * @brief sort routines using C qsort
 *
 * @copyright Copyright  (C)  2021 Jörg Behrens <behrens@dkrz.de>
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

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>

#include "xt/xt_sort.h"
#include "xt/quicksort.h"
#include "xt/sort_common.h"
#include "core/ppm_xfuncs.h"

void xt_qsort_index (Xt_int * v_idx, int n, int * v_pos, int reset_pos) {

  /*  Initialization */
  int *v_pos_orig =  v_pos;

  if (!v_pos) v_pos = xmalloc((size_t)n  * sizeof(v_pos[0]));

  if (v_pos != v_pos_orig || reset_pos)
    xt_assign_id_map_int((size_t)n, v_pos, 0);

  xt_qsort_permutation(v_idx, (size_t)n, v_pos);
  if (v_pos != v_pos_orig) free(v_pos);
}

static int compare_idxpos_type(void const * a_, void const * b_) {

  idxpos_type const * a = (idxpos_type const *)a_;
  idxpos_type const * b = (idxpos_type const *)b_;

  int ret = (a->idx > b->idx) - (a->idx < b->idx);
  if (!ret) ret = a->pos - b->pos;

  return ret;
}

void xt_qsort_permutation(Xt_int *idx, size_t n, int *pos) {

  idxpos_type * idxpos = xmalloc(n * sizeof(*idxpos));

  for (size_t i = 0; i < n; ++i) {
    idxpos[i].idx = idx[i];
    idxpos[i].pos = pos[i];
  }

  xt_qsort_idxpos(idxpos, n);

  for (size_t i = 0; i < n; ++i) {
    idx[i] = idxpos[i].idx;
    pos[i] = idxpos[i].pos;
  }

  free(idxpos);
}

void xt_qsort_idxpos(idxpos_type *idxpos, size_t n) {

  qsort(idxpos, n, sizeof(*idxpos), compare_idxpos_type);
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
