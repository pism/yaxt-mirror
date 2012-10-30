#ifndef MERGESORT_H
#define MERGESORT_H

#include "sort_common.h"

/** mergesort changing structured values
  *
  * @param[in,out]  v            data to be sorted
  * @param[in]      n            number of elements in v
 **/
void xt_mergesort_idxpos(idxpos_type *v, int n);

/** mergesort changing values and indices
  *
  * @param[in,out]  a            data to be sorted
  * @param[in]      n            length of data
  * @param[in,out]  idx          old index of sorted returned a
  * @param[in]      reset_index  override given idx by identity idx
  */
void xt_mergesort_index (Xt_idx * a, int n, int * idx, int reset_index);

#endif
