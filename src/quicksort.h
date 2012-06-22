#ifndef QUICKSORT_H
#define QUICKSORT_H

#include "sort_common.h"

/** quicksort changing structured values
  *
  * @param[in,out]  v            data to be sorted
  * @param[in]      n            length of data
 **/
void xt_quicksort_idxpos(idxpos_type *v, int n);

/** quicksort changing values and indices
  *
  * @param[in,out]  a            data to be sorted
  * @param[in]      n            length of data
  * @param[in,out]  idx          old index of sorted returned a
  * @param[in]      reset_index  override given idx by identity idx
  */
void xt_quicksort_index (Xt_idx * a, int n, int * idx, int reset_index);

#endif
