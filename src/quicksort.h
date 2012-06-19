#ifndef QUICKSORT_H
#define QUICKSORT_H

typedef struct sort_idxpos_struct {
  Xt_idx idx;    //!< first criterion
  int pos;       //!< second criterion
} sort_idxpos_type;

/** quicksort changing structured values
  *
  * @param[in,out]  v            data to be sorted
  * @param[in]      n            length of data
 **/
void quicksort_idxpos(sort_idxpos_type *v, int n);

/** quicksort changing values and indices
  *
  * @param[in,out]  a            data to be sorted
  * @param[in]      n            length of data
  * @param[in,out]  idx          old index of sorted returned a
  * @param[in]      reset_index  override given idx by identity idx
  */
void quicksort_index (Xt_idx * a, int n, int * idx, int reset_index);

#endif
