#include "xt/xt_sort.h"

void (*xt_sort_index)(Xt_idx * a, int n, int * idx, int reset_index)
= xt_quicksort_index;
void (*xt_sort_idxpos)(idxpos_type *v, int n) = xt_quicksort_idxpos;

void xt_sort_init() {

  //! \todo getenv: key = XT_SORT_METHOD, possible values: { QUICKSORT, MERGESORT }; we might accept uppercase and lowercase values

  /* xt_sort_index  = ; */
  /* xt_sort_idxpos = ; */
}
