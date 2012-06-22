#ifndef XT_SORT_H
#define XT_SORT_H

#include "quicksort.h"
#include "mergesort.h"

extern void (*xt_sort_index)(Xt_idx * a, int n, int * idx, int reset_index);
extern void (*xt_sort_idxpos)(idxpos_type *v, int n);

#endif // XT_SORT_H
