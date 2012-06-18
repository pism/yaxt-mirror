#ifndef QUICKSORT_H
#define QUICKSORT_H

typedef struct sort_idxpos_struct {
  Xt_idx idx; // first crit
  int pos; // second crit
} sort_idxpos_type;

void quicksort_idxpos(sort_idxpos_type  *v, int n);

void quicksort_index (Xt_idx * a, int n, int * idx, int reset_index);

#endif
