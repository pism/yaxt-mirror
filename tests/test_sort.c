#include <assert.h>

//#define VERBOSE

#include "tests.h"
#include "xt_sort.h"

void test1(void (*anysort_index)(Xt_idx * a, int n, int * idx, int reset_index)) {
  int tsize = 8;
  Xt_idx ivec[]        = {7,5,3,1,2,2,3,3};
  Xt_idx sorted_ivec[] = {1,2,2,3,3,3,5,7};
  int pvec[]           = {1,2,3,4,5,6,7,8};
  int sorted_pvec[]    = {3,4,5,2,6,7,1,0};

  tsize = sizeof(ivec) /sizeof(ivec[0]);
  assert(tsize == sizeof(pvec) /sizeof(pvec[0]));

  xt_sort_index(ivec, tsize, pvec, 1);
  for(int i=0; i<tsize; i++) {
    if ( ivec[i] != sorted_ivec[i] )
      PUT_ERR("wrong sorting values\n");
    if ( pvec[i] != sorted_pvec[i] )
      PUT_ERR("wrong sorting positions\n");
  }
}

int main(void) {

  test1(xt_quicksort_index);
  test1(xt_mergesort_index);
  test1(xt_sort_index);

  return TEST_EXIT_CODE;
}
