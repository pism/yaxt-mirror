#include "yaxt.h"

#include <stdlib.h>

/** \file quicksort.c
  * \brief Non-recursive stack version of Quicksort
  *
  *   ... from N. Wirth's Pascal Book, 'Algorithms + Data Structures = Programms'.
  *       by Alan Miller ( 19 Jul 1995 )
  *
  *  taken from:
  * - http://www.nag.com/nagware/examples.asp
  * - http://www.nag.com/nagware/Examples/nur.f90
  *
  *  see also:
  * - http://en.wikipedia.org/wiki/Quicksort
  *
  * @param[in,out]  a            data to be sorted
  * @param[in]      n            length of data
  * @param[out]     idx          old index of sorted returned a
  *
  * \todo consider changing comparision to include position in order
  *       to get a conservative sorting
 **/

void quicksort_index (Xt_idx * a, int n, int * idx) {

  int pointer_inc = 64;
  int pointer_size;

  int * stackl;
  int * stackr;

  int i, j, k, l, r, s;
  int ww;
  Xt_idx w, x;

  /*  Initialization */

  pointer_size = pointer_inc;

  stackl = malloc ( pointer_size * sizeof(stackl[0]) );
  stackr = malloc ( pointer_size * sizeof(stackr[0]) );

  s = 0;
  stackl[0] = 0;
  stackr[0] = n-1;

  /*  Start sorting
            
  ... keep taking the top request from the stack until s = 0.

  */

  while  ( s > -1 ) {

    l = stackl[s];
    r = stackr[s];
    s = s - 1;

    /* ... keep splitting a[l], ... ,a[r] until l>= r. */

    while ( l < r ) {

      i = l;
      j = r;
      k = (l+r) / 2;
      x = a[k];

      /* Search from lower end */

      for (;;) {
        for (;;) {
          if ( a[i] < x ) {
            i = i + 1;
          }
          else break;
        }

        /* Search from upper end */

        for (;;) {
          if ( x < a[j] ) {
            j = j - 1;
          }
          else break;
        }

        /* Swap positions i & j */

        if ( i <= j ) {

        w      = a[i];
        a[i]   = a[j];        
        a[j]   = w;
	//! \todo do not allow idx == NULL, remove cond here
        if (idx != NULL) {
           ww     = idx[i];
           idx[i] = idx[j];
           idx[j] = ww;
        }

        i = i + 1;
        j = j - 1;

        if ( i > j ) break;
      }
      else break;
      }

      if ( j-l >= r-i ) {

      if ( l < j ) {

      if ( s > pointer_size ) {
        pointer_size = pointer_size + pointer_inc;
        stackl = realloc (stackl, pointer_size * sizeof(stackl[0]) );
        stackr = realloc (stackr, pointer_size * sizeof(stackr[0]) );
      }

      s = s + 1;
      stackl[s] = l;
      stackr[s] = j;

      }

      l = i;

      }

      else {

      if ( i < r ) {

      if ( s > pointer_size ) {
        pointer_size = pointer_size + pointer_inc;
        stackl = realloc (stackl, pointer_size * sizeof(stackl[0]) );
        stackr = realloc (stackr, pointer_size * sizeof(stackr[0]) );
      }

      s = s + 1;
      stackl[s] = i;
      stackr[s] = r;

      }

      r = j;

      }

    } /* ( l < r ) */

  } /* ( s /= 0 ) */

  free(stackl);
  free(stackr);
}
