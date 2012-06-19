#include "yaxt.h"

#include <stdlib.h>

#include "quicksort.h"

/** \file quicksort.c
  * \brief Non-recursive stack version of Quicksort
  *
  *       based on N. Wirth's Pascal Book, 'Algorithms + Data Structures = Programms'.
  *       by Alan Miller ( 19 Jul 1995 )
  *
  *  based on:
  * - http://www.nag.com/nagware/examples.asp
  * - http://www.nag.com/nagware/Examples/nur.f90
  *
  *  see also:
  * - http://en.wikipedia.org/wiki/Quicksort
  */

void quicksort_idxpos(sort_idxpos_type  *v, int n) {

  int pointer_inc = 64;
  int pointer_size;

  int * stackl;
  int * stackr;

  int i, j, k, l, r, s;
  sort_idxpos_type w, x;
 
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
      x = v[k];

      /* Search from lower end */

      while(1) {
        while(1)
          if ((v[i].idx < x.idx) || (v[i].idx == x.idx && v[i].pos < x.pos))
            i = i + 1;
          else
            break;

        /* Search from upper end */

        while(1)
          if ((x.idx < v[j].idx) || (x.idx == v[j].idx && x.pos < v[j].pos))
            j = j - 1;
          else
            break;

        /* Swap positions i & j */

        if ( i <= j ) {

          w      = v[i];
          v[i]   = v[j];
          v[j]   = w;

          i = i + 1;
          j = j - 1;

          if ( i > j ) break;
        } else
          break;
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
      } else {

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

void quicksort_index (Xt_idx * a, int n, int * idx, int reset_index) {

  //! \todo: consider using long long instead of sort_idxpos_type
  //!        I believe we should check for available sorting libs within configure, too

  int i;
  sort_idxpos_type *v;
  v = malloc ( n  * sizeof(v[0]) );
 
  if (idx && !reset_index) {
    for(i=0; i<n; i++) {
      v[i].idx = a[i];
      v[i].pos = idx[i];
    }
  } else {
    for(i=0; i<n; i++) {
      v[i].idx = a[i];
      v[i].pos = i;
    }
  }

  quicksort_idxpos(v,n);

  if (idx) {
    for(i=0; i<n; i++) {
      a[i] = v[i].idx;
      idx[i] = v[i].pos;
    }
  } else {
    for(i=0; i<n; i++) {
      a[i] = v[i].idx;
    }
  }

  free(v);
}
