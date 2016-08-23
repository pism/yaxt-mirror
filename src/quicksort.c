/**
 * @file quicksort.c
 * @brief Non-recursive stack version of Quicksort
 *
 *  based on N. Wirth's Pascal Book, 'Algorithms + Data Structures = Programms'.
 *  by Alan Miller ( 19 Jul 1995 )
 *
 *  based on:
 * - http://www.nag.com/nagware/examples.asp
 * - http://www.nag.com/nagware/Examples/nur.f90
 *
 *  see also:
 * - http://en.wikipedia.org/wiki/Quicksort
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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
 * URL: https://doc.redmine.dkrz.de/yaxt/html/
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

#include "xt/quicksort.h"
#include "core/ppm_xfuncs.h"

void xt_quicksort_idxpos(idxpos_type  *v, int n) {

  enum { pointer_inc = 64 };
  unsigned pointer_size;

  int * stackl;
  int * stackr;

  int i, j, k, l, r;
  unsigned s;
  idxpos_type w, x;

  /*  Initialization */

  pointer_size = pointer_inc;

  stackl = xmalloc((size_t)pointer_size * sizeof(stackl[0]));
  stackr = xmalloc((size_t)pointer_size * sizeof(stackr[0]));

  s = 1;
  stackl[0] = 0;
  stackr[0] = n-1;

  /*  Start sorting

  ... keep taking the top request from the stack until s = 0.

  */

  while  ( s > 0 ) {

    --s;
    l = stackl[s];
    r = stackr[s];

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

          if ( s >= pointer_size ) {
            pointer_size = pointer_size + pointer_inc;
            stackl = xrealloc(stackl, (size_t)pointer_size * sizeof(stackl[0]));
            stackr = xrealloc(stackr, (size_t)pointer_size * sizeof(stackr[0]));
          }

          stackl[s] = l;
          stackr[s] = j;
          ++s;
        }
        l = i;
      } else {

        if ( i < r ) {

          if ( s >= pointer_size ) {
            pointer_size = pointer_size + pointer_inc;
            stackl = xrealloc(stackl, (size_t)pointer_size * sizeof(stackl[0]));
            stackr = xrealloc(stackr, (size_t)pointer_size * sizeof(stackr[0]));
          }

          stackl[s] = i;
          stackr[s] = r;
          ++s;
        }
        r = j;
      }
    } /* ( l < r ) */
  } /* ( s /= 0 ) */

  free(stackl);
  free(stackr);
}

void xt_quicksort_index(Xt_int * v_idx, int n, int * v_pos, int reset_pos) {

  enum { pointer_inc = 64 };
  unsigned pointer_size;

  int * stackl;
  int * stackr;

  int i, j, k, l, r;
  unsigned s;
  Xt_int w_idx, x_idx;
  int w_pos, x_pos;

  /*  Initialization */
  int *v_pos_orig =  v_pos;

  if (!v_pos) v_pos = malloc((size_t)n  * sizeof(v_pos[0]));

  if (v_pos != v_pos_orig || reset_pos) {
    for(i=0; i<n; i++) {
      v_pos[i] = i;
    }
  }

  pointer_size = pointer_inc;

  stackl = xmalloc((size_t)pointer_size * sizeof(stackl[0]));
  stackr = xmalloc((size_t)pointer_size * sizeof(stackr[0]));

  s = 1;
  stackl[0] = 0;
  stackr[0] = n-1;

  /*  Start sorting

  ... keep taking the top request from the stack until s = 0.

  */

  while  ( s > 0 ) {

    --s;
    l = stackl[s];
    r = stackr[s];

    /* ... keep splitting a[l], ... ,a[r] until l>= r. */

    while ( l < r ) {

      i = l;
      j = r;
      k = (l+r) / 2;
      x_idx = v_idx[k];
      x_pos = v_pos[k];

      /* Search from lower end */

      while(1) {
        while(1)
          if ((v_idx[i] < x_idx) || (v_idx[i] == x_idx && v_pos[i] < x_pos))
            i = i + 1;
          else
            break;

        /* Search from upper end */

        while(1)
          if ((x_idx < v_idx[j]) || (x_idx == v_idx[j] && x_pos < v_pos[j]))
            j = j - 1;
          else
            break;

        /* Swap positions i & j */

        if ( i <= j ) {

          w_idx      = v_idx[i];
          v_idx[i]   = v_idx[j];
          v_idx[j]   = w_idx;

          w_pos      = v_pos[i];
          v_pos[i]   = v_pos[j];
          v_pos[j]   = w_pos;

          i = i + 1;
          j = j - 1;

          if ( i > j ) break;
        } else
          break;
      }

      if ( j-l >= r-i ) {

        if ( l < j ) {

          if ( s >= pointer_size ) {
            pointer_size = pointer_size + pointer_inc;
            stackl = xrealloc(stackl, (size_t)pointer_size * sizeof(stackl[0]));
            stackr = xrealloc(stackr, (size_t)pointer_size * sizeof(stackr[0]));
          }

          stackl[s] = l;
          stackr[s] = j;
          ++s;
        }
        l = i;
      } else {

        if ( i < r ) {

          if ( s >= pointer_size ) {
            pointer_size = pointer_size + pointer_inc;
            stackl = xrealloc(stackl, (size_t)pointer_size * sizeof(stackl[0]));
            stackr = xrealloc(stackr, (size_t)pointer_size * sizeof(stackr[0]));
          }

          stackl[s] = i;
          stackr[s] = r;
          ++s;
        }
        r = j;
      }
    } /* ( l < r ) */
  } /* ( s /= 0 ) */

  free(stackl);
  free(stackr);
  if (v_pos != v_pos_orig) free(v_pos);
}



static inline int *med3_int(int *a, int *b, int *c)
{
  return *a < *b ?
      (*b < *c ? b : (*a < *c ? a : c ))
    : (*c < *b ? b : (*a < *c ? a : c ));
}

static inline void
vecswap_int(int *restrict a, int *restrict b, size_t n)
{
  for (size_t i = 0; i < n; ++i) {
    int t = a[i];
    a[i] = b[i];
    b[i] = t;
  }
}

/*
 * xt_quicksort_int was derived from the ScalES-PPM library code,
 * which is derived from genometools, which in turn is derived from
 * the FreeBSD libc.
 */
/*
  Modifications for integration with genometools
  2008 Thomas Jahns <Thomas.Jahns@gmx.net>

  The advertising clause 3. was removed due to the corresponding
  revoke by William Hoskins on July 22, 1999.
  <ftp://ftp.cs.berkeley.edu/pub/4bsd/README.Impt.License.Change>
*/
/*-
 * Copyright (c) 1992, 1993
 *        The Regents of the University of California.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

void xt_quicksort_int(int *a, size_t n)
{
#define MIN(a,b) ((a) < (b) ? (a) : (b))

  while (1) {
    bool swap_cnt = false;
    if (n < 7) {
      for (int *pm = a + 1; pm < a + n; ++pm)
        for (int *pl = pm; pl > a && pl[-1] > *pl; --pl) {
          int t = *pl; *pl = pl[-1]; pl[-1] = t;
        }
      return;
    }
    {
      int *pm = a + (n / 2);
      if (n > 7) {
        int *pl = a;
        int *pn = a + (n - 1);
        if (n > 40) {
          size_t d = n / 8;
          pl = med3_int(pl,  pl + d,  pl + 2 * d);
          pm = med3_int(pm - d,  pm,  pm + d);
          pn = med3_int(pn - 2 * d,  pn - d,  pn);
        }
        pm = med3_int(pl,  pm,  pn);
      }
      int t = *a; *a = *pm; *pm = t;
    }
    int *pa = a + 1, *pb = pa;
    int *pc = a + n - 1, *pd = pc;
    int pivot = *a;
    for (;;) {
      while (pb <= pc && *pb <= pivot) {
        if (*pb == pivot) {
          swap_cnt = true;
          int t = *pa; *pa = *pb; *pb = t;
          ++pa;
        }
        ++pb;
      }
      while (pb <= pc && *pc >= pivot) {
        if (*pc == pivot) {
          swap_cnt = true;
          int t = *pc; *pc = *pd; *pd = t;
          --pd;
        }
        --pc;
      }
      if (pb > pc)
        break;
      int t = *pb; *pb = *pc; *pc = t;
      swap_cnt = true;
      ++pb;
      --pc;
    }
    if (!swap_cnt) {  /* Switch to insertion sort */
      for (int *pm = a + 1; pm < a + n; ++pm)
        for (int *pl = pm; pl > a && pl[-1] > *pl; --pl) {
          int t = *pl; *pl = pl[-1]; pl[-1] = t;
        }
      return;
    }

    int *pn = a + n;
    ptrdiff_t pdiff = MIN(pa - a, pb - pa);
    vecswap_int(a, pb - pdiff, (size_t)pdiff);
    pdiff = MIN(pd - pc, pn - pd - 1);
    vecswap_int(pb, pn - pdiff, (size_t)pdiff);
    if ((pdiff = pb - pa) > (ptrdiff_t)1)
      xt_quicksort_int(a, (size_t)pdiff);
    if ((pdiff = pd - pc) > (ptrdiff_t)1) {
      /* Iterate rather than recurse to save stack space */
      a = pn - pdiff;
      n = (size_t)pdiff;
    }
    else
      break;
  }
#undef MIN
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
