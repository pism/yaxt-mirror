/**
 * @file xt_quicksort_base.h
 * @brief macros to create quicksort implementations
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
#ifndef XT_QUICKSORT_BASE_H
#define XT_QUICKSORT_BASE_H
#define TOKEN_PASTE(a,b) a##_##b
#define NAME_COMPOSE(a,b) TOKEN_PASTE(a,b)
#endif
#ifndef SORT_TYPE
#error "must define type to sort on"
#endif
#ifndef SORT_TYPE_SUFFIX
#error "must define suffix for type to name functions"
#endif
#ifndef SORT_TYPE_CMP_LT
#error "must define suffix for type to name functions"
#endif

#define MED3 NAME_COMPOSE(med3,SORT_TYPE_SUFFIX)
#define VECSWAP NAME_COMPOSE(vecswap,SORT_TYPE_SUFFIX)
#define XT_QUICKSORT NAME_COMPOSE(xt_quicksort,SORT_TYPE_SUFFIX)

static inline SORT_TYPE *
MED3(SORT_TYPE *a, SORT_TYPE *b, SORT_TYPE *c)
{
  return SORT_TYPE_CMP_LT(*a,*b) ?
    (SORT_TYPE_CMP_LT(*b,*c) ? b : (SORT_TYPE_CMP_LT(*a,*c) ? a : c ))
    : (SORT_TYPE_CMP_LT(*c,*b) ? b : (SORT_TYPE_CMP_LT(*a,*c) ? a : c ));
}

#ifndef SWAP
#define SWAP(a,b) do {                          \
    SORT_TYPE t = a; a = b; b = t;              \
  } while (0)
#endif

static inline void
VECSWAP(SORT_TYPE *restrict a, SORT_TYPE *restrict b, size_t n)
{
  for (size_t i = 0; i < n; ++i) {
    SORT_TYPE t = a[i];
    a[i] = b[i];
    b[i] = t;
  }
}

void XT_QUICKSORT(SORT_TYPE *a, size_t n)
{
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

  while (1) {
    bool swap_cnt = false;
    if (n < 7) {
      for (SORT_TYPE *pm = a + 1; pm < a + n; ++pm)
        for (SORT_TYPE *pl = pm; pl > a && SORT_TYPE_CMP_LT(*pl, pl[-1]); --pl)
          SWAP(*pl, pl[-1]);
      return;
    }
    {
      SORT_TYPE *pm = a + (n / 2);
      if (n > 7) {
        SORT_TYPE *pl = a;
        SORT_TYPE *pn = a + (n - 1);
        if (n > 40) {
          size_t d = n / 8;
          pl = MED3(pl,  pl + d,  pl + 2 * d);
          pm = MED3(pm - d,  pm,  pm + d);
          pn = MED3(pn - 2 * d,  pn - d,  pn);
        }
        pm = MED3(pl,  pm,  pn);
      }
      SWAP(*a, *pm);
    }
    SORT_TYPE *pa = a + 1, *pb = pa;
    SORT_TYPE *pc = a + n - 1, *pd = pc;
    SORT_TYPE pivot = *a;
    for (;;) {
      while (pb <= pc && SORT_TYPE_CMP_LE(*pb, pivot)) {
        if (SORT_TYPE_CMP_EQ(*pb, pivot)) {
          swap_cnt = true;
          SWAP(*pa, *pb);
          ++pa;
        }
        ++pb;
      }
      while (pb <= pc && SORT_TYPE_CMP_LE(pivot, *pc)) {
        if (SORT_TYPE_CMP_EQ(*pc, pivot)) {
          swap_cnt = true;
          SWAP(*pc, *pd);
          --pd;
        }
        --pc;
      }
      if (pb > pc)
        break;
      SWAP(*pb, *pc);
      swap_cnt = true;
      ++pb;
      --pc;
    }
    if (!swap_cnt) {  /* Switch to insertion sort */
      for (SORT_TYPE *pm = a + 1; pm < a + n; ++pm)
        for (SORT_TYPE *pl = pm; pl > a && SORT_TYPE_CMP_LT(*pl, pl[-1]); --pl)
          SWAP(*pl, pl[-1]);
      return;
    }

    SORT_TYPE *pn = a + n;
    ptrdiff_t pdiff = MIN(pa - a, pb - pa);
    VECSWAP(a, pb - pdiff, (size_t)pdiff);
    pdiff = MIN(pd - pc, pn - pd - 1);
    VECSWAP(pb, pn - pdiff, (size_t)pdiff);
    if ((pdiff = pb - pa) > (ptrdiff_t)1)
      XT_QUICKSORT(a, (size_t)pdiff);
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
