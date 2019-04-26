/**
 * @file xt_memcmp.c
 *
 * @copyright Copyright  (C)  2019 Jörg Behrens <behrens@dkrz.de>
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

#include <stdlib.h>
#include <string.h>
#include "core/core.h"

#include "xt_memcmp.h"

#define DEF_CMP_BYTE(X) \
static int cmp_byte_ ## X (const void * a, const void * b) { \
  for (size_t i = 0; i < ( X ); ++i) { \
    int ret = (int)(((const unsigned char *)a)[i] - \
                    ((const unsigned char *)b)[i]); \
    if (ret) return ret; \
  } \
  return 0; \
}

DEF_CMP_BYTE(1)
DEF_CMP_BYTE(2)
DEF_CMP_BYTE(3)
DEF_CMP_BYTE(4)
DEF_CMP_BYTE(5)
DEF_CMP_BYTE(6)
DEF_CMP_BYTE(7)
DEF_CMP_BYTE(8)
DEF_CMP_BYTE(9)
DEF_CMP_BYTE(10)
DEF_CMP_BYTE(11)
DEF_CMP_BYTE(12)
DEF_CMP_BYTE(13)
DEF_CMP_BYTE(14)
DEF_CMP_BYTE(15)
DEF_CMP_BYTE(16)
DEF_CMP_BYTE(17)
DEF_CMP_BYTE(18)
DEF_CMP_BYTE(19)
DEF_CMP_BYTE(20)
DEF_CMP_BYTE(21)
DEF_CMP_BYTE(22)
DEF_CMP_BYTE(23)
DEF_CMP_BYTE(24)
DEF_CMP_BYTE(25)
DEF_CMP_BYTE(26)
DEF_CMP_BYTE(27)
DEF_CMP_BYTE(28)
DEF_CMP_BYTE(29)
DEF_CMP_BYTE(30)
DEF_CMP_BYTE(31)
DEF_CMP_BYTE(32)
DEF_CMP_BYTE(33)
DEF_CMP_BYTE(34)
DEF_CMP_BYTE(35)
DEF_CMP_BYTE(36)
DEF_CMP_BYTE(37)
DEF_CMP_BYTE(38)
DEF_CMP_BYTE(39)
DEF_CMP_BYTE(40)
DEF_CMP_BYTE(41)
DEF_CMP_BYTE(42)
DEF_CMP_BYTE(43)
DEF_CMP_BYTE(44)
DEF_CMP_BYTE(45)
DEF_CMP_BYTE(46)
DEF_CMP_BYTE(47)
DEF_CMP_BYTE(48)
DEF_CMP_BYTE(49)
DEF_CMP_BYTE(50)
DEF_CMP_BYTE(51)
DEF_CMP_BYTE(52)
DEF_CMP_BYTE(53)
DEF_CMP_BYTE(54)
DEF_CMP_BYTE(55)
DEF_CMP_BYTE(56)
DEF_CMP_BYTE(57)
DEF_CMP_BYTE(58)
DEF_CMP_BYTE(59)
DEF_CMP_BYTE(60)
DEF_CMP_BYTE(61)
DEF_CMP_BYTE(62)
DEF_CMP_BYTE(63)
DEF_CMP_BYTE(64)

#undef DEF_CMP_BYTE

Xt_memcmp_func xt_memcmp_get_func(size_t data_size) {

#define CASE_BYTE_CMP(X) \
  case( X ): { \
    return cmp_byte_ ## X ;\
  }

  switch (data_size) {
    CASE_BYTE_CMP(1)
    CASE_BYTE_CMP(2)
    CASE_BYTE_CMP(3)
    CASE_BYTE_CMP(4)
    CASE_BYTE_CMP(5)
    CASE_BYTE_CMP(6)
    CASE_BYTE_CMP(7)
    CASE_BYTE_CMP(8)
    CASE_BYTE_CMP(9)
    CASE_BYTE_CMP(10)
    CASE_BYTE_CMP(11)
    CASE_BYTE_CMP(12)
    CASE_BYTE_CMP(13)
    CASE_BYTE_CMP(14)
    CASE_BYTE_CMP(15)
    CASE_BYTE_CMP(16)
    CASE_BYTE_CMP(17)
    CASE_BYTE_CMP(18)
    CASE_BYTE_CMP(19)
    CASE_BYTE_CMP(20)
    CASE_BYTE_CMP(21)
    CASE_BYTE_CMP(22)
    CASE_BYTE_CMP(23)
    CASE_BYTE_CMP(24)
    CASE_BYTE_CMP(25)
    CASE_BYTE_CMP(26)
    CASE_BYTE_CMP(27)
    CASE_BYTE_CMP(28)
    CASE_BYTE_CMP(29)
    CASE_BYTE_CMP(30)
    CASE_BYTE_CMP(31)
    CASE_BYTE_CMP(32)
    CASE_BYTE_CMP(33)
    CASE_BYTE_CMP(34)
    CASE_BYTE_CMP(35)
    CASE_BYTE_CMP(36)
    CASE_BYTE_CMP(37)
    CASE_BYTE_CMP(38)
    CASE_BYTE_CMP(39)
    CASE_BYTE_CMP(40)
    CASE_BYTE_CMP(41)
    CASE_BYTE_CMP(42)
    CASE_BYTE_CMP(43)
    CASE_BYTE_CMP(44)
    CASE_BYTE_CMP(45)
    CASE_BYTE_CMP(46)
    CASE_BYTE_CMP(47)
    CASE_BYTE_CMP(48)
    CASE_BYTE_CMP(49)
    CASE_BYTE_CMP(50)
    CASE_BYTE_CMP(51)
    CASE_BYTE_CMP(52)
    CASE_BYTE_CMP(53)
    CASE_BYTE_CMP(54)
    CASE_BYTE_CMP(55)
    CASE_BYTE_CMP(56)
    CASE_BYTE_CMP(57)
    CASE_BYTE_CMP(58)
    CASE_BYTE_CMP(59)
    CASE_BYTE_CMP(60)
    CASE_BYTE_CMP(61)
    CASE_BYTE_CMP(62)
    CASE_BYTE_CMP(63)
    CASE_BYTE_CMP(64)
    default:
      Xt_abort(Xt_default_comm, "unsupported element size", __FILE__, __LINE__);
      return NULL;
  }
#undef CASE_BYTE_CMP
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
