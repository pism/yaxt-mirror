/**
 * @file xt_ddt_internal.h
 * @brief internal utility routines for manual handling of MPI DDT's
 *
 * contains internal utility routines for handling manual MPI DDT's
 *
 * @copyright Copyright  (C)  2022 Jörg Behrens <behrens@dkrz.de>
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
 *
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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

#ifndef XT_DDT_INTERNAL_H
#define XT_DDT_INTERNAL_H

#include "xt/xt_ddt.h"



/**
 * determines the type of memory associated with the provided pointer
 * @param[in] ptr pointer returned by a allocation routine
 * @returns type of memory
 */
enum xt_ddt_memtype xt_ddt_get_memtype(const void *ptr);

/**
 * packs the data from the source buffer into destination buffer
 * @param[in]  ddt     xt_ddt object
 * @param[in]  src     source buffer
 * @param[out] dst     destination buffer
 * @param[in]  memtype type of source and destination memory
 * @remark both src and dst have to point to the same type of memory
 */
void xt_ddt_pack_internal(
  Xt_ddt ddt, void const * src, void * dst, enum xt_ddt_memtype memtype);

/**
 * inpacks the data from the source buffer into destination buffer
 * @param[in]  ddt     xt_ddt object
 * @param[in]  src     source buffer
 * @param[out] dst     destination buffer
 * @param[in]  memtype type of source and destination memory
 * @remark both src and dst have to point to the same type of memory
 */
void xt_ddt_unpack_internal(
  Xt_ddt ddt, void const * src, void * dst, enum xt_ddt_memtype memtype);

/**
 * allocates memory of the specified type
 * @param[in] alloc_size size of memory to be allocated
 * @param[in] memtype    type of memory to be allocated
 * @return allocated memory
 */
void * xt_ddt_malloc(size_t alloc_size, enum xt_ddt_memtype memtype);

/**
 * frees memory that was allocated using \xt_ddt_malloc
 * @param[in] ptr     pointer to memory
 * @param[in] memtype type of memory to be freed
 */
void xt_ddt_free(void * ptr, enum xt_ddt_memtype memtype);

/**
 * copies memory; src and dst can have different memory types
 * @param[in] dst         pointer to destination
 * @param[in] src         pointer to source
 * @param[in] buffer_size size of memory to be copied
 * @param[in] src_memtype type of source memory
 * @param[in] tgt_memtype type of target memory
 */
void xt_ddt_memcpy(
  void * dst, void const * src, size_t buffer_size,
  enum xt_ddt_memtype dst_memtype, enum xt_ddt_memtype src_memtype);

#endif // XT_DDT_INTERNAL_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
