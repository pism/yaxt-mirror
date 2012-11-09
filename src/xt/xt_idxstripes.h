/**
 * @file xt_idxstripes.h
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://redmine.dkrz.de/doc/yaxt/html/index.html
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

#ifndef XT_IDXSTRIPES_H
#define XT_IDXSTRIPES_H

#include "xt/xt_idxlist.h"

/** \example test_idxstripes.c
 */

/**
 * Generates an index list that is built up of stripes of indices.
 *
 * @param[in] stripes     array defining the stripes
 * @param[in] num_stripes number of stripes
 */
Xt_idxlist
xt_idxstripes_new(struct Xt_stripe const * stripes, Xt_count num_stripes);

Xt_idxlist xt_idxstripes_unpack(void *buffer, int buffer_size, int *position,
                                MPI_Comm comm);

/**
 * converts an array with indices into an array of stripes that represent
 * the same indices (with stride 1)
 *
 * @param[in]  indices     array containing the indices that are to be converted
 * @param[in]  num_indices number of indices in array "indices"
 * @param[out] stripes     array containing the stripes (user is responsible
 *                         for freeing the memory associated with the stripes)
 * @param[out] num_stripes number of stripes in array "stripes"
 */
void xt_idxstripes_convert_to_stripes(Xt_idx const * indices,
                                      Xt_count num_indices,
                                      struct Xt_stripe ** stripes,
                                      Xt_count * num_stripes);

#endif // XT_IDXSTRIPES_H
