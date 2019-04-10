/**
 * @file xt_xmap_dist_radix.h
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

#ifndef XT_XMAP_DIST_RADIX_H
#define XT_XMAP_DIST_RADIX_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xt/xt_core.h"

/** \example test_xmap_dist_radix.c
 */
/** \example test_xmap_dist_radix_f.f90
 */
/** \example test_xmap_dist_radix_parallel.c
 */
/** \example test_xmap_dist_radix_parallel_f.f90
 */
/** @example test_xmap_dist_radix_fail.c
 */
/** @example test_xmap_dist_radix_fail_f.f90
 */

/**
 * constructor for an exchange map \n
 * this operation is collective over all processes in comm \n
 * it uses all-to-all communication for the initialisation
 *
 * @param[in] src_ids     source identifier array
 * @param[in] num_src_ids number of entries in src_ids
 * @param[in] dst_ids     destination identifier array
 * @param[in] num_dst_ids number of entries in dst_ids
 * @param[in] id_size     size of single identifier
 * @param[in] comm        MPI communicator that contains all processes
 *                        that take part in the exchange (xt_xmap_dist_radix_new
 *                        will make its own copy of comm)
 * @remark This constructor will do a binary matching between source and
 *         destination identifiers, therefore the user has to ensure, that no
 *         identifer contains undefined bits (e.g. due to alignment issues).
 * @remark This constructor supports intra communicators as well as inter
 *         communications.
 */
Xt_xmap xt_xmap_dist_radix_new(void const * src_ids, int num_src_ids,
                               void const * dst_ids, int num_dst_ids,
                               size_t id_size, MPI_Comm comm);

#endif // XT_XMAP_DIST_RADIX_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
