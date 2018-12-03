/**
 * @file xt_dmap.h
 * @brief dmap
 *
 * contains declaration the dmap data structure, which
 * is derived from a xt_xmap
 *
 * @copyright Copyright  (C)  2018 Jörg Behrens <behrens@dkrz.de>
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

#ifndef XT_DMAP_H
#define XT_DMAP_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "xt/xt_core.h"

/** \example test_dmap_parallel.c
 */

enum xt_reorder_type {
  XT_REORDER_NONE,    //! no reordering
  XT_REORDER_SEND_UP, //! optimise data access on sender side
  XT_REORDER_RECV_UP, //! optimise data access on receiver side
};


/**
 * repeat the pattern stored in the data map
 * @param[in,out] dmap              data map
 * @param[in]     num_repetitions   number of repetitions of the given data map
 * @param[in]     src_displacements displacements for source repetitions
 * @param[in]     dst_displacements displacements for destination repetitions
 * @remark This routine needs to be called for all processes associated with the
 *         data map. The argument num_repetitions needs to be identical on all
 *         processes.
 * @remarkt The displacements are in reference to the extent for the
 *          base_datatype provided to \ref xt_dmap_generate_redist.
 */
void xt_dmap_repeat(Xt_dmap dmap, int num_repetitions,
                    const int src_displacements[num_repetitions],
                    const int dst_displacements[num_repetitions]);

/**
 * reorder offsets stored in data map
 * @param[in,out] dmap data map
 * @param[in]     type reorder pattern to be used
 * @remark this routine is collective
 */
void xt_dmap_reorder(Xt_dmap dmap, enum xt_reorder_type type);

/**
 * generate a redistributation from a data map
 * @param[in] dmap          data map
 * @param[in] base_datatype MPI datatype of single element
 *                          in the data to be exchanged
 */
Xt_redist xt_dmap_generate_redist(Xt_dmap dmap, MPI_Datatype base_datatype);
void xt_dmap_delete(Xt_dmap dmap);

#endif // XT_DMAP_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
