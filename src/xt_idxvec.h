#ifndef XT_IDXVEC_H
#define XT_IDXVEC_H

#include "xt/xt_core.h"
#include "xt/xt_idxlist.h"

/** \example test_idxvec.c
 */

/**
 * generates a new index list based on an index vector
 *
 * @param[in] idxlist     array containing the global indices
 * @param[in] num_indices number of indices in idxlist
 * @return returns a Xt_idxlist$
 */
Xt_idxlist xt_idxvec_new(Xt_idx *idxlist, Xt_count num_indices);

Xt_idxlist
xt_idxvec_from_stripes_new(struct Xt_stripe * stripes, Xt_count num_stripes);

Xt_idxlist xt_idxvec_unpack(void *buffer, int buffer_size, int *position,
                            MPI_Comm comm);

#endif // XT_IDXVEC_H
