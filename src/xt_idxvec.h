#ifndef XT_IDXVEC_H
#define XT_IDXVEC_H

#include "yaxt.h"
#include "xt_idxlist.h"

/**
 * generates a new index list based on an index vector
 *
 * @param[in] idxlist     array containing the global indices
 * @param[in] num_indices number of indices in idxlist
 * @return returns a xt_idxlist
 */
xt_idxlist xt_idxvec_new(Xt_idx *idxlist, Xt_count num_indices);

xt_idxlist xt_idxvec_unpack(void *buffer, int buffer_size, int *position,
                            MPI_Comm comm);

#endif // XT_IDXVEC_H
