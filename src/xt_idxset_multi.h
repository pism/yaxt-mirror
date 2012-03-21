#ifndef XT_IDXSET_MULTI_H
#define XT_IDXSET_MULTI_H

#include "xt_idxset.h"

/**
 * generates a new index set based on an index list
 *
 * @param[in] idxsets     array containing a number of index sets
 * @param[in] num_idxsets number of index sets in idxsets
 * @return returns a xt_idxset
 */
xt_idxset xt_idxset_multi_new(xt_idxset *idxsets, unsigned num_idxsets);

xt_idxset xt_idxset_multi_unpack(void *buffer, int buffer_size, int *position,
                                 MPI_Comm comm);

#endif // XT_IDXSET_MULTI_H
