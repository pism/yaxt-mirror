#ifndef XT_IDXLIST_H
#define XT_IDXLIST_H

#include "xt_idxset.h"

/**
 * generates a new index set based on an index list
 *
 * @param[in] idxlist     array containing the global indices
 * @param[in] num_indices number of indices in idxlist
 * @return returns a xt_idxset
 */
xt_idxset xt_idxlist_new(int *idxlist, unsigned num_indices);

xt_idxset xt_idxlist_unpack(void *buffer, int buffer_size, int *position,
                            MPI_Comm comm);

#endif // XT_IDXLIST_H
