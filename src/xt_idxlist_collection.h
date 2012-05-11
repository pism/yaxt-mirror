#ifndef XT_IDXLIST_COLLECTION_H
#define XT_IDXLIST_COLLECTION_H

#include "xt_idxlist.h"

/**
 * generates a new index vector based on a collection of index vectors
 *
 * @param[in] idxlists     array containing a number of index vectors
 * @param[in] num_idxlists number of index vectors in idxlists
 * @return returns a xt_idxlist
 */
xt_idxlist xt_idxlist_collection_new(xt_idxlist *idxlists, unsigned num_idxlists);

xt_idxlist xt_idxlist_collection_unpack(void *buffer, int buffer_size, int *position,
                                        MPI_Comm comm);

#endif // XT_IDXLIST_COLLECTION_H
