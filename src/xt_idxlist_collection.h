#ifndef XT_IDXLIST_COLLECTION_H
#define XT_IDXLIST_COLLECTION_H

#include "xt_idxlist.h"

/** \example test_idxlist_collection.c
 */

/**
 * generates a new index list based on a collection of index lists
 *
 * @param[in] idxlists     array containing a number of index lists
 * @param[in] num_idxlists number of index lists in idxlists
 * @return returns a Xt_idxlist$
 */
Xt_idxlist
xt_idxlist_collection_new(Xt_idxlist *idxlists, unsigned num_idxlists);

Xt_idxlist
xt_idxlist_collection_unpack(void *buffer, int buffer_size, int *position,
                             MPI_Comm comm);

#endif // XT_IDXLIST_COLLECTION_H
