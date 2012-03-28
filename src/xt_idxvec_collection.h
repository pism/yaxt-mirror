#ifndef XT_IDXVEC_COLLECTION_H
#define XT_IDXVEC_COLLECTION_H

#include "xt_idxvec.h"

/**
 * generates a new index vector based on a collection of index vectors
 *
 * @param[in] idxvecs     array containing a number of index vectors
 * @param[in] num_idxvecs number of index vectors in idxvecs
 * @return returns a xt_idxvec
 */
xt_idxvec xt_idxvec_collection_new(xt_idxvec *idxvecs, unsigned num_idxvecs);

xt_idxvec xt_idxvec_collection_unpack(void *buffer, int buffer_size, int *position,
                                      MPI_Comm comm);

#endif // XT_IDXVEC_COLLECTION_H
