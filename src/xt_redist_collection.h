#ifndef XT_REDIST_COLLECTION_H
#define XT_REDIST_COLLECTION_H

#include "xt_redist.h"

/** \example test_redist_collection.c
 */

/** \example test_redist_collection_parallel.c
 */

/**
 * constructor for a redistribution collection that is comprised
 * of multiple other redistributions
 * @param[in] redists     redistributions
 * @param[in] num_redists number of redistributions
 * @param[in] cache_size  number of cached datatypes required for
 *                        different input array displacements \n
 *                        n = -1 => default cache size \n
 *                        n =  0 => no cache \n
 *                        n >  0 => cache size = n
 * @param[in] comm        MPI communicator
 *
 * \remarks all redistributions need to be based on the same
 *          MPI communicator
 */
Xt_redist xt_redist_collection_new(Xt_redist * redists, unsigned num_redists,
                                   int cache_size, MPI_Comm comm);

#endif // XT_REDIST_COLLECTION_H
