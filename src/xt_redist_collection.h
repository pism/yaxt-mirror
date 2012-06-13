#ifndef XT_REDIST_COLLECTION_H
#define XT_REDIST_COLLECTION_H

#include "xt_redist.h"

/**
 * constructor for a redistribution collection that is comprised
 * of multiple other redistributions
 * @param[in] redists     redistributions
 * @param[in] num_redists number of redistributions
 * @param[in] comm        MPI communicator
 *
 * \remarks all redistributions need to be based on the same
 *          MPI communicator
 */
Xt_redist xt_redist_collection_new(Xt_redist * redists, unsigned num_redists,
                                   MPI_Comm comm);

#endif // XT_REDIST_COLLECTION_H
