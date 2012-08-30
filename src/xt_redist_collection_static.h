#ifndef XT_REDIST_COLLECTION_STATIC_H
#define XT_REDIST_COLLECTION_STATIC_H

#include "xt_redist.h"

/**
 * constructor for a redistribution collection that is comprised
 * of multiple other redistributions
 * @param[in] redists           redistributions
 * @param[in] num_redists       number of redistributions
 * @param[in] src_displacements array of displacements of the source
 *                              input arrays for the exchange
 * @param[in] dst_displacements array of displacements of the destination
 *                              input arrays for the exchange
 * @param[in] comm              MPI communicator
 *
 * \remarks all redistributions need to be based on the same
 *          MPI communicator
 */
Xt_redist xt_redist_collection_static_new(Xt_redist * redists, unsigned num_redists,
                                          MPI_Aint * src_displacements,
                                          MPI_Aint * dst_displacements, MPI_Comm comm);

#endif // XT_REDIST_COLLECTION_STATIC_H
