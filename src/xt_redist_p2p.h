#ifndef XT_REDIST_P2P_H
#define XT_REDIST_P2P_H

#include <mpi.h>

#include "xt_redist.h"
#include "xt_xmap.h"

/**
 * constructor for a redistribution using point to point
 * communication for the exchange
 * @param[in] xmap     exchange map
 * @param[in] datatype MPI datatype of single element
 *                     in the data to be exchanged
 */
Xt_redist xt_redist_p2p_new(Xt_xmap xmap, MPI_Datatype datatype);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange
 * @param[in] xmap     exchange map
 * @param[in] offsets  array containing for all elements in
 *                     the source index list passed to the
 *                     exchange map the position of the respective
 *                     element in the input array passed to the
 *                     exchange routine
 * @param[in] datatype MPI datatype of single element
 *                     in the data to be exchanged
 */
Xt_redist xt_redist_p2p_off_new(Xt_xmap xmap, ssize_t *offsets,
                                MPI_Datatype datatype);

#endif // XT_REDIST_P2P_H
