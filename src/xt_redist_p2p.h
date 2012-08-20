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
 * @param[in] xmap         exchange map
 * @param[in] src_offsets  array containing for all elements in
 *                         the source index list passed to the
 *                         exchange map the position of the respective
 *                         element in the input array passed to the
 *                         exchange routine
 * @param[in] dst_offsets  array containing for all elements in
 *                         the destination index list passed to the
 *                         exchange map the position of the respective
 *                         element in the output array passed to the
 *                         exchange routine
 * @param[in] datatype     MPI datatype of single element
 *                         in the data to be exchanged
 */
Xt_redist xt_redist_p2p_off_new(Xt_xmap xmap, int *src_offsets,
                                int *dst_offsets, MPI_Datatype datatype);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange, special case: blocks with offsets
 * @param[in] xmap         exchange map
 * @param[in] src_block_offsets  array containing for all source index space positions
 *                               of xmap the offsets for blocks in data space
 * @param[in] src_block_sizes  source block lengths in unit of elements
 * @param[in] src_block_num  number of src blocks
 * @param[in] dst_offsets  array containing for all destination index space positions
 *                         of xmap the offsets for blocks in data space
 * @param[in] dst_block_sizes  destination block lengths in unit of elements
 * @param[in] dst_block_num  number of dst blocks
 * @param[in] datatype     MPI datatype of a single element in data space,
 *                         all elements have the same
 * \remarks NULL offsets arguments mean: use zero based prefix sum of block sizes as effective offsets
 */
Xt_redist xt_redist_p2p_blocks_off_new(Xt_xmap xmap,
                                       int *src_block_offsets, int *src_block_sizes, int src_block_num,
                                       int *dst_block_offsets, int *dst_block_sizes, int dst_block_num,
                                       MPI_Datatype datatype);

/**
 * constructor for a redistribution using point to point
 * communication for the exchange, special case: blocks without explicitoffsets
 * @param[in] xmap         exchange map
 * @param[in] src_block_sizes  source block lengths in unit of elements
 * @param[in] src_block_num  number of src blocks
 * @param[in] dst_block_sizes  destination block lengths in unit of elements
 * @param[in] dst_block_num  number of dst blocks
 * @param[in] datatype     MPI datatype of a single element in data space,
 *                         all elements have the same
 * \remarks calls \ref xt_redist_p2p_blocks_off_new with NULL offsets
 */
Xt_redist xt_redist_p2p_blocks_new(Xt_xmap xmap,
                                   int *src_block_sizes, int src_block_num,
                                   int *dst_block_sizes, int dst_block_num,
                                   MPI_Datatype datatype);

#endif // XT_REDIST_P2P_H
