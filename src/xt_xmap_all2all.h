#ifndef XT_XMAP_ALL2ALL_H
#define XT_XMAP_ALL2ALL_H

#include "xt/xt_core.h"

/** \example test_xmap_all2all.c
 */

/** \example test_xmap_all2all_parallel.c
 */

/**
 * constructor for an exchange map \n
 * this operation is collective over all processes in comm \n
 * it uses a all to all communication for the initialisation
 *
 * @param[in] src_idxlist source index list
 * @param[in] dst_idxlist destination index list
 * @param[in] comm        MPI communicator that contains all processes
 *                        that take part in the exchange (xt_xmap_all2all_new
 *                        will make its own copy of comm)
 */
Xt_xmap xt_xmap_all2all_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                            MPI_Comm comm);

#endif // XT_XMAP_ALL2ALL_H
