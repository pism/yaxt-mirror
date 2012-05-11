#ifndef XT_IDXLIST_T_H
#define XT_IDXLIST_T_H

#include "yaxt.h"
#include "xt_idxlist.h"

/**
 * generates a new index vector based on an index list
 *
 * @param[in] idxlist     array containing the global indices
 * @param[in] num_indices number of indices in idxlist
 * @return returns a xt_idxlist
 */
xt_idxlist xt_idxlist_t_new(Xt_idx *idxlist, size_t num_indices);

xt_idxlist xt_idxlist_t_unpack(void *buffer, int buffer_size, int *position,
                               MPI_Comm comm);

#endif // XT_IDXLIST_T_H
