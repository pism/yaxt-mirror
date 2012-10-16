#ifndef XT_IDXSECTION_H
#define XT_IDXSECTION_H

#include "xt_idxlist.h"

/** \example test_idxsection.c
 */

/**
 * generates an index list that is comprised of a section of a set of
 * indices that are arranged in an n-dimensional....don't know how to
 * describe it...
 * @param[in] start          lowest index of the global array (typically 0 or 1)
 * @param[in] num_dimensions number of dimensions
 * @param[in] global_size    global size of each dimension
 * @param[in] local_size     size of the local section in each dimension
 * @param[in] local_start    index of the first index of the local part
 *                           within the global array for each dimension (C-style)
 * 
 */
Xt_idxlist xt_idxsection_new(Xt_idx start, unsigned num_dimensions,
                             unsigned * global_size, unsigned * local_size,
                             unsigned * local_start);

Xt_idxlist xt_idxsection_unpack(void *buffer, int buffer_size, int *position,
                                MPI_Comm comm);

#endif // XT_IDXSECTION_H
