#ifndef XT_IDXSTRIPES_H
#define XT_IDXSTRIPES_H

#include "xt_idxlist.h"

/** \example test_idxstripes.c
 */

/**
 * Generates an index list that is built up of stripes of indices.
 *
 * @param[in] stripes     array defining the stripes
 * @param[in] num_stripes number of stripes
 */
Xt_idxlist
xt_idxstripes_new(struct Xt_stripe const * stripes, Xt_count num_stripes);

Xt_idxlist xt_idxstripes_unpack(void *buffer, int buffer_size, int *position,
                                MPI_Comm comm);

/**
 * converts an array with indices into an array of stripes that represent
 * the same indices (with stride 1)
 *
 * @param[in]  indices     array containing the indices that are to be converted
 * @param[in]  num_indices number of indices in array "indices"
 * @param[out] stripes     array containing the stripes (user is responsible
 *                         for freeing the memory associated with the stripes)
 * @param[out] num_stripes number of stripes in array "stripes"
 */
void xt_idxstripes_convert_to_stripes(Xt_idx const * indices,
                                      Xt_count num_indices,
                                      struct Xt_stripe ** stripes,
                                      Xt_count * num_stripes);

#endif // XT_IDXSTRIPES_H
