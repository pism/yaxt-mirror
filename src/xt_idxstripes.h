#ifndef XT_IDXSTRIPES_H
#define XT_IDXSTRIPES_H

#include "xt_idxlist.h"

struct Xt_stripe {

  Xt_idx start;
  unsigned extent;
  int stride;
};

/**
 * Generates an index list that is build up of stripes of indices.
 *
 * @param[in] stripes     array defining the stripes
 * @param[in] num_stripes number of stripes
 *
 * @remarks this routine does not check for overlaps between stripes
 *          and assumes that there are no
 */
Xt_idxlist xt_idxstripes_new(struct Xt_stripe * stripes, unsigned num_stripes);

xt_idxlist xt_idxstripes_unpack(void *buffer, int buffer_size, int *position,
                                MPI_Comm comm);

#endif // XT_IDXSTRIPES_H
