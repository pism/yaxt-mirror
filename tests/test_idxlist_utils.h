#ifndef TEST_IDXLIST_UTILS_H
#define TEST_IDXLIST_UTILS_H

#include "tests.h"

#include "yaxt.h"
#include "xt_idxlist.h"

void check_idxlist(xt_idxlist idxlist, Xt_idx * ref_indices,
                   unsigned ref_num_indices) {

  unsigned i;

  if (xt_idxlist_get_num_indices(idxlist) != ref_num_indices)
    PUT_ERR("wrong number of indices in idxlist\n");

  Xt_idx indices[ref_num_indices];

  xt_idxlist_get_indices(idxlist, indices);

  Xt_idx index;

  Xt_idx reverse_indices[ref_num_indices];
  int position, positions[ref_num_indices], reverse_positions[ref_num_indices];

  for (i = 0; i < ref_num_indices; ++i)
    reverse_indices[i] = ref_indices[ref_num_indices - i - 1];

  if (xt_idxlist_get_position_of_indices(idxlist, ref_indices,
                                         ref_num_indices, positions, 0))
    PUT_ERR("get_position_of_indices failed\n");

  if (xt_idxlist_get_position_of_indices(idxlist, reverse_indices,
                                         ref_num_indices, reverse_positions, 0))
    PUT_ERR("get_position_of_indices failed for reverse indices\n");

  for (i = 0; i < ref_num_indices; ++i) {

    if (indices[i] != ref_indices[i])
      PUT_ERR("index list contains wrong indices\n");

    if (xt_idxlist_get_index(idxlist, i, &index) || index != ref_indices[i])
      PUT_ERR("error in xt_idxlist_get_index\n");

    if (xt_idxlist_get_position_of_index(idxlist, ref_indices[i], &position) ||
        position != i)
      PUT_ERR("error in xt_idxlist_get_position_of_index\n");

    if (!xt_idxlist_get_position_of_index_off(idxlist, ref_indices[i], &position,
                                              position+1))
      PUT_ERR("error in xt_idxlist_get_position_of_index_off\n");

    if (positions[i] != i)
      PUT_ERR("get_position_of_indices returned wrong position\n");

    if (reverse_positions[i] != ref_num_indices - i - 1)
      PUT_ERR("get_position_of_indices (with reverse indices) returned wrong position\n");
  }

  if (!xt_idxlist_get_index(idxlist, -1, &index))
    PUT_ERR("error in xt_idxlist_get_index\n");

  if (!xt_idxlist_get_index(idxlist, ref_num_indices, &index))
    PUT_ERR("error in xt_idxlist_get_index\n");
}
#endif // TEST_IDXLIST_UTILS_H
