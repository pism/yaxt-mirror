#ifndef TEST_IDXLIST_UTILS_H
#define TEST_IDXLIST_UTILS_H

#include "tests.h"

#include "yaxt.h"
#include "xt_idxlist.h"

void ref_get_position_of_indices(Xt_idx * indices, int num_indices, Xt_idx * selection, int selection_size,
                                 int * positions, int repeat_first_match) {
  int used[num_indices];
  int found;

  for (int i = 0; i < selection_size; i++) {
    used[i] = 0;
  }

  for (int i = 0; i < selection_size; i++) {

    found = 0;
    for (int j = 0; j < num_indices; j++) {
      if (selection[i] == indices[j] && !used[j]) {
        positions[i]=j;
        if (!repeat_first_match) used[j] = 1;
        found = 1;
        break;
      }
    }
    if (!found)
      PUT_ERR("ref_get_position_of_indices: internal test error\n");
  }
}


void check_idxlist(xt_idxlist idxlist, Xt_idx * ref_indices,
                   unsigned ref_num_indices) {

  unsigned i;

  int num_indices;
  num_indices = xt_idxlist_get_num_indices(idxlist);

  if (num_indices != ref_num_indices)
    PUT_ERR("wrong number of indices in idxlist\n");

  Xt_idx indices[ref_num_indices];
  Xt_idx index;

  xt_idxlist_get_indices(idxlist, indices);

  {

    int position;
    int ref_positions[ref_num_indices];

    ref_get_position_of_indices(indices, num_indices, ref_indices, ref_num_indices,
                                ref_positions, 1);

    for (i = 0; i < ref_num_indices; ++i) {

      if (indices[i] != ref_indices[i])
        PUT_ERR("index list contains wrong indices\n");

      if (xt_idxlist_get_index(idxlist, i, &index) || index != ref_indices[i])
        PUT_ERR("error in xt_idxlist_get_index\n");

      if (xt_idxlist_get_position_of_index(idxlist, ref_indices[i], &position) ||
          position != ref_positions[i])
        PUT_ERR("error in xt_idxlist_get_position_of_index\n");
    }

    ref_get_position_of_indices(indices, num_indices, ref_indices, ref_num_indices,
                                ref_positions, 0);
    position = -1;
    for (i = 0; i < ref_num_indices; ++i) {
      if (xt_idxlist_get_position_of_index_off(idxlist, ref_indices[i], &position,
                                                position+1))
        PUT_ERR("error in xt_idxlist_get_position_of_index_off\n");

      if (position != ref_positions[i])
        PUT_ERR("get_position_of_indices returned wrong position\n");
    }
  }

  {
    int ref_positions[ref_num_indices], positions[ref_num_indices];

    ref_get_position_of_indices(indices, num_indices, ref_indices, ref_num_indices,
                                ref_positions, 0);

    if (xt_idxlist_get_position_of_indices(idxlist, ref_indices,
                                           ref_num_indices, positions, 1))
      PUT_ERR("get_position_of_indices failed\n");

    for (i = 0; i < ref_num_indices; ++i)
      if (positions[i] != ref_positions[i])
        PUT_ERR("check_idxlist: get_position_of_indices failed (case 1)\n");
  }

  {
    int ref_positions[ref_num_indices], positions[ref_num_indices];

    ref_get_position_of_indices(indices, num_indices, ref_indices, ref_num_indices,
                                ref_positions, 1);

    if (xt_idxlist_get_position_of_indices(idxlist, ref_indices,
                                           ref_num_indices, positions, 0))
      PUT_ERR("get_position_of_indices failed\n");

    for (i = 0; i < ref_num_indices; ++i)
      if (positions[i] != ref_positions[i])
        PUT_ERR("check_idxlist: get_position_of_indices failed (case 2)\n");
  }

  {
    Xt_idx reverse_ref_indices[ref_num_indices];

    for (i = 0; i < ref_num_indices; ++i)
      reverse_ref_indices[i] = ref_indices[ref_num_indices - i - 1];

    {
      int ref_positions[ref_num_indices], positions[ref_num_indices];

      ref_get_position_of_indices(indices, num_indices, reverse_ref_indices, ref_num_indices,
                                  ref_positions, 0);

      if (xt_idxlist_get_position_of_indices(idxlist, reverse_ref_indices,
                                             ref_num_indices, positions, 1))
        PUT_ERR("get_position_of_indices failed\n");

      for (i = 0; i < ref_num_indices; ++i)
        if (positions[i] != ref_positions[i])
          PUT_ERR("check_idxlist: get_position_of_indices failed (case 3)\n");
    }

    {
      int ref_positions[ref_num_indices], positions[ref_num_indices];

      ref_get_position_of_indices(indices, num_indices, reverse_ref_indices, ref_num_indices,
                                  ref_positions, 1);

      if (xt_idxlist_get_position_of_indices(idxlist, reverse_ref_indices,
                                             ref_num_indices, positions, 0))
        PUT_ERR("get_position_of_indices failed\n");

      for (i = 0; i < ref_num_indices; ++i)
        if (positions[i] != ref_positions[i])
          PUT_ERR("check_idxlist: get_position_of_indices failed (case 4)\n");
    }
  }

  if (!xt_idxlist_get_index(idxlist, -1, &index))
    PUT_ERR("error in xt_idxlist_get_index\n");

  if (!xt_idxlist_get_index(idxlist, ref_num_indices, &index))
    PUT_ERR("error in xt_idxlist_get_index\n");
}
#endif // TEST_IDXLIST_UTILS_H
