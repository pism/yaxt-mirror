#ifndef TEST_IDXLIST_UTILS_H
#define TEST_IDXLIST_UTILS_H

#include "tests.h"

#include "yaxt.h"
#include "xt_idxlist.h"

void ref_get_positions_of_indices(Xt_idx * indices, Xt_count num_indices,
                                 Xt_idx * selection, int selection_size,
                                 int * positions, int repeat_first_match) {
  int used[num_indices];
  int found;

  for (int i = 0; i < selection_size; i++) {
    used[i] = 0;
  }

  for (int i = 0; i < selection_size; i++) {

    found = 0;
    for (Xt_count j = 0; j < num_indices; j++) {
      if (selection[i] == indices[j] && !used[j]) {
        positions[i]=j;
        if (!repeat_first_match) used[j] = 1;
        found = 1;
        break;
      }
    }
    if (!found)
      PUT_ERR("ref_get_positions_of_indices: internal test error\n");
  }
}

Xt_idx get_min_index(Xt_idx * indices, Xt_count num_indices) {

  Xt_idx index;
  Xt_count i;

  if (num_indices > 0)
    index = indices[0];
  else
    index = 0;

  for (i = 1; i < num_indices; ++i)
    if (index > indices[i])
      index = indices[i];

  return index;
}

Xt_idx get_max_index(Xt_idx * indices, Xt_count num_indices) {

  Xt_idx index;
  Xt_count i;

  if (num_indices > 0)
    index = indices[0];
  else
    index = 0;

  for (i = 1; i < num_indices; ++i)
    if (index < indices[i])
      index = indices[i];

  return index;
}

void check_idxlist(Xt_idxlist idxlist, Xt_idx * ref_indices,
                   Xt_count ref_num_indices) {

  Xt_count i;

  Xt_count num_indices;
  num_indices = xt_idxlist_get_num_indices(idxlist);

  if (num_indices != ref_num_indices)
    PUT_ERR("wrong number of indices in idxlist\n");

  Xt_idx indices[ref_num_indices];
  Xt_idx const *indices_const;
  Xt_idx index;

  xt_idxlist_get_indices(idxlist, indices);
  indices_const = xt_idxlist_get_indices_const(idxlist);
  {
    for(int i=0; i<ref_num_indices; i++)
      if (ref_indices[i] != indices_const[i])
        PUT_ERR("check_idxlist: ref_indices[i] != indices_const[i]\n");
  }

  {

    int position;
    int ref_positions[ref_num_indices];

    ref_get_positions_of_indices(indices, num_indices, ref_indices, ref_num_indices,
                                ref_positions, 1);

    for (i = 0; i < ref_num_indices; ++i) {

      if (indices[i] != ref_indices[i])
        PUT_ERR("index list contains wrong indices\n");

      if (xt_idxlist_get_index_at_position(idxlist, i, &index) || index != ref_indices[i])
        PUT_ERR("error in xt_idxlist_get_index_at_position\n");

      if (xt_idxlist_get_position_of_index(idxlist, ref_indices[i], &position) ||
          position != ref_positions[i])
        PUT_ERR("error in xt_idxlist_get_position_of_index\n");
    }

    ref_get_positions_of_indices(indices, num_indices, ref_indices, ref_num_indices,
                                ref_positions, 0);
    position = -1;
    for (i = 0; i < ref_num_indices; ++i) {
      if (xt_idxlist_get_position_of_index_off(idxlist, ref_indices[i], &position,
                                                position+1))
        PUT_ERR("error in xt_idxlist_get_position_of_index_off\n");

      if (position != ref_positions[i])
        PUT_ERR("get_positions_of_indices returned wrong position\n");
    }
  }

  {
    int ref_positions[ref_num_indices], positions[ref_num_indices];

    ref_get_positions_of_indices(indices, num_indices, ref_indices, ref_num_indices,
                                ref_positions, 0);

    if (xt_idxlist_get_positions_of_indices(idxlist, ref_indices,
                                           ref_num_indices, positions, 1))
      PUT_ERR("get_positions_of_indices failed\n");

    for (i = 0; i < ref_num_indices; ++i)
      if (positions[i] != ref_positions[i])
        PUT_ERR("check_idxlist: get_positions_of_indices failed (case 1)\n");
  }

  {
    int ref_positions[ref_num_indices], positions[ref_num_indices];

    ref_get_positions_of_indices(indices, num_indices, ref_indices, ref_num_indices,
                                ref_positions, 1);

    if (xt_idxlist_get_positions_of_indices(idxlist, ref_indices,
                                           ref_num_indices, positions, 0))
      PUT_ERR("get_positions_of_indices failed\n");

    for (i = 0; i < ref_num_indices; ++i)
      if (positions[i] != ref_positions[i])
        PUT_ERR("check_idxlist: get_positions_of_indices failed (case 2)\n");
  }

  {
    Xt_idx reverse_ref_indices[ref_num_indices];

    for (i = 0; i < ref_num_indices; ++i)
      reverse_ref_indices[i] = ref_indices[ref_num_indices - i - 1];

    {
      int ref_positions[ref_num_indices], positions[ref_num_indices];

      ref_get_positions_of_indices(indices, num_indices, reverse_ref_indices, ref_num_indices,
                                  ref_positions, 0);

      if (xt_idxlist_get_positions_of_indices(idxlist, reverse_ref_indices,
                                             ref_num_indices, positions, 1))
        PUT_ERR("get_positions_of_indices failed\n");

      for (i = 0; i < ref_num_indices; ++i)
        if (positions[i] != ref_positions[i])
          PUT_ERR("check_idxlist: get_positions_of_indices failed (case 3)\n");
    }

    {
      int ref_positions[ref_num_indices], positions[ref_num_indices];

      ref_get_positions_of_indices(indices, num_indices, reverse_ref_indices, ref_num_indices,
                                  ref_positions, 1);

      if (xt_idxlist_get_positions_of_indices(idxlist, reverse_ref_indices,
                                             ref_num_indices, positions, 0))
        PUT_ERR("get_positions_of_indices failed\n");

      for (i = 0; i < ref_num_indices; ++i)
        if (positions[i] != ref_positions[i])
          PUT_ERR("check_idxlist: get_positions_of_indices failed (case 4)\n");
    }
  }

  if (ref_num_indices > 0) {

    if (get_min_index(ref_indices, ref_num_indices) !=
        xt_idxlist_get_min_index(idxlist))
      PUT_ERR("check_idxlist: get_min_index failed\n")

    if (get_max_index(ref_indices, ref_num_indices) !=
        xt_idxlist_get_max_index(idxlist))
      PUT_ERR("check_idxlist: get_max_index failed\n")
  }

  if (!xt_idxlist_get_index_at_position(idxlist, -1, &index))
    PUT_ERR("error in xt_idxlist_get_index_at_position\n");

  if (!xt_idxlist_get_index_at_position(idxlist, ref_num_indices, &index))
    PUT_ERR("error in xt_idxlist_get_index_at_position\n");
}
#endif // TEST_IDXLIST_UTILS_H
