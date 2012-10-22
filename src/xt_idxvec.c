#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "yaxt.h"
#include "xt_idxlist.h"
#include "xt_idxvec.h"
#include "xt_idxstripes.h"
#include "xt_mpi.h"
#include "core/ppm_xfuncs.h"
#include "quicksort.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

static void
idxvec_delete(Xt_idxlist data);

static size_t
idxvec_get_pack_size(Xt_idxlist data, MPI_Comm comm);

static void
idxvec_pack(Xt_idxlist data, void *buffer, int buffer_size,
            int *position, MPI_Comm comm);

static Xt_idxlist
idxvec_get_intersection(Xt_idxlist idxlist_a, Xt_idxlist idxlist_b);

static Xt_idxlist
idxvec_copy(Xt_idxlist idxlist);

static Xt_count
idxvec_get_num_indices(Xt_idxlist idxlist);

static void
idxvec_get_indices(Xt_idxlist idxlist, Xt_idx *indices);

static void
idxvec_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                         Xt_count * num_stripes);

static int
idxvec_get_index_at_position(Xt_idxlist idxlist, int position, Xt_idx * index);

static int
idxvec_get_indices_at_positions(Xt_idxlist idxlist, int *positions, Xt_count num, Xt_idx *index, Xt_idx undef_idx);

static int
idxvec_get_position_of_index(Xt_idxlist idxlist, Xt_idx index, int * position);

static int
idxvec_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                 int * position, int offset);

static int
idxvec_get_positions_of_indices(Xt_idxlist idxlist, Xt_idx const * indices,
                               Xt_count num_indices, int * positions,
                               int single_match_only);

static Xt_idx
idxvec_get_min_index(Xt_idxlist idxlist);

static Xt_idx
idxvec_get_max_index(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxvec_vtable = {
  .delete                      = idxvec_delete,
  .get_pack_size               = idxvec_get_pack_size,
  .pack                        = idxvec_pack,
  .get_intersection            = idxvec_get_intersection,
  .copy                        = idxvec_copy,
  .get_num_indices             = idxvec_get_num_indices,
  .get_indices                 = idxvec_get_indices,
  .get_index_stripes           = idxvec_get_index_stripes,
  .get_index_at_position       = idxvec_get_index_at_position,
  .get_indices_at_positions    = idxvec_get_indices_at_positions,
  .get_position_of_index       = idxvec_get_position_of_index,
  .get_positions_of_indices    = idxvec_get_positions_of_indices,
  .get_position_of_index_off   = idxvec_get_position_of_index_off,
  .get_positions_of_indices_off = NULL,
  .get_min_index               = idxvec_get_min_index,
  .get_max_index               = idxvec_get_max_index
};

// index vector data structure
struct xt_idxvec {

   const struct xt_idxlist_vtable *vtable;

   Xt_idx *vector;
   Xt_count num_indices;

   // internal array used to optimise access to vector data
   Xt_idx *sorted_vector;        // sorted version of vector
   int    *sorted_vec_positions; // original positions of the
                                 // indices in sorted_vector
  /*
    we have the following relations:
    sorted_vector[i-1] <= sorted_vector[i],
    vector[sorted_vec_positions[i]] = sorted_vector[i]
   */
};

Xt_idxlist xt_idxvec_new(Xt_idx *idxvec, Xt_count num_indices) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = xmalloc(1 * sizeof(*idxvec_data));

   idxvec_data->vtable = &idxvec_vtable;

   if (num_indices != 0) {

      Xt_count vector_size;

      vector_size = num_indices * sizeof(*(idxvec_data->vector));

      idxvec_data->vector = xmalloc(vector_size);

      memcpy(idxvec_data->vector, idxvec, vector_size);

   } else {

      idxvec_data->vector = NULL;
   }

   idxvec_data->num_indices = num_indices;
   idxvec_data->sorted_vector = NULL;
   idxvec_data->sorted_vec_positions = NULL;

   return (void*)idxvec_data;
}

static int stripes_overlap(struct Xt_stripe a, struct Xt_stripe b) {

  Xt_idx a_min, a_max, b_min, b_max;

  if (a.stride >= 0)
    a_min = a.start, a_max = a.start + (a.nstrides-1) * a.stride;
  else
    a_max = a.start, a_min = a.start + (a.nstrides-1) * a.stride;

  if (b.stride >= 0)
    b_min = b.start, b_max = b.start + (b.nstrides-1) * b.stride;
  else
    b_max = b.start, b_min = b.start + (b.nstrides-1) * b.stride;

  return MIN(a_max, b_max) >= MAX(a_min, b_min);

}

static int
decode_stripe(struct Xt_stripe stripe, Xt_idx * sorted_vector,
              int * sorted_vec_pos, int pos_offset) {

  Xt_count i;

  if (stripe.stride >= 0) {
    for (i = 0; i < stripe.nstrides; ++i) {
      sorted_vector[i] = stripe.start + i * stripe.stride;
      sorted_vec_pos[i] = pos_offset + i;
    }
  } else {
    for (i = 0; i < stripe.nstrides; ++i) {
        Xt_count j = stripe.nstrides-i-1;
        sorted_vector[j] = stripe.start + i * stripe.stride;
        sorted_vec_pos[j] = pos_offset + i;
      }
  }

  return stripe.nstrides;
}

static void
generate_sorted_vector_from_stripes(struct Xt_stripe * stripes, Xt_count num_stripes,
                                    struct xt_idxvec * idxvec) {

  if (num_stripes == 0) {
    idxvec->sorted_vector = NULL;
    idxvec->sorted_vec_positions = NULL;
    return;
  }

  Xt_idx * stripe_min;
  int * sorted_stripe_min_pos;

  stripe_min = xmalloc(num_stripes * sizeof(*stripe_min));
  sorted_stripe_min_pos = xmalloc(num_stripes * sizeof(*sorted_stripe_min_pos));

  Xt_count i, j;

  for(i = 0; i < num_stripes; ++i)
    if (stripes[i].stride >= 0)
      stripe_min[i] = stripes[i].start;
    else
      stripe_min[i] = stripes[i].start + stripes[i].stride * (stripes[i].nstrides - 1);

  xt_quicksort_index(stripe_min, num_stripes, sorted_stripe_min_pos, 1);

  int * sorted_pos_prefix_sum, * orig_pos_prefix_sum;

  sorted_pos_prefix_sum = xmalloc(num_stripes * sizeof(*sorted_pos_prefix_sum));
  orig_pos_prefix_sum = xmalloc(num_stripes * sizeof(*orig_pos_prefix_sum));

  orig_pos_prefix_sum[0] = 0;
  for (i = 1; i < num_stripes; ++i)
    orig_pos_prefix_sum[i] = orig_pos_prefix_sum[i-1] + stripes[i-1].nstrides;

  for (i = 0; i < num_stripes; ++i)
    sorted_pos_prefix_sum[i] = orig_pos_prefix_sum[sorted_stripe_min_pos[i]];

  free(orig_pos_prefix_sum);

  int * overlap_flag; // does i'th stripe overlap with (i+1)'th stripe

  overlap_flag = xmalloc(num_stripes * sizeof(*overlap_flag));

  for (i = 0; i < num_stripes - 1; ++i)
    overlap_flag[i] = stripes_overlap(stripes[sorted_stripe_min_pos[i]],
                                      stripes[sorted_stripe_min_pos[i+1]]);
  overlap_flag[num_stripes - 1] = 0;

  idxvec->sorted_vector = xmalloc(idxvec->num_indices *
                                  sizeof(*(idxvec->sorted_vector)));
  idxvec->sorted_vec_positions = xmalloc(idxvec->num_indices *
                                         sizeof(*(idxvec->sorted_vec_positions)));;

  Xt_count offset = 0;

  i = 0;
  while (i < num_stripes) {

    int do_overlap;

    do_overlap = overlap_flag[i];

    Xt_count num_selection = 1;

    while(i + num_selection < num_stripes &&
          overlap_flag[i + num_selection] == do_overlap) ++num_selection;

    num_selection += do_overlap;

    Xt_count curr_offset = 0;

    for (j = 0; j < num_selection; ++j)
      curr_offset += decode_stripe(stripes[sorted_stripe_min_pos[i+j]],
                                   idxvec->sorted_vector + offset + curr_offset,
                                   idxvec->sorted_vec_positions + offset + curr_offset,
                                   sorted_pos_prefix_sum[i+j]);

    if (do_overlap)
      xt_quicksort_index(idxvec->sorted_vector + offset, curr_offset,
                         idxvec->sorted_vec_positions + offset, 0);
      

    offset += curr_offset;
    i += num_selection;
  }

  free(sorted_pos_prefix_sum);
  free(overlap_flag);
  free(sorted_stripe_min_pos);
  free(stripe_min);
}

Xt_idxlist 
xt_idxvec_from_stripes_new(struct Xt_stripe * stripes, Xt_count num_stripes) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = xmalloc(1 * sizeof(*idxvec_data));

   idxvec_data->vtable = &idxvec_vtable;

   Xt_count i, j, k;
   Xt_count num_indices = 0;

   for (i = 0; i < num_stripes; ++i)
      num_indices += stripes[i].nstrides;

   Xt_idx * indices;

   indices = xmalloc(num_indices * sizeof(*indices));

   k = 0;
   for (i = 0; i < num_stripes; ++i)
      for (j = 0; j < stripes[i].nstrides; ++j)
         indices[k++] = stripes[i].start + j * stripes[i].stride;

   idxvec_data->vector = indices;
   idxvec_data->num_indices = num_indices;

   generate_sorted_vector_from_stripes(stripes, num_stripes, idxvec_data);

   return (void*)idxvec_data;
}

static void idxvec_delete(Xt_idxlist data) {

   if (((struct xt_idxvec*)data)->sorted_vector !=
       ((struct xt_idxvec*)data)->vector)
      free(((struct xt_idxvec*)data)->sorted_vector);
   free(((struct xt_idxvec*)data)->sorted_vec_positions);
   free(((struct xt_idxvec*)data)->vector);
   free(data);
}

static size_t idxvec_get_pack_size(Xt_idxlist data, MPI_Comm comm) {

   struct xt_idxvec *idxvec = (struct xt_idxvec*)data;

   int size_xt_idx, size_xt_count, size_int_type;

   mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);
   mpi_err_handler(MPI_Pack_size(1, Xt_count_dt, comm, &size_xt_count), comm);
   mpi_err_handler(MPI_Pack_size(idxvec->num_indices, Xt_idx_dt, comm,
                                 &size_xt_idx), comm);

   return size_xt_count + size_xt_idx + size_int_type;
}

void idxvec_pack(Xt_idxlist data, void *buffer, int buffer_size,
                 int *position, MPI_Comm comm) {

  assert(data);
  struct xt_idxvec *idxvec = (struct xt_idxvec*)data;
  int type = VECTOR;

  assert(idxvec->num_indices <= INT_MAX);

  mpi_err_handler(MPI_Pack(&(type), 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(&(idxvec->num_indices), 1, Xt_count_dt, buffer,
                           buffer_size, position, comm), comm);
  if (idxvec->num_indices != 0)
    mpi_err_handler(MPI_Pack(idxvec->vector, idxvec->num_indices,
                             Xt_idx_dt, buffer,
                             buffer_size, position, comm), comm);
}

Xt_idxlist xt_idxvec_unpack(void *buffer, int buffer_size, int *position,
                            MPI_Comm comm) {

   struct xt_idxvec *idxvec;
   Xt_count num_indices;

   idxvec = xmalloc(1 * sizeof(*idxvec));

   idxvec->vtable = &idxvec_vtable;

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              &num_indices, 1, Xt_count_dt, comm), comm);
   assert(num_indices >= 0);
   idxvec->num_indices = num_indices;

   if (num_indices != 0) {

      idxvec->vector = xmalloc(idxvec->num_indices * sizeof(*(idxvec->vector)));
      mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                                 idxvec->vector, num_indices,
                                 Xt_idx_dt, comm),comm);
   } else {

      idxvec->vector = NULL;
   }

   idxvec->sorted_vector = NULL;
   idxvec->sorted_vec_positions = NULL;

   return (Xt_idxlist)idxvec;
}

static Xt_idx *
get_sorted_vector(struct xt_idxvec * idxvec) {

  if (idxvec->sorted_vector != NULL)
    return idxvec->sorted_vector;

  idxvec->sorted_vec_positions = xmalloc(idxvec->num_indices *
                                         sizeof(*(idxvec->sorted_vec_positions)));

  Xt_count i;

  // check if we are already sorted:
  for (i = 1; i < idxvec->num_indices; ++i)
    if (idxvec->vector[i-1] > idxvec->vector[i])
      break;

  // we are done if we are already sorted
  if (i == idxvec->num_indices) {
    // gen id-map:
    for (i = 0; i < idxvec->num_indices; ++i)
      idxvec->sorted_vec_positions[i] = i;    
    return idxvec->sorted_vector = idxvec->vector;
  }
  Xt_idx *sorted_vector;

  sorted_vector = xmalloc(idxvec->num_indices * sizeof(*sorted_vector));

  memcpy(sorted_vector, idxvec->vector,
         idxvec->num_indices * sizeof(*sorted_vector));

  xt_quicksort_index(sorted_vector, idxvec->num_indices,
                  idxvec->sorted_vec_positions, 1);

  idxvec->sorted_vector = sorted_vector;

  return sorted_vector;
}

static Xt_idxlist
idxvec_get_intersection(Xt_idxlist idxlist_src, Xt_idxlist idxlist_dst) {

  if (idxlist_dst->vtable != &idxvec_vtable) return NULL;

  // both lists are index vectors:

  struct xt_idxvec *idxvec_src, *idxvec_dst;

  idxvec_src = (struct xt_idxvec*)idxlist_src;
  idxvec_dst = (struct xt_idxvec*)idxlist_dst;

  struct xt_idxvec *inter_vector;

  inter_vector = xmalloc(1 * sizeof(*inter_vector));

  inter_vector->num_indices = 0;
  inter_vector->vector =
    xmalloc(idxvec_dst->num_indices * sizeof(*(inter_vector->vector)));

  Xt_idx *sorted_src_vector, *sorted_dst_vector;

  // get sorted indices of source and destination

  sorted_src_vector = get_sorted_vector(idxvec_src);
  sorted_dst_vector = get_sorted_vector(idxvec_dst);

  Xt_count i, j;

  // compute the intersection

  for (i = 0, j = 0; i < idxvec_dst->num_indices; ++i) {

    while (j < idxvec_src->num_indices &&
           sorted_src_vector[j] < sorted_dst_vector[i]) ++j;

    if (j >= idxvec_src->num_indices) break;

    if (sorted_src_vector[j] == sorted_dst_vector[i])
      inter_vector->vector[inter_vector->num_indices++] = sorted_dst_vector[i];
  }

  if (inter_vector->num_indices) {

    inter_vector->vector = xrealloc(inter_vector->vector,
                                    inter_vector->num_indices *
                                    sizeof(*(inter_vector->vector)));
  } else {
    free(inter_vector->vector);
    inter_vector->vector = NULL;
  }

  inter_vector->vtable = &idxvec_vtable;
  inter_vector->sorted_vector = NULL;
  inter_vector->sorted_vec_positions = NULL;

  return (Xt_idxlist)inter_vector;
}

static Xt_idxlist
idxvec_copy(Xt_idxlist idxlist) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = (struct xt_idxvec *)idxlist;

   return xt_idxvec_new(idxvec_data->vector, idxvec_data->num_indices);
}

static Xt_count
idxvec_get_num_indices(Xt_idxlist idxlist) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = (struct xt_idxvec *)idxlist;

   return idxvec_data->num_indices;
}

static void
idxvec_get_indices(Xt_idxlist idxlist, Xt_idx *indices) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = (struct xt_idxvec *)idxlist;

   memcpy(indices, idxvec_data->vector,
          idxvec_data->num_indices * sizeof(*indices));
}

static void
idxvec_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                         Xt_count * num_stripes) {

  struct xt_idxvec *idxvec_data;

  idxvec_data = (struct xt_idxvec *)idxlist;

  xt_idxstripes_convert_to_stripes(idxvec_data->vector,
                                   idxvec_data->num_indices,
                                   stripes, num_stripes);
}

static int
idxvec_get_index_at_position(Xt_idxlist idxlist, int position, Xt_idx * index) {

  struct xt_idxvec *idxvec_data;

  idxvec_data = (struct xt_idxvec *)idxlist;

  if (position < 0 || position >= idxvec_data->num_indices)
    return 1;

  *index = idxvec_data->vector[position];

  return 0;
}

static int
idxvec_get_indices_at_positions(Xt_idxlist idxlist, int *positions, Xt_count num_pos, Xt_idx *index, Xt_idx undef_idx) {

  struct xt_idxvec *idxvec;

  idxvec = (struct xt_idxvec *)idxlist;
  int p;
  int max_p = idxvec->num_indices-1;
  Xt_idx *v = idxvec->vector;
  int undef_count = 0;
  for (int ip = 0; ip < num_pos; ip++) {
    p = positions[ip];
    if (p < 0 || p > max_p) {
      index[ip] = undef_idx;
      undef_count++;
      continue;
    }
    index[ip] = v[p];
  }

  return undef_count;
}

static int
idxvec_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                 int * position, int offset) {

  struct xt_idxvec *idxvec_data;

  idxvec_data = (struct xt_idxvec *)idxlist;

  if ((offset < 0) || (offset >= idxvec_data->num_indices))
    return 1;

  Xt_idx * sorted_vector;

  sorted_vector = get_sorted_vector(idxvec_data);

  if ((index < sorted_vector[0]) ||
      (index > sorted_vector[idxvec_data->num_indices-1]))
    return 1;

  Xt_count lb, ub;

  // bisection to find one matching position:
  lb = 0;
  ub = idxvec_data->num_indices - 1;

  while (sorted_vector[lb] < index) {

    Xt_count middle = (ub + lb + 1) / 2;

    if (sorted_vector[middle] <= index)
      lb = middle;
    else
      if (ub == middle)
        return 1;
      else
        ub = middle;
  }

  // find left most match:
  while (lb > 0 && idxvec_data->sorted_vector[lb-1] == index) --lb;

  // go forward until offset condition is satisfied:
  while (lb < idxvec_data->num_indices - 1 &&  // boundary condition
         idxvec_data->sorted_vec_positions[lb] < offset && // ignore postions left of offset
         idxvec_data->sorted_vector[lb] == index) ++lb;  // check if index is valid

  // check if position is invalid:
  if (lb >= idxvec_data->num_indices ||
      sorted_vector[lb] != index)
    return 1; // failure

  // result:
  *position = idxvec_data->sorted_vec_positions[lb];
  return 0;
}

static int
idxvec_get_position_of_index(Xt_idxlist idxlist, Xt_idx index, int * position) {

  return idxvec_get_position_of_index_off(idxlist, index, position, 0);
}

static int idx_vec_is_sorted(Xt_idx const *idx, unsigned n) {

  if (n<2) return 1;

  for (unsigned i = 1; i < n; i++)
    if (idx[i] < idx[i-1]) return 0;

  return 1;
}

int idxvec_get_positions_of_indices(Xt_idxlist body_idxlist, Xt_idx const * selection_idx,
                                   Xt_count num_selection, int * positions,
                                   int single_match_only) {

  if (num_selection == 0) return 0;

  int selection_is_ordered = idx_vec_is_sorted(selection_idx, num_selection);
  /// \todo try linear scan of sorted data instead (requires performance test first)

  Xt_idx const *sorted_selection;
  int *sorted_selection_pos = NULL;
  Xt_idx *tmp_idx = NULL;

  if (selection_is_ordered) {
    sorted_selection = selection_idx;
    sorted_selection_pos = NULL; // not used
  } else {
    size_t memsize;
    memsize = num_selection * sizeof(*sorted_selection);
    tmp_idx = xmalloc(memsize);
    memcpy(tmp_idx, selection_idx, memsize);

    sorted_selection_pos = xmalloc(num_selection * sizeof(*sorted_selection_pos));
    xt_quicksort_index(tmp_idx, num_selection, sorted_selection_pos, 1);
    sorted_selection = tmp_idx;
  }

  // on the target side we want single_match_only,
  // on the source side we don't
  struct xt_idxvec *body_idxvec;
  body_idxvec = (struct xt_idxvec *)body_idxlist;
  Xt_idx *sorted_body = get_sorted_vector(body_idxvec);
  int *sorted_body_pos = body_idxvec->sorted_vec_positions;
  int body_num = body_idxvec->num_indices;

  Xt_idx isel, last_isel;
  int lb_start, post_match_step, lb, ub;

  if (single_match_only)
    post_match_step = 1;
  else
    post_match_step = 0;

  lb_start = 0;
  last_isel = sorted_selection[0];
  for (int i=0; i<num_selection; i++) {
    isel = sorted_selection[i];
    assert(isel >= last_isel);

    // bisection to find one matching position:
    lb = lb_start;
    ub = body_num - 1;
    while (ub-lb>1) {
      int middle = (ub + lb + 1) / 2;
      if (sorted_body[middle] <= isel)
        lb = middle;
      else
        ub = middle;
    }

    // search is now narrowed to two positions, select one of them:
    if (isel == sorted_body[ub]) {
      lb = ub;
    } else {
      if (isel != sorted_body[lb]) return 1; // no match
    }

    // find left most match >= lb_start
    while (lb > lb_start && sorted_body[lb-1] == isel) lb--;

    // check if position is invalid:
    if (sorted_body[lb] != isel) return 1; // failure

    // result:
    int j;
    if (selection_is_ordered)
      j = i;
    else
      j = sorted_selection_pos[i];
    positions[j] = sorted_body_pos[lb];

    // prepare next search:
    lb_start = lb + post_match_step;
  }

  if (tmp_idx) free(tmp_idx);
  if (sorted_selection_pos) free(sorted_selection_pos);

  return 0;
}

static Xt_idx
idxvec_get_min_index(Xt_idxlist idxlist) {

  struct xt_idxvec *idxvec_data;

  idxvec_data = (struct xt_idxvec *)idxlist;

  return get_sorted_vector(idxvec_data)[0];
}

static Xt_idx
idxvec_get_max_index(Xt_idxlist idxlist) {

  struct xt_idxvec *idxvec_data;

  idxvec_data = (struct xt_idxvec *)idxlist;

  return get_sorted_vector(idxvec_data)[idxvec_data->num_indices-1];
}
