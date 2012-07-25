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

static void
idxvec_delete(xt_idxlist data);

static size_t
idxvec_get_pack_size(xt_idxlist data, MPI_Comm comm);

static void
idxvec_pack(xt_idxlist data, void *buffer, int buffer_size,
            int *position, MPI_Comm comm);

static xt_idxlist
idxvec_get_intersection(xt_idxlist idxlist_a, xt_idxlist idxlist_b);

static xt_idxlist
idxvec_copy(xt_idxlist idxlist);

static size_t
idxvec_get_num_indices(xt_idxlist idxlist);

static void
idxvec_get_indices(xt_idxlist idxlist, Xt_idx *indices);

static void
idxvec_get_index_stripes(xt_idxlist idxlist, struct Xt_stripe ** stripes,
                         unsigned * num_stripes);

static int
idxvec_get_index(Xt_idxlist idxlist, int position, Xt_idx * index);

static int
idxvec_get_position_of_index(Xt_idxlist idxlist, Xt_idx index, int * position);

static int
idxvec_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                 int * position, int offset);

static int
idxvec_get_position_of_indices(Xt_idxlist idxlist, Xt_idx * indices,
                               unsigned num_indices, int * positions,
                               int single_match_only);

static const struct xt_idxlist_vtable idxvec_vtable = {
              .delete = idxvec_delete,
              .get_pack_size =  idxvec_get_pack_size,
              .pack = idxvec_pack,
              .get_intersection = idxvec_get_intersection,
              .copy = idxvec_copy,
              .get_num_indices = idxvec_get_num_indices,
              .get_indices = idxvec_get_indices,
              .get_index_stripes = idxvec_get_index_stripes,
              .get_index =  idxvec_get_index,
              .get_position_of_index = idxvec_get_position_of_index,
              .get_position_of_indices = idxvec_get_position_of_indices,
              .get_position_of_index_off = idxvec_get_position_of_index_off,
              .get_position_of_indices_off = NULL};

// index vector data structure
struct xt_idxvec {

   const struct xt_idxlist_vtable *vtable;

   Xt_idx *vector;
   size_t num_indices;

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

struct xt_idxlist *xt_idxvec_new(Xt_idx *idxvec, size_t num_indices) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = xmalloc(1 * sizeof(*idxvec_data));

   idxvec_data->vtable = &idxvec_vtable;

   if (num_indices != 0) {

      size_t vector_size;

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

static void idxvec_delete(struct xt_idxlist *data) {

   if (((struct xt_idxvec*)data)->sorted_vector !=
       ((struct xt_idxvec*)data)->vector)
      free(((struct xt_idxvec*)data)->sorted_vector);
   free(((struct xt_idxvec*)data)->sorted_vec_positions);
   free(((struct xt_idxvec*)data)->vector);
   free(data);
}

static size_t idxvec_get_pack_size(struct xt_idxlist *data, MPI_Comm comm) {

   struct xt_idxvec *idxvec = (struct xt_idxvec*)data;

   int size_unsigned, size_int, size_int_type;

   mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);
   mpi_err_handler(MPI_Pack_size(idxvec->num_indices, Xt_idx_dt, comm,
                                 &size_int), comm);
   mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_unsigned), comm);

   return size_unsigned + size_int + size_int_type;
}

void idxvec_pack(struct xt_idxlist *data, void *buffer, int buffer_size,
                 int *position, MPI_Comm comm) {

  assert(data);
  struct xt_idxvec *idxvec = (struct xt_idxvec*)data;
  int type = VECTOR;
  int num_indices = (int)(idxvec->num_indices);

  assert(idxvec->num_indices <= INT_MAX);

  mpi_err_handler(MPI_Pack(&(type), 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(&num_indices, 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  if (num_indices != 0)
    mpi_err_handler(MPI_Pack(idxvec->vector, (int)idxvec->num_indices,
                             Xt_idx_dt, buffer,
                             buffer_size, position, comm), comm);
}

xt_idxlist xt_idxvec_unpack(void *buffer, int buffer_size, int *position,
                            MPI_Comm comm) {

   struct xt_idxvec *idxvec;
   int num_indices;

   idxvec = xmalloc(1 * sizeof(*idxvec));

   idxvec->vtable = &idxvec_vtable;

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              &num_indices, 1, MPI_INT, comm), comm);
   assert(num_indices >= 0);
   idxvec->num_indices = (size_t)num_indices;

   if (num_indices != 0) {

      idxvec->vector = xmalloc(idxvec->num_indices * sizeof(*(idxvec->vector)));
      mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                                 idxvec->vector, idxvec->num_indices,
                                 Xt_idx_dt, comm),comm);
   } else {

      idxvec->vector = NULL;
   }

   idxvec->sorted_vector = NULL;
   idxvec->sorted_vec_positions = NULL;

   return (xt_idxlist)idxvec;
}

static Xt_idx *
get_sorted_vector(struct xt_idxvec * idxvec) {

  if (idxvec->sorted_vector != NULL)
    return idxvec->sorted_vector;

  idxvec->sorted_vec_positions = xmalloc(idxvec->num_indices *
                                         sizeof(*(idxvec->sorted_vec_positions)));

  unsigned i;

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

static xt_idxlist
idxvec_get_intersection(xt_idxlist idxlist_src, xt_idxlist idxlist_dst) {


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

  size_t i, j;

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

  return (xt_idxlist)inter_vector;
}

static xt_idxlist
idxvec_copy(xt_idxlist idxlist) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = (struct xt_idxvec *)idxlist;

   return xt_idxvec_new(idxvec_data->vector, idxvec_data->num_indices);
}

static size_t
idxvec_get_num_indices(xt_idxlist idxlist) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = (struct xt_idxvec *)idxlist;

   return idxvec_data->num_indices;
}

static void
idxvec_get_indices(xt_idxlist idxlist, Xt_idx *indices) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = (struct xt_idxvec *)idxlist;

   memcpy(indices, idxvec_data->vector,
          idxvec_data->num_indices * sizeof(*indices));
}

static void
idxvec_get_index_stripes(xt_idxlist idxlist, struct Xt_stripe ** stripes,
                         unsigned * num_stripes) {

  struct xt_idxvec *idxvec_data;

  idxvec_data = (struct xt_idxvec *)idxlist;

  xt_idxstripes_convert_to_stripes(idxvec_data->vector,
                                   idxvec_data->num_indices,
                                   stripes, num_stripes);
}

static int
idxvec_get_index(Xt_idxlist idxlist, int position, Xt_idx * index) {

  struct xt_idxvec *idxvec_data;

  idxvec_data = (struct xt_idxvec *)idxlist;

  if (position < 0 || position >= idxvec_data->num_indices)
    return 1;

  *index = idxvec_data->vector[position];

  return 0;
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

  int lb, ub;

  // bisection to find one matching position:
  lb = 0;
  ub = idxvec_data->num_indices - 1;

  while (sorted_vector[lb] < index) {

    int middle = (ub + lb + 1) / 2;

    if (sorted_vector[middle] <= index)
      lb = middle;
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

static int idx_vec_is_sorted(Xt_idx *idx, int n) {
  if (n<2) return 1;

  for (int i=1; i<n; i++) {
    if (idx[i]<idx[i-1]) return 0;
  }
  return 1;
}

int idxvec_get_position_of_indices(Xt_idxlist body_idxlist, Xt_idx * selection_idx,
                                   unsigned num_selection, int * positions,
                                   int single_match_only) {

  if (num_selection == 0) return 0;

  int selection_is_ordered = idx_vec_is_sorted(selection_idx, num_selection);

  Xt_idx *sorted_selection;
  int *sorted_selection_pos;

  if (selection_is_ordered) {
    sorted_selection = selection_idx;
    sorted_selection_pos = NULL; // not used
  } else {
    int memsize;
    memsize = num_selection * sizeof(*sorted_selection);
    sorted_selection = xmalloc(memsize);
    memcpy(sorted_selection, selection_idx, memsize);

    sorted_selection_pos = xmalloc(num_selection * sizeof(*sorted_selection_pos));
    xt_quicksort_index(sorted_selection, num_selection, sorted_selection_pos, 1);
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
    while (sorted_body[lb] < isel) {
      int middle = (ub + lb + 1) / 2;
      if (sorted_body[middle] <= isel)
        lb = middle;
      else
        ub = middle;
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

  if (sorted_selection != selection_idx) free(sorted_selection);
  if (sorted_selection_pos) free(sorted_selection_pos);

  return 0;
}
