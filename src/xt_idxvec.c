/**
 * @file xt_idxvec.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://redmine.dkrz.de/doc/yaxt/html/index.html
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xt/xt_core.h"
#include "xt/xt_idxlist.h"
#include "xt_idxlist_internal.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_mpi.h"
#include "xt_idxlist_unpack.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"
#include "xt/quicksort.h"
#include "instr.h"

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

static Xt_int
idxvec_get_num_indices(Xt_idxlist idxlist);

static void
idxvec_get_indices(Xt_idxlist idxlist, Xt_int *indices);

static Xt_int const*
idxvec_get_indices_const(Xt_idxlist idxlist);

static void
idxvec_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                         Xt_int * num_stripes);

static int
idxvec_get_index_at_position(Xt_idxlist idxlist, int position, Xt_int * index);

static int
idxvec_get_indices_at_positions(Xt_idxlist idxlist, int *positions, Xt_int num, Xt_int *index, Xt_int undef_idx);

static int
idxvec_get_position_of_index(Xt_idxlist idxlist, Xt_int index, int * position);

static int
idxvec_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                 int * position, int offset);

static int
idxvec_get_positions_of_indices(Xt_idxlist idxlist, const Xt_int *indices,
                                Xt_int num_indices, int *positions,
                                int single_match_only);

static Xt_int
idxvec_get_min_index(Xt_idxlist idxlist);

static Xt_int
idxvec_get_max_index(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxvec_vtable = {
  .delete                      = idxvec_delete,
  .get_pack_size               = idxvec_get_pack_size,
  .pack                        = idxvec_pack,
  .get_intersection            = idxvec_get_intersection,
  .copy                        = idxvec_copy,
  .get_num_indices             = idxvec_get_num_indices,
  .get_indices                 = idxvec_get_indices,
  .get_indices_const           = idxvec_get_indices_const,
  .get_index_stripes           = idxvec_get_index_stripes,
  .get_index_at_position       = idxvec_get_index_at_position,
  .get_indices_at_positions    = idxvec_get_indices_at_positions,
  .get_position_of_index       = idxvec_get_position_of_index,
  .get_positions_of_indices    = idxvec_get_positions_of_indices,
  .get_position_of_index_off   = idxvec_get_position_of_index_off,
  .get_positions_of_indices_off = NULL,
  .get_min_index               = idxvec_get_min_index,
  .get_max_index               = idxvec_get_max_index,
  .get_bounding_box             = NULL
};

typedef struct Xt_idxvec_ *Xt_idxvec;

// index vector data structure
struct Xt_idxvec_ {

  struct Xt_idxlist_ parent;

   Xt_int *vector;
   Xt_uint num_indices;

   // internal array used to optimise access to vector data
   Xt_int *sorted_vector;        // sorted version of vector
   int    *sorted_vec_positions; // original positions of the
                                 // indices in sorted_vector
  /*
    we have the following relations:
    sorted_vector[i-1] <= sorted_vector[i],
    vector[sorted_vec_positions[i]] = sorted_vector[i]
   */
};


Xt_idxlist xt_idxvec_new(const Xt_int *idxvec, Xt_int num_indices) {
  INSTR_DEF(t_idxvec_new,"xt_idxvec_new");
  INSTR_START(t_idxvec_new);
  if (num_indices == 0)
    return xt_idxempty_new();
  Xt_idxvec restrict idxvec_data = xmalloc(sizeof (*idxvec_data));

  Xt_idxlist_init(&idxvec_data->parent, &idxvec_vtable);

  if (num_indices > 0) {

    size_t vector_size = (size_t)num_indices * sizeof(*(idxvec_data->vector));

    idxvec_data->vector = xmalloc(vector_size);

    memcpy(idxvec_data->vector, idxvec, vector_size);
  } else {
    die("number of indices passed to xt_idxvec_new must not be negative!");
  }
  idxvec_data->num_indices = (Xt_uint)num_indices;
  idxvec_data->sorted_vector = NULL;
  idxvec_data->sorted_vec_positions = NULL;

  INSTR_STOP(t_idxvec_new);
  return (void*)idxvec_data;
}

static int stripes_overlap(struct Xt_stripe a, struct Xt_stripe b) {

  Xt_int a_min, a_max, b_min, b_max;

  if (a.stride >= 0)
    a_min = a.start, a_max = (Xt_int)(a.start + (a.nstrides-1) * a.stride);
  else
    a_max = a.start, a_min = (Xt_int)(a.start + (a.nstrides-1) * a.stride);

  if (b.stride >= 0)
    b_min = b.start, b_max = (Xt_int)(b.start + (b.nstrides-1) * b.stride);
  else
    b_max = b.start, b_min = (Xt_int)(b.start + (b.nstrides-1) * b.stride);

  return MIN(a_max, b_max) >= MAX(a_min, b_min);

}

static int
decode_stripe(struct Xt_stripe stripe, Xt_int * sorted_vector,
              int * sorted_vec_pos, int pos_offset) {

  Xt_int i;

  if (stripe.stride >= 0) {
    for (i = 0; i < stripe.nstrides; ++i) {
      sorted_vector[i] = (Xt_int)(stripe.start + i * stripe.stride);
      sorted_vec_pos[i] = pos_offset + (int)i;
    }
  } else {
    for (i = 0; i < stripe.nstrides; ++i) {
      Xt_int j = (Xt_int)(stripe.nstrides - i - 1);
      sorted_vector[j] = (Xt_int)(stripe.start + i * stripe.stride);
      sorted_vec_pos[j] = pos_offset + (int)i;
    }
  }

  return (int)stripe.nstrides;
}

static void
generate_sorted_vector_from_stripes(const struct Xt_stripe stripes[],
                                    Xt_int num_stripes,
                                    Xt_idxvec idxvec) {

  if (num_stripes == 0) {
    idxvec->sorted_vector = NULL;
    idxvec->sorted_vec_positions = NULL;
    return;
  }

  Xt_int * stripe_min;
  int * sorted_stripe_min_pos;

  stripe_min = xmalloc((size_t)num_stripes * sizeof(*stripe_min));
  sorted_stripe_min_pos = xmalloc((size_t)num_stripes * sizeof(*sorted_stripe_min_pos));

  Xt_int i, j;

  for(i = 0; i < num_stripes; ++i)
    if (stripes[i].stride >= 0)
      stripe_min[i] = stripes[i].start;
    else
      stripe_min[i] = (Xt_int)(stripes[i].start
                               + stripes[i].stride
                               * (stripes[i].nstrides - 1));

  xt_quicksort_index(stripe_min, (int)num_stripes, sorted_stripe_min_pos, 1);

  int * sorted_pos_prefix_sum, * orig_pos_prefix_sum;

  sorted_pos_prefix_sum
    = xmalloc((size_t)num_stripes * sizeof(*sorted_pos_prefix_sum));
  orig_pos_prefix_sum
    = xmalloc((size_t)num_stripes * sizeof(*orig_pos_prefix_sum));

  orig_pos_prefix_sum[0] = 0;
  for (i = 1; i < num_stripes; ++i)
    orig_pos_prefix_sum[i] = orig_pos_prefix_sum[i-1] + (int)stripes[i-1].nstrides;

  for (i = 0; i < num_stripes; ++i)
    sorted_pos_prefix_sum[i] = orig_pos_prefix_sum[sorted_stripe_min_pos[i]];

  free(orig_pos_prefix_sum);

  int * overlap_flag; // does i'th stripe overlap with (i+1)'th stripe

  overlap_flag = xmalloc((size_t)num_stripes * sizeof(*overlap_flag));

  for (i = 0; i < num_stripes - 1; ++i)
    overlap_flag[i] = stripes_overlap(stripes[sorted_stripe_min_pos[i]],
                                      stripes[sorted_stripe_min_pos[i+1]]);
  overlap_flag[num_stripes - 1] = 0;

  idxvec->sorted_vector = xmalloc((size_t)idxvec->num_indices *
                                  sizeof(*(idxvec->sorted_vector)));
  idxvec->sorted_vec_positions = xmalloc((size_t)idxvec->num_indices *
                                         sizeof(*(idxvec->sorted_vec_positions)));;

  Xt_int offset = 0;

  i = 0;
  while (i < num_stripes) {

    int do_overlap;

    do_overlap = overlap_flag[i];

    Xt_int num_selection = 1;

    while(i + num_selection < num_stripes &&
          overlap_flag[i + num_selection] == do_overlap) ++num_selection;

    num_selection = (Xt_int)(num_selection + do_overlap);

    Xt_int curr_offset = 0;

    for (j = 0; j < num_selection; ++j)
      curr_offset
        = (Xt_int)(curr_offset
                   + decode_stripe(stripes[sorted_stripe_min_pos[i+j]],
                                   idxvec->sorted_vector + offset
                                   + curr_offset,
                                   idxvec->sorted_vec_positions + offset
                                   + curr_offset,
                                   sorted_pos_prefix_sum[i+j]));

    if (do_overlap)
      xt_quicksort_index(idxvec->sorted_vector + offset, (int)curr_offset,
                         idxvec->sorted_vec_positions + offset, 0);

    offset = (Xt_int)(offset + curr_offset);
    i = (Xt_int)(i + num_selection);
  }

  free(sorted_pos_prefix_sum);
  free(overlap_flag);
  free(sorted_stripe_min_pos);
  free(stripe_min);
}

Xt_idxlist
xt_idxvec_from_stripes_new(const struct Xt_stripe stripes[],
                           Xt_int num_stripes) {

   Xt_idxvec idxvec_data = xmalloc(sizeof (*idxvec_data));

   Xt_idxlist_init(&idxvec_data->parent, &idxvec_vtable);

   Xt_uint num_indices = 0;

   /* FIXME: guard against overflow here */
   for (Xt_int i = 0; i < num_stripes; ++i)
     num_indices = (Xt_uint)(num_indices + (Xt_uint)stripes[i].nstrides);

   Xt_int *indices = xmalloc((size_t)num_indices * sizeof(*indices));

   Xt_int k = -1;
   for (Xt_int i = 0; i < num_stripes; ++i)
      for (Xt_int j = 0; j < stripes[i].nstrides; ++j)
        indices[++k] = (Xt_int)(stripes[i].start + j * stripes[i].stride);

   idxvec_data->vector = indices;
   idxvec_data->num_indices = num_indices;

   generate_sorted_vector_from_stripes(stripes, num_stripes, idxvec_data);

   return (void*)idxvec_data;
}

static void idxvec_delete(Xt_idxlist data) {

   if (((Xt_idxvec)data)->sorted_vector !=
       ((Xt_idxvec)data)->vector)
      free(((Xt_idxvec)data)->sorted_vector);
   free(((Xt_idxvec)data)->sorted_vec_positions);
   free(((Xt_idxvec)data)->vector);
   free(data);
}

static size_t idxvec_get_pack_size(Xt_idxlist data, MPI_Comm comm) {

  Xt_idxvec idxvec = (Xt_idxvec)data;
  int size_xt_idx, size_xt_count, size_int_type;

  assert(idxvec->num_indices <= INT_MAX);

  xt_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);
  xt_mpi_call(MPI_Pack_size(1, Xt_uint_dt, comm, &size_xt_count), comm);
  xt_mpi_call(MPI_Pack_size((int)idxvec->num_indices, Xt_int_dt, comm,
                            &size_xt_idx), comm);

  return (size_t)size_xt_count + (size_t)size_xt_idx + (size_t)size_int_type;
}

void idxvec_pack(Xt_idxlist data, void *buffer, int buffer_size,
                 int *position, MPI_Comm comm) {

  assert(data);
  Xt_idxvec idxvec = (Xt_idxvec)data;
  int type = VECTOR;

  assert(idxvec->num_indices <= INT_MAX);

  xt_mpi_call(MPI_Pack(&(type), 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  xt_mpi_call(MPI_Pack(&(idxvec->num_indices), 1, Xt_uint_dt, buffer,
                           buffer_size, position, comm), comm);
  if (idxvec->num_indices != 0)
    xt_mpi_call(MPI_Pack(idxvec->vector, (int)idxvec->num_indices,
                         Xt_int_dt, buffer,
                         buffer_size, position, comm), comm);
}

Xt_idxlist xt_idxvec_unpack(void *buffer, int buffer_size, int *position,
                            MPI_Comm comm) {

  Xt_idxvec idxvec = xmalloc(sizeof (*idxvec));
  Xt_uint num_indices;


  Xt_idxlist_init(&idxvec->parent, &idxvec_vtable);

  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         &num_indices, 1, Xt_uint_dt, comm), comm);
  assert(num_indices <= INT_MAX);
  idxvec->num_indices = num_indices;

  if (num_indices != 0) {

    idxvec->vector
      = xmalloc((size_t)idxvec->num_indices * sizeof(*(idxvec->vector)));
    xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                           idxvec->vector, (int)num_indices,
                           Xt_int_dt, comm),comm);
  } else {

    idxvec->vector = NULL;
  }

  idxvec->sorted_vector = NULL;
  idxvec->sorted_vec_positions = NULL;

  return (Xt_idxlist)idxvec;
}

static Xt_int *
get_sorted_vector(Xt_idxvec idxvec) {

  if (idxvec->sorted_vector != NULL)
    return idxvec->sorted_vector;

  idxvec->sorted_vec_positions = xmalloc((size_t)idxvec->num_indices *
                                         sizeof(*(idxvec->sorted_vec_positions)));

  Xt_uint i;

  // check if we are already sorted:
  for (i = 1; i < idxvec->num_indices; ++i)
    if (idxvec->vector[i-1] > idxvec->vector[i])
      break;

  // we are done if we are already sorted
  if (i == idxvec->num_indices) {
    // gen id-map:
    for (i = 0; i < idxvec->num_indices; ++i)
      idxvec->sorted_vec_positions[i] = (int)i;
    return idxvec->sorted_vector = idxvec->vector;
  }
  Xt_int *sorted_vector
    = xmalloc((size_t)idxvec->num_indices * sizeof(*sorted_vector));

  memcpy(sorted_vector, idxvec->vector,
         (size_t)idxvec->num_indices * sizeof(*sorted_vector));

  xt_quicksort_index(sorted_vector, (int)idxvec->num_indices,
                     idxvec->sorted_vec_positions, 1);

  idxvec->sorted_vector = sorted_vector;

  return sorted_vector;
}

static Xt_idxlist
idxvec_get_intersection(Xt_idxlist idxlist_src, Xt_idxlist idxlist_dst) {

  if (idxlist_dst->vtable != &idxvec_vtable) return NULL;

  // both lists are index vectors:

  Xt_idxvec idxvec_src = (Xt_idxvec)idxlist_src,
    idxvec_dst = (Xt_idxvec)idxlist_dst;

  Xt_idxvec inter_vector = xmalloc(sizeof (*inter_vector));

  inter_vector->num_indices = 0;
  inter_vector->vector =
    xmalloc((size_t)idxvec_dst->num_indices * sizeof(*(inter_vector->vector)));

  Xt_int *sorted_src_vector, *sorted_dst_vector;

  // get sorted indices of source and destination

  sorted_src_vector = get_sorted_vector(idxvec_src);
  sorted_dst_vector = get_sorted_vector(idxvec_dst);

  // compute the intersection

  for (Xt_uint i = 0, j = 0; i < idxvec_dst->num_indices; ++i) {

    while (j < idxvec_src->num_indices &&
           sorted_src_vector[j] < sorted_dst_vector[i]) ++j;

    if (j >= idxvec_src->num_indices) break;

    if (sorted_src_vector[j] == sorted_dst_vector[i])
      inter_vector->vector[inter_vector->num_indices++] = sorted_dst_vector[i];
  }

  if (inter_vector->num_indices) {

    inter_vector->vector = xrealloc(inter_vector->vector,
                                    (size_t)inter_vector->num_indices *
                                    sizeof(*(inter_vector->vector)));
  } else {
    free(inter_vector->vector);
    inter_vector->vector = NULL;
  }

  Xt_idxlist_init(&inter_vector->parent, &idxvec_vtable);
  inter_vector->sorted_vector = NULL;
  inter_vector->sorted_vec_positions = NULL;

  return (Xt_idxlist)inter_vector;
}

static Xt_idxlist
idxvec_copy(Xt_idxlist idxlist) {

   Xt_idxvec idxvec_data = (Xt_idxvec)idxlist;

   return xt_idxvec_new(idxvec_data->vector, (Xt_int)idxvec_data->num_indices);
}

static Xt_int
idxvec_get_num_indices(Xt_idxlist idxlist) {

   Xt_idxvec idxvec_data = (Xt_idxvec)idxlist;

   return (Xt_int)idxvec_data->num_indices;
}

static void
idxvec_get_indices(Xt_idxlist idxlist, Xt_int *indices) {

   Xt_idxvec idxvec_data = (Xt_idxvec)idxlist;

   memcpy(indices, idxvec_data->vector,
          (size_t)idxvec_data->num_indices * sizeof(*indices));
}

static Xt_int const*
idxvec_get_indices_const(Xt_idxlist idxlist) {
  Xt_idxvec idxvec = (Xt_idxvec)idxlist;

  return idxvec->vector;
}


static void
idxvec_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                         Xt_int * num_stripes) {

  Xt_idxvec idxvec_data = (Xt_idxvec)idxlist;

  xt_idxstripes_convert_to_stripes(idxvec_data->vector,
                                   (Xt_int)idxvec_data->num_indices,
                                   stripes, num_stripes);
}

static int
idxvec_get_index_at_position(Xt_idxlist idxlist, int position, Xt_int * index) {

  Xt_idxvec idxvec_data = (Xt_idxvec)idxlist;

  if (position < 0 || (Xt_uint)position >= idxvec_data->num_indices)
    return 1;

  *index = idxvec_data->vector[position];

  return 0;
}

static int
idxvec_get_indices_at_positions(Xt_idxlist idxlist, int *positions,
                                Xt_int num_pos, Xt_int *index,
                                Xt_int undef_idx) {

  Xt_idxvec idxvec = (Xt_idxvec)idxlist;
  Xt_uint num_indices = idxvec->num_indices;
  Xt_int *v = idxvec->vector;

  int undef_count = 0;
  for (int ip = 0; ip < num_pos; ip++) {
    int p = positions[ip];
    if (p >= 0 && (unsigned)p < num_indices) {
      index[ip] = v[p];
    } else {
      index[ip] = undef_idx;
      undef_count++;
    }
  }

  return undef_count;
}

static int
idxvec_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                 int * position, int offset) {

  Xt_idxvec idxvec_data = (Xt_idxvec)idxlist;

  *position = -1;

  if ((offset < 0) || ((unsigned)offset >= idxvec_data->num_indices))
    return 1;

  Xt_int * sorted_vector;

  sorted_vector = get_sorted_vector(idxvec_data);

  if ((index < sorted_vector[0]) ||
      (index > sorted_vector[idxvec_data->num_indices-1]))
    return 1;

  Xt_uint lb, ub;

  // bisection to find one matching position:
  lb = 0;
  ub = (Xt_uint)(idxvec_data->num_indices - 1);

  while (sorted_vector[lb] < index) {

    Xt_uint middle = (Xt_uint)((Xt_uint)((Xt_uint)(ub + lb) + 1U) / 2U);

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
idxvec_get_position_of_index(Xt_idxlist idxlist, Xt_int index, int * position) {

  return idxvec_get_position_of_index_off(idxlist, index, position, 0);
}

static int idx_vec_is_sorted(Xt_int const *idx, Xt_int n) {

  if (n<2) return 1;

  for (Xt_int i = 1; i < n; i++)
    if (idx[i] < idx[i-1]) return 0;

  return 1;
}

static int
idxvec_get_positions_of_indices(Xt_idxlist body_idxlist,
                                const Xt_int *selection_idx,
                                Xt_int num_selection, int *positions,
                                int single_match_only) {

  if (num_selection == 0) return 0;

  int selection_is_ordered = idx_vec_is_sorted(selection_idx, num_selection);
  /// \todo try linear scan of sorted data instead (requires performance test first)

  Xt_int const *sorted_selection;
  int *sorted_selection_pos = NULL;
  Xt_int *tmp_idx = NULL;

  if (selection_is_ordered) {
    sorted_selection = selection_idx;
    sorted_selection_pos = NULL; // not used
  } else {
    size_t memsize = (size_t)num_selection * sizeof(*sorted_selection);
    tmp_idx = xmalloc(memsize);
    memcpy(tmp_idx, selection_idx, memsize);

    sorted_selection_pos
      = xmalloc((size_t)num_selection * sizeof(*sorted_selection_pos));
    xt_quicksort_index(tmp_idx, (int)num_selection, sorted_selection_pos, 1);
    sorted_selection = tmp_idx;
  }

  /* motivation for usage of single_match_only:
   *  on the target side we want single_match_only,
   *  on the source side we don't
   */
  Xt_idxvec body_idxvec = (Xt_idxvec)body_idxlist;
  Xt_int *sorted_body = get_sorted_vector(body_idxvec);
  int *sorted_body_pos = body_idxvec->sorted_vec_positions;
  int search_end = (int)body_idxvec->num_indices -1;
  int num_unmatched = 0;
  int match_pos;
  Xt_int isel, last_isel;
  int search_start, post_match_step, lb, ub;

  if (single_match_only)
    // after the match we will move on one step in order to avoid matching the same position again
    post_match_step = 1;
  else
    post_match_step = 0;

  search_start = 0;
  last_isel = sorted_selection[0];
  for (int i=0; i<num_selection; i++) {
    if (search_start<=search_end) {

      isel = sorted_selection[i];
      assert(isel >= last_isel);
      // bisection to find one matching position:
      lb = search_start;
      ub = search_end;
      while (ub-lb>1) {
        int middle = (ub + lb + 1) / 2;
        if (sorted_body[middle] <= isel)
          lb = middle;
        else
          ub = middle;
      }

      // search is now narrowed to two positions, select one of them:
      if (isel == sorted_body[lb]) {
        match_pos = lb;
      } else if (isel == sorted_body[ub]){
        match_pos = ub;
      } else {
        match_pos = -1;
        num_unmatched++;
      }

      if (match_pos>0) {
        // find left most match >= search_start (bisection can lead to any match >= search_start)
        while (match_pos > search_start && sorted_body[match_pos-1] == isel) match_pos--;
      }

      if (match_pos>=0) {
        /// \todo Do we really need this check?
        // check (again?) if position is invalid:
        if (sorted_body[match_pos] != isel) {
          match_pos = -1;
          num_unmatched++;
        }
      }

    } else {
      match_pos = -1;
      num_unmatched++;
    }

    // result:
    int j;
    if (selection_is_ordered)
      j = i;
    else
      j = sorted_selection_pos[i];

    // update positions and prepare next search:
    if (match_pos >= 0) {
      positions[j] = sorted_body_pos[match_pos];
      search_start = match_pos + post_match_step;
    } else {
      positions[j] = -1;
    }

  }

  if (tmp_idx) free(tmp_idx);
  if (sorted_selection_pos) free(sorted_selection_pos);

  return num_unmatched;
}

static Xt_int
idxvec_get_min_index(Xt_idxlist idxlist) {

  Xt_idxvec idxvec_data = (Xt_idxvec)idxlist;

  if (!idxvec_data->num_indices)
    die("idxvec_get_min_index: empty index vector");

  return get_sorted_vector(idxvec_data)[0];
}

static Xt_int
idxvec_get_max_index(Xt_idxlist idxlist) {

  Xt_idxvec idxvec_data = (Xt_idxvec)idxlist;

  if (!idxvec_data->num_indices)
    die("idxvec_get_max_index: empty index vector");

  return get_sorted_vector(idxvec_data)[idxvec_data->num_indices-1];
}
