/**
 * @file xt_idxvec.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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
 * URL: https://doc.redmine.dkrz.de/yaxt/html/
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xt/xt_core.h"
#include "xt/xt_idxlist.h"
#include "xt_idxlist_internal.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_idxvec.h"
#include "xt_idxvec_internal.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_mpi.h"
#include "xt_idxlist_unpack.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"
#include "xt_stripe_util.h"
#include "xt/quicksort.h"
#include "instr.h"

static void
idxvec_delete(Xt_idxlist data);

static size_t
idxvec_get_pack_size(Xt_idxlist data, MPI_Comm comm);

static void
idxvec_pack(Xt_idxlist data, void *buffer, int buffer_size,
            int *position, MPI_Comm comm);

static Xt_idxlist
idxvec_copy(Xt_idxlist idxlist);

static void
idxvec_get_indices(Xt_idxlist idxlist, Xt_int *indices);

static Xt_int const*
idxvec_get_indices_const(Xt_idxlist idxlist);

static void
idxvec_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                         int * num_stripes);

static int
idxvec_get_index_at_position(Xt_idxlist idxlist, int position, Xt_int * index);

static int
idxvec_get_indices_at_positions(Xt_idxlist idxlist, const int *positions,
                                int num, Xt_int *index, Xt_int undef_idx);

static int
idxvec_get_position_of_index(Xt_idxlist idxlist, Xt_int index, int * position);

static int
idxvec_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                 int * position, int offset);

static int
idxvec_get_positions_of_indices(Xt_idxlist idxlist, const Xt_int *indices,
                                int num_indices, int *positions,
                                int single_match_only);

static Xt_int
idxvec_get_min_index(Xt_idxlist idxlist);

static Xt_int
idxvec_get_max_index(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxvec_vtable = {
  .delete                      = idxvec_delete,
  .get_pack_size               = idxvec_get_pack_size,
  .pack                        = idxvec_pack,
  .copy                        = idxvec_copy,
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
  .get_bounding_box            = NULL,
  .idxlist_pack_code           = VECTOR,
};

typedef struct Xt_idxvec_ *Xt_idxvec;

// index vector data structure
struct Xt_idxvec_ {

  struct Xt_idxlist_ parent;

  const Xt_int *vector;

   // internal array used to optimise access to vector data
  const Xt_int *sorted_vector;        // sorted version of vector
   int    *sorted_vec_positions; // original positions of the
                                 // indices in sorted_vector
  /*
    we have the following relations:
    sorted_vector[i-1] <= sorted_vector[i],
    vector[sorted_vec_positions[i]] = sorted_vector[i]
   */
};


Xt_idxlist xt_idxvec_new(const Xt_int *idxvec, int num_indices) {
  INSTR_DEF(t_idxvec_new,"xt_idxvec_new")
  if (num_indices > 0)
    ;
  else if (num_indices == 0)
    return xt_idxempty_new();
  else
    die("number of indices passed to xt_idxvec_new must not be negative!");

  INSTR_START(t_idxvec_new);
  size_t vector_size = (size_t)num_indices * sizeof (idxvec[0]),
    header_size = ((sizeof (struct Xt_idxvec_) + sizeof (Xt_int) - 1)
                   /sizeof (Xt_int)) * sizeof (Xt_int);
  struct Xt_idxvec_ *restrict idxvec_obj = xmalloc(header_size + vector_size);
  Xt_idxlist_init(&idxvec_obj->parent, &idxvec_vtable, num_indices);

  Xt_int *vector_assign = (Xt_int *)(void *)((unsigned char *)idxvec_obj + header_size);
  idxvec_obj->vector = vector_assign;
  memcpy(vector_assign, idxvec, vector_size);
  idxvec_obj->sorted_vector = NULL;
  idxvec_obj->sorted_vec_positions = NULL;

  INSTR_STOP(t_idxvec_new);
  return (void *)idxvec_obj;
}

Xt_idxlist xt_idxvec_prealloc_new(const Xt_int *idxvec, int num_indices)
{
  if (num_indices > 0)
    ;
  else if (num_indices == 0)
    return xt_idxempty_new();
  else
    die("number of indices passed to xt_idxvec_new must not be negative!");
  struct Xt_idxvec_ *restrict idxvec_obj = xmalloc(sizeof (*idxvec_obj));
  Xt_idxlist_init(&idxvec_obj->parent, &idxvec_vtable, num_indices);
  idxvec_obj->vector = idxvec;
  idxvec_obj->sorted_vector = NULL;
  idxvec_obj->sorted_vec_positions = NULL;
  return (void *)idxvec_obj;
}

static int
decode_stripe(struct Xt_stripe stripe, Xt_int * sorted_vector,
              int * sorted_vec_pos, int pos_offset) {

  if (stripe.stride >= 0) {
    for (int i = 0; i < stripe.nstrides; ++i) {
      sorted_vector[i] = (Xt_int)(stripe.start + i * stripe.stride);
      sorted_vec_pos[i] = pos_offset + i;
    }
  } else {
    for (int i = 0; i < stripe.nstrides; ++i) {
      int j = stripe.nstrides - i - 1;
      sorted_vector[j] = (Xt_int)(stripe.start + i * stripe.stride);
      sorted_vec_pos[j] = pos_offset + i;
    }
  }

  return stripe.nstrides;
}

#define MAX(a,b) ((a) >= (b) ? (a) : (b))

static void
generate_sorted_vector_from_stripes(const struct Xt_stripe stripes[],
                                    int num_stripes,
                                    Xt_idxvec idxvec) {

  if (num_stripes == 0) {
    idxvec->sorted_vector = NULL;
    idxvec->sorted_vec_positions = NULL;
    return;
  }

  Xt_int *restrict sorted_vector_assign
    = xmalloc((size_t)idxvec->parent.num_indices
              * sizeof(*(idxvec->sorted_vector)));
  idxvec->sorted_vector = sorted_vector_assign;
  idxvec->sorted_vec_positions
    = xmalloc((size_t)idxvec->parent.num_indices *
              sizeof(*(idxvec->sorted_vec_positions)));

  /* stripe_minmax[0][i] is the minimal index in stripe i at first,
   * later of sorted stripe i, stripe_minmax[1][i] is the
   * corresponding maximal index */
  Xt_int (*restrict stripe_minmax)[num_stripes]
    = xmalloc(2 * sizeof(*stripe_minmax));
  int *restrict sorted_stripe_min_pos
    = xmalloc((size_t)num_stripes * 3 * sizeof(*sorted_stripe_min_pos));

  for(int i = 0; i < num_stripes; ++i) {
    Xt_int ofs = (Xt_int)(stripes[i].stride * (stripes[i].nstrides - 1)),
      mask = Xt_isign_mask(ofs);
    stripe_minmax[0][i] = (Xt_int)(stripes[i].start + (ofs & mask));
  }

  xt_quicksort_index(stripe_minmax[0], (int)num_stripes,
                     sorted_stripe_min_pos, 1);

  int *restrict sorted_pos_prefix_sum, *restrict orig_pos_prefix_sum;

  sorted_pos_prefix_sum
    = sorted_stripe_min_pos + num_stripes;
  orig_pos_prefix_sum
    = xmalloc((size_t)num_stripes * sizeof(*orig_pos_prefix_sum));

  orig_pos_prefix_sum[0] = 0;
  for (int i = 1; i < num_stripes; ++i)
    orig_pos_prefix_sum[i] = orig_pos_prefix_sum[i-1] + stripes[i-1].nstrides;

  for (int i = 0; i < num_stripes; ++i) {
    sorted_pos_prefix_sum[i] = orig_pos_prefix_sum[sorted_stripe_min_pos[i]];
    int sorted_pos = sorted_stripe_min_pos[i];
    Xt_int ofs = (Xt_int)(stripes[sorted_pos].stride
                          * (stripes[sorted_pos].nstrides - 1)),
      mask = Xt_isign_mask(ofs);
    stripe_minmax[1][i] = (Xt_int)(stripes[sorted_pos].start + (ofs & ~mask));
  }

  free(orig_pos_prefix_sum);

  /* i'th stripe overlaps with overlap_count[i] following stripes, or
   * is part of a stretch of this many overlapping stripes, if
   * overlap_count[i] is > 0, in case overlap_count[i] <= 0, this many
   * non-overlapping stripes follow after negation + 1 */
  int *restrict overlap_count
    = sorted_stripe_min_pos + 2 * num_stripes;
  for (int i = 0; i < num_stripes - 1; ++i) {
    bool do_overlap = stripe_minmax[1][i] >= stripe_minmax[0][i + 1];
    if (do_overlap) {
      int j = i + 1;
      /* range_max_idx is the maximal index encountered in a rage of
       * overlapping stripes, only stop when a stripe starting at
       * index larger than this is encountered */
      Xt_int range_max_idx = MAX(stripe_minmax[1][i], stripe_minmax[1][i+1]);
      while (j + 1 < num_stripes
             && stripe_minmax[0][j + 1] <= range_max_idx) {
        range_max_idx = MAX(range_max_idx, stripe_minmax[1][j+1]);
        ++j;
      }
      overlap_count[i] = j - i;
      i = j;
    } else {
      int j = i + 1;
      while (j + 1 < num_stripes
             && stripe_minmax[0][j + 1] > stripe_minmax[1][j])
        ++j;
      overlap_count[i] = -(j - i - 1);
      i = j - 1;
    }
  }
  overlap_count[num_stripes - 1] = 0;

  Xt_int offset = 0;

  size_t i = 0;
  while (i < (size_t)num_stripes) {

    bool do_overlap = overlap_count[i] > 0;
    size_t num_selection = (size_t)(abs(overlap_count[i])) + 1;
    Xt_int curr_offset = 0;

    for (size_t j = 0; j < num_selection; ++j)
      curr_offset
        = (Xt_int)(curr_offset
                   + decode_stripe(stripes[sorted_stripe_min_pos[i+j]],
                                   sorted_vector_assign + offset
                                   + curr_offset,
                                   idxvec->sorted_vec_positions + offset
                                   + curr_offset,
                                   sorted_pos_prefix_sum[i+j]));

    if (do_overlap)
      xt_quicksort_index(sorted_vector_assign + offset, (int)curr_offset,
                         idxvec->sorted_vec_positions + offset, 0);

    offset = (Xt_int)(offset + curr_offset);
    i += num_selection;
  }

  free(sorted_stripe_min_pos);
  free(stripe_minmax);
}

Xt_idxlist
xt_idxvec_from_stripes_new(const struct Xt_stripe stripes[],
                           int num_stripes) {

  long long num_indices = 0;

  for (int i = 0; i < num_stripes; ++i)
    num_indices += stripes[i].nstrides;
  assert((sizeof (long long) > sizeof (int)) & (num_indices <= INT_MAX)
         & (num_indices >= 0));

  size_t vector_size = (size_t)num_indices * sizeof (Xt_int),
    header_size = ((sizeof (struct Xt_idxvec_) + sizeof (Xt_int) - 1)
                   /sizeof (Xt_int)) * sizeof (Xt_int);
  Xt_idxvec idxvec_obj = xmalloc(header_size + vector_size);

  Xt_int *restrict indices
    = (Xt_int *)(void *)((unsigned char *)idxvec_obj + header_size);
  idxvec_obj->vector = indices;

  size_t k = (size_t)-1;
  for (int i = 0; i < num_stripes; ++i)
    for (int j = 0; j < stripes[i].nstrides; ++j)
      indices[++k] = (Xt_int)(stripes[i].start + j * stripes[i].stride);

  Xt_idxlist_init(&idxvec_obj->parent, &idxvec_vtable, (int)num_indices);

  generate_sorted_vector_from_stripes(stripes, num_stripes, idxvec_obj);

  return (Xt_idxlist)idxvec_obj;
}

static void idxvec_delete(Xt_idxlist obj) {

  if (((Xt_idxvec)obj)->sorted_vector !=
      ((Xt_idxvec)obj)->vector)
    free((void *)(((Xt_idxvec)obj)->sorted_vector));
   free(((Xt_idxvec)obj)->sorted_vec_positions);
   free(obj);
}

static size_t idxvec_get_pack_size(Xt_idxlist obj, MPI_Comm comm) {

  Xt_idxvec idxvec = (Xt_idxvec)obj;
  int size_xt_idx, size_int_type;

  xt_mpi_call(MPI_Pack_size(2, MPI_INT, comm, &size_int_type), comm);
  xt_mpi_call(MPI_Pack_size(idxvec->parent.num_indices, Xt_int_dt, comm,
                            &size_xt_idx), comm);

  return (size_t)size_xt_idx + (size_t)size_int_type;
}

void idxvec_pack(Xt_idxlist obj, void *buffer, int buffer_size,
                 int *position, MPI_Comm comm) {

  assert(obj);
  Xt_idxvec idxvec = (Xt_idxvec)obj;
  int type = VECTOR;

  xt_mpi_call(MPI_Pack(&(type), 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  xt_mpi_call(MPI_Pack(&(idxvec->parent.num_indices), 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  if (idxvec->parent.num_indices != 0)
    xt_mpi_call(MPI_Pack((Xt_int *)idxvec->vector, idxvec->parent.num_indices,
                         Xt_int_dt, buffer,
                         buffer_size, position, comm), comm);
}

Xt_idxlist xt_idxvec_unpack(void *buffer, int buffer_size, int *position,
                            MPI_Comm comm) {

  int num_indices;

  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         &num_indices, 1, MPI_INT, comm), comm);

  size_t vector_size = (size_t)num_indices * sizeof (Xt_int),
    header_size = ((sizeof (struct Xt_idxvec_) + sizeof (Xt_int) - 1)
                   /sizeof (Xt_int)) * sizeof (Xt_int);
  Xt_idxvec idxvec = xmalloc(header_size + vector_size);
  Xt_idxlist_init(&idxvec->parent, &idxvec_vtable, num_indices);

  Xt_int *vector_assign = (Xt_int *)(void *)((unsigned char *)idxvec + header_size);
  idxvec->vector = vector_assign;
  if (num_indices != 0) {
    xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                           vector_assign, num_indices,
                           Xt_int_dt, comm), comm);
  } else {
    fputs("warning: implementation generated empty vector!\n", stderr);
    idxvec->vector = NULL;
  }

  idxvec->sorted_vector = NULL;
  idxvec->sorted_vec_positions = NULL;

  return (Xt_idxlist)idxvec;
}

static const Xt_int *
get_sorted_vector(Xt_idxvec idxvec) {

  if (idxvec->sorted_vector != NULL)
    return idxvec->sorted_vector;

  unsigned num_indices = (unsigned)idxvec->parent.num_indices;
  idxvec->sorted_vec_positions = xmalloc((size_t)num_indices *
                                         sizeof(*(idxvec->sorted_vec_positions)));


  int sorted = 1;
  // check if we are already sorted:
  for (unsigned i = 1; i < num_indices; ++i)
    sorted &= (idxvec->vector[i-1] <= idxvec->vector[i]);

  // we are done if we are already sorted
  if (sorted) {
    // gen id-map:
    for (unsigned i = 0; i < num_indices; ++i)
      idxvec->sorted_vec_positions[i] = (int)i;
    return idxvec->sorted_vector = idxvec->vector;
  }
  Xt_int *sorted_vector
    = xmalloc((size_t)num_indices * sizeof(*sorted_vector));

  memcpy(sorted_vector, idxvec->vector,
         (size_t)num_indices * sizeof(*sorted_vector));

  xt_quicksort_index(sorted_vector, (int)num_indices,
                     idxvec->sorted_vec_positions, 1);

  idxvec->sorted_vector = sorted_vector;

  return sorted_vector;
}

Xt_idxlist
xt_idxvec_get_intersection(Xt_idxlist idxlist_src, Xt_idxlist idxlist_dst) {

  // both lists are index vectors:

  Xt_idxvec idxvec_src = (Xt_idxvec)idxlist_src,
    idxvec_dst = (Xt_idxvec)idxlist_dst;


  unsigned num_indices_inter = 0,
    num_indices_src = (unsigned)idxvec_src->parent.num_indices,
    num_indices_dst = (unsigned)idxvec_dst->parent.num_indices;

  size_t vector_size = num_indices_dst * sizeof (idxvec_dst->vector[0]),
    header_size = ((sizeof (struct Xt_idxvec_) + sizeof (Xt_int) - 1)
                   /sizeof (Xt_int)) * sizeof (Xt_int);

  Xt_idxvec inter_vector = xmalloc(header_size + vector_size);

  Xt_int *vector_assign
    = (Xt_int *)(void *)((unsigned char *)inter_vector + header_size);
  inter_vector->vector = vector_assign;

  const Xt_int *restrict sorted_src_vector, *restrict sorted_dst_vector;

  // get sorted indices of source and destination

  sorted_src_vector = get_sorted_vector(idxvec_src);
  sorted_dst_vector = get_sorted_vector(idxvec_dst);

  // compute the intersection

  for (unsigned i = 0, j = 0; i < num_indices_dst; ++i) {

    while (j < num_indices_src &&
           sorted_src_vector[j] < sorted_dst_vector[i]) ++j;
    if (j >= num_indices_src) break;
    if (sorted_src_vector[j] == sorted_dst_vector[i])
      vector_assign[num_indices_inter++] = sorted_dst_vector[i];
  }

  if (num_indices_inter) {
    vector_size = (size_t)num_indices_inter * sizeof (idxvec_dst->vector[0]);
    inter_vector = xrealloc(inter_vector, header_size + vector_size);
    inter_vector->vector
      = (Xt_int *)(void *)((unsigned char *)inter_vector + header_size);
  } else {
    free(inter_vector);
    return xt_idxempty_new();
  }

  Xt_idxlist_init(&inter_vector->parent, &idxvec_vtable, (int)num_indices_inter);
  inter_vector->sorted_vector = NULL;
  inter_vector->sorted_vec_positions = NULL;

  return (Xt_idxlist)inter_vector;
}

static Xt_idxlist
idxvec_copy(Xt_idxlist idxlist) {

   Xt_idxvec idxvec_obj = (Xt_idxvec)idxlist;

   return xt_idxvec_new(idxvec_obj->vector, idxvec_obj->parent.num_indices);
}

static void
idxvec_get_indices(Xt_idxlist idxlist, Xt_int *indices) {

   Xt_idxvec idxvec_obj = (Xt_idxvec)idxlist;

   memcpy(indices, idxvec_obj->vector,
          (size_t)idxvec_obj->parent.num_indices * sizeof(*indices));
}

static Xt_int const*
idxvec_get_indices_const(Xt_idxlist idxlist) {
  Xt_idxvec idxvec = (Xt_idxvec)idxlist;

  return idxvec->vector;
}


static void
idxvec_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                         int * num_stripes) {

  Xt_idxvec idxvec_obj = (Xt_idxvec)idxlist;

  xt_convert_indices_to_stripes(idxvec_obj->vector,
                                idxvec_obj->parent.num_indices,
                                stripes, num_stripes);
}

static int
idxvec_get_index_at_position(Xt_idxlist idxlist, int position, Xt_int * index) {

  Xt_idxvec idxvec_obj = (Xt_idxvec)idxlist;

  if (position < 0 || position >= idxvec_obj->parent.num_indices)
    return 1;

  *index = idxvec_obj->vector[position];

  return 0;
}

static int
idxvec_get_indices_at_positions(Xt_idxlist idxlist, const int *restrict positions,
                                int num_pos, Xt_int *index,
                                Xt_int undef_idx) {

  Xt_idxvec idxvec = (Xt_idxvec)idxlist;
  int num_indices = idxvec->parent.num_indices;
  const Xt_int *restrict v = idxvec->vector;

  int undef_count = 0;
  for (int ip = 0; ip < num_pos; ip++) {
    int p = positions[ip];
    if (p >= 0 && p < num_indices) {
      index[ip] = v[p];
    } else {
      index[ip] = undef_idx;
      undef_count++;
    }
  }

  return undef_count;
}

/**
 * \todo check datatype of variables lb, ub and middle
 */
static int
idxvec_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                 int * position, int offset) {

  Xt_idxvec idxvec_obj = (Xt_idxvec)idxlist;

  *position = -1;

  int num_indices = idxvec_obj->parent.num_indices;
  if ((offset < 0) || (offset >= num_indices))
    return 1;

  const Xt_int *sorted_vector = get_sorted_vector(idxvec_obj);

  if ((index < sorted_vector[0]) ||
      (index > sorted_vector[num_indices-1]))
    return 1;

  // bisection to find one matching position:
  int lb = 0;
  int ub = num_indices - 1;

  while (sorted_vector[lb] < index) {

    int middle = (ub + lb + 1)/2;

    if (sorted_vector[middle] <= index)
      lb = middle;
    else
      if (ub == middle)
        return 1;
      else
        ub = middle;
  }

  // find left most match:
  while (lb > 0 && idxvec_obj->sorted_vector[lb-1] == index) --lb;

  // go forward until offset condition is satisfied:
  while (lb < num_indices - 1 &&  // boundary condition
         idxvec_obj->sorted_vec_positions[lb] < offset && // ignore postions left of offset
         idxvec_obj->sorted_vector[lb] == index) ++lb;  // check if index is valid

  // check if position is invalid:
  if (lb >= num_indices ||
      sorted_vector[lb] != index)
    return 1; // failure

  // result:
  *position = idxvec_obj->sorted_vec_positions[lb];
  return 0;
}

static int
idxvec_get_position_of_index(Xt_idxlist idxlist, Xt_int index, int * position) {

  return idxvec_get_position_of_index_off(idxlist, index, position, 0);
}

static int idx_vec_is_sorted(Xt_int const *idx, int n) {

  if (n>=2)
    for (int i = 1; i < n; i++)
      if (idx[i] < idx[i-1]) return 0;

  return 1;
}

static int
idxvec_get_positions_of_indices(Xt_idxlist body_idxlist,
                                const Xt_int *selection_idx,
                                int num_selection, int *positions,
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
    size_t idx_memsize = (size_t)num_selection * sizeof(*sorted_selection),
      pos_memsize = (size_t)num_selection * sizeof(*sorted_selection_pos),
      /* round pos_memsize up to next multiple of sizeof (int) */
      pos_ofs = ((idx_memsize + sizeof (int) - 1)
                 & ((size_t)-(ssize_t)(sizeof(int)))),
      /* compute size of merged allocation */
      alloc_size = pos_ofs + pos_memsize;

    tmp_idx = xmalloc(alloc_size);
    memcpy(tmp_idx, selection_idx, idx_memsize);

    sorted_selection_pos
      = (void *)((unsigned char *)tmp_idx + pos_ofs);
    xt_quicksort_index(tmp_idx, num_selection, sorted_selection_pos, 1);
    sorted_selection = tmp_idx;
  }

  /* motivation for usage of single_match_only:
   *  on the target side we want single_match_only,
   *  on the source side we don't
   */
  Xt_idxvec body_idxvec = (Xt_idxvec)body_idxlist;
  const Xt_int *sorted_body = get_sorted_vector(body_idxvec);
  int *sorted_body_pos = body_idxvec->sorted_vec_positions;
  int search_end = body_idxvec->parent.num_indices - 1;
  int num_unmatched = 0;

  // after the match we will move on one step in order to avoid matching the same position again
  int post_match_step = single_match_only != 0;

  int i=0;
  for (int search_start = 0, ub_guess_ofs = 1;
       i<num_selection && search_start<=search_end;
       ++i) {
    int selection_pos = selection_is_ordered ? i : sorted_selection_pos[i];
    Xt_int isel = sorted_selection[i];
    // bisection to find one matching position:
    int ub = search_start + ub_guess_ofs;
    ub = ub <= search_end ? ub : search_end;
    int lb = search_start;
     /* guess too low? */
    if (sorted_body[ub] < isel) {
      lb = ub + 1;
      lb = lb <= search_end ? lb : search_end;
      ub = search_end;
    }
    while (ub-lb>1) {
      int middle = (ub + lb + 1) / 2;
      /* todo: make branch free with mask/inv mask by predicate */
      if (sorted_body[middle] <= isel)
        lb = middle;
      else
        ub = middle;
    }

    int match_pos;
    // search is now narrowed to two positions, select one of them:
    if (isel == sorted_body[lb]) {
      match_pos = lb;
    } else if (isel == sorted_body[ub]){
      match_pos = ub;
    } else {
      num_unmatched++;
      positions[selection_pos] = -1;
      continue;
    }

    // find left most match >= search_start (bisection can lead to any match >= search_start)
    while (match_pos > search_start && sorted_body[match_pos-1] == isel)
      --match_pos;

    // result:
    // update positions and prepare next search:
    positions[selection_pos] = sorted_body_pos[match_pos];
    ub_guess_ofs = match_pos - search_start;
    search_start = match_pos + post_match_step;
  }
  if (i < num_selection) {
    if (selection_is_ordered)
      do {
        num_unmatched++;
        positions[i] = -1;
      } while (++i < num_selection);
    else
      do {
        num_unmatched++;
        positions[sorted_selection_pos[i]] = -1;
      } while (++i < num_selection);
  }
  if (tmp_idx) free(tmp_idx);

  return num_unmatched;
}

static Xt_int
idxvec_get_min_index(Xt_idxlist idxlist) {

  Xt_idxvec idxvec_obj = (Xt_idxvec)idxlist;

  if (!idxvec_obj->parent.num_indices)
    die("idxvec_get_min_index: empty index vector");

  return get_sorted_vector(idxvec_obj)[0];
}

static Xt_int
idxvec_get_max_index(Xt_idxlist idxlist) {

  Xt_idxvec idxvec_obj = (Xt_idxvec)idxlist;

  if (!idxvec_obj->parent.num_indices)
    die("idxvec_get_max_index: empty index vector");

  return get_sorted_vector(idxvec_obj)[idxvec_obj->parent.num_indices-1];
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
