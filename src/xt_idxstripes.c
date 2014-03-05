/**
 * @file xt_idxstripes.c
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

#include "xt/xt_idxlist.h"
#include "xt_idxlist_internal.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxstripes.h"
#include "xt_idxstripes_internal.h"
#include "xt/xt_mpi.h"
#include "xt_idxlist_unpack.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "ensure_array_size.h"
#include "instr.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

static void
idxstripes_delete(Xt_idxlist data);

static size_t
idxstripes_get_pack_size(Xt_idxlist data, MPI_Comm comm);

static void
idxstripes_pack(Xt_idxlist data, void *buffer, int buffer_size,
                int *position, MPI_Comm comm);

static Xt_idxlist
idxstripes_copy(Xt_idxlist idxlist);

static int
idxstripes_get_num_indices(Xt_idxlist idxlist);

static void
idxstripes_get_indices(Xt_idxlist idxlist, Xt_int *indices);

static const Xt_int *
idxstripes_get_indices_const(Xt_idxlist idxlist);

static void
idxstripes_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                             int * num_stripes);

static int
idxstripes_get_index_at_position(Xt_idxlist idxlist, int position,
                                 Xt_int * index);

static int
idxstripes_get_indices_at_positions(Xt_idxlist idxlist, int *positions,
                                    int num, Xt_int *index,
                                    Xt_int undef_idx);

static int
idxstripes_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                 int * position);

static int
idxstripes_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                     int * position, int offset);

static Xt_int
idxstripes_get_min_index(Xt_idxlist idxlist);

static Xt_int
idxstripes_get_max_index(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxstripes_vtable = {
  .delete                       = idxstripes_delete,
  .get_pack_size                = idxstripes_get_pack_size,
  .pack                         = idxstripes_pack,
  .copy                         = idxstripes_copy,
  .get_num_indices              = idxstripes_get_num_indices,
  .get_indices                  = idxstripes_get_indices,
  .get_indices_const            = idxstripes_get_indices_const,
  .get_index_stripes            = idxstripes_get_index_stripes,
  .get_index_at_position        = idxstripes_get_index_at_position,
  .get_indices_at_positions     = idxstripes_get_indices_at_positions,
  .get_position_of_index        = idxstripes_get_position_of_index,
  .get_positions_of_indices     = NULL,
  .get_position_of_index_off    = idxstripes_get_position_of_index_off,
  .get_positions_of_indices_off = NULL,
  .get_min_index                = idxstripes_get_min_index,
  .get_max_index                = idxstripes_get_max_index,
  .get_bounding_box             = NULL,
  .idxlist_pack_code            = STRIPES,
};

static MPI_Datatype stripe_dt;

void
xt_idxstripes_initialize()
{
  struct Xt_stripe stripe;

  MPI_Aint base_address, start_address, nstrides_address, stride_address;

  MPI_Get_address(&stripe, &base_address);
  MPI_Get_address(&stripe.start, &start_address);
  MPI_Get_address(&stripe.stride, &stride_address);
  MPI_Get_address(&stripe.nstrides, &nstrides_address);

  int block_lengths[5] = {1,1,1,1,1};
  MPI_Aint displacements[5] = {0, start_address - base_address,
                               stride_address - base_address,
                               nstrides_address - base_address,
                               sizeof(stripe)};
  MPI_Datatype types[5] = {MPI_LB, Xt_int_dt, Xt_int_dt, MPI_INT, MPI_UB};

  xt_mpi_call(MPI_Type_create_struct(5, block_lengths, displacements,
                                     types, &stripe_dt), Xt_default_comm);
  xt_mpi_call(MPI_Type_commit(&stripe_dt), Xt_default_comm);
}

void
xt_idxstripes_finalize()
{
  xt_mpi_call(MPI_Type_free(&stripe_dt), Xt_default_comm);
}

typedef struct Xt_idxstripes_ *Xt_idxstripes;

struct Xt_idxstripes_ {

  struct Xt_idxlist_ parent;

  struct Xt_stripe * stripes;
  int num_stripes;

  struct Xt_stripe * stripes_sorted;
  int stripes_overlap;

  int index_num_cache;
  Xt_int *index_array_cache;
};

static int compare_xtstripes(const void * a, const void * b)
{
  return (((struct Xt_stripe*)a)->start > ((struct Xt_stripe*)b)->start) -
         (((struct Xt_stripe*)a)->start < ((struct Xt_stripe*)b)->start);
}

static void generate_stripes_sorted(Xt_idxstripes stripes) {
  INSTR_DEF(instr,"generate_stripes_sorted");
  INSTR_START(instr);

  struct Xt_stripe *sp = stripes->stripes_sorted
    = xmalloc((size_t)stripes->num_stripes * sizeof (*sp));

  memcpy(sp, stripes->stripes, (size_t)stripes->num_stripes * sizeof (*sp));

  for (int i = 0; i < stripes->num_stripes; ++i)
    if (sp[i].stride < 0) {
      sp[i].start = (Xt_int)(sp[i].start + (sp[i].nstrides-1) * sp[i].stride);
      sp[i].stride = (Xt_int)(-sp[i].stride);
    }

  qsort(sp, (size_t)stripes->num_stripes, sizeof (*sp), compare_xtstripes);
  INSTR_STOP(instr);
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

Xt_idxlist
xt_idxstripes_new(struct Xt_stripe const * stripes, int num_stripes) {
  INSTR_DEF(instr,"xt_idxstripes_new");
  INSTR_START(instr);

  Xt_idxlist result;

  if (num_stripes > 0) {
    size_t header_size = ((sizeof (struct Xt_idxstripes_)
                           + sizeof (struct Xt_stripe) - 1)
                          / sizeof (struct Xt_stripe))
      * sizeof (struct Xt_stripe),
      body_size = sizeof (struct Xt_stripe) * (size_t)num_stripes;
    Xt_idxstripes idxstripes = xmalloc(header_size + body_size);
    Xt_idxlist_init(&idxstripes->parent, &idxstripes_vtable);
    idxstripes->num_stripes = num_stripes;
    idxstripes->index_num_cache = 0;
    idxstripes->index_array_cache = NULL;

    idxstripes->stripes
      = (struct Xt_stripe *)((unsigned char *)idxstripes + header_size);
    memcpy(idxstripes->stripes, stripes,
           (size_t)num_stripes * sizeof(*(idxstripes->stripes)));
    int sign_err = 0;
    for(int i = 0; i < num_stripes; i++) {
      if (idxstripes->stripes[i].nstrides<0) {
        idxstripes->stripes[i].nstrides = 0;
        sign_err++;
      }
      idxstripes->index_num_cache
        = idxstripes->index_num_cache + (int)stripes[i].nstrides;
    }
    if (sign_err)
      fputs("WARNING: xt_idxstripes_new called with invalid stripes\n", stderr);

    generate_stripes_sorted(idxstripes);

    idxstripes->stripes_overlap = 0;
    for (int i = 0; i < num_stripes - 1; ++i)
      if (stripes_overlap(idxstripes->stripes_sorted[i],
                          idxstripes->stripes_sorted[i+1])) {
        idxstripes->stripes_overlap = 1;
        break;
      }
    result = (Xt_idxlist)idxstripes;
  } else
    result = xt_idxempty_new();
  INSTR_STOP(instr);
  return result;
}

static void
idxstripes_delete(Xt_idxlist data) {

  if (data == NULL) return;

  Xt_idxstripes stripes = (Xt_idxstripes)data;

  free(stripes->stripes_sorted);
  free(stripes->index_array_cache);
  free(stripes);
}

static size_t
idxstripes_get_pack_size(Xt_idxlist data, MPI_Comm comm) {

  Xt_idxstripes stripes = (Xt_idxstripes)data;

  int size_int_type, size_num_stripes, size_stripes =  0;

  xt_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);
  xt_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &size_num_stripes), comm);
  if (stripes->num_stripes)
    xt_mpi_call(MPI_Pack_size(stripes->num_stripes, stripe_dt, comm,
                              &size_stripes), comm);

  return (size_t)size_int_type + (size_t)size_num_stripes
    + (size_t)size_stripes;
}

static void
idxstripes_pack(Xt_idxlist data, void *buffer, int buffer_size,
                int *position, MPI_Comm comm) {
  INSTR_DEF(instr,"idxstripes_pack");
  INSTR_START(instr);

  assert(data);
  Xt_idxstripes stripes = (Xt_idxstripes)data;
  int type = STRIPES;

  xt_mpi_call(MPI_Pack(&(type), 1, MPI_INT, buffer,
                       buffer_size, position, comm), comm);
  xt_mpi_call(MPI_Pack(&stripes->num_stripes, 1, MPI_INT, buffer,
                       buffer_size, position, comm), comm);
  if (stripes->num_stripes)
    xt_mpi_call(MPI_Pack(stripes->stripes, stripes->num_stripes, stripe_dt, buffer,
                         buffer_size, position, comm), comm);
  INSTR_STOP(instr);
}

Xt_idxlist xt_idxstripes_unpack(void *buffer, int buffer_size, int *position,
                                MPI_Comm comm) {

  INSTR_DEF(instr,"xt_idxstripes_unpack");
  INSTR_START(instr);

  int num_stripes;
  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         &num_stripes, 1, MPI_INT, comm), comm);

  Xt_idxlist result;
  if (num_stripes) {
    size_t header_size = ((sizeof (struct Xt_idxstripes_)
                           + sizeof (struct Xt_stripe) - 1)
                          / sizeof (struct Xt_stripe))
      * sizeof (struct Xt_stripe),
      body_size = sizeof (struct Xt_stripe) * (size_t)num_stripes;
    Xt_idxstripes idxstripes = xmalloc(header_size + body_size);
    Xt_idxlist_init(&idxstripes->parent, &idxstripes_vtable);
    idxstripes->index_num_cache = -1; // negative num means: undef
    idxstripes->index_array_cache = NULL;
    idxstripes->num_stripes = num_stripes;

    idxstripes->stripes
      = (struct Xt_stripe *)((unsigned char *)idxstripes + header_size);
    xt_mpi_call(MPI_Unpack(buffer, buffer_size, position, idxstripes->stripes,
                           idxstripes->num_stripes, stripe_dt, comm),comm);

    generate_stripes_sorted(idxstripes);

    idxstripes->stripes_overlap = 0;
    for (int i = 0; i < num_stripes - 1; ++i)
      if (stripes_overlap(idxstripes->stripes_sorted[i],
                          idxstripes->stripes_sorted[i+1])) {
        idxstripes->stripes_overlap = 1;
        break;
      }
    result = (Xt_idxlist)idxstripes;
  } else
    result = xt_idxempty_new();

  INSTR_STOP(instr);
  return result;
}

void xt_idxstripes_convert_to_stripes(Xt_int const * indices,
                                      int num_indices,
                                      struct Xt_stripe ** stripes,
                                      int * num_stripes) {

  INSTR_DEF(instr,"xt_idxstripes_convert_to_stripes");
  INSTR_START(instr);

  struct Xt_stripe * temp_stripes = NULL;
  size_t temp_stripes_array_size = 0;
  size_t num_temp_stripes = 0;

  int i, j;

  i = 0;
  while(i < num_indices) {

    ++num_temp_stripes;

    ENSURE_ARRAY_SIZE(temp_stripes, temp_stripes_array_size, num_temp_stripes);

    j = 1;

    while ((i + j) < num_indices && indices[i] + (Xt_int)j == indices[i + j]) ++j;

    temp_stripes[num_temp_stripes-1].start  = indices[i];
    temp_stripes[num_temp_stripes-1].nstrides = (Xt_int)j;
    temp_stripes[num_temp_stripes-1].stride = 1;

    i = i + j;
  }

  *stripes = xrealloc(temp_stripes, num_temp_stripes * sizeof(*temp_stripes));
  *num_stripes = (int)num_temp_stripes;
  INSTR_STOP(instr);
}

static Xt_idxlist
compute_intersection_fallback(Xt_idxstripes idxstripes_src,
                              Xt_idxstripes idxstripes_dst) {
  INSTR_DEF(instr,"compute_intersection_fallback");
  INSTR_START(instr);

  Xt_idxlist idxvec_from_stripes_src;
  Xt_idxlist idxvec_from_stripes_dst;

  idxvec_from_stripes_src
    = xt_idxvec_from_stripes_new(idxstripes_src->stripes,
                                 idxstripes_src->num_stripes);
  idxvec_from_stripes_dst
    = xt_idxvec_from_stripes_new(idxstripes_dst->stripes,
                                 idxstripes_dst->num_stripes);

  Xt_idxlist intersection;

  intersection = xt_idxlist_get_intersection(idxvec_from_stripes_src,
                                             idxvec_from_stripes_dst);

  xt_idxlist_delete(idxvec_from_stripes_src);
  xt_idxlist_delete(idxvec_from_stripes_dst);
  INSTR_STOP(instr);
  return intersection;
}

// least common multiple
static Xt_int lcm(Xt_int a, Xt_int b) {

  Xt_int temp = (Xt_int)(a * b);

  while (a != b)
    if (a > b) a = (Xt_int)(a - b);
    else       b = (Xt_int)(b - a);

  return (Xt_int)(temp / a);
}

static struct Xt_stripe
get_stripe_intersection(struct Xt_stripe * stripe_a,
                        struct Xt_stripe * stripe_b) {

  INSTR_DEF(instr,"get_stripe_intersection");
  INSTR_START(instr);

  unsigned i;

  if (stripe_a->start > stripe_b->start) {
    struct Xt_stripe * temp = stripe_a;
    stripe_a = stripe_b;
    stripe_b = temp;
  }

  Xt_int stripe_lcm = lcm(stripe_a->stride, stripe_b->stride);

  Xt_int norm_start = (Xt_int)(stripe_b->start - stripe_a->start);

  for (i = 0;
       i < (unsigned)(stripe_lcm / stripe_a->stride);
       ++i, norm_start = (Xt_int)(norm_start + stripe_b->stride))
    if (!(norm_start % stripe_a->stride)) break;

  struct Xt_stripe intersection;

  if (i == (unsigned)(stripe_lcm / stripe_a->stride)) {

    intersection.start = 0;
    intersection.stride = 0;
    intersection.nstrides = 0;

  } else {

    Xt_int stripe_a_end, stripe_b_end;

    intersection.start = (Xt_int)(norm_start + stripe_a->start);
    intersection.stride = stripe_lcm;
    stripe_a_end = (Xt_int)(stripe_a->start
                            + stripe_a->stride * (stripe_a->nstrides-1));
    stripe_b_end = (Xt_int)(stripe_b->start
                            + stripe_b->stride * (stripe_b->nstrides-1));
    intersection.nstrides
      = (int)((MIN(stripe_a_end, stripe_b_end) - intersection.start + 1 +
              (stripe_lcm - 1)) / stripe_lcm);
  }
  INSTR_STOP(instr);
  return intersection;
}

// this routine only works for idxstripes that have stripes_overlap==0
static Xt_idxlist
idxstripes_compute_intersection(Xt_idxstripes idxstripes_src,
                                Xt_idxstripes idxstripes_dst) {

  INSTR_DEF(instr,"idxstripes_compute_intersection");
  INSTR_START(instr);

  struct Xt_stripe * inter_stripes = NULL;
  size_t num_inter_stripes = 0;
  size_t inter_stripes_array_size = 0;

  int i_src, i_dst;

  i_src = 0;
  i_dst = 0;
  while (i_src < idxstripes_src->num_stripes &&
         i_dst < idxstripes_dst->num_stripes) {

    while (i_src < idxstripes_src->num_stripes &&
           idxstripes_src->stripes_sorted[i_src].start +
           idxstripes_src->stripes_sorted[i_src].stride *
           (idxstripes_src->stripes_sorted[i_src].nstrides-1) <
           idxstripes_dst->stripes_sorted[i_dst].start) ++i_src;

    if ( i_src >= idxstripes_src->num_stripes ) break;

    while (i_dst < idxstripes_dst->num_stripes &&
           idxstripes_dst->stripes_sorted[i_dst].start +
           idxstripes_dst->stripes_sorted[i_dst].stride *
           (idxstripes_dst->stripes_sorted[i_dst].nstrides-1) <
           idxstripes_src->stripes_sorted[i_src].start) ++i_dst;

    if ( i_dst >= idxstripes_dst->num_stripes ) break;

    if (stripes_overlap(idxstripes_src->stripes_sorted[i_src],
                        idxstripes_dst->stripes_sorted[i_dst])) {

      ENSURE_ARRAY_SIZE(inter_stripes, inter_stripes_array_size,
                        num_inter_stripes+1);

      inter_stripes[num_inter_stripes] =
        get_stripe_intersection(idxstripes_src->stripes_sorted+i_src,
                                idxstripes_dst->stripes_sorted+i_dst);

      if (inter_stripes[num_inter_stripes].nstrides > 0) ++num_inter_stripes;
    }

    if (idxstripes_dst->stripes_sorted[i_dst].start +
        idxstripes_dst->stripes_sorted[i_dst].stride *
        (idxstripes_dst->stripes_sorted[i_dst].nstrides+1) <
        idxstripes_src->stripes_sorted[i_src].start +
        idxstripes_src->stripes_sorted[i_src].stride *
        (idxstripes_src->stripes_sorted[i_src].nstrides+1))
      i_dst++;
    else
      i_src++;
  }

  Xt_idxlist inter;

  inter = xt_idxstripes_new(inter_stripes, (int)num_inter_stripes);

  free(inter_stripes);

  INSTR_STOP(instr);
  return inter;
}

Xt_idxlist
xt_idxstripes_get_intersection(Xt_idxlist idxlist_src, Xt_idxlist idxlist_dst)
{

  if (idxlist_dst->vtable != &idxstripes_vtable) return NULL;

  // both lists are index stripes:

  Xt_idxstripes idxstripes_src, idxstripes_dst;

  idxstripes_src = (Xt_idxstripes)idxlist_src;
  idxstripes_dst = (Xt_idxstripes)idxlist_dst;

  if (idxstripes_src->stripes_overlap ||
      idxstripes_dst->stripes_overlap) {

    return compute_intersection_fallback(idxstripes_src, idxstripes_dst);
  } else
    return idxstripes_compute_intersection(idxstripes_src, idxstripes_dst);
}

static Xt_idxlist
idxstripes_copy(Xt_idxlist idxlist) {

  Xt_idxstripes stripes = (Xt_idxstripes)idxlist;

  return xt_idxstripes_new(stripes->stripes, stripes->num_stripes);
}

static int
idxstripes_get_num_indices(Xt_idxlist idxlist) {   ///< \todo use cached value

  Xt_idxstripes stripes;

  stripes = (Xt_idxstripes)idxlist;

  if (stripes->index_num_cache >= 0) return stripes->index_num_cache;

  int size = 0;

  for (int i = 0; i < stripes->num_stripes; ++i)
    size = (int)(size + stripes->stripes[i].nstrides);

  stripes->index_num_cache = size;
  return size;
}

static void
idxstripes_get_indices(Xt_idxlist idxlist, Xt_int *indices) {
  INSTR_DEF(instr,"idxstripes_get_indices");
  INSTR_START(instr);

  /// \todo use memcpy with index_array_cache if available
  Xt_idxstripes stripes = (Xt_idxstripes)idxlist;

  --indices;
  for (int i = 0; i < stripes->num_stripes; ++i)
    for (Xt_int j = 0; j < stripes->stripes[i].nstrides; ++j)
      *(++indices)
        = (Xt_int)(stripes->stripes[i].start + j * stripes->stripes[i].stride);

  INSTR_STOP(instr);
}

static Xt_int const*
idxstripes_get_indices_const(Xt_idxlist idxlist) {

  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;

  if (idxstripes->index_array_cache) return idxstripes->index_array_cache;

  int index_num = idxstripes_get_num_indices(idxlist);

  Xt_int *tmp_index_array
    = xmalloc((size_t)index_num * sizeof( *(idxstripes->index_array_cache) ) );

  idxstripes_get_indices(idxlist, tmp_index_array);

  idxstripes->index_array_cache = tmp_index_array;

  return idxstripes->index_array_cache;
}

static void
idxstripes_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                             int * num_stripes) {

  INSTR_DEF(instr,"idxstripes_get_index_stripes");
  INSTR_START(instr);

  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;

  struct Xt_stripe * temp_stripes = NULL;
  size_t temp_stripes_array_size = 0;
  size_t num_temp_stripes = 0;

  Xt_int j;

  for (int i = 0; i < idxstripes->num_stripes; ++i) {

    if (idxstripes->stripes[i].stride == 1) {

      ++num_temp_stripes;

      ENSURE_ARRAY_SIZE(temp_stripes, temp_stripes_array_size,
                        num_temp_stripes);

      temp_stripes[num_temp_stripes-1] = idxstripes->stripes[i];

    } else {

      ENSURE_ARRAY_SIZE(temp_stripes, temp_stripes_array_size,
                        num_temp_stripes
                        + (size_t)idxstripes->stripes[i].nstrides);

      for (j = 0; j < idxstripes->stripes[i].nstrides; ++j)  {

        temp_stripes[num_temp_stripes].start
          = (Xt_int)(idxstripes->stripes[i].start
                     + j * idxstripes->stripes[i].stride);
        temp_stripes[num_temp_stripes].nstrides = 1;
        temp_stripes[num_temp_stripes].stride = 1;

        ++num_temp_stripes;
      }
    }
  }

  *stripes = xrealloc(temp_stripes, num_temp_stripes * sizeof(*temp_stripes));
  *num_stripes = (int)num_temp_stripes;

  INSTR_STOP(instr);
}

static int
idxstripes_get_index_at_position(Xt_idxlist idxlist, int position,
                                 Xt_int * index) {

  INSTR_DEF(instr,"idxstripes_get_index_at_position");
  INSTR_START(instr);

  int retval = 1;

  Xt_idxstripes stripes = (Xt_idxstripes)idxlist;

  if (position < 0) goto fun_exit;

  for (int i = 0; i < stripes->num_stripes; ++i)
    if (position >= stripes->stripes[i].nstrides)
      position-= (int)stripes->stripes[i].nstrides;
    else {
      *index = (Xt_int)(stripes->stripes[i].start
                        + position * stripes->stripes[i].stride);
      retval = 0;
      break;
    }

 fun_exit: ;
  INSTR_STOP(instr);
  return retval;
}


static int
idxstripes_get_indices_at_positions(Xt_idxlist idxlist, int *positions,
                                    int num_pos, Xt_int *index,
                                    Xt_int undef_idx) {

  INSTR_DEF(instr,"idxstripes_get_indices_at_positions");
  INSTR_START(instr);

  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
  struct Xt_stripe *stripes = idxstripes->stripes;

  int max_pos = idxstripes_get_num_indices(idxlist)-1;
  int seek_pos;
  int sub_pos = 0;
  int stripe_start_pos = 0;
  int istripe = 0;
  int undef_count = 0;

  for (int ipos = 0; ipos < num_pos; ipos++) {

    seek_pos = positions[ipos];

    if (seek_pos < 0 || seek_pos > max_pos) {
      index[ipos] = undef_idx;
      undef_count++;
      continue;
    }

    while (seek_pos < stripe_start_pos) {
      istripe--;
      if (istripe < 0)
        die("idxstripes_get_indices_at_positions: internal error:"
            " crossed 0-boundary");
      stripe_start_pos -= (int)stripes[istripe].nstrides;
    }

    while (seek_pos > stripe_start_pos + stripes[istripe].nstrides - 1) {
      stripe_start_pos += (int)stripes[istripe].nstrides;
      istripe++;
      if (istripe >= idxstripes->num_stripes)
        die("idxstripes_get_indices_at_positions: internal error:"
            " crossed boundary");
    }

    sub_pos = seek_pos - stripe_start_pos;
    index[ipos]
      = (Xt_int)(stripes[istripe].start + sub_pos * stripes[istripe].stride);
  }

  INSTR_STOP(instr);

  return undef_count;
}

static int
idxstripes_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                 int * position) {

  return idxstripes_get_position_of_index_off(idxlist, index, position, 0);
}

static int
idxstripes_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                     int * position, int offset) {

  INSTR_DEF(instr,"idxstripes_get_position_of_index_off");
  INSTR_START(instr);

  int retval = 1;

  Xt_idxstripes stripes = (Xt_idxstripes)idxlist;

  int i = 0;
  Xt_int position_offset = 0;

  while(i < stripes->num_stripes &&
        position_offset + stripes->stripes[i].nstrides <= offset)
    position_offset
      = (Xt_int)(position_offset + stripes->stripes[i++].nstrides);

  for (; i < stripes->num_stripes;
       position_offset
         = (Xt_int)(position_offset + stripes->stripes[i++].nstrides)) {

    if ((stripes->stripes[i].stride > 0 && index < stripes->stripes[i].start)
        || (stripes->stripes[i].stride < 0
            && index > stripes->stripes[i].start))
      continue;

    Xt_int rel_start
      = (Xt_int)(index - stripes->stripes[i].start);

    if (rel_start%stripes->stripes[i].stride) continue;

    if (rel_start/stripes->stripes[i].stride >= stripes->stripes[i].nstrides)
      continue;

    *position = (int)(rel_start/stripes->stripes[i].stride + position_offset);

    retval = 0;
    goto fun_exit;
  }

  *position = -1;

 fun_exit: ;
  INSTR_STOP(instr);
  return retval;
}

static Xt_int
idxstripes_get_min_index(Xt_idxlist idxlist) {

  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;

  if (idxstripes->num_stripes > 0)
    return idxstripes->stripes_sorted[0].start;
  else
    return 0;
}

static Xt_int
idxstripes_get_max_index(Xt_idxlist idxlist) {

  INSTR_DEF(instr,"idxstripes_get_max_index");
  INSTR_START(instr);

  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;

  Xt_int index = idxstripes_get_min_index(idxlist);

  for (int i = 0; i < idxstripes->num_stripes; ++i)
    index = (Xt_int)(MAX(index, idxstripes->stripes_sorted[i].start +
                         idxstripes->stripes_sorted[i].stride *
                         (idxstripes->stripes_sorted[i].nstrides-1)));
  INSTR_STOP(instr);
  return index;
}
