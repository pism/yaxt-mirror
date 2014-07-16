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
#include "xt_arithmetic_util.h"
#include "xt/xt_idxlist.h"
#include "xt_idxlist_internal.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxstripes.h"
#include "xt_idxstripes_internal.h"
#include "xt/xt_mpi.h"
#include "xt_idxlist_unpack.h"
#include "xt_cover.h"
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
idxstripes_get_indices_at_positions(Xt_idxlist idxlist, const int *positions,
                                    int num, Xt_int *index,
                                    Xt_int undef_idx);
static int
idxstripes_get_pos_exts_of_index_stripes(Xt_idxlist idxlist,
                                         int num_stripes,
                                         const struct Xt_stripe *stripes,
                                         int *num_ext,
                                         struct Xt_pos_ext **pos_ext,
                                         int single_match_only);

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
  .get_indices                  = idxstripes_get_indices,
  .get_indices_const            = idxstripes_get_indices_const,
  .get_index_stripes            = idxstripes_get_index_stripes,
  .get_index_at_position        = idxstripes_get_index_at_position,
  .get_indices_at_positions     = idxstripes_get_indices_at_positions,
  .get_position_of_index        = idxstripes_get_position_of_index,
  .get_positions_of_indices     = NULL,
  .get_pos_exts_of_index_stripes = idxstripes_get_pos_exts_of_index_stripes,
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
  bool stripes_do_overlap;

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

  size_t num_stripes = (size_t)stripes->num_stripes;
  struct Xt_stripe *restrict sp = stripes->stripes_sorted;
  memcpy(sp, stripes->stripes, num_stripes * sizeof (*sp));

  for (size_t i = 0; i < num_stripes; ++i)
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

  return (a_max >= b_min) & (a_min <= b_max);
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
    Xt_idxstripes idxstripes = xmalloc(header_size + 2 * body_size);
    idxstripes->num_stripes = num_stripes;
    idxstripes->index_array_cache = NULL;

    idxstripes->stripes
      = (struct Xt_stripe *)((unsigned char *)idxstripes + header_size);
    memcpy(idxstripes->stripes, stripes,
           (size_t)num_stripes * sizeof(*(idxstripes->stripes)));
    struct Xt_stripe *restrict stripes = idxstripes->stripes;
    int sign_err = (stripes[0].nstrides < 0), sorted = 1,
      have_negative_stride_sign = (stripes[0].stride < 0);
    if (sign_err)
      stripes[0].nstrides = 0;
    long long num_indices = (long long)stripes[0].nstrides;
    for(size_t i = 1; i < (size_t)num_stripes; i++) {
      if (stripes[i].nstrides<0) {
        stripes[i].nstrides = 0;
        sign_err++;
      }
      have_negative_stride_sign |= (stripes[i].stride < 0);
      sorted &= (stripes[i].start >= stripes[i - 1].start);
      num_indices += (long long)stripes[i].nstrides;
    }
    if (sign_err)
      fputs("WARNING: xt_idxstripes_new called with invalid stripes\n", stderr);

    assert(num_indices <= INT_MAX);
    Xt_idxlist_init(&idxstripes->parent, &idxstripes_vtable, (int)num_indices);

    if (sorted && !have_negative_stride_sign)
    {
      idxstripes = xrealloc(idxstripes, header_size + body_size);
      idxstripes->stripes_sorted = idxstripes->stripes
        = (struct Xt_stripe *)((unsigned char *)idxstripes + header_size);
    }
    else
    {
      idxstripes->stripes_sorted = idxstripes->stripes + num_stripes;
      generate_stripes_sorted(idxstripes);
    }

    idxstripes->stripes_do_overlap = false;
    for (int i = 0; i < num_stripes - 1; ++i)
      if (stripes_overlap(idxstripes->stripes_sorted[i],
                          idxstripes->stripes_sorted[i+1])) {
        idxstripes->stripes_do_overlap = true;
        break;
      }
    result = (Xt_idxlist)idxstripes;
  } else
    result = xt_idxempty_new();
  INSTR_STOP(instr);
  return result;
}

Xt_idxlist
xt_idxstripes_prealloc_new(struct Xt_stripe const * stripes, int num_stripes)
{
  Xt_idxlist result;

  if (num_stripes > 0) {
    size_t header_size = ((sizeof (struct Xt_idxstripes_)
                           + sizeof (struct Xt_stripe) - 1)
                          / sizeof (struct Xt_stripe))
      * sizeof (struct Xt_stripe),
      body_size = sizeof (struct Xt_stripe) * (size_t)num_stripes;
    Xt_idxstripes idxstripes = xmalloc(header_size + body_size);
    idxstripes->num_stripes = num_stripes;
    idxstripes->index_array_cache = NULL;

    idxstripes->stripes = (struct Xt_stripe *)stripes;
    int sign_err = (stripes[0].nstrides < 0), sorted = 1,
      have_negative_stride_sign = (stripes[0].stride < 0);
    long long num_indices = (long long)stripes[0].nstrides;
    for(int i = 1; i < num_stripes; i++) {
      sign_err |= (stripes[i].nstrides < 0);
      have_negative_stride_sign |= (stripes[i].stride < 0);
      num_indices += (long long)stripes[i].nstrides;
      sorted &= (stripes[i].start >=stripes[i - 1].start);
    }
    assert(num_indices <= INT_MAX);
    if (sign_err)
      die("ERROR: xt_idxstripes_prealloc_new called with invalid stripes");

    Xt_idxlist_init(&idxstripes->parent, &idxstripes_vtable, (int)num_indices);

    if (sorted && !have_negative_stride_sign)
    {
      idxstripes = xrealloc(idxstripes, sizeof (struct Xt_idxstripes_));
      idxstripes->stripes_sorted = (struct Xt_stripe *)stripes;
    }
    else
    {
      idxstripes->stripes_sorted
        = (struct Xt_stripe *)((unsigned char *)idxstripes + header_size);
      generate_stripes_sorted(idxstripes);
    }
    idxstripes->stripes_do_overlap = false;
    for (int i = 0; i < num_stripes - 1; ++i)
      if (stripes_overlap(idxstripes->stripes_sorted[i],
                          idxstripes->stripes_sorted[i+1])) {
        idxstripes->stripes_do_overlap = true;
        break;
      }
    result = (Xt_idxlist)idxstripes;
  } else
    result = xt_idxempty_new();
  return result;
}


static void
idxstripes_delete(Xt_idxlist data) {

  if (data == NULL) return;

  Xt_idxstripes stripes = (Xt_idxstripes)data;

  free(stripes->index_array_cache);
  free(stripes);
}

static size_t
idxstripes_get_pack_size(Xt_idxlist data, MPI_Comm comm) {

  Xt_idxstripes stripes = (Xt_idxstripes)data;

  int size_int_type, size_stripes =  0;

  xt_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);
  if (stripes->num_stripes)
    xt_mpi_call(MPI_Pack_size(stripes->num_stripes, stripe_dt, comm,
                              &size_stripes), comm);

  return (size_t)(2 + (size_stripes>0)) * (size_t)size_int_type + (size_t)size_stripes;
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
  {
    int sorted = (stripes->stripes_sorted == stripes->stripes);
    xt_mpi_call(MPI_Pack(&sorted, 1, MPI_INT,
                         buffer, buffer_size, position, comm), comm);
    xt_mpi_call(MPI_Pack(stripes->stripes, stripes->num_stripes, stripe_dt,
                         buffer, buffer_size, position, comm), comm);
  }
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
    int sorted;
    xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                           &sorted, 1, MPI_INT, comm), comm);
    size_t header_size = ((sizeof (struct Xt_idxstripes_)
                           + sizeof (struct Xt_stripe) - 1)
                          / sizeof (struct Xt_stripe))
      * sizeof (struct Xt_stripe),
      body_size = sizeof (struct Xt_stripe) * (size_t)num_stripes;
    Xt_idxstripes idxstripes = xmalloc(header_size
                                       + (size_t)(!sorted + 1) * body_size);
    idxstripes->index_array_cache = NULL;
    idxstripes->num_stripes = num_stripes;

    idxstripes->stripes
      = (struct Xt_stripe *)((unsigned char *)idxstripes + header_size);
    xt_mpi_call(MPI_Unpack(buffer, buffer_size, position, idxstripes->stripes,
                           idxstripes->num_stripes, stripe_dt, comm),comm);
    if (sorted)
    {
      idxstripes->stripes_sorted = idxstripes->stripes;
    }
    else
    {
      idxstripes->stripes_sorted = idxstripes->stripes + num_stripes;
      generate_stripes_sorted(idxstripes);
    }

    long long num_indices = idxstripes->stripes[0].nstrides;
    idxstripes->stripes_do_overlap = false;
    int i;
    for (i = 1; i < num_stripes; ++i) {
      if (stripes_overlap(idxstripes->stripes_sorted[i - 1],
                          idxstripes->stripes_sorted[i])) {
        idxstripes->stripes_do_overlap = true;
        break;
      }
      num_indices += idxstripes->stripes[i].nstrides;
    }
    for (; i < num_stripes; ++i)
      num_indices += idxstripes->stripes[i].nstrides;
    assert(num_indices <= INT_MAX);
    result = (Xt_idxlist)idxstripes;
    Xt_idxlist_init(&idxstripes->parent, &idxstripes_vtable, (int)num_indices);
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


struct extended_gcd {
  Xt_int gcd, u, v;
};

/* computes egcd of two positive integers a and b such that
   egcd.gcd == egcd.u * a + egcd.v * b */
static inline struct extended_gcd extended_gcd(Xt_int a, Xt_int b) {
  Xt_int t = 1, u = 1, v = 0, s = 0;
  while (b>0)
  {
    Xt_int q = (Xt_int)(a / b);
    Xt_int prev_a = a;
    a = b;
    b = (Xt_int)(prev_a - q * b);
    Xt_int prev_u = u;
    u = s;
    s = (Xt_int)(prev_u - q * s);
    Xt_int prev_v = v;
    v = t;
    t = (Xt_int)(prev_v - q * t);
  }
  return (struct extended_gcd){ .gcd = a, .u = u, .v = v };
}

/* This implementation uses the method outlined in
 * James M. Stichnoth,
 * Efficient Compilation of Array Statements for Private Memory Multicomputers
 * February, 1993  CMU-CS-93-109
 */
static struct Xt_stripe
get_stripe_intersection(struct Xt_stripe stripe_a,
                        struct Xt_stripe stripe_b) {

  INSTR_DEF(instr,"get_stripe_intersection");
  INSTR_START(instr);

  struct Xt_bounded_stripe {
    Xt_int min, max, stride, representative;
  };

  struct Xt_bounded_stripe bsa, bsb, bsi;

  bsa.stride = stripe_a.stride;
  if (stripe_a.stride >= 0)
    bsa.min = stripe_a.start,
      bsa.max = (Xt_int)(stripe_a.start
                         + stripe_a.stride * (stripe_a.nstrides - 1));
  else
    bsa.min = (Xt_int)(stripe_a.start
                       + stripe_a.stride * (stripe_a.nstrides - 1)),
      bsa.max = stripe_a.start;
  bsa.representative = stripe_a.start;

  bsb.stride = stripe_b.stride;
  if (stripe_b.stride >= 0)
    bsb.min = stripe_b.start,
      bsb.max = (Xt_int)(stripe_b.start
                         + stripe_b.stride * (stripe_b.nstrides - 1));
  else
    bsb.min = (Xt_int)(stripe_b.start
                       + stripe_b.stride * (stripe_b.nstrides - 1)),
      bsb.max = stripe_b.start;
  bsb.representative = stripe_b.start;

  struct extended_gcd eg = extended_gcd(XT_INT_ABS(bsa.stride),
                                        XT_INT_ABS(bsb.stride));
  bsi.min = MAX(bsa.min, bsb.min);
  bsi.max = MIN(bsa.max, bsb.max);
  bsi.stride = (Xt_int)((XT_INT_ABS(stripe_a.stride) * stripe_b.stride)/eg.gcd);
  int even_divide
    = ((bsb.representative - bsa.representative) % eg.gcd) == 0;
  bsi.representative
    = (Xt_int)(bsa.representative
               + ((bsb.representative - bsa.representative) * eg.u * bsa.stride
                  / eg.gcd));

  struct Xt_stripe intersection;
  /* compute minimal bsi representative >= bsi.min */
  Xt_int min_rep
    = (Xt_int)(bsi.representative
               - (bsi.representative - bsi.min) / bsi.stride * bsi.stride);
  /* requires two's complement integers */
  int strides_mask = ~(((even_divide) & (bsi.min <= bsi.max)
                        & (min_rep <= bsi.max)) - 1);
  int nstrides = (int)((bsi.max - min_rep)/bsi.stride + Xt_isign(bsi.stride));
  Xt_int max_rep
    = (Xt_int)(min_rep + (nstrides - Xt_isign(bsi.stride)) * bsi.stride);
  intersection.start = (Xt_int)((bsa.stride >= 0) ? min_rep : max_rep);
  intersection.stride
    = (Xt_int)(Xt_isign(bsa.stride) * XT_INT_ABS(bsi.stride));
  intersection.nstrides = abs(nstrides) & strides_mask;

  INSTR_STOP(instr);
  return intersection;
}

// this routine only works for idxstripes where !stripes_do_overlap
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
        get_stripe_intersection(idxstripes_src->stripes_sorted[i_src],
                                idxstripes_dst->stripes_sorted[i_dst]);

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
  // both lists are index stripes:
  Xt_idxstripes idxstripes_src = (Xt_idxstripes)idxlist_src,
    idxstripes_dst = (Xt_idxstripes)idxlist_dst;

  if (idxstripes_src->stripes_do_overlap ||
      idxstripes_dst->stripes_do_overlap) {

    return compute_intersection_fallback(idxstripes_src, idxstripes_dst);
  } else
    return idxstripes_compute_intersection(idxstripes_src, idxstripes_dst);
}

static Xt_idxlist
idxstripes_copy(Xt_idxlist idxlist) {

  Xt_idxstripes stripes = (Xt_idxstripes)idxlist;

  return xt_idxstripes_new(stripes->stripes, stripes->num_stripes);
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

  int num_indices = idxlist->num_indices;

  Xt_int *tmp_index_array
    = xmalloc((size_t)num_indices * sizeof( *(idxstripes->index_array_cache) ) );

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
idxstripes_get_indices_at_positions(Xt_idxlist idxlist, const int *positions,
                                    int num_pos, Xt_int *index,
                                    Xt_int undef_idx) {

  INSTR_DEF(instr,"idxstripes_get_indices_at_positions");
  INSTR_START(instr);

  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
  struct Xt_stripe *stripes = idxstripes->stripes;

  int max_pos = idxlist->num_indices - 1;
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

static inline void
append_ext(struct Xt_pos_ext pos_ext, size_t *num_pos_exts,
           size_t *size_pos_exts,
           struct Xt_pos_ext **pos_exts)
{
  size_t num_pos_exts_ = *num_pos_exts, size_pos_exts_ = *size_pos_exts;
  struct Xt_pos_ext *restrict pos_exts_ = *pos_exts;
  if (xt_pos_ext_is_appendable(pos_exts_[num_pos_exts_ - 1], pos_ext))
  {
    if (num_pos_exts_ + 1 == size_pos_exts_)
    {
      size_pos_exts_ += 16;
      /* offsetting by 1 necessary to keep the terminator in place */
      *pos_exts = pos_exts_ = (struct Xt_pos_ext *)
        xrealloc(pos_exts_ - 1, (size_pos_exts_ + 1) * sizeof (*pos_exts_)) + 1;
      *size_pos_exts = size_pos_exts_;
    }
    pos_exts_[num_pos_exts_] = pos_ext;
    *num_pos_exts = num_pos_exts_ + 1;
  }
  else {
    /* merge new ext with previous */
    pos_exts_[num_pos_exts_ - 1].size
      =  isign(pos_ext.start - pos_exts_[num_pos_exts_ - 1].start)
      * (abs(pos_exts_[num_pos_exts_ - 1].size) + abs(pos_ext.size));
  }
}

static inline intptr_t
roundAddrUpToMultiple(intptr_t addr, intptr_t factor)
{
  return (addr + factor - 1)/factor * factor;
}

static inline int
ext_overlap(struct Xt_pos_ext a, struct Xt_pos_ext b)
{
  return ((b.start >= a.start) & (b.start < a.start + a.size))
    | ((a.start >= b.start) & (a.start < b.start + b.size));
}

struct Xt_stripe_minmax {
  Xt_int min, max;
};

static inline struct Xt_stripe_minmax
stripe2minmax(struct Xt_stripe stripe)
{
  Xt_int min = stripe.start, max = min;
  Xt_int stride = stripe.stride;
  int nstrides = stripe.nstrides;
  Xt_int mask = Xt_isign_mask(stride);
  stride = (Xt_int)((Xt_int)(nstrides - 1) * stride);
  min = (Xt_int)(min + (Xt_int)(mask & stride));
  max = (Xt_int)(max + (Xt_int)((Xt_int)~mask & stride));
  return (struct Xt_stripe_minmax){ .min = min, .max = max};
}

struct Xt_stripes_lookup {
  size_t num_stripes;
  const struct Xt_stripe *stripes;
  const struct Xt_stripe_minmax *minmax;
  const int *stripes_nstrides_psum;
};

static inline void
create_stripes_lookup(struct Xt_stripes_lookup *restrict db,
                      size_t num_db_stripes,
                      const struct Xt_stripe *restrict stripes) {
  db->stripes = stripes;
  /* using num_db_stripes + 1 ensures re-aligning is always possible */
  struct Xt_stripe_minmax *restrict db_minmax
    = xmalloc((num_db_stripes + 1)
              * (sizeof (db->minmax[0])
                 + sizeof (db->stripes_nstrides_psum[0])));
  int *restrict db_stripes_nstrides_psum
    = (int *)roundAddrUpToMultiple(
      (intptr_t)(db_minmax + num_db_stripes),
      (intptr_t)sizeof (*db_stripes_nstrides_psum));
  db_stripes_nstrides_psum[0] = 0;
  for (size_t j = 0; j < num_db_stripes; ++j) {
    db_minmax[j] = stripe2minmax(stripes[j]);
    db_stripes_nstrides_psum[j + 1]
      = db_stripes_nstrides_psum[j] + stripes[j].nstrides;
  }
  db->minmax = db_minmax;
  db->num_stripes = num_db_stripes;
  db->stripes_nstrides_psum = db_stripes_nstrides_psum;
}

static inline void
destroy_stripes_lookup(struct Xt_stripes_lookup *restrict db) {
  free((void *)db->minmax);
}



static size_t
idxstripes_get_pos_exts_of_index_stripe(
  struct Xt_stripe query,
  const struct Xt_stripes_lookup *restrict db,
  size_t *num_ext_,
  size_t *size_ext,
  struct Xt_pos_ext **pos_exts,
  bool single_match_only,
  size_t db_start);

int
idxstripes_get_pos_exts_of_index_stripes(Xt_idxlist idxlist,
                                         int num_stripes,
                                         const struct Xt_stripe stripes[num_stripes],
                                         int *num_ext,
                                         struct Xt_pos_ext **pos_exts,
                                         int single_match_only)
{
  size_t unmatched = 0;
  struct Xt_pos_ext_vec result;
  result.num_pos_ext = 0;
  struct Xt_idxstripes_ *idxstripes = (struct Xt_idxstripes_ *)idxlist;
  if (num_stripes > 0)
  {
    struct Xt_stripes_lookup stripes_db;
    result.size_pos_ext = (size_t)MAX(MIN(idxstripes->num_stripes,
                                          num_stripes), 8);
    result.pos_ext = xmalloc((result.size_pos_ext + 1)
                             * sizeof (*result.pos_ext));
    /* put non-concatenable *terminator* at offset -1 */
    result.pos_ext[0] = (struct Xt_pos_ext){.start = INT_MIN, .size = -1 };
    ++result.pos_ext;
    create_stripes_lookup(&stripes_db,
                          (size_t)idxstripes->num_stripes, idxstripes->stripes);
    for (size_t i = 0; i < (size_t)num_stripes; ++i) {
      struct Xt_stripe query = stripes[i];
      unmatched += idxstripes_get_pos_exts_of_index_stripe(
        query, &stripes_db,
        &result.num_pos_ext, &result.size_pos_ext, &result.pos_ext,
        single_match_only != 0, 0);
    }
    --(result.pos_ext);
    memmove(result.pos_ext, result.pos_ext + 1,
            sizeof (*result.pos_ext) * result.num_pos_ext);
    *pos_exts = xrealloc(result.pos_ext,
                         sizeof (*result.pos_ext) * result.num_pos_ext);
    destroy_stripes_lookup(&stripes_db);
  }
  *num_ext = (int)result.num_pos_ext;
  return (int)unmatched;
}

static size_t
conditional_pos_ext_insert(struct Xt_stripe query,
                           struct Xt_pos_ext pos_ext2add,
                           const struct Xt_stripes_lookup *restrict db,
                           size_t *num_ext,
                           size_t *size_ext,
                           struct Xt_pos_ext **pos_exts,
                           size_t db_start);

struct unmatched_tail
{
  size_t unmatched;
  struct Xt_stripe query_tail;
};

static struct unmatched_tail
idxstripes_complex_get_pos_exts_of_index_stripe(
  struct Xt_stripe query,
  const struct Xt_stripes_lookup *restrict stripes_lookup,
  size_t *num_ext_,
  size_t *size_ext,
  struct Xt_pos_ext **pos_exts,
  bool single_match_only,
  size_t db_start);

static size_t
pos_ext_insert(struct Xt_stripe query,
               struct Xt_pos_ext pos_ext2add,
               const struct Xt_stripes_lookup *stripes_lookup,
               size_t *num_ext_,
               size_t *size_ext,
               struct Xt_pos_ext **pos_exts,
               bool single_match_only,
               size_t db_pos)
{
  size_t unmatched = 0;
  if (single_match_only)
    unmatched += conditional_pos_ext_insert(
      query, pos_ext2add, stripes_lookup,
      num_ext_, size_ext, pos_exts, db_pos);
  else
    append_ext(pos_ext2add, num_ext_, size_ext, pos_exts);
  return unmatched;
}


size_t
idxstripes_get_pos_exts_of_index_stripe(
  struct Xt_stripe query,
  const struct Xt_stripes_lookup *restrict db,
  size_t *num_ext_,
  size_t *size_ext,
  struct Xt_pos_ext **pos_exts,
  bool single_match_only,
  size_t db_start)
{
  size_t unmatched = 0;
  struct Xt_stripe_minmax query_minmax = stripe2minmax(query);
  size_t num_db_stripes = db->num_stripes;
  const struct Xt_stripe *restrict db_stripes = db->stripes;
  const struct Xt_stripe_minmax *restrict db_minmax = db->minmax;
  const int *restrict db_stripes_nstrides_psum = db->stripes_nstrides_psum;
  for (size_t j = db_start; j < num_db_stripes; ++j) {
    if ((query_minmax.max >= db_minmax[j].min)
        & (query_minmax.min <= db_minmax[j].max)) {
      struct Xt_stripe db_stripe = db_stripes[j];
      Xt_int stride = query.stride;
      /* determine if db_stripe and query can be easily aligned */
      if ((stride > 0)
          & (stride == db_stripe.stride)
          & ((query.start - db_stripe.start) % stride == 0)) {
        /* divide query into skipped, matching and left-over parts,
         * where skipped and left-over are non-matching */
        Xt_int overlap_start = MAX(query.start, db_stripe.start);
        /* handle skipped query part */
        int skipLen = (int)((overlap_start - query.start) / stride);
        if (skipLen)
        {
          struct Xt_stripe query_head = {
            .start = query.start,
            .stride = stride,
            .nstrides = skipLen };
          unmatched
            += idxstripes_get_pos_exts_of_index_stripe(
              query_head, db,
              num_ext_, size_ext, pos_exts, single_match_only, j + 1);
          query.start = (Xt_int)(query.start + stride * (Xt_int)skipLen);
          query.nstrides -= skipLen;
        }
        int db_stripe_skip
          = (int)((overlap_start - db_stripe.start) / stride);
        int overlap_nstrides
          = imin(query.nstrides, db_stripe.nstrides - db_stripe_skip);

        unmatched +=
          pos_ext_insert(
            (struct Xt_stripe){ .start = query.start, .stride = query.stride,
                .nstrides = overlap_nstrides},
            (struct Xt_pos_ext){ .start
                = db_stripes_nstrides_psum[j] + db_stripe_skip,
                .size = overlap_nstrides
                }, db, num_ext_, size_ext, pos_exts, single_match_only, j);

        int query_leftover = query.nstrides - overlap_nstrides;
        if (query_leftover)
        {
          query = (struct Xt_stripe){
            .start = (Xt_int)(overlap_start + stride * overlap_nstrides),
            .stride = stride,
            .nstrides = query_leftover };
          query_minmax = stripe2minmax(query);
          continue;
        }
        else
          goto search_finished;
      }
      else
      {
        /* Handle complex overlap */
        struct unmatched_tail search_result =
          idxstripes_complex_get_pos_exts_of_index_stripe(
            query, db, num_ext_, size_ext, pos_exts, single_match_only, j);
        unmatched += search_result.unmatched;
        if (search_result.query_tail.nstrides)
        {
          query = search_result.query_tail;
          continue;
        }
        else
          goto search_finished;
      }
    }
  }
  /* query wasn't found, add all indices to unmatched */
  unmatched += (size_t)query.nstrides;
search_finished:
  return unmatched;
}

/* todo: add indices into pos_exts which are sorted by end/first of ranges */
static size_t
conditional_pos_ext_insert(struct Xt_stripe query,
                           struct Xt_pos_ext pos_ext2add,
                           const struct Xt_stripes_lookup *restrict db,
                           size_t *num_ext,
                           size_t *size_ext,
                           struct Xt_pos_ext **pos_exts,
                           size_t db_start)
{
  size_t unmatched = 0;
  struct Xt_pos_ext *restrict pos_exts_ = *pos_exts;
  size_t num_ext_ = *num_ext;
  Xt_int stride = query.stride;
  if (pos_ext2add.size == -1)
    pos_ext2add.size = 1;
  int querySizeMaskNeg = isign_mask(pos_ext2add.size);
  int pos_ext2add_s = pos_ext2add.start
    + (querySizeMaskNeg & (pos_ext2add.size + 1)),
    pos_ext2add_e = pos_ext2add.start
    + (~querySizeMaskNeg & (pos_ext2add.size - 1));
  /* single_match_only, never re-match positions */
  for (size_t k = 0; k < num_ext_; ++k) {
    int dbSizeMaskNeg = isign_mask(pos_exts_[k].size),
      db_s = pos_exts_[k].start
      + (dbSizeMaskNeg & (pos_exts_[k].size + 1)),
      db_e = pos_exts_[k].start
      + (~dbSizeMaskNeg & (pos_exts_[k].size - 1));
    /* does overlap exist? */
    if ((db_s > pos_ext2add_e) | (pos_ext2add_s > db_e))
      continue;
    else {
      /* determine length of overlap parts */
      int lowQuerySkip = db_s - pos_ext2add_s;
      int lowDbSkip = -lowQuerySkip;
      lowQuerySkip = (int)((unsigned)(lowQuerySkip + abs(lowQuerySkip))/2);
      lowDbSkip = (int)((unsigned)(lowDbSkip + abs(lowDbSkip))/2);
      int overlapLen = MIN(db_e - db_s - lowDbSkip + 1,
                           abs(pos_ext2add.size) - lowQuerySkip);
      int highQuerySkip = abs(pos_ext2add.size) - lowQuerySkip - overlapLen;
      /* then adjust lengths to direction of overlap (from
       * perspective of pos_ext2add */
      int querySkipLen = (~querySizeMaskNeg & lowQuerySkip)
        | (querySizeMaskNeg & -highQuerySkip),
        queryTailLen = (querySizeMaskNeg & -lowQuerySkip)
        | (~querySizeMaskNeg & highQuerySkip);
      if (querySkipLen)
      {
        struct Xt_stripe query_skip = {
          .start = query.start,
          .stride = stride,
          .nstrides = abs(querySkipLen),
        };
        struct Xt_pos_ext pos_ext2add_skip = {
          .start = pos_ext2add.start,
          .size = querySkipLen
        };
        unmatched
          += conditional_pos_ext_insert(
            query_skip, pos_ext2add_skip, db,
            num_ext,size_ext, pos_exts, db_start + 1);
        num_ext_ = *num_ext;
        pos_exts_ = *pos_exts;
        query.start = (Xt_int)(query.start
                               + stride * (Xt_int)abs(querySkipLen));
        query.nstrides -= abs(querySkipLen);
        pos_ext2add.start += querySkipLen;
        pos_ext2add.size -= querySkipLen;
      }
      /* head part of (remaining) query matches part of already inserted index
       * range */
      struct Xt_stripe query_head = {
        .start = query.start,
        .stride = stride,
        .nstrides = abs(overlapLen) };
      /* restart search for overlapping part within following ranges */
      unmatched
        += idxstripes_get_pos_exts_of_index_stripe(
          query_head, db, num_ext, size_ext, pos_exts, true, db_start + 1);
      num_ext_ = *num_ext;
      pos_exts_ = *pos_exts;
      if (queryTailLen) {
        /* shorten query accordingly */
        query.nstrides -= abs(overlapLen);
        query.start = (Xt_int)(query.start + stride * (Xt_int)abs(overlapLen));
        pos_ext2add.start += overlapLen;
        pos_ext2add.size -= overlapLen;
      }
      else
        /* whole range handled, return */
        return unmatched;
    }
  }
  /* remaining extent does not overlap any existing one */
  append_ext(pos_ext2add, num_ext, size_ext, pos_exts);
  return unmatched;
}

static struct unmatched_tail
idxstripes_complex_get_pos_exts_of_index_stripe(
  struct Xt_stripe query,
  const struct Xt_stripes_lookup *restrict db,
  size_t *num_ext_,
  size_t *size_ext,
  struct Xt_pos_ext **pos_exts,
  bool single_match_only,
  size_t db_stripe_pos)
{
  size_t unmatched = 0;
  const struct Xt_stripe *restrict db_stripes = db->stripes;
  struct Xt_stripe overlap = get_stripe_intersection(query,
                                                     db_stripes[db_stripe_pos]);
  if (overlap.nstrides == 0)
    return (struct unmatched_tail){ .unmatched = 0, .query_tail = query};

  int skipped = (int)((overlap.start - query.start)/query.stride);
  if (skipped)
  {
    unmatched += idxstripes_get_pos_exts_of_index_stripe(
      (struct Xt_stripe){ .start = query.start, .stride = query.stride,
          .nstrides = skipped}, db, num_ext_, size_ext, pos_exts,
      single_match_only, db_stripe_pos + 1);
    query.start = (Xt_int)(query.start + skipped * query.stride);
    query.nstrides -= skipped;
  }
  /* Since overlap.nstrides > 0, overlap.start always matches, but
   * depending on the stride the remaining overlapping indices might or might
   * not be consecutive in the query, in the latter case intervening
   * parts need to be searched for too.
   * Stripes of length 1 are naturally consecutive.
   */
  int db_stripe_skip
    = (int)((overlap.start - db_stripes[db_stripe_pos].start)
            / db_stripes[db_stripe_pos].stride);
  int db_pos = db->stripes_nstrides_psum[db_stripe_pos] + db_stripe_skip;
  if (((overlap.stride == query.stride)
       & (overlap.stride == db_stripes[db_stripe_pos].stride))
      | (overlap.nstrides == 1))
  {
    unmatched += pos_ext_insert(overlap, (struct Xt_pos_ext){
        .start = db_pos,
          .size = overlap.nstrides },
      db, num_ext_, size_ext, pos_exts, single_match_only, db_stripe_pos);
    query.nstrides -= overlap.nstrides;
    query.start = (Xt_int)(query.start + overlap.nstrides * query.stride);
  }
  else if ((overlap.stride == query.stride)
           & (overlap.stride == -db_stripes[db_stripe_pos].stride))
  {
    /* all parts of the overlap can be used directly,
       but are inversely sorted in db_stripe */
    unmatched += pos_ext_insert(overlap, (struct Xt_pos_ext){
        .start = db_pos, .size = -overlap.nstrides },
      db, num_ext_, size_ext, pos_exts, single_match_only, db_stripe_pos);
    query.nstrides -= overlap.nstrides;
    query.start = (Xt_int)(query.start + overlap.nstrides * query.stride);
  }
  else if (overlap.stride == query.stride)
  {
    /* all parts of the overlap can be used but are non-consecutive
     * in db_stripe */
    int db_step = (int)(overlap.stride/db_stripes[db_stripe_pos].stride);
    for (int i = 0; i < overlap.nstrides; ++i, db_pos += db_step)
      unmatched += pos_ext_insert(overlap, (struct Xt_pos_ext){
          .start = db_pos, .size = 1 },
        db, num_ext_, size_ext, pos_exts, single_match_only, db_stripe_pos);
    query.nstrides -= overlap.nstrides;
    query.start = (Xt_int)(query.start + overlap.nstrides * query.stride);
  }
  else /* overlap.stride != query.stride => overlap.stride > query.stride */
  {
    /*
     * query.start = (Xt_int)(query.start
     *                        + overlap.stride * overlap.nstrides);
     */
    int stride_step = (int)(overlap.stride / query.stride);
    int db_step = (int)(overlap.stride/db_stripes[db_stripe_pos].stride);
    do {
      struct Xt_stripe consecutive_overlap = { .start = query.start,
                                               .stride = 1,
                                               .nstrides = 1 },
        intervening = { .start = (Xt_int)(query.start + query.stride),
                        .stride = query.stride,
                        .nstrides = stride_step - 1 };
      /* split off start index, then handle intervening */
      unmatched +=
        pos_ext_insert(consecutive_overlap, (struct Xt_pos_ext){
            .start = db_pos, .size = 1 },
          db, num_ext_, size_ext, pos_exts, single_match_only,
          db_stripe_pos + 1);
      db_pos += db_step;
      unmatched +=
        idxstripes_get_pos_exts_of_index_stripe(
          intervening,
          db, num_ext_, size_ext, pos_exts,
          single_match_only, db_stripe_pos + 1);
      query.nstrides -= stride_step;
      query.start = (Xt_int)(query.start + query.stride * stride_step);
      overlap.start = (Xt_int)(overlap.start + overlap.stride);
    } while (--overlap.nstrides);
  }
  return (struct unmatched_tail){ .unmatched = unmatched, .query_tail = query};
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
