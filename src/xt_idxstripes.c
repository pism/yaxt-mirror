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

typedef struct Xt_idxstripes_ *Xt_idxstripes;

struct Xt_stripes_sort {
  struct Xt_stripe_minmax range; /* minimal and maximal position of stripe */
  int position;                  /* position of stripe this range
                                  * corresponds to (permutation
                                  * obtained from sorting) */
  int inv_position;              /* position in stripes_sort when
                                  * indexed with unsorted index */
};

struct Xt_idxstripes_ {

  struct Xt_idxlist_ parent;

  struct Xt_stripe * stripes;
  Xt_int min_index, max_index;
  int num_stripes;
  bool stripes_do_overlap;

  Xt_int *index_array_cache;
  struct Xt_stripes_sort stripes_sort[];
};

static int compare_xtstripes(const void * a_, const void * b_)
{
  const struct Xt_stripes_sort *restrict a = a_, *restrict b = b_;
  return ((a->range.min > b->range.min) -
          (a->range.min < b->range.min));
}

static void
idxstripes_aggregate(Xt_idxstripes idxstripes,
                     const char *caller)
{
  struct Xt_stripe *restrict stripes = idxstripes->stripes;
  struct Xt_stripes_sort *restrict stripes_sort = idxstripes->stripes_sort;
  Xt_int min, max;
  {
    struct Xt_stripe_minmax stripe_range = stripe2minmax(stripes[0]);
    stripes_sort[0].range = stripe_range;
    stripes_sort[0].position = 0;
    min = stripe_range.min;
    max = stripe_range.max;
  }
  size_t num_stripes = (size_t)idxstripes->num_stripes;
  long long num_indices = (long long)stripes[0].nstrides;
  int sign_err = stripes[0].nstrides;
  for (size_t i = 1; i < num_stripes; ++i) {
    struct Xt_stripe_minmax stripe_range = stripe2minmax(stripes[i]);
    stripes_sort[i].range = stripe_range;
    stripes_sort[i].position = (int)i;
    min = MIN(stripe_range.min, min);
    max = MAX(stripe_range.max, max);
    num_indices += (long long)stripes[i].nstrides;
    sign_err |= stripes[i].nstrides;
  }
  /* shift sign bit into position 0 */
  sign_err >>= (sizeof (sign_err) * CHAR_BIT - 1);
  if (sign_err) {
    static const char template[] = "ERROR: %s called with invalid stripes";
    size_t buf_size = sizeof (template) - 2 + strlen(caller);
    char *msg = xmalloc(buf_size);
    snprintf(msg, buf_size, template, caller);
    die(msg);
  }
  assert(num_indices <= INT_MAX);
  qsort(stripes_sort, num_stripes, sizeof (*stripes_sort), compare_xtstripes);

  stripes_sort[stripes_sort[0].position].inv_position = 0;
  bool stripes_do_overlap = false;
  for (size_t i = 1; i < num_stripes; ++i) {
    stripes_do_overlap
      |= stripes_sort[i - 1].range.max >= stripes_sort[i].range.min;
    stripes_sort[stripes_sort[i].position].inv_position = (int)i;
  }

  idxstripes->stripes_do_overlap = stripes_do_overlap;
  idxstripes->min_index = min;
  idxstripes->max_index = max;
  idxstripes->index_array_cache = NULL;
  Xt_idxlist_init(&idxstripes->parent, &idxstripes_vtable, (int)num_indices);
}


Xt_idxlist
xt_idxstripes_new(struct Xt_stripe const * stripes, int num_stripes) {
  INSTR_DEF(instr,"xt_idxstripes_new");
  INSTR_START(instr);

  Xt_idxlist result;

  if (num_stripes > 0) {
    size_t header_size = ((sizeof (struct Xt_idxstripes_)
                           + (sizeof (struct Xt_stripes_sort)
                              * (size_t)num_stripes)
                           + sizeof (struct Xt_stripe) - 1)
                          / sizeof (struct Xt_stripe))
      * sizeof (struct Xt_stripe),
      body_size = sizeof (struct Xt_stripe) * (size_t)num_stripes;
    Xt_idxstripes idxstripes = xmalloc(header_size + body_size);
    idxstripes->num_stripes = num_stripes;
    idxstripes->stripes
      = (struct Xt_stripe *)((unsigned char *)idxstripes + header_size);
    memcpy(idxstripes->stripes, stripes,
           (size_t)num_stripes * sizeof(*(idxstripes->stripes)));

    idxstripes_aggregate(idxstripes, __func__);
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
                           + ((size_t)num_stripes
                              * sizeof (struct Xt_stripes_sort))
                           + sizeof (struct Xt_stripe) - 1)
                          / sizeof (struct Xt_stripe))
      * sizeof (struct Xt_stripe);
    Xt_idxstripes idxstripes = xmalloc(header_size);
    idxstripes->num_stripes = num_stripes;
    idxstripes->stripes = (struct Xt_stripe *)stripes;
    idxstripes_aggregate(idxstripes, __func__);
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

  int size_header, size_stripes =  0;

  xt_mpi_call(MPI_Pack_size(2, MPI_INT, comm, &size_header), comm);
  if (stripes->num_stripes)
    xt_mpi_call(MPI_Pack_size(stripes->num_stripes, stripe_dt, comm,
                              &size_stripes), comm);

  return (size_t)size_header + (size_t)size_stripes;
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
  int num_stripes = stripes->num_stripes;
  xt_mpi_call(MPI_Pack(&stripes->num_stripes, 1, MPI_INT, buffer,
                       buffer_size, position, comm), comm);
  if (num_stripes)
    xt_mpi_call(MPI_Pack(stripes->stripes, num_stripes, stripe_dt,
                         buffer, buffer_size, position, comm), comm);
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
                           + ((size_t)num_stripes
                              * sizeof (struct Xt_stripes_sort))
                           + sizeof (struct Xt_stripe) - 1)
                          / sizeof (struct Xt_stripe))
      * sizeof (struct Xt_stripe),
      body_size = sizeof (struct Xt_stripe) * (size_t)num_stripes;
    Xt_idxstripes idxstripes = xmalloc(header_size + body_size);
    idxstripes->num_stripes = num_stripes;
    idxstripes->stripes
      = (struct Xt_stripe *)((unsigned char *)idxstripes + header_size);
    xt_mpi_call(MPI_Unpack(buffer, buffer_size, position, idxstripes->stripes,
                           num_stripes, stripe_dt, comm),comm);
    idxstripes_aggregate(idxstripes, __func__);
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

  struct Xt_stripe *restrict inter_stripes = NULL;
  size_t num_inter_stripes = 0;
  size_t inter_stripes_array_size = 0;

  const struct Xt_stripes_sort *restrict src_stripes_sort
    = idxstripes_src->stripes_sort,
    *restrict dst_stripes_sort = idxstripes_dst->stripes_sort;
  const struct Xt_stripe *restrict stripes_src = idxstripes_src->stripes,
    *restrict stripes_dst = idxstripes_dst->stripes;

  size_t i_src = 0, i_dst = 0;
  size_t num_stripes_src = (size_t)idxstripes_src->num_stripes,
    num_stripes_dst = (size_t)idxstripes_dst->num_stripes;

  while ((i_src < num_stripes_src) &
         (i_dst < num_stripes_dst)) {

    while (i_src < num_stripes_src &&
           src_stripes_sort[i_src].range.max
           < dst_stripes_sort[i_dst].range.min) ++i_src;

    if ( i_src >= num_stripes_src ) break;

    while (i_dst < num_stripes_dst &&
           dst_stripes_sort[i_dst].range.max
           < src_stripes_sort[i_src].range.min) ++i_dst;

    if ( i_dst >= num_stripes_dst ) break;

    if ((src_stripes_sort[i_src].range.min
         <= dst_stripes_sort[i_dst].range.max)
        & (src_stripes_sort[i_src].range.max
           >= dst_stripes_sort[i_dst].range.min)) {
      ENSURE_ARRAY_SIZE(inter_stripes, inter_stripes_array_size,
                        num_inter_stripes+1);

      struct Xt_stripe intersection_stripe;
      inter_stripes[num_inter_stripes] = intersection_stripe =
        get_stripe_intersection(stripes_src[src_stripes_sort[i_src].position],
                                stripes_dst[dst_stripes_sort[i_dst].position]);
      num_inter_stripes += intersection_stripe.nstrides > 0;
    }

    if (dst_stripes_sort[i_dst].range.max
        < src_stripes_sort[i_src].range.max)
      i_dst++;
    else
      i_src++;
  }

  if (num_inter_stripes) {
    /* invert negative and merge consecutive stripes */
    struct Xt_stripe prev_stripe = inter_stripes[0];
    if (prev_stripe.stride < 0) {
      prev_stripe.start
        = (Xt_int)(prev_stripe.start
                   + (prev_stripe.stride * (Xt_int)(prev_stripe.nstrides - 1)));
      prev_stripe.stride = (Xt_int)-prev_stripe.stride;
    }
    inter_stripes[0] = prev_stripe;
    size_t j = 0;
    for (size_t i = 1; i < num_inter_stripes; ++i) {
      struct Xt_stripe stripe = inter_stripes[i];
      if (stripe.stride < 0) {
        stripe.start = (Xt_int)(stripe.start
                                + stripe.stride * (Xt_int)(stripe.nstrides - 1));
        stripe.stride = (Xt_int)-stripe.stride;
      }
      if ((stripe.stride == prev_stripe.stride)
          & (stripe.start
             == (prev_stripe.start
                 + prev_stripe.stride * (Xt_int)prev_stripe.nstrides)))
      {
        prev_stripe.nstrides += stripe.nstrides;
        inter_stripes[j].nstrides = prev_stripe.nstrides;
      }
      else
      {
        inter_stripes[++j] = stripe;
        prev_stripe = stripe;
      }
    }
    num_inter_stripes = j + 1;
  }

  Xt_idxlist inter = xt_idxstripes_new(inter_stripes, (int)num_inter_stripes);

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
append_ext(struct Xt_pos_ext pos_ext, struct Xt_pos_ext_vec *restrict result)
{
  size_t num_pos_exts_ = result->num_pos_ext,
    size_pos_exts_ = result->size_pos_ext;
  struct Xt_pos_ext *restrict pos_exts_ = result->pos_ext;
  if (xt_pos_ext_is_appendable(pos_exts_[num_pos_exts_ - 1], pos_ext))
  {
    if (num_pos_exts_ + 1 == size_pos_exts_)
    {
      size_pos_exts_ += 16;
      /* offsetting by 1 necessary to keep the terminator in place */
      result->pos_ext = pos_exts_ = (struct Xt_pos_ext *)
        xrealloc(pos_exts_ - 1, (size_pos_exts_ + 1) * sizeof (*pos_exts_)) + 1;
      result->size_pos_ext = size_pos_exts_;
    }
    pos_exts_[num_pos_exts_] = pos_ext;
    result->num_pos_ext = num_pos_exts_ + 1;
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
  struct Xt_pos_ext_vec *restrict result,
  struct Xt_pos_ext_vec *restrict cover,
  bool single_match_only,
  size_t db_start);

int
idxstripes_get_pos_exts_of_index_stripes(
  Xt_idxlist idxlist,
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
    result.size_pos_ext = (size_t)MAX(MIN(idxstripes->num_stripes,
                                          num_stripes), 8);
    result.pos_ext = xmalloc((result.size_pos_ext + 1)
                             * sizeof (*result.pos_ext));
    /* put non-concatenable *terminator* at offset -1 */
    result.pos_ext[0] = (struct Xt_pos_ext){.start = INT_MIN, .size = -1 };
    ++result.pos_ext;
    struct Xt_pos_ext_vec cover;
    xt_cover_start(&cover, result.size_pos_ext);
    struct Xt_stripes_lookup stripes_db;
    create_stripes_lookup(&stripes_db,
                          (size_t)idxstripes->num_stripes, idxstripes->stripes);
    for (size_t i = 0; i < (size_t)num_stripes; ++i) {
      struct Xt_stripe query = stripes[i];
      unmatched += idxstripes_get_pos_exts_of_index_stripe(
        query, &stripes_db, &result, &cover, single_match_only != 0, 0);
    }
    --(result.pos_ext);
    memmove(result.pos_ext, result.pos_ext + 1,
            sizeof (*result.pos_ext) * result.num_pos_ext);
    *pos_exts = xrealloc(result.pos_ext,
                         sizeof (*result.pos_ext) * result.num_pos_ext);
    destroy_stripes_lookup(&stripes_db);
    xt_cover_finish(&cover);
  }
  *num_ext = (int)result.num_pos_ext;
  return (int)unmatched;
}

static size_t
conditional_pos_ext_insert(struct Xt_stripe query,
                           struct Xt_pos_ext pos_ext2add,
                           const struct Xt_stripes_lookup *restrict db,
                           size_t db_start,
                           struct Xt_pos_ext_vec *restrict result,
                           struct Xt_pos_ext_vec *restrict cover,
                           size_t cover_start);

struct unmatched_tail
{
  size_t unmatched;
  struct Xt_stripe query_tail;
};

static struct unmatched_tail
idxstripes_complex_get_pos_exts_of_index_stripe(
  struct Xt_stripe query,
  const struct Xt_stripes_lookup *restrict stripes_lookup,
  struct Xt_pos_ext_vec *restrict result,
  struct Xt_pos_ext_vec *restrict cover,
  bool single_match_only,
  size_t db_start);

static size_t
pos_ext_insert(struct Xt_stripe query,
               struct Xt_pos_ext pos_ext2add,
               const struct Xt_stripes_lookup *stripes_lookup,
               struct Xt_pos_ext_vec *restrict result,
               struct Xt_pos_ext_vec *restrict cover,
               bool single_match_only,
               size_t db_pos)
{
  size_t unmatched = 0;
  if (single_match_only)
    unmatched += conditional_pos_ext_insert(
      query, pos_ext2add, stripes_lookup, db_pos, result, cover, 0);
  else
    append_ext(pos_ext2add, result);
  return unmatched;
}


size_t
idxstripes_get_pos_exts_of_index_stripe(
  struct Xt_stripe query,
  const struct Xt_stripes_lookup *restrict db,
  struct Xt_pos_ext_vec *restrict result,
  struct Xt_pos_ext_vec *restrict cover,
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
              query_head, db, result, cover, single_match_only, j + 1);
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
                }, db, result, cover, single_match_only, j);

        if (!(query.nstrides -= overlap_nstrides))
          goto search_finished;
        else {
          query.start = (Xt_int)(overlap_start + stride * overlap_nstrides);
          query_minmax = stripe2minmax(query);
          continue;
        }
      }
      else
      {
        /* Handle complex overlap */
        struct unmatched_tail search_result =
          idxstripes_complex_get_pos_exts_of_index_stripe(
            query, db, result, cover, single_match_only, j);
        unmatched += search_result.unmatched;
        if (!search_result.query_tail.nstrides)
          goto search_finished;
        else {
          query = search_result.query_tail;
          continue;
        }
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
                           size_t db_start,
                           struct Xt_pos_ext_vec *restrict result,
                           struct Xt_pos_ext_vec *restrict cover,
                           size_t cover_start)
{
  /* single_match_only is true => never re-match positions */
  size_t unmatched = 0;
  Xt_int stride = query.stride;
  if (pos_ext2add.size == -1)
    pos_ext2add.size = 1;
  int querySizeMaskNeg = isign_mask(pos_ext2add.size);
tail_search:
  ;
  int pos_ext2add_s = pos_ext2add.start
    + (querySizeMaskNeg & (pos_ext2add.size + 1)),
    pos_ext2add_e = pos_ext2add.start
    + (~querySizeMaskNeg & (pos_ext2add.size - 1));
  struct Xt_pos_range query_range
    = { .start = pos_ext2add_s, .end = pos_ext2add_e };
  /* does overlap exist? */
  size_t overlap_pos =
    xt_cover_insert_or_overlap(cover, query_range, true, cover_start);
  if (overlap_pos == SIZE_MAX) {
    /* remaining extent does not overlap any existing one */
    append_ext(pos_ext2add, result);
  } else {
    struct Xt_pos_ext *restrict pos_exts_ = cover->pos_ext;
    int dbSizeMaskNeg = isign_mask(pos_exts_[overlap_pos].size),
      db_s = pos_exts_[overlap_pos].start
      + (dbSizeMaskNeg & (pos_exts_[overlap_pos].size + 1)),
      db_e = pos_exts_[overlap_pos].start
      + (~dbSizeMaskNeg & (pos_exts_[overlap_pos].size - 1));
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
          query_skip, pos_ext2add_skip, db, db_start + 1,
          result, cover, overlap_pos + 1);
      pos_exts_ = result->pos_ext;
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
        query_head, db, result, cover, true, db_start + 1);
    pos_exts_ = result->pos_ext;
    if (queryTailLen) {
      /* shorten query accordingly */
      query.nstrides -= abs(overlapLen);
      query.start = (Xt_int)(query.start + stride * (Xt_int)abs(overlapLen));
      pos_ext2add.start += overlapLen;
      pos_ext2add.size -= overlapLen;
      cover_start = overlap_pos + 1;
      goto tail_search;
    }
    /* whole range handled, return */
  }
  return unmatched;
}

static struct unmatched_tail
idxstripes_complex_get_pos_exts_of_index_stripe(
  struct Xt_stripe query,
  const struct Xt_stripes_lookup *restrict db,
  struct Xt_pos_ext_vec *restrict result,
  struct Xt_pos_ext_vec *restrict cover,
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
          .nstrides = skipped}, db, result, cover,
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
      db, result, cover, single_match_only, db_stripe_pos);
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
      db, result, cover, single_match_only, db_stripe_pos);
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
        db, result, cover, single_match_only, db_stripe_pos);
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
          db, result, cover, single_match_only,
          db_stripe_pos + 1);
      db_pos += db_step;
      unmatched +=
        idxstripes_get_pos_exts_of_index_stripe(
          intervening, db, result, cover, single_match_only, db_stripe_pos + 1);
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
  return idxstripes->min_index;
}

static Xt_int
idxstripes_get_max_index(Xt_idxlist idxlist) {
  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
  return idxstripes->max_index;
}
