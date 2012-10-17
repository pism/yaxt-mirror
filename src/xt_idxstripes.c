#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xt_idxlist.h"
#include "xt_idxvec.h"
#include "xt_idxstripes.h"
#include "xt_mpi.h"
#include "core/ppm_xfuncs.h"
#include "ensure_array_size.h"

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
idxstripes_get_intersection(Xt_idxlist idxlist_a, Xt_idxlist idxlist_b);

static Xt_idxlist
idxstripes_copy(Xt_idxlist idxlist);

static Xt_count
idxstripes_get_num_indices(Xt_idxlist idxlist);

static void
idxstripes_get_indices(Xt_idxlist idxlist, Xt_idx *indices);

static void
idxsection_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                             Xt_count * num_stripes);

static int
idxstripes_get_index_at_position(Xt_idxlist idxlist, int position, Xt_idx * index);

static int
idxstripes_get_position_of_index(Xt_idxlist idxlist, Xt_idx index, int * position);

static int
idxstripes_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                     int * position, int offset);

static Xt_idx
idxstripes_get_min_index(Xt_idxlist idxlist);

static Xt_idx
idxstripes_get_max_index(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxstripes_vtable = {
  .delete                      = idxstripes_delete,
  .get_pack_size               = idxstripes_get_pack_size,
  .pack                        = idxstripes_pack,
  .get_intersection            = idxstripes_get_intersection,
  .copy                        = idxstripes_copy,
  .get_num_indices             = idxstripes_get_num_indices,
  .get_indices                 = idxstripes_get_indices,
  .get_index_stripes           = idxsection_get_index_stripes,
  .get_index_at_position       = idxstripes_get_index_at_position,
  .get_indices_at_positions    = NULL,
  .get_position_of_index       = idxstripes_get_position_of_index,
  .get_positions_of_indices    = NULL,
  .get_position_of_index_off   = idxstripes_get_position_of_index_off,
  .get_positions_of_indices_off = NULL,
  .get_min_index               = idxstripes_get_min_index,
  .get_max_index               = idxstripes_get_max_index};

static MPI_Datatype xt_get_stripes_datatype() {

  static MPI_Datatype type = MPI_DATATYPE_NULL;


  if (type == MPI_DATATYPE_NULL) {

  struct Xt_stripe stripe;

    MPI_Aint base_address, start_address, nstrides_address, stride_address;

    MPI_Get_address(&stripe, &base_address);
    MPI_Get_address(&stripe.start, &start_address);
    MPI_Get_address(&stripe.nstrides, &nstrides_address);
    MPI_Get_address(&stripe.stride, &stride_address);

    int block_lengths[5] = {1,1,1,1,1};
    MPI_Aint displacements[5] = {0, start_address - base_address,
                                 nstrides_address - base_address,
                                 stride_address - base_address,
                                 sizeof(stripe)};
    MPI_Datatype types[5] = {MPI_LB, Xt_idx_dt, Xt_count_dt, Xt_idx_dt, MPI_UB};

    mpi_err_handler(MPI_Type_create_struct(5, block_lengths, displacements,
                                           types, &type), MPI_COMM_WORLD);
    mpi_err_handler(MPI_Type_commit(&type), MPI_COMM_WORLD);
  }

  return type;
}

struct xt_idxstripes {

  const struct xt_idxlist_vtable *vtable;

  struct Xt_stripe * stripes;
  Xt_count num_stripes;

  struct Xt_stripe * stripes_sorted;
  int stripes_overlap;
};

static int compare_xtstripes(const void * a, const void * b)
{
  return (((struct Xt_stripe*)a)->start > ((struct Xt_stripe*)b)->start) -
         (((struct Xt_stripe*)a)->start < ((struct Xt_stripe*)b)->start);
}

static void generate_stripes_sorted(struct xt_idxstripes * stripes) {

  stripes->stripes_sorted = malloc(stripes->num_stripes * sizeof(*(stripes->stripes_sorted)));

  memcpy(stripes->stripes_sorted, stripes->stripes,
         stripes->num_stripes * sizeof(*(stripes->stripes_sorted)));

  Xt_count i;

  for (i = 0; i < stripes->num_stripes; ++i)
    if (stripes->stripes_sorted[i].stride < 0) {

      stripes->stripes_sorted[i].start = stripes->stripes_sorted[i].start +
                                         (stripes->stripes_sorted[i].nstrides-1) *
                                         stripes->stripes_sorted[i].stride;
      stripes->stripes_sorted[i].stride = -stripes->stripes_sorted[i].stride;
    }

  qsort(stripes->stripes_sorted, stripes->num_stripes,
        sizeof(*(stripes->stripes_sorted)), compare_xtstripes);
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

Xt_idxlist xt_idxstripes_new(struct Xt_stripe const * stripes, Xt_count num_stripes) {

  struct xt_idxstripes *idxstripes;

  idxstripes = xmalloc(1 * sizeof(*idxstripes));

  idxstripes->vtable      = &idxstripes_vtable;
  idxstripes->num_stripes = num_stripes;

  if (num_stripes > 0) {
    idxstripes->stripes = xmalloc(num_stripes * sizeof(*(idxstripes->stripes)));
    memcpy(idxstripes->stripes, stripes, num_stripes * sizeof(*(idxstripes->stripes)));
    int sign_err = 0;
    for(Xt_count i=0; i<num_stripes; i++) {
      if (idxstripes->stripes[i].nstrides<0) {
        idxstripes->stripes[i].nstrides = 0;
        sign_err++;
      }
    }
    if (sign_err) fprintf(stderr, "WARNING: xt_idxstripes_new called with invalid stripes\n");

    generate_stripes_sorted(idxstripes);

    idxstripes->stripes_overlap = 0;
    for (Xt_count i = 0; i < idxstripes->num_stripes - 1; ++i)
      if (stripes_overlap(idxstripes->stripes_sorted[i],
                          idxstripes->stripes_sorted[i+1])) {
        idxstripes->stripes_overlap = 1;
        break;
      }

  } else {
    idxstripes->stripes = NULL;
    idxstripes->stripes_sorted = NULL;
    idxstripes->stripes_overlap = 0;
  }

  return (Xt_idxlist)idxstripes;
}

static void
idxstripes_delete(Xt_idxlist data) {

  struct xt_idxstripes *stripes;

  if (data == NULL) return;

  stripes = (struct xt_idxstripes*)data;

  free(stripes->stripes_sorted);
  free(stripes->stripes);
  free(stripes);
}

static size_t
idxstripes_get_pack_size(Xt_idxlist data, MPI_Comm comm) {

  struct xt_idxstripes *stripes = (struct xt_idxstripes*)data;

  int size_int_type, size_unsigned, size_stripes;

  mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);
  mpi_err_handler(MPI_Pack_size(1, Xt_count_dt, comm, &size_unsigned), comm);
  if (stripes->num_stripes)
    mpi_err_handler(MPI_Pack_size(stripes->num_stripes, xt_get_stripes_datatype(),
                                  comm, &size_stripes), comm);

  return size_int_type + size_unsigned + size_stripes;
}

static void
idxstripes_pack(Xt_idxlist data, void *buffer, int buffer_size,
                int *position, MPI_Comm comm) {

  assert(data);
  struct xt_idxstripes *stripes = (struct xt_idxstripes*)data;
  int type = STRIPES;

  mpi_err_handler(MPI_Pack(&(type), 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(&(stripes->num_stripes), 1, MPI_UNSIGNED, buffer,
                           buffer_size, position, comm), comm);
  if (stripes->num_stripes)
    mpi_err_handler(MPI_Pack(stripes->stripes, (int)stripes->num_stripes,
                             xt_get_stripes_datatype(), buffer,
                             buffer_size, position, comm), comm);
}

Xt_idxlist xt_idxstripes_unpack(void *buffer, int buffer_size, int *position,
                                MPI_Comm comm) {

  struct xt_idxstripes *stripes;

  stripes = xmalloc(1 * sizeof(*stripes));

  stripes->vtable = &idxstripes_vtable;

  mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                             &(stripes->num_stripes), 1, MPI_UNSIGNED,
                             comm), comm);

  if (stripes->num_stripes) {
    stripes->stripes = xmalloc(stripes->num_stripes * sizeof(*(stripes->stripes)));
    mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                            stripes->stripes, stripes->num_stripes,
                            xt_get_stripes_datatype(), comm),comm);

    generate_stripes_sorted(stripes);

    stripes->stripes_overlap = 0;
    for (Xt_count i = 0; i < stripes->num_stripes - 1; ++i)
      if (stripes_overlap(stripes->stripes_sorted[i],
                          stripes->stripes_sorted[i+1])) {
        stripes->stripes_overlap = 1;
        break;
      }
  } else {
    stripes->stripes = NULL;
    stripes->stripes_sorted = NULL;
    stripes->stripes_overlap = 0;
  }

  return (Xt_idxlist)stripes;
}

void xt_idxstripes_convert_to_stripes(Xt_idx const * indices,
                                      Xt_count num_indices,
                                      struct Xt_stripe ** stripes,
                                      Xt_count * num_stripes) {

  struct Xt_stripe * temp_stripes = NULL;
  unsigned temp_stripes_array_size = 0;
  Xt_count num_temp_stripes = 0;

  Xt_count i, j;

  i = 0;
  while(i < num_indices) {

    ++num_temp_stripes;

    ENSURE_ARRAY_SIZE(temp_stripes, temp_stripes_array_size, num_temp_stripes);

    j = 1;

    while ((i + j) < num_indices && indices[i] + j == indices[i + j]) ++j;

    temp_stripes[num_temp_stripes-1].start  = indices[i];
    temp_stripes[num_temp_stripes-1].nstrides = j;
    temp_stripes[num_temp_stripes-1].stride = 1;

    i += j;
  }

  *stripes = xrealloc(temp_stripes, num_temp_stripes * sizeof(*temp_stripes));
  *num_stripes = num_temp_stripes;
}

static Xt_idxlist
compute_intersection_fallback(struct xt_idxstripes * idxstripes_src,
                              struct xt_idxstripes * idxstripes_dst) {

  Xt_idxlist idxvec_from_stripes_src;
  Xt_idxlist idxvec_from_stripes_dst;

  idxvec_from_stripes_src = xt_idxvec_from_stripes_new(idxstripes_src->stripes,
                                                       idxstripes_src->num_stripes);
  idxvec_from_stripes_dst = xt_idxvec_from_stripes_new(idxstripes_dst->stripes,
                                                       idxstripes_dst->num_stripes);

  Xt_idxlist intersection;

  intersection = xt_idxlist_get_intersection(idxvec_from_stripes_src,
                                             idxvec_from_stripes_dst);

  xt_idxlist_delete(idxvec_from_stripes_src);
  xt_idxlist_delete(idxvec_from_stripes_dst);

  return intersection;
}

// least common multiple
Xt_idx lcm(Xt_idx a, Xt_idx b) {

  Xt_idx temp;

  temp = a * b;

  while (a != b)
    if (a > b) a -= b;
    else       b -= a;

  return temp / a;
}

static struct Xt_stripe
get_stripe_intersection(struct Xt_stripe * stripe_a,
                        struct Xt_stripe * stripe_b) {

  Xt_count i;

  if (stripe_a->start > stripe_b->start) {
    struct Xt_stripe * temp = stripe_a;
    stripe_a = stripe_b;
    stripe_b = temp;
  }

  Xt_idx stripe_lcm = lcm(stripe_a->stride, stripe_b->stride);

  Xt_idx norm_start = stripe_b->start - stripe_a->start;

  for (i = 0; i < stripe_lcm / stripe_a->stride; ++i, norm_start += stripe_b->stride)
    if (!(norm_start % stripe_a->stride)) break;

  struct Xt_stripe intersection;

  if (i == stripe_lcm / stripe_a->stride) {

    intersection.start = 0;
    intersection.stride = 0;
    intersection.nstrides = 0;

  } else {

    Xt_idx stripe_a_end, stripe_b_end;

    intersection.start = norm_start + stripe_a->start;
    intersection.stride = stripe_lcm;
    stripe_a_end = stripe_a->start + stripe_a->stride * (stripe_a->nstrides-1);
    stripe_b_end = stripe_b->start + stripe_b->stride * (stripe_b->nstrides-1);
    intersection.nstrides = (MIN(stripe_a_end, stripe_b_end) - intersection.start + 1 +
                             (stripe_lcm - 1)) /
                            stripe_lcm;
  }

  return intersection;
}

// this routine only works for idxstripes that have stripes_overlap==0
static Xt_idxlist
idxstripes_compute_intersection(struct xt_idxstripes *idxstripes_src,
                                struct xt_idxstripes *idxstripes_dst) {

  struct Xt_stripe * inter_stripes = NULL;
  Xt_count num_inter_stripes = 0;
  Xt_count inter_stripes_array_size = 0;

  Xt_count i_src, i_dst;

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

      ENSURE_ARRAY_SIZE(inter_stripes, inter_stripes_array_size, num_inter_stripes+1);

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

  inter = xt_idxstripes_new(inter_stripes, num_inter_stripes);

  free(inter_stripes);

  return inter;
}

static Xt_idxlist
idxstripes_get_intersection(Xt_idxlist idxlist_src, Xt_idxlist idxlist_dst) {

  if (idxlist_dst->vtable != &idxstripes_vtable) return NULL;

  // both lists are index stripes:

  struct xt_idxstripes *idxstripes_src, *idxstripes_dst;

  idxstripes_src = (struct xt_idxstripes*)idxlist_src;
  idxstripes_dst = (struct xt_idxstripes*)idxlist_dst;

  if (idxstripes_src->stripes_overlap ||
      idxstripes_dst->stripes_overlap) {

    return compute_intersection_fallback(idxstripes_src, idxstripes_dst);
  } else
    return idxstripes_compute_intersection(idxstripes_src, idxstripes_dst);
}

static Xt_idxlist
idxstripes_copy(Xt_idxlist idxlist) {

  struct xt_idxstripes *stripes;

  stripes = (struct xt_idxstripes*)idxlist;

  return xt_idxstripes_new(stripes->stripes, stripes->num_stripes);
}

static Xt_count
idxstripes_get_num_indices(Xt_idxlist idxlist) {

  struct xt_idxstripes * stripes;

  stripes = (struct xt_idxstripes *)idxlist;

  Xt_count i;
  Xt_count size = 0;

  for (i = 0; i < stripes->num_stripes; ++i)
    size += stripes->stripes[i].nstrides;

  return size;
}

static void
idxstripes_get_indices(Xt_idxlist idxlist, Xt_idx *indices) {

  struct xt_idxstripes *stripes;

  stripes = (struct xt_idxstripes*)idxlist;

  Xt_count i, j;

  for (i = 0; i < stripes->num_stripes; ++i)
    for (j = 0; j < stripes->stripes[i].nstrides; ++j)
      *(indices++) = stripes->stripes[i].start + j * stripes->stripes[i].stride;
}

static void
idxsection_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                             Xt_count * num_stripes) {

  struct xt_idxstripes *idxstripes;

  idxstripes = (struct xt_idxstripes*)idxlist;

  struct Xt_stripe * temp_stripes = NULL;
  unsigned temp_stripes_array_size = 0;
  Xt_count num_temp_stripes = 0;

  Xt_count i, j;

  for (i = 0; i < idxstripes->num_stripes; ++i) {

    if (idxstripes->stripes[i].stride == 1) {

      ++num_temp_stripes;

      ENSURE_ARRAY_SIZE(temp_stripes, temp_stripes_array_size, num_temp_stripes);

      temp_stripes[num_temp_stripes-1] = idxstripes->stripes[i];

    } else {

      ENSURE_ARRAY_SIZE(temp_stripes, temp_stripes_array_size,
                        num_temp_stripes + idxstripes->stripes[i].nstrides);

      for (j = 0; j < idxstripes->stripes[i].nstrides; ++j)  {

        temp_stripes[num_temp_stripes].start = idxstripes->stripes[i].start +
                                               j * idxstripes->stripes[i].stride;
        temp_stripes[num_temp_stripes].nstrides = 1;
        temp_stripes[num_temp_stripes].stride = 1;

        ++num_temp_stripes;
      }
    }
  }

  *stripes = xrealloc(temp_stripes, num_temp_stripes * sizeof(*temp_stripes));
  *num_stripes = num_temp_stripes;
}

static int
idxstripes_get_index_at_position(Xt_idxlist idxlist, int position, Xt_idx * index) {

  struct xt_idxstripes *stripes;

  stripes = (struct xt_idxstripes*)idxlist;

  if (position < 0) return 1;

  Xt_count i;

  for (i = 0; i < stripes->num_stripes; ++i)
    if (position >= stripes->stripes[i].nstrides)
      position-= stripes->stripes[i].nstrides;
    else {
      *index = stripes->stripes[i].start + position * stripes->stripes[i].stride;
      return 0;
    }

  return 1;
}

static int
idxstripes_get_position_of_index(Xt_idxlist idxlist, Xt_idx index, int * position) {

  struct xt_idxstripes *stripes;

  stripes = (struct xt_idxstripes*)idxlist;

  Xt_count i;
  Xt_count position_offset = 0;

  for (i = 0; i < stripes->num_stripes;
       position_offset += stripes->stripes[i++].nstrides) {

    if ((stripes->stripes[i].stride > 0 && index < stripes->stripes[i].start) ||
        (stripes->stripes[i].stride < 0 && index > stripes->stripes[i].start))continue;

    Xt_idx rel_start;

    rel_start = index - stripes->stripes[i].start;

    if (rel_start%stripes->stripes[i].stride) continue;

    if (rel_start/stripes->stripes[i].stride > stripes->stripes[i].nstrides) continue;

    *position = rel_start/stripes->stripes[i].stride + position_offset;

    return 0;
  }

  return 1;
}

static int
idxstripes_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                     int * position, int offset) {

  int temp_position;

  if (idxstripes_get_position_of_index(idxlist, index, &temp_position))
    return 1;

  if (temp_position < offset)
    return 1;

  *position = temp_position;

  return 0;
}

static Xt_idx
idxstripes_get_min_index(Xt_idxlist idxlist) {

  struct xt_idxstripes * idxstripes;

  idxstripes = (struct xt_idxstripes *)idxlist;

  if (idxstripes->num_stripes > 0)
    return idxstripes->stripes_sorted[0].start;
  else
    return 0;
}

static Xt_idx
idxstripes_get_max_index(Xt_idxlist idxlist) {

  struct xt_idxstripes * idxstripes;

  idxstripes = (struct xt_idxstripes *)idxlist;

  Xt_count i;
  Xt_idx index = idxstripes_get_min_index(idxlist);

  for (i = 0; i < idxstripes->num_stripes; ++i)
    index = MAX(index, idxstripes->stripes_sorted[i].start +
                       idxstripes->stripes_sorted[i].stride *
                       (idxstripes->stripes_sorted[i].nstrides-1));

  return index;
}
