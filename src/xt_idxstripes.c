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

static void
idxstripes_delete(xt_idxlist data);

static size_t
idxstripes_get_pack_size(xt_idxlist data, MPI_Comm comm);

static void
idxstripes_pack(xt_idxlist data, void *buffer, int buffer_size,
                int *position, MPI_Comm comm);

static xt_idxlist
idxstripes_get_intersection(xt_idxlist idxlist_a, xt_idxlist idxlist_b);

static xt_idxlist
idxstripes_copy(xt_idxlist idxlist);

static size_t
idxstripes_get_num_indices(xt_idxlist idxlist);

static void
idxstripes_get_indices(xt_idxlist idxlist, Xt_idx *indices);

static void
idxsection_get_index_stripes(xt_idxlist idxlist, struct Xt_stripe ** stripes,
                             unsigned * num_stripes);

static int
idxstripes_get_index(Xt_idxlist idxlist, int position, Xt_idx * index);

static int
idxstripes_get_position_of_index(Xt_idxlist idxlist, Xt_idx index, int * position);

static int
idxstripes_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                     int * position, int offset);

static const struct xt_idxlist_vtable idxstripes_vtable = {
  .delete                      = idxstripes_delete,
  .get_pack_size               = idxstripes_get_pack_size,
  .pack                        = idxstripes_pack,
  .get_intersection            = idxstripes_get_intersection,
  .copy                        = idxstripes_copy,
  .get_num_indices             = idxstripes_get_num_indices,
  .get_indices                 = idxstripes_get_indices,
  .get_index_stripes           = idxsection_get_index_stripes,
  .get_index                   = idxstripes_get_index,
  .get_position_of_index       = idxstripes_get_position_of_index,
  .get_position_of_indices     = NULL,
  .get_position_of_index_off   = idxstripes_get_position_of_index_off,
  .get_position_of_indices_off = NULL};

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
    MPI_Datatype types[5] = {MPI_LB, Xt_idx_dt, MPI_INT, Xt_idx_dt, MPI_UB};

    mpi_err_handler(MPI_Type_create_struct(5, block_lengths, displacements,
                                           types, &type), MPI_COMM_WORLD);
    mpi_err_handler(MPI_Type_commit(&type), MPI_COMM_WORLD);
  }

  return type;
}

struct xt_idxstripes {

  const struct xt_idxlist_vtable *vtable;

  struct Xt_stripe * stripes;
  unsigned num_stripes;
};

Xt_idxlist xt_idxstripes_new(struct Xt_stripe const * stripes, unsigned num_stripes) {

  struct xt_idxstripes *idxstripes;

  idxstripes = xmalloc(1 * sizeof(*idxstripes));

  idxstripes->vtable      = &idxstripes_vtable;
  idxstripes->num_stripes = num_stripes;

  if (num_stripes > 0) {
    idxstripes->stripes = xmalloc(num_stripes * sizeof(*(idxstripes->stripes)));
    memcpy(idxstripes->stripes, stripes, num_stripes * sizeof(*(idxstripes->stripes)));
    int sign_err = 0;
    for(int i=0; i<num_stripes; i++) {
      if (idxstripes->stripes[i].nstrides<0) {
        idxstripes->stripes[i].nstrides = 0;
        sign_err++;
      }
    }
    if (sign_err) fprintf(stderr, "WARNING: xt_idxstripes_new called with invalid stripes\n");
  } else
    idxstripes->stripes = NULL;

  return (Xt_idxlist)idxstripes;
}

static void
idxstripes_delete(xt_idxlist data) {

  struct xt_idxstripes *stripes;

  if (data == NULL) return;

  stripes = (struct xt_idxstripes*)data;

  free(stripes->stripes);
  free(stripes);
}

static size_t
idxstripes_get_pack_size(xt_idxlist data, MPI_Comm comm) {

  struct xt_idxstripes *stripes = (struct xt_idxstripes*)data;

  int size_int_type, size_unsigned, size_stripes;

  mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);
  mpi_err_handler(MPI_Pack_size(1, MPI_UNSIGNED, comm, &size_unsigned), comm);
  if (stripes->num_stripes)
    mpi_err_handler(MPI_Pack_size(stripes->num_stripes, xt_get_stripes_datatype(),
                                  comm, &size_stripes), comm);

  return size_int_type + size_unsigned + size_stripes;
}

static void
idxstripes_pack(xt_idxlist data, void *buffer, int buffer_size,
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

xt_idxlist xt_idxstripes_unpack(void *buffer, int buffer_size, int *position,
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
  } else
    stripes->stripes = NULL;

  return (xt_idxlist)stripes;
}

static int compare_xtidx(const void * a, const void * b)
{
  return (*(Xt_idx*)a > *(Xt_idx*)b) - (*(Xt_idx*)a < *(Xt_idx*)b);
}

void xt_idxstripes_convert_to_stripes(Xt_idx const * indices,
                                      unsigned num_indices,
                                      struct Xt_stripe ** stripes,
                                      unsigned * num_stripes) {

  struct Xt_stripe * temp_stripes = NULL;
  unsigned temp_stripes_array_size = 0;
  unsigned num_temp_stripes = 0;

  unsigned i, j;

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
compute_intersection_brute_force(xt_idxlist idxlist_src,
                                 xt_idxlist idxlist_dst) {

  size_t src_size, dst_size;

  src_size = idxstripes_get_num_indices(idxlist_src);
  dst_size = idxstripes_get_num_indices(idxlist_dst);

  Xt_idx * src_indices, * dst_indices;

  src_indices = xmalloc(src_size * sizeof(*src_indices));
  dst_indices = xmalloc(dst_size * sizeof(*dst_indices));

  // get array of indices of source and target

  idxstripes_get_indices(idxlist_src, src_indices);
  idxstripes_get_indices(idxlist_dst, dst_indices);

  // sort indices
  qsort(src_indices, src_size, sizeof (*src_indices), compare_xtidx);
  qsort(dst_indices, dst_size, sizeof (*dst_indices), compare_xtidx);

  size_t inter_size = dst_size;

  Xt_idx * inter_indices;

  inter_indices = xmalloc(inter_size * sizeof(*inter_indices));

  size_t i, j;

  // compute intersection
  for (i = 0, j = 0, inter_size = 0; i < dst_size && j < src_size; ++i) {

    while (j < src_size && src_indices[j] < dst_indices[i]) ++j;

    if (j < src_size && src_indices[j] == dst_indices[i])
      inter_indices[inter_size++] = src_indices[j];
  }

  struct Xt_stripe * inter_stripes;
  unsigned num_inter_stripes;

  // generate new idx stripes
  xt_idxstripes_convert_to_stripes(inter_indices, inter_size,
                                   &inter_stripes, &num_inter_stripes);

  Xt_idxlist intersection;

  intersection = xt_idxstripes_new(inter_stripes, num_inter_stripes);

  free(inter_stripes);
  free(inter_indices);
  free(src_indices);
  free(dst_indices);

  return intersection;
}

static int check_for_simple_strides(struct xt_idxstripes * stripe) {

  if (stripe->num_stripes == 0)
    return 0;

  unsigned i;

  if (stripe->stripes[0].stride != 1)
    return 0;

  for (i = 1; i < stripe->num_stripes; ++i)
    if ((stripe->stripes[i].stride != 1) &&
        ((stripe->stripes[i-1].start + stripe->stripes[i-1].nstrides) <=
         stripe->stripes[i].start))
      return 0;

  return 1;
}

static Xt_idxlist
compute_intersection_simple(struct xt_idxstripes *idxstripes_src,
                            struct xt_idxstripes *idxstripes_dst) {

  unsigned src_pos, dst_pos;

  src_pos = 0;
  dst_pos = 0;

  struct Xt_stripe * inter_stripes = NULL;
  unsigned inter_stripes_array_size = 0;
  unsigned num_inter_stripes = 0;

  while (src_pos < idxstripes_src->num_stripes &&
         dst_pos < idxstripes_dst->num_stripes) {

    while (src_pos < idxstripes_src->num_stripes &&
           idxstripes_src->stripes[src_pos].start +
           idxstripes_src->stripes[src_pos].nstrides <=
           idxstripes_dst->stripes[dst_pos].start) ++src_pos;

      if (src_pos >= idxstripes_src->num_stripes) break;

    while (dst_pos < idxstripes_dst->num_stripes &&
           idxstripes_dst->stripes[dst_pos].start +
           idxstripes_dst->stripes[dst_pos].nstrides <=
           idxstripes_src->stripes[src_pos].start) ++dst_pos;

      if (dst_pos >= idxstripes_dst->num_stripes) break;

    Xt_idx start, end;

    start = (idxstripes_src->stripes[src_pos].start >
             idxstripes_dst->stripes[dst_pos].start)?
            (idxstripes_src->stripes[src_pos].start):
            (idxstripes_dst->stripes[dst_pos].start);

    if ((idxstripes_src->stripes[src_pos].start +
         idxstripes_src->stripes[src_pos].nstrides) <
        (idxstripes_dst->stripes[dst_pos].start +
         idxstripes_dst->stripes[dst_pos].nstrides)) {

      end = idxstripes_src->stripes[src_pos].start +
            idxstripes_src->stripes[src_pos].nstrides;
      src_pos++;
    } else {
      end = idxstripes_dst->stripes[dst_pos].start +
            idxstripes_dst->stripes[dst_pos].nstrides;
      dst_pos++;
    }

    if (end > start) {

      num_inter_stripes++;

      ENSURE_ARRAY_SIZE(inter_stripes, inter_stripes_array_size, num_inter_stripes);

      inter_stripes[num_inter_stripes-1].start  = start;
      inter_stripes[num_inter_stripes-1].nstrides = end - start;
      inter_stripes[num_inter_stripes-1].stride = 1;
    }
  }

  Xt_idxlist intersection;

  intersection = xt_idxstripes_new(inter_stripes, num_inter_stripes);

  free(inter_stripes);

  return intersection;
}

static Xt_idxlist
idxstripes_get_intersection(xt_idxlist idxlist_src, xt_idxlist idxlist_dst) {

  if (idxlist_dst->vtable != &idxstripes_vtable) return NULL;

  // both lists are index stripes:

  struct xt_idxstripes *idxstripes_src, *idxstripes_dst;

  idxstripes_src = (struct xt_idxstripes*)idxlist_src;
  idxstripes_dst = (struct xt_idxstripes*)idxlist_dst;

  // check for requirements of a simple intersection solution
  if (check_for_simple_strides(idxstripes_src) &&
      check_for_simple_strides(idxstripes_dst)) {

    return compute_intersection_simple(idxstripes_src, idxstripes_dst);
  } else // compute the intersection using brute force
    return compute_intersection_brute_force(idxlist_src, idxlist_dst);
}

static xt_idxlist
idxstripes_copy(xt_idxlist idxlist) {

  struct xt_idxstripes *stripes;

  stripes = (struct xt_idxstripes*)idxlist;

  return xt_idxstripes_new(stripes->stripes, stripes->num_stripes);
}

static size_t
idxstripes_get_num_indices(xt_idxlist idxlist) {

  struct xt_idxstripes * stripes;

  stripes = (struct xt_idxstripes *)idxlist;

  unsigned i;
  size_t size = 0;

  for (i = 0; i < stripes->num_stripes; ++i)
    size += stripes->stripes[i].nstrides;

  return size;
}

static void
idxstripes_get_indices(xt_idxlist idxlist, Xt_idx *indices) {

  struct xt_idxstripes *stripes;

  stripes = (struct xt_idxstripes*)idxlist;

  unsigned i, j;

  for (i = 0; i < stripes->num_stripes; ++i)
    for (j = 0; j < stripes->stripes[i].nstrides; ++j)
      *(indices++) = stripes->stripes[i].start + j * stripes->stripes[i].stride;
}

static void
idxsection_get_index_stripes(xt_idxlist idxlist, struct Xt_stripe ** stripes,
                             unsigned * num_stripes) {

  struct xt_idxstripes *idxstripes;

  idxstripes = (struct xt_idxstripes*)idxlist;

  struct Xt_stripe * temp_stripes = NULL;
  unsigned temp_stripes_array_size = 0;
  unsigned num_temp_stripes = 0;

  unsigned i, j;

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
idxstripes_get_index(Xt_idxlist idxlist, int position, Xt_idx * index) {

  struct xt_idxstripes *stripes;

  stripes = (struct xt_idxstripes*)idxlist;

  if (position < 0) return 1;

  unsigned i;

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

  unsigned i;
  int position_offset = 0;

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
