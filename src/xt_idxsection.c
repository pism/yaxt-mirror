#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "xt_idxlist.h"
#include "xt_idxvec.h"
#include "xt_idxsection.h"
#include "xt_mpi.h"
#include "core/ppm_xfuncs.h"

static void
idxsection_delete(xt_idxlist data);

static size_t
idxsection_get_pack_size(xt_idxlist data, MPI_Comm comm);

static void
idxsection_pack(xt_idxlist data, void *buffer, int buffer_size,
                int *position, MPI_Comm comm);

static xt_idxlist
idxsection_get_intersection(xt_idxlist idxlist_a, xt_idxlist idxlist_b);

static xt_idxlist
idxsection_copy(xt_idxlist idxlist);

static size_t
idxsection_get_num_indices(xt_idxlist idxlist);

static void
idxsection_get_indices(xt_idxlist idxlist, Xt_idx *indices);

static void
idxsection_get_index_stripes(xt_idxlist idxlist, struct Xt_stripe ** stripes,
                             unsigned * num_stripes);

static int
idxsection_get_index(Xt_idxlist idxlist, int position, Xt_idx * index);

static int
idxsection_get_position_of_index(Xt_idxlist idxlist, Xt_idx index, int * position);

static int
idxsection_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                     int * position, int offset);

static const struct xt_idxlist_vtable idxsection_vtable = {
              idxsection_delete,
              idxsection_get_pack_size,
              idxsection_pack,
              idxsection_get_intersection,
              idxsection_copy,
              idxsection_get_num_indices,
              idxsection_get_indices,
              idxsection_get_index_stripes,
              idxsection_get_index,
              idxsection_get_position_of_index,
              idxsection_get_position_of_index_off};

struct xt_idxsection {

  const struct xt_idxlist_vtable *vtable;

  Xt_idx global_start_index;
  Xt_idx local_start_index;
  unsigned ndim;
  unsigned * global_size;
  unsigned * local_size;
  unsigned * local_start;
  unsigned * global_strides;
  unsigned * local_strides;
};

Xt_idxlist xt_idxsection_new(Xt_idx start, unsigned num_dimensions,
                             unsigned * global_size, unsigned * local_size,
                             unsigned * local_start) {

  struct xt_idxsection *idxsection;

  idxsection = xmalloc(1 * sizeof(*idxsection));

  idxsection->vtable             = &idxsection_vtable;
  idxsection->global_start_index = start;
  idxsection->ndim               = num_dimensions;
  idxsection->global_size        = xmalloc(num_dimensions * sizeof(*(idxsection->global_size)));
  idxsection->local_size         = xmalloc(num_dimensions * sizeof(*(idxsection->local_size)));
  idxsection->local_start        = xmalloc(num_dimensions * sizeof(*(idxsection->local_start)));
  idxsection->global_strides     = xmalloc(num_dimensions * sizeof(*(idxsection->global_strides)));
  idxsection->local_strides      = xmalloc(num_dimensions * sizeof(*(idxsection->local_strides)));


  unsigned i;

  for (i = 0; i < num_dimensions; ++i) {
    idxsection->global_size[i] = global_size[i];
    idxsection->local_size[i]  = local_size[i];
    idxsection->local_start[i] = local_start[i];
  }

  idxsection->local_start_index = start;
  idxsection->global_strides[num_dimensions-1] = 1;
  idxsection->local_strides[num_dimensions-1] = 1;

  for (i = num_dimensions - 2; i < num_dimensions; --i) {
    idxsection->global_strides[i] = idxsection->global_strides[i+1] * global_size[i+1];
    idxsection->local_strides[i] = idxsection->local_strides[i+1] * local_size[i+1];
  }

  for (i = num_dimensions - 1; i < num_dimensions; --i)
    idxsection->local_start_index += idxsection->global_strides[i] * local_start[i];

  return (Xt_idxlist)idxsection;
}

static void
idxsection_delete(xt_idxlist data) {

  struct xt_idxsection *section;

  if (data == NULL) return;

  section = (struct xt_idxsection*)data;

  free(section->local_size);
  free(section->local_start);
  free(section->global_size);
  free(section->global_strides);
  free(section->local_strides);
  free(section);
}

static size_t
idxsection_get_pack_size(xt_idxlist data, MPI_Comm comm) {

  struct xt_idxsection *section = (struct xt_idxsection*)data;

  int size_int_type, size_unsigned, size_int;

  mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);
  mpi_err_handler(MPI_Pack_size(2, Xt_idx_dt, comm, &size_int), comm);
  mpi_err_handler(MPI_Pack_size(1 + section->ndim * 5, MPI_UNSIGNED,
                                comm, &size_unsigned), comm);

  return size_int_type + size_unsigned + size_int;
}

static void
idxsection_pack(xt_idxlist data, void *buffer, int buffer_size,
                int *position, MPI_Comm comm) {

  assert(data);
  struct xt_idxsection *section = (struct xt_idxsection*)data;
  int type = SECTION;

  mpi_err_handler(MPI_Pack(&(type), 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(&(section->global_start_index), 1,
                           Xt_idx_dt, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(&(section->local_start_index), 1,
                           Xt_idx_dt, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(&(section->ndim), 1, MPI_UNSIGNED, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(section->global_size, section->ndim, MPI_UNSIGNED, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(section->local_size, section->ndim, MPI_UNSIGNED, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(section->local_start, section->ndim, MPI_UNSIGNED, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(section->global_strides, section->ndim, MPI_UNSIGNED, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(section->local_strides, section->ndim, MPI_UNSIGNED, buffer,
                           buffer_size, position, comm), comm);
}

xt_idxlist xt_idxsection_unpack(void *buffer, int buffer_size, int *position,
                                MPI_Comm comm) {

  struct xt_idxsection *section;

  section = xmalloc(1 * sizeof(*section));

  section->vtable = &idxsection_vtable;

  mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                             &(section->global_start_index), 1, Xt_idx_dt,
                             comm), comm);
  mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                             &(section->local_start_index), 1, Xt_idx_dt,
                             comm), comm);
  mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                             &(section->ndim), 1, MPI_UNSIGNED,
                             comm), comm);
  assert(section->ndim > 0);

  section->global_size = xmalloc(section->ndim * sizeof(*(section->global_size)));
  section->local_size = xmalloc(section->ndim * sizeof(*(section->global_size)));
  section->local_start = xmalloc(section->ndim * sizeof(*(section->global_size)));
  section->global_strides = xmalloc(section->ndim * sizeof(*(section->global_strides)));
  section->local_strides = xmalloc(section->ndim * sizeof(*(section->local_strides)));

  mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                          section->global_size, section->ndim,
                          MPI_UNSIGNED, comm),comm);
  mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                          section->local_size, section->ndim,
                          MPI_UNSIGNED, comm),comm);
  mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                          section->local_start, section->ndim,
                          MPI_UNSIGNED, comm),comm);
  mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                          section->global_strides, section->ndim,
                          MPI_UNSIGNED, comm),comm);
  mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                          section->local_strides, section->ndim,
                          MPI_UNSIGNED, comm),comm);

  return (xt_idxlist)section;
}

static Xt_idxlist
idxsection_get_intersection(xt_idxlist idxlist_src, xt_idxlist idxlist_dst) {

  if (idxlist_dst->vtable != &idxsection_vtable) return NULL;

  // both lists are index section:

  struct xt_idxsection *idxsection_src, *idxsection_dst;

  idxsection_src = (struct xt_idxsection*)idxlist_src;
  idxsection_dst = (struct xt_idxsection*)idxlist_dst;

  if (idxsection_src->ndim != idxsection_dst->ndim ||
      idxsection_src->global_start_index != idxsection_dst->global_start_index)
    return NULL;

  unsigned i;

  // the size of first global dimension is irrelevant,
  // the others have to be identically
  for (i = 1; i < idxsection_src->ndim; ++i)
    if (idxsection_src->global_size[i] != idxsection_dst->global_size[i])
      return NULL;

  unsigned *local_start, *local_size;

  local_start = xmalloc(idxsection_src->ndim * sizeof(*local_start));
  local_size = xmalloc(idxsection_src->ndim * sizeof(*local_size));

  for (i = 0; i < idxsection_src->ndim; ++i) {
  
    unsigned src_start, src_end, dst_start, dst_end, local_end;

    src_start = idxsection_src->local_start[i];
    dst_start = idxsection_dst->local_start[i];
    src_end = idxsection_src->local_start[i] + idxsection_src->local_size[i];
    dst_end = idxsection_dst->local_start[i] + idxsection_dst->local_size[i];
  
    local_start[i] = (src_start > dst_start)?src_start:dst_start;
    local_end = (src_end > dst_end)?dst_end:src_end;

    if (local_end <= local_start[i]) {
      free(local_start);
      free(local_size);
      return xt_idxvec_new(NULL, 0);
    }

    local_size[i] = local_end - local_start[i];
  }

  struct xt_idxlist *intersection;

  intersection = xt_idxsection_new(idxsection_src->global_start_index,
                                   idxsection_src->ndim,
                                   idxsection_src->global_size,
                                   local_size, local_start);

  free(local_start);
  free(local_size);

  return intersection;
}

static xt_idxlist
idxsection_copy(xt_idxlist idxlist) {

  struct xt_idxsection *section;

  section = (struct xt_idxsection*)idxlist;

  return xt_idxsection_new(section->global_start_index, section->ndim,
                           section->global_size, section->local_size,
                           section->local_start);
}

static size_t
idxsection_get_num_indices(xt_idxlist idxlist) {

  struct xt_idxsection * section;

  section = (struct xt_idxsection *)idxlist;

  unsigned i;
  size_t size = 1;
  
  for (i = 0; i < section->ndim; ++i)
    size *= section->local_size[i];

  return size;
}

static void
idxsection_get_indices(xt_idxlist idxlist, Xt_idx *indices) {

  struct xt_idxsection *section;

  section = (struct xt_idxsection*)idxlist;

  unsigned * curr_local_position;

  curr_local_position = xcalloc(section->ndim, sizeof(*curr_local_position));

  unsigned i, j, num_indices;

  num_indices = idxsection_get_num_indices(idxlist);
  
  for (i = 0; i < num_indices; ++i) {

    for (j = 0, *indices = section->local_start_index; j < section->ndim; ++j)
      *indices += curr_local_position[j] * section->global_strides[j];

    ++indices;

    for (j = section->ndim - 1; j < section->ndim; --j)
      if (curr_local_position[j] < section->local_size[j] - 1) {
        curr_local_position[j]++;
        break;
      } else
        curr_local_position[j] = 0;
  }

  free(curr_local_position);
}

static void
idxsection_get_index_stripes(xt_idxlist idxlist, struct Xt_stripe ** stripes,
                             unsigned * num_stripes) {

  struct xt_idxsection *section;

  section = (struct xt_idxsection*)idxlist;

  unsigned i, j;

  if (section->local_size[section->ndim-1] == 0)
    *num_stripes = 0;
  else
    *num_stripes = 1;

  for (i = 0; i < section->ndim-1; ++i)
    *num_stripes *= section->local_size[i];

  if (num_stripes == 0) {

    *stripes = NULL;
    *num_stripes = 0;
    return;
  }

  *stripes = xmalloc(*num_stripes * sizeof(**stripes));

  unsigned * curr_local_position;

  curr_local_position = xcalloc(section->ndim-1, sizeof(*curr_local_position));

  for (i = 0; i < *num_stripes; ++i) {

    (*stripes)[i].start  = section->local_start_index;
    (*stripes)[i].extent = section->local_size[section->ndim-1];
    (*stripes)[i].stride = 1;

    for (j = 0; j < section->ndim - 1; ++j)
      (*stripes)[i].start += curr_local_position[j] * section->global_strides[j];

    for (j = section->ndim - 2; j < section->ndim; --j)
      if (curr_local_position[j] < section->local_size[j] - 1) {
        curr_local_position[j]++;
        break;
      } else
        curr_local_position[j] = 0;
  }

  free(curr_local_position);
}

static int
idxsection_get_index(Xt_idxlist idxlist, int position, Xt_idx * index) {

  struct xt_idxsection *section;

  section = (struct xt_idxsection*)idxlist;

  if (position < 0) return 1;
  
  Xt_idx temp_index;

  temp_index = section->local_start_index;

  unsigned dim;
  unsigned curr_local_position;

  for (dim = 0; dim < section->ndim; ++dim) {

    curr_local_position = position / section->local_strides[dim];

    if (curr_local_position >= section->local_size[dim]) return 1;

    temp_index += curr_local_position * section->global_strides[dim];

    position %= section->local_strides[dim];
  }

  *index = temp_index;

  return 0;
}

static int
idxsection_get_position_of_index(Xt_idxlist idxlist, Xt_idx index, int * position) {

  struct xt_idxsection *section;

  section = (struct xt_idxsection*)idxlist;

  index -= section->global_start_index;

  unsigned i;
  int temp_position = 0;
  unsigned curr_global_position, curr_local_position;

  for (i = 0; i < section->ndim; ++i) {

    curr_global_position = index / section->global_strides[i];
    index %= section->global_strides[i];

    if (curr_global_position < section->local_start[i]) return 1;
    
    curr_local_position = curr_global_position - section->local_start[i];

    if (curr_local_position >= section->local_size[i]) return 1;

    temp_position += curr_local_position * section->local_strides[i];
  }

  *position = temp_position;

  return 0;
}

static int
idxsection_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                 int * position, int offset) {

  struct xt_idxsection *section;

  section = (struct xt_idxsection*)idxlist;

  int temp_position;

  if (idxsection_get_position_of_index(idxlist, index, &temp_position))
    return 1;

  if (temp_position < offset)
    return 1;

  *position = temp_position;

  return 0;
}
