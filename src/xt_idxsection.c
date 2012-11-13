/**
 * @file xt_idxsection.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xt/xt_idxlist.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxsection.h"
#include "xt/xt_mpi.h"
#include "xt/mergesort.h"
#include "xt/quicksort.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"


static void
idxsection_delete(Xt_idxlist data);

static size_t
idxsection_get_pack_size(Xt_idxlist data, MPI_Comm comm);

static void
idxsection_pack(Xt_idxlist data, void *buffer, int buffer_size,
                int *position, MPI_Comm comm);

static Xt_idxlist
idxsection_get_intersection(Xt_idxlist idxlist_a, Xt_idxlist idxlist_b);

static Xt_idxlist
idxsection_copy(Xt_idxlist idxlist);

static Xt_count
idxsection_get_num_indices(Xt_idxlist idxlist);

static void
idxsection_get_indices(Xt_idxlist idxlist, Xt_int *indices);

static const Xt_int *
idxsection_get_indices_const(Xt_idxlist idxlist);

static void
idxsection_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                             Xt_count * num_stripes);

static int
idxsection_get_index_at_position(Xt_idxlist idxlist, int position,
                                 Xt_int * index);

static int
idxsection_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                 int * position);

static int
idxsection_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                     int * position, int offset);
static int
idxsection_get_positions_of_indices(Xt_idxlist body_idxlist, Xt_int const * selection_idx,
                                    Xt_count num_selection, int * positions,
                                    int single_match_only);

static Xt_int
idxsection_get_min_index(Xt_idxlist idxlist);

static Xt_int
idxsection_get_max_index(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxsection_vtable = {
              .delete                      = idxsection_delete,
              .get_pack_size               = idxsection_get_pack_size,
              .pack                        = idxsection_pack,
              .get_intersection            = idxsection_get_intersection,
              .copy                        = idxsection_copy,
              .get_num_indices             = idxsection_get_num_indices,
              .get_indices                 = idxsection_get_indices,
              .get_indices_const           = idxsection_get_indices_const,
              .get_index_stripes           = idxsection_get_index_stripes,
              .get_index_at_position       = idxsection_get_index_at_position,
              .get_indices_at_positions    = NULL,
              .get_position_of_index       = idxsection_get_position_of_index,
              .get_positions_of_indices    = idxsection_get_positions_of_indices,
              .get_position_of_index_off   = idxsection_get_position_of_index_off,
              .get_positions_of_indices_off = NULL,
              .get_min_index               = idxsection_get_min_index,
              .get_max_index               = idxsection_get_max_index};

/* descriptor for per-dimension extent and stride */
struct dim_desc
{
  Xt_int global_size, global_stride, local_start, local_size, local_stride;
};

struct xt_idxsection {

  const struct xt_idxlist_vtable *vtable;

  Xt_int global_start_index;
  Xt_int local_start_index;
  unsigned ndim;
  struct dim_desc *dims;

  Xt_count index_num_cache;
  Xt_int *index_array_cache;

  Xt_int min_index_cache;
  Xt_int max_index_cache;

};

Xt_idxlist xt_idxsection_new(Xt_int start, unsigned num_dimensions,
                             Xt_int global_size[num_dimensions],
                             Xt_int local_size[num_dimensions],
                             Xt_int local_start[num_dimensions]) {

  struct xt_idxsection *idxsection;

  idxsection = xmalloc(1 * sizeof(*idxsection));

  idxsection->vtable             = &idxsection_vtable;
  idxsection->global_start_index = start;
  idxsection->ndim               = num_dimensions;
  idxsection->dims = xmalloc(num_dimensions * sizeof (idxsection->dims[0]));

  idxsection->index_array_cache  = NULL;

  unsigned i;

  if (num_dimensions<=0)
    idxsection->index_num_cache    = 0;
  else
    idxsection->index_num_cache    = 1;

  for (i = 0; i < num_dimensions; ++i) {
    idxsection->dims[i].global_size = global_size[i];
    idxsection->dims[i].local_size  = local_size[i];
    idxsection->dims[i].local_start = local_start[i];
    idxsection->index_num_cache *= local_size[i];
  }

  idxsection->local_start_index = start;
  idxsection->dims[num_dimensions - 1].global_stride = 1;
  idxsection->dims[num_dimensions - 1].local_stride = 1;

  for (i = num_dimensions - 2; i < num_dimensions; --i) {
    idxsection->dims[i].global_stride
      = idxsection->dims[i+1].global_stride * global_size[i + 1];
    idxsection->dims[i].local_stride
      = idxsection->dims[i + 1].local_stride * local_size[i + 1];
  }

  for (i = num_dimensions - 1; i < num_dimensions; --i)
    idxsection->local_start_index
      += idxsection->dims[i].global_stride * local_start[i];

  idxsection_get_index_at_position((Xt_idxlist)idxsection, 0,
                                   &(idxsection->min_index_cache));
  idxsection_get_index_at_position(
    (Xt_idxlist)idxsection,
    idxsection_get_num_indices((Xt_idxlist)idxsection) - 1,
    &(idxsection->max_index_cache));

  return (Xt_idxlist)idxsection;
}

Xt_idxlist
xt_idxsection_new2(Xt_int start, unsigned num_dimensions,
                   const struct dim_desc dims[num_dimensions])
{
}

static void
idxsection_delete(Xt_idxlist data) {

  struct xt_idxsection *section;

  if (data == NULL) return;

  section = (struct xt_idxsection*)data;

  free(section->dims);
  free(section->index_array_cache);
  free(section);
}

static size_t
idxsection_get_pack_size(Xt_idxlist data, MPI_Comm comm) {

  struct xt_idxsection *section = (struct xt_idxsection*)data;

  int size_int_type, size_ndim, size_dims, size_int;

  xt_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);
  xt_mpi_call(MPI_Pack_size(2, Xt_int_dt, comm, &size_int), comm);
  xt_mpi_call(MPI_Pack_size(1, MPI_UNSIGNED, comm, &size_ndim), comm);
  xt_mpi_call(MPI_Pack_size(section->ndim * 5, Xt_int_dt,
                            comm, &size_dims), comm);

  return size_int_type + size_ndim + size_dims + size_int;
}

static void
idxsection_pack(Xt_idxlist data, void *buffer, int buffer_size,
                int *position, MPI_Comm comm) {

  assert(data);
  struct xt_idxsection *section = (struct xt_idxsection*)data;
  int type = SECTION;

  xt_mpi_call(MPI_Pack(&(type), 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  xt_mpi_call(MPI_Pack(&(section->global_start_index), 1,
                           Xt_int_dt, buffer,
                           buffer_size, position, comm), comm);
  xt_mpi_call(MPI_Pack(&(section->local_start_index), 1,
                           Xt_int_dt, buffer,
                           buffer_size, position, comm), comm);
  xt_mpi_call(MPI_Pack(&(section->ndim), 1, MPI_UNSIGNED, buffer,
                           buffer_size, position, comm), comm);
  xt_mpi_call(MPI_Pack(section->dims, section->ndim * 5, Xt_int_dt,
                       buffer, buffer_size, position, comm), comm);
}

Xt_idxlist xt_idxsection_unpack(void *buffer, int buffer_size, int *position,
                                MPI_Comm comm) {

  struct xt_idxsection *section;

  section = xmalloc(1 * sizeof(*section));

  section->vtable = &idxsection_vtable;
  section->index_array_cache = NULL;
  section-> index_num_cache = -1; // == undef

  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                             &(section->global_start_index), 1, Xt_int_dt,
                             comm), comm);
  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                             &(section->local_start_index), 1, Xt_int_dt,
                             comm), comm);
  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                             &(section->ndim), 1, MPI_UNSIGNED,
                             comm), comm);
  assert(section->ndim > 0);

  section->dims = xmalloc(section->ndim * sizeof(section->dims[0]));

  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         section->dims, section->ndim * 5,
                         Xt_int_dt, comm),comm);

  // fill other cache values:
  idxsection_get_index_at_position((Xt_idxlist)section, 0, &(section->min_index_cache));
  idxsection_get_index_at_position((Xt_idxlist)section, idxsection_get_num_indices((Xt_idxlist)section)-1,
                                   &(section->max_index_cache));

  return (Xt_idxlist)section;
}

static Xt_idxlist
idxsection_get_intersection_with_other_idxlist(struct xt_idxsection *src_idxsection,
                                               Xt_idxlist dst_idxlist) {

  // intersection between an idxsection and a general idxlist:
  //
  // performance picture:
  //  - src_idxsection is treated as too big for elemental transforms/access
  //  - dst_idxlist is considered to be small enough (subdomain like) for elemental usage

  Xt_count num_dst_idx = xt_idxlist_get_num_indices(dst_idxlist);

  Xt_int const* dst_idx = xt_idxlist_get_indices_const(dst_idxlist);
  int * pos;
  int single_match_only = 0;
  Xt_count i, j;

  for (i = 1; i < num_dst_idx; ++i)
    if (dst_idx[i] < dst_idx[i-1])
      break;

  Xt_int const * sorted_dst_idx;
  Xt_int * temp_dst_idx = NULL;
      
  // if the destination indices are not sorted
  if (num_dst_idx > 1 && i != num_dst_idx) {

    temp_dst_idx = xmalloc(num_dst_idx * sizeof(*temp_dst_idx));
    memcpy(temp_dst_idx, dst_idx, num_dst_idx * sizeof(*temp_dst_idx));

    xt_mergesort_index(temp_dst_idx, num_dst_idx, NULL, 0);

    sorted_dst_idx = temp_dst_idx;
  } else
    sorted_dst_idx = dst_idx;

  pos = xmalloc(num_dst_idx * sizeof(*pos));
  int num_unmatched = idxsection_get_positions_of_indices((Xt_idxlist)src_idxsection,
                                                          sorted_dst_idx,
                                                          num_dst_idx, pos,
                                                          single_match_only);
  Xt_count num_inter_idx = num_dst_idx - num_unmatched;
  Xt_int *intersection = xmalloc(num_inter_idx * sizeof(*intersection));

  for(i = 0, j = 0; i < num_dst_idx && j < num_inter_idx; i++)
    if (pos[i] >= 0)
      intersection[j++] = sorted_dst_idx[i];

  free(temp_dst_idx);
  free(pos);

  //! \todo use xt_idxvec_pinned_new here to avoid copy
  Xt_idxlist result = xt_idxvec_new(intersection, num_inter_idx);

  free(intersection);
  return result;
  // return xt_idxvec_pinned_new(intersection, num_inter_idx);
}

static Xt_idxlist
idxsection_get_intersection(Xt_idxlist idxlist_src, Xt_idxlist idxlist_dst) {

  if (idxlist_dst->vtable != &idxsection_vtable)
    return idxsection_get_intersection_with_other_idxlist((struct xt_idxsection *)idxlist_src,
                                                          idxlist_dst);

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
    if (idxsection_src->dims[i].global_size
        != idxsection_dst->dims[i].global_size)
      return idxsection_get_intersection_with_other_idxlist(
        (struct xt_idxsection *)idxlist_src, idxlist_dst);

  Xt_int *local_start, *local_size, *global_size;

  local_start = xmalloc(idxsection_src->ndim * sizeof(*local_start));
  local_size = xmalloc(idxsection_src->ndim * sizeof(*local_size));
  global_size = xmalloc(idxsection_src->ndim * sizeof(*global_size));

  for (i = 0; i < idxsection_src->ndim; ++i) {

    Xt_int src_start, src_end, dst_start, dst_end, local_end;

    src_start = idxsection_src->dims[i].local_start;
    dst_start = idxsection_dst->dims[i].local_start;
    src_end = idxsection_src->dims[i].local_start
      + idxsection_src->dims[i].local_size;
    dst_end = idxsection_dst->dims[i].local_start
      + idxsection_dst->dims[i].local_size;

    local_start[i] = (src_start > dst_start)?src_start:dst_start;
    local_end = (src_end > dst_end)?dst_end:src_end;

    if (local_end <= local_start[i]) {
      free(global_size);
      free(local_size);
      free(local_start);
      return xt_idxvec_new(NULL, 0);
    }

    local_size[i] = local_end - local_start[i];
    global_size[i] = idxsection_src->dims[i].global_size;
  }

  struct Xt_idxlist *intersection;

  intersection = xt_idxsection_new(idxsection_src->global_start_index,
                                   idxsection_src->ndim,
                                   global_size, local_size, local_start);

  free(global_size);
  free(local_size);
  free(local_start);

  return intersection;
}

static Xt_idxlist
idxsection_copy(Xt_idxlist idxlist) {

  struct xt_idxsection *src = (struct xt_idxsection*)idxlist;
  struct xt_idxsection *idxsection = xmalloc(1 * sizeof(*idxsection));

  *idxsection = *src;

  unsigned num_dimensions = src->ndim;

  idxsection->dims = xmalloc(num_dimensions * sizeof (idxsection->dims[0]));
  idxsection->index_array_cache  = NULL;

  memcpy(idxsection->dims, src->dims, num_dimensions * sizeof (src->dims[0]));

  return (Xt_idxlist)idxsection;
}

static Xt_count
idxsection_get_num_indices(Xt_idxlist idxlist) {

  struct xt_idxsection * section;

  section = (struct xt_idxsection *)idxlist;
  if (section->index_num_cache >= 0) return section->index_num_cache;

  unsigned i;
  size_t size = 1;

  for (i = 0; i < section->ndim; ++i)
    size *= section->dims[i].local_size;

  section->index_num_cache = size;
  return size;
}


static Xt_count
idxsection_get_indices_any(Xt_int start_index, Xt_int *indices,
                           unsigned ndim, struct dim_desc dims[ndim])
{
  if (ndim == 1)
  {
    for (Xt_count i = 0; i < dims[0].local_size; ++i)
      indices[i] = start_index + i;
    return dims[0].local_size;
  }
  else
  {
    Xt_count indices_written = 0;
    assert(ndim > 1);
    for (Xt_int dim_ofs = 0; dim_ofs < dims[0].local_size; ++dim_ofs)
    {
      indices_written
        += idxsection_get_indices_any(start_index
                                      + dim_ofs * dims[0].global_stride,
                                      indices + indices_written,
                                      ndim - 1, dims + 1);
    }
    return indices_written;
  }
}

static void
idxsection_get_indices(Xt_idxlist idxlist, Xt_int *indices) {

  struct xt_idxsection *section;

  section = (struct xt_idxsection*)idxlist;

  Xt_count num_indices = idxsection_get_num_indices(idxlist);

  if (num_indices == 0) return;

  if (section->index_array_cache != NULL) {

    memcpy(indices, section->index_array_cache, num_indices * sizeof(*indices));
    return;
  } else
    section->index_array_cache = xmalloc(num_indices * sizeof(*(section->index_array_cache)));

  // if this routine was called by idxsection_get_indices_const
  if (indices == NULL)
    indices = section->index_array_cache;

  idxsection_get_indices_any(section->local_start_index, indices,
                             section->ndim, section->dims);

  if (section->index_array_cache != indices)
    memcpy(section->index_array_cache, indices, num_indices * sizeof(*indices));
}

static Xt_int const*
idxsection_get_indices_const(Xt_idxlist idxlist) {

  struct xt_idxsection *idxsection;

  idxsection = (struct xt_idxsection*)idxlist;

  if (idxsection->index_array_cache == NULL)
    idxsection_get_indices(idxlist, NULL);

  return idxsection->index_array_cache;
}


static void
idxsection_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                             Xt_count * num_stripes) {

  struct xt_idxsection *section;

  section = (struct xt_idxsection*)idxlist;

  unsigned i, j;

  if (section->dims[section->ndim-1].local_size == 0)
    *num_stripes = 0;
  else
    *num_stripes = 1;

  for (i = 0; i < section->ndim-1; ++i)
    *num_stripes *= section->dims[i].local_size;

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
    (*stripes)[i].nstrides = section->dims[section->ndim-1].local_size;
    (*stripes)[i].stride = 1;

    for (j = 0; j < section->ndim - 1; ++j)
      (*stripes)[i].start
        += curr_local_position[j] * section->dims[j].global_stride;

    for (j = section->ndim - 2; j < section->ndim; --j)
      if (curr_local_position[j] < section->dims[j].local_size - 1) {
        curr_local_position[j]++;
        break;
      } else
        curr_local_position[j] = 0;
  }

  free(curr_local_position);
}

static int
idxsection_get_index_at_position(Xt_idxlist idxlist, int position,
                                 Xt_int * index) {

  struct xt_idxsection *section;

  section = (struct xt_idxsection*)idxlist;

  if (position < 0) return 1;

  Xt_int temp_index;

  temp_index = section->local_start_index;

  unsigned dim;
  unsigned curr_local_position;

  for (dim = 0; dim < section->ndim; ++dim) {

    curr_local_position = position / section->dims[dim].local_stride;

    if (curr_local_position >= section->dims[dim].local_size) return 1;

    temp_index += curr_local_position * section->dims[dim].global_stride;

    position %= section->dims[dim].local_stride;
  }

  *index = temp_index;

  return 0;
}

static int
idxsection_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                 int * position) {

  struct xt_idxsection *section = (struct xt_idxsection*)idxlist;
  *position = -1;

  if (index < section->min_index_cache || index > section->max_index_cache)
    return 1;

  index -= section->global_start_index;

  unsigned i;
  int temp_position = 0;
  unsigned curr_global_position, curr_local_position;

  for (i = 0; i < section->ndim; ++i) {

    curr_global_position = index / section->dims[i].global_stride;
    index %= section->dims[i].global_stride;

    if (curr_global_position < section->dims[i].local_start) return 1;

    curr_local_position = curr_global_position - section->dims[i].local_start;

    if (curr_local_position >= section->dims[i].local_size) return 1;

    temp_position += curr_local_position * section->dims[i].local_stride;
  }

  *position = temp_position;

  return 0;
}

int idxsection_get_positions_of_indices_v1(Xt_idxlist body_idxlist, Xt_int const * selection_idx,
                                           Xt_count num_selection, int * positions,
                                           int single_match_only) {

  if (num_selection < 1) return 0; //nothing to do

  if (num_selection == 1)
    return idxsection_get_position_of_index(body_idxlist, *selection_idx, positions);

  int num_unmatched = 0;

  if (!single_match_only) {
    // this is the easy case, we don't care about multiple uses of the same position
    for (Xt_count i = 0; i < num_selection; ++i)
      num_unmatched += idxsection_get_position_of_index(body_idxlist, selection_idx[i],
                                                        &positions[i]);
    return num_unmatched;
  }

  Xt_count i;

  for (i = 1; i < num_selection; ++i)
    if (selection_idx[i] < selection_idx[i-1])
      break;

  // indices are sorted
  if (i == num_selection) {

    Xt_int prev_index = selection_idx[0] - 1; // we need an index that is
                                              // different from the current one
  
    for (i = 0; i < num_selection; i++) {

      Xt_int curr_index = selection_idx[i];

      if (prev_index != curr_index) {

        num_unmatched += idxsection_get_position_of_index(body_idxlist, curr_index,
                                                          positions + i);
        prev_index = curr_index;

      } else {
        // for an idxsection there is a unique map from indices to positions,
        // we got the same index again, so there is no match left:
        positions[i] = -1;
        num_unmatched++;
      }
    }
  // indices are not sorted
  } else {
    // the remaining (single_match_only) case follows:
    idxpos_type *v;

    v = xmalloc(num_selection * sizeof(*v) );
    for (int i=0; i<num_selection; i++) {
      v[i].idx = selection_idx[i];
      v[i].pos = i;
    }
    xt_mergesort_idxpos(v, num_selection);
    Xt_int last_jx = v[0].idx - 1; // any index that is not equal v[0].idx will do
    for (i=0; i<num_selection; i++) {
      int j = v[i].pos;
      Xt_int jx = v[i].idx;
      if (jx != last_jx) {
        num_unmatched += idxsection_get_position_of_index(body_idxlist, jx, &positions[j]);
        last_jx = jx;
      } else {
        // for an idxsection there is a unique map from indices to positions,
        // we got the same index again, so there is no match left:
        positions[j] = -1;
        num_unmatched++;
      }
    }
    free(v);
  }

  return num_unmatched;
}

static
int idxsection_get_positions_of_indices_v2(Xt_idxlist body_idxlist,
                                           Xt_int const * selection_idx,
                                           Xt_count num_selection, int * positions,
                                           int single_match_only) {

  if (num_selection < 1) return 0; //nothing to do

  if (num_selection == 1)
    return idxsection_get_position_of_index(body_idxlist, *selection_idx, positions);

  Xt_count i, j;
  Xt_int * temp_selection_idx = NULL;
  Xt_int const * sorted_selection_idx;
  int * selection_pos = NULL;

  for (i = 1; i < num_selection; ++i)
    if (selection_idx[i] < selection_idx[i-1])
      break;

  if (i != num_selection) {

    temp_selection_idx = xmalloc(num_selection * sizeof(*temp_selection_idx));
    memcpy(temp_selection_idx, selection_idx,
           num_selection * sizeof(*temp_selection_idx));
    selection_pos = xmalloc(num_selection * sizeof(*selection_pos));

    xt_quicksort_index(temp_selection_idx, num_selection, selection_pos, 1);

    sorted_selection_idx = temp_selection_idx;

  } else
    sorted_selection_idx = selection_idx;

  Xt_int const * body_indices = idxsection_get_indices_const(body_idxlist);
  Xt_count num_body_indices = idxsection_get_num_indices(body_idxlist);

  if (!single_match_only) {

    if (selection_pos == NULL) {
      for (i = 0, j = 0; i < num_selection && j < num_body_indices; ++i) {

        while(j < num_body_indices && body_indices[j] < sorted_selection_idx[i]) ++j;

        if (j >= num_body_indices) break;

        positions[i] = (body_indices[j] == sorted_selection_idx[i])?j:-1;
      }
    } else {
      for (i = 0, j = 0; i < num_selection && j < num_body_indices; ++i) {

        while(j < num_body_indices && body_indices[j] < sorted_selection_idx[i]) ++j;

        if (j >= num_body_indices) break;

        positions[selection_pos[i]] = (body_indices[j] == sorted_selection_idx[i])?j:-1;
      }
    }
  } else {

    Xt_int last_idx = sorted_selection_idx[0] - 1;

    if (selection_pos == NULL) {
      for (i = 0, j = 0; i < num_selection && j < num_body_indices; ++i) {

        while(j < num_body_indices && body_indices[j] < sorted_selection_idx[i]) ++j;

        if (j >= num_body_indices) break;

        positions[i] = ((last_idx == sorted_selection_idx[i]) ||
                        (body_indices[j] != sorted_selection_idx[i]))?-1:j;

        last_idx = sorted_selection_idx[i];
      }
    } else {
      for (i = 0, j = 0; i < num_selection && j < num_body_indices; ++i) {

        while(j < num_body_indices && body_indices[j] < sorted_selection_idx[i]) ++j;

        if (j >= num_body_indices) break;

        positions[selection_pos[i]] = ((last_idx == sorted_selection_idx[i]) ||
                                       (body_indices[j] != sorted_selection_idx[i]))?-1:j;

        last_idx = sorted_selection_idx[i];
      }
    }
  }

  if (selection_pos == NULL)
    for (; i < num_selection; ++i)
      positions[i] = -1;
  else
    for (; i < num_selection; ++i)
      positions[selection_pos[i]] = -1;

  free(temp_selection_idx);
  free(selection_pos);

  int num_unmatched = 0;

  for (Xt_count i = 0; i < num_selection; ++i)
    if (positions[i] == -1)
      ++num_unmatched;

  return num_unmatched;
}

static
int idxsection_get_positions_of_indices(Xt_idxlist body_idxlist,
                                        Xt_int const * selection_idx,
                                        Xt_count num_selection, int * positions,
                                        int single_match_only) {

  struct xt_idxsection *section;
  section = (struct xt_idxsection*)body_idxlist;

  Xt_count num_section_indices = idxsection_get_num_indices(body_idxlist);

  if ((section->index_array_cache != NULL) ||
      ((num_section_indices * sizeof(Xt_int) <= 128 * 1024 * 1024) &&
       (num_section_indices <= 1000 * num_selection)))
    return idxsection_get_positions_of_indices_v2(body_idxlist, selection_idx,
                                                  num_selection, positions,
                                                  single_match_only);
  else
    return idxsection_get_positions_of_indices_v1(body_idxlist, selection_idx,
                                                  num_selection, positions,
                                                  single_match_only);
}

static int
idxsection_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                 int * position, int offset) {

  struct xt_idxsection *section;

  section = (struct xt_idxsection*)idxlist;

  int temp_position;
  // we make use of the uniqueness of the index-to-position relation:
  if (idxsection_get_position_of_index(idxlist, index, &temp_position))
    return 1;

  if (temp_position < offset)
    return 1;

  *position = temp_position;

  return 0;
}

static Xt_int
idxsection_get_min_index(Xt_idxlist idxlist) {

  struct xt_idxsection *section = (struct xt_idxsection*)idxlist;
  return section->min_index_cache;
}

static Xt_int
idxsection_get_max_index(Xt_idxlist idxlist) {

  struct xt_idxsection *section = (struct xt_idxsection*)idxlist;
  return section->max_index_cache;
}
