#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "mpi.h"

#include "yaxt.h"
#include "xt_idxlist.h"
#include "xt_idxlist_collection.h"
#include "xt_mpi.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "ensure_array_size.h"

static void
idxlist_collection_delete(Xt_idxlist data);

static size_t
idxlist_collection_get_pack_size(Xt_idxlist data, MPI_Comm comm);

static void
idxlist_collection_pack(Xt_idxlist data, void *buffer, int buffer_size,
                        int *position, MPI_Comm comm);

static Xt_idxlist
idxlist_collection_get_intersection(Xt_idxlist idxlist_src,
                                    Xt_idxlist idxlist_dst);

static Xt_idxlist
idxlist_collection_copy(Xt_idxlist idxlist);

static Xt_count
idxlist_collection_get_num_indices(Xt_idxlist idxlist);

static void
idxlist_collection_get_indices(Xt_idxlist idxlist, Xt_idx *indices);

static void
idxlist_collection_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                                     Xt_count * num_stripes);

static int
idxlist_collection_get_index_at_position(Xt_idxlist idxlist, int position, Xt_idx * index);

static int
idxlist_collection_get_position_of_index(Xt_idxlist idxlist, Xt_idx index, int * position);

static int
idxlist_collection_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                             int * position, int offset);

static Xt_idx
idxlist_collection_get_min_index(Xt_idxlist idxlist);

static Xt_idx
idxlist_collection_get_max_index(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxlist_collection_vtable = {
  .delete                      = idxlist_collection_delete,
  .get_pack_size               = idxlist_collection_get_pack_size,
  .pack                        = idxlist_collection_pack,
  .get_intersection            = idxlist_collection_get_intersection,
  .copy                        = idxlist_collection_copy,
  .get_num_indices             = idxlist_collection_get_num_indices,
  .get_indices                 = idxlist_collection_get_indices,
  .get_index_stripes           = idxlist_collection_get_index_stripes,
  .get_index_at_position       = idxlist_collection_get_index_at_position,
  .get_indices_at_positions    = NULL,
  .get_position_of_index       = idxlist_collection_get_position_of_index,
  .get_positions_of_indices    = NULL,
  .get_position_of_index_off   = idxlist_collection_get_position_of_index_off,
  .get_positions_of_indices_off = NULL,
  .get_min_index               = idxlist_collection_get_min_index,
  .get_max_index               = idxlist_collection_get_max_index
};

struct xt_idxlist_collection {

   const struct xt_idxlist_vtable *vtable;

   Xt_idxlist * idxlists;
   unsigned num_idxlists;
};

Xt_idxlist xt_idxlist_collection_new(Xt_idxlist *idxlists, unsigned num_idxlists) {

   struct xt_idxlist_collection *collectionlist;

   collectionlist = xmalloc(1 * sizeof(*collectionlist));

   collectionlist->vtable = &idxlist_collection_vtable;
   collectionlist->num_idxlists = num_idxlists;
   collectionlist->idxlists = xmalloc(num_idxlists * sizeof(*(collectionlist->idxlists)));

   unsigned i;

   for (i = 0; i < num_idxlists; ++i)
      collectionlist->idxlists[i] = xt_idxlist_copy(idxlists[i]);

   return (Xt_idxlist)collectionlist;
}

void idxlist_collection_delete(Xt_idxlist data) {

   struct xt_idxlist_collection *collectionlist = (struct xt_idxlist_collection*)data;

   unsigned i;

   for (i = 0; i < collectionlist->num_idxlists; ++i)
      xt_idxlist_delete(collectionlist->idxlists[i]);

   free(collectionlist->idxlists);
   free(collectionlist);
}

size_t idxlist_collection_get_pack_size(Xt_idxlist data, MPI_Comm comm) {

   struct xt_idxlist_collection *collectionlist = (struct xt_idxlist_collection*)data;

   int size_int, size_unsigned;
   size_t size_idxlists = 0;
   unsigned i;

   mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_int), comm);
   mpi_err_handler(MPI_Pack_size(1, MPI_UNSIGNED, comm, &size_unsigned), comm);

   for (i = 0; i < collectionlist->num_idxlists; ++i)
      size_idxlists += xt_idxlist_get_pack_size(collectionlist->idxlists[i], comm);

   return size_int + size_unsigned + size_idxlists;
}

void idxlist_collection_pack(Xt_idxlist data, void *buffer, int buffer_size, int *position,
                             MPI_Comm comm) {

   struct xt_idxlist_collection *collectionlist = (struct xt_idxlist_collection*)data;
   int type = COLLECTION;

   mpi_err_handler(MPI_Pack(&(type), 1, MPI_INT, buffer,
                            buffer_size, position, comm), comm);
   mpi_err_handler(MPI_Pack(&(collectionlist->num_idxlists), 1, MPI_UNSIGNED, buffer,
                            buffer_size, position, comm), comm);

   unsigned i;

   for (i = 0; i < collectionlist->num_idxlists; ++i)
      xt_idxlist_pack(collectionlist->idxlists[i], buffer, buffer_size, position, comm);
}

Xt_idxlist xt_idxlist_collection_unpack(void *buffer, int buffer_size, int *position,
                                        MPI_Comm comm) {

   struct xt_idxlist_collection *collectionlist;

   collectionlist = xmalloc(1 * sizeof(*collectionlist));

   collectionlist->vtable = &idxlist_collection_vtable;

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              &(collectionlist->num_idxlists), 1, MPI_UNSIGNED, comm),
                   comm);
   
   collectionlist->idxlists = xmalloc(collectionlist->num_idxlists *
                                     sizeof(*(collectionlist->idxlists)));

   unsigned i;

   for (i = 0; i < collectionlist->num_idxlists; ++i)
      collectionlist->idxlists[i] = xt_idxlist_unpack(buffer, buffer_size,
                                                      position, comm);

   return (Xt_idxlist)collectionlist;
}

Xt_idxlist idxlist_collection_get_intersection(Xt_idxlist idxlist_src,
                                               Xt_idxlist idxlist_dst) {

   return NULL;
}

Xt_idxlist idxlist_collection_copy(Xt_idxlist idxlist) {

   struct xt_idxlist_collection *collectionlist = (struct xt_idxlist_collection*)idxlist;

   return xt_idxlist_collection_new(collectionlist->idxlists, collectionlist->num_idxlists);
}

static Xt_count
idxlist_collection_get_num_indices(Xt_idxlist idxlist) {

   struct xt_idxlist_collection *collectionlist = (struct xt_idxlist_collection*)idxlist;

   unsigned i;
   Xt_count num_indices = 0;

   for (i = 0; i < collectionlist->num_idxlists; ++i)
      num_indices += xt_idxlist_get_num_indices(collectionlist->idxlists[i]);

   return num_indices;
}

static void
idxlist_collection_get_indices(Xt_idxlist idxlist, Xt_idx *indices) {

   struct xt_idxlist_collection *collectionlist = (struct xt_idxlist_collection*)idxlist;

   unsigned i;
   Xt_count offlist = 0;

   for (i = 0; i < collectionlist->num_idxlists; ++i) {

      xt_idxlist_get_indices(collectionlist->idxlists[i], indices+offlist);
      offlist += xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
   }
}

static void
idxlist_collection_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe ** stripes,
                                     Xt_count * num_stripes) {

  struct xt_idxlist_collection *collectionlist = (struct xt_idxlist_collection*)idxlist;

  struct Xt_stripe * temp_stripes = NULL;
  unsigned temp_stripes_array_size = 0;
  Xt_count num_temp_stripes = 0;

  unsigned i;

  if (collectionlist->num_idxlists > 0) {
    xt_idxlist_get_index_stripes(collectionlist->idxlists[0], &temp_stripes,
                                 &num_temp_stripes);
    temp_stripes_array_size = num_temp_stripes;
  }

  for (i = 1; i < collectionlist->num_idxlists; ++i) {

    struct Xt_stripe * curr_stripes;
    Xt_count curr_num_stripes;

    xt_idxlist_get_index_stripes(collectionlist->idxlists[i], &curr_stripes,
                                 &curr_num_stripes);

    ENSURE_ARRAY_SIZE(temp_stripes, temp_stripes_array_size,
                      num_temp_stripes + curr_num_stripes);

    memcpy(temp_stripes + num_temp_stripes, curr_stripes,
           curr_num_stripes * sizeof(*curr_stripes));

    free(curr_stripes);

    num_temp_stripes += curr_num_stripes;
  }

  *stripes = xrealloc(temp_stripes, num_temp_stripes * sizeof(*temp_stripes));
  *num_stripes = num_temp_stripes;
}

static int
idxlist_collection_get_index_at_position(Xt_idxlist idxlist, int position, Xt_idx * index) {

  struct xt_idxlist_collection *collectionlist = (struct xt_idxlist_collection*)idxlist;

  for (Xt_count i = 0; i < collectionlist->num_idxlists; ++i) {
    Xt_count n = xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
    if (position >= n)
      position -= n;
    else {
      return xt_idxlist_get_index_at_position(collectionlist->idxlists[i], position, index);
    }
  }
  return 1;

}

static int
idxlist_collection_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                             int * position, int offset) {

  struct xt_idxlist_collection *collectionlist = (struct xt_idxlist_collection*)idxlist;

  unsigned i;
  Xt_count curr_num_indices = 0;

  int idxlist_offsets = 0;

  assert(offset >= 0);

  i = 0;

  do {
      idxlist_offsets += curr_num_indices;
      curr_num_indices = xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
  } while (idxlist_offsets + curr_num_indices <= offset &&
           ++i < collectionlist->num_idxlists);

  offset -= idxlist_offsets;

  for (;i < collectionlist->num_idxlists; ++i)
    if (!xt_idxlist_get_position_of_index_off(collectionlist->idxlists[i], index,
                                              position, offset)) {
      *position += idxlist_offsets;
      return 0;
    } else {
      idxlist_offsets += xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
      offset = 0;
    }

  return 1;
}

static int
idxlist_collection_get_position_of_index(Xt_idxlist idxlist, Xt_idx index,
                                         int * position) {

  return idxlist_collection_get_position_of_index_off(idxlist, index, position, 0);
}

static Xt_idx
idxlist_collection_get_min_index(Xt_idxlist idxlist) {

  struct xt_idxlist_collection *collectionlist = (struct xt_idxlist_collection*)idxlist;

  Xt_idx tmp_min, min = 0;

  if (collectionlist->num_idxlists == 0) die("idxlist_collection_get_min_index empty collection");

  min = xt_idxlist_get_min_index(collectionlist->idxlists[0]);

  for (unsigned i = 1; i < collectionlist->num_idxlists; ++i)
    if ((tmp_min = xt_idxlist_get_min_index(collectionlist->idxlists[i])) > min)
      min = tmp_min;

  return min;
}

static Xt_idx
idxlist_collection_get_max_index(Xt_idxlist idxlist) {

  struct xt_idxlist_collection *collectionlist = (struct xt_idxlist_collection*)idxlist;

  Xt_idx tmp_max, max = 0;

  if (collectionlist->num_idxlists == 0) die("idxlist_collection_get_max_index empty collection");

  max = xt_idxlist_get_max_index(collectionlist->idxlists[0]);

  for (unsigned i = 1; i < collectionlist->num_idxlists; ++i)
    if ((tmp_max = xt_idxlist_get_max_index(collectionlist->idxlists[i])) > max)
      max = tmp_max;

  return max;
}
