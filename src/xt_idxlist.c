#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xt/xt_core.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_mpi.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"

void xt_idxlist_delete(Xt_idxlist idxlist) {

   idxlist->vtable->delete(idxlist);
}

size_t xt_idxlist_get_pack_size(Xt_idxlist idxlist,
                                MPI_Comm comm) {

   return idxlist->vtable->get_pack_size(idxlist, comm);
}

void xt_idxlist_pack(Xt_idxlist idxlist, void *buffer,
                     int buffer_size, int *position,
                     MPI_Comm comm) {

   idxlist->vtable->pack(idxlist, buffer, buffer_size,
                         position, comm);
}

Xt_idxlist xt_idxlist_unpack(void *buffer, int buffer_size,
                             int *position, MPI_Comm comm) {

   int type;

   // unpack the type of the index list to be unpacked
   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                   &type, 1, MPI_INT, comm), comm);

   unsigned i;
   extern struct xt_idxlist_unpack_routine_table xt_idxlist_unpack_routines[];
   extern unsigned xt_num_unpack_routines;

   for (i = 0; i < xt_num_unpack_routines; ++i) {

      if (type == xt_idxlist_unpack_routines[i].type) {
         return xt_idxlist_unpack_routines[i].xt_idxlist_unpack(buffer,
                                                                buffer_size,
                                                                position, comm);
         break;
      }
   }

  Xt_abort(comm, "xt_idxlist_unpack: unknown index list type",
           __FILE__, __LINE__);

  return NULL;
}

Xt_idxlist xt_idxlist_get_intersection(Xt_idxlist idxlist_src,
                                       Xt_idxlist idxlist_dst) {

   Xt_idxlist intersection;

   intersection
     = idxlist_src->vtable->get_intersection(idxlist_src, idxlist_dst);

   // if the get_intersection routine was not able to compute the intersection
   if (intersection == NULL) {

      Xt_count num_indices_src, num_indices_dst;

      num_indices_src = xt_idxlist_get_num_indices(idxlist_src);
      num_indices_dst = xt_idxlist_get_num_indices(idxlist_dst);

      Xt_idx *indices_src, *indices_dst;
      indices_src = xmalloc(num_indices_src * sizeof (indices_src[0]));
      indices_dst = xmalloc(num_indices_dst * sizeof (indices_dst[0]));

      xt_idxlist_get_indices(idxlist_src, indices_src);
      xt_idxlist_get_indices(idxlist_dst, indices_dst);

      Xt_idxlist idxvec_src, idxvec_dst;

      idxvec_src = xt_idxvec_new(indices_src, num_indices_src);
      idxvec_dst = xt_idxvec_new(indices_dst, num_indices_dst);

      intersection = xt_idxlist_get_intersection(idxvec_src, idxvec_dst);

      xt_idxlist_delete(idxvec_src);
      xt_idxlist_delete(idxvec_dst);
      free(indices_dst);
      free(indices_src);
   }

   return intersection;
}

Xt_idxlist xt_idxlist_copy(Xt_idxlist idxlist) {

   return idxlist->vtable->copy(idxlist);
}

Xt_count xt_idxlist_get_num_indices(Xt_idxlist idxlist) {

   return idxlist->vtable->get_num_indices(idxlist);
}

void xt_idxlist_get_indices(Xt_idxlist idxlist, Xt_idx *indices) {

   idxlist->vtable->get_indices(idxlist, indices);
}


const Xt_idx *xt_idxlist_get_indices_const(Xt_idxlist idxlist) {
  if (idxlist->vtable->get_indices != NULL)
    return idxlist->vtable->get_indices_const (idxlist);

  die("xt_idxlist_get_indices_const: fatal error: "
      "get_indices_const not implemented");
  return NULL;
}


void xt_idxlist_get_index_stripes(Xt_idxlist idxlist,
                                  struct Xt_stripe ** stripes,
                                  Xt_count * num_stripes) {

  if (idxlist->vtable->get_index_stripes != NULL) {

    idxlist->vtable->get_index_stripes(idxlist, stripes,
                                       num_stripes);

  } else { // fall-back solution

    Xt_count num_indices;

    num_indices = xt_idxlist_get_num_indices(idxlist);

    Xt_idx *indices;
    indices = xmalloc(num_indices * sizeof (indices[0]));

    xt_idxlist_get_indices(idxlist, indices);

    xt_idxstripes_convert_to_stripes(indices, num_indices,
                                     stripes, num_stripes);

    free(indices);
  }
}

int xt_idxlist_get_index_at_position(Xt_idxlist idxlist, int position,
                                     Xt_idx * index) {

  return idxlist->vtable->get_index_at_position(idxlist, position, index);
}

int
xt_idxlist_get_indices_at_positions(Xt_idxlist idxlist, int *positions,
                                    Xt_count num_pos, Xt_idx *indices,
                                    Xt_idx undef_idx) {

  if ( idxlist->vtable->get_indices_at_positions ) {

    return idxlist->vtable->get_indices_at_positions(idxlist, positions,
                                                     num_pos, indices,
                                                     undef_idx);

  } else {

    // fallback solution using xt_idxlist_get_index_at_position:
    int undef_count = 0;
    for (int ip=0; ip<num_pos; ip++) {
      int p = positions[ip];
      if (xt_idxlist_get_index_at_position(idxlist, p, &indices[ip]) != 0) {
        indices[ip] = undef_idx;
        undef_count++;
      }
    }
    return undef_count;

  }
}

int xt_idxlist_get_position_of_index(Xt_idxlist idxlist, Xt_idx index,
                                     int * position) {

  return idxlist->vtable->get_position_of_index(idxlist, index, position);
}

int
xt_idxlist_get_positions_of_indices(Xt_idxlist idxlist, Xt_idx const * indices,
                                    Xt_count num_indices, int * positions,
                                    int single_match_only) {

  if (idxlist->vtable->get_positions_of_indices != NULL)
    return idxlist->vtable->get_positions_of_indices(idxlist, indices,
                                                     num_indices,
                                                     positions,
                                                     single_match_only);
  else {
    Xt_count num_tmp_indices;
    num_tmp_indices = xt_idxlist_get_num_indices(idxlist);
    Xt_idx *tmp_indices;
    tmp_indices = xmalloc(num_tmp_indices * sizeof (tmp_indices[0]));
    xt_idxlist_get_indices(idxlist, tmp_indices);
    Xt_idxlist tmp_idxvec;
    tmp_idxvec = xt_idxvec_new(tmp_indices, num_tmp_indices);
    int retval
      = tmp_idxvec->vtable->get_positions_of_indices(tmp_idxvec, indices,
                                                     num_indices, positions,
                                                     single_match_only);
    xt_idxlist_delete(tmp_idxvec);
    free(tmp_indices);
    return retval;
  }
}

int xt_idxlist_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                         int * position, int offset) {

  return idxlist->vtable->get_position_of_index_off(idxlist, index,
                                                    position, offset);
}

int
xt_idxlist_get_positions_of_indices_off(Xt_idxlist idxlist,
                                        const Xt_idx *indices,
                                        Xt_count num_indices, int * positions,
                                        int * offsets) {

  if (idxlist->vtable->get_positions_of_indices_off != NULL)
    return idxlist->vtable->get_positions_of_indices_off(idxlist, indices,
                                                         num_indices,
                                                         positions, offsets);

  Xt_count i;

  for (i = 0; i < num_indices; ++i) {

    int ret_val, temp_position, offset;

    if (offsets == NULL)
      offset = 0;
    else
      offset = offsets[i];

    ret_val
      = idxlist->vtable->get_position_of_index_off(idxlist, indices[i],
                                                   &temp_position, offset);

    if (ret_val)
      return ret_val;

    positions[i] = temp_position;
  }

  return 0;
}

Xt_idx xt_idxlist_get_min_index(Xt_idxlist idxlist) {

  return idxlist->vtable->get_min_index(idxlist);
}

Xt_idx xt_idxlist_get_max_index(Xt_idxlist idxlist) {

  return idxlist->vtable->get_max_index(idxlist);
}
