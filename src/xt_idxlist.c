#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "yaxt.h"
#include "xt_idxlist.h"
#include "xt_idxvec.h"
#include "xt_mpi.h"

void xt_idxlist_delete(struct xt_idxlist *idxlist) {

   idxlist->vtable->delete(idxlist);
}

size_t xt_idxlist_get_pack_size(struct xt_idxlist *idxlist,
                                MPI_Comm comm) {

   return idxlist->vtable->get_pack_size(idxlist, comm);
}

void xt_idxlist_pack(struct xt_idxlist *idxlist, void *buffer,
                     int buffer_size, int *position,
                     MPI_Comm comm) {

   idxlist->vtable->pack(idxlist, buffer, buffer_size,
                         position, comm);
}

xt_idxlist xt_idxlist_unpack(void *buffer, int buffer_size,
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
         return xt_idxlist_unpack_routines[i].xt_idxlist_unpack(buffer, buffer_size,
                                                                position, comm);
         break;
      }
   }

   fprintf(stderr, "xt_idxlist_unpack: unknown index list type\n");
   exit(EXIT_FAILURE);
}

xt_idxlist xt_idxlist_get_intersection(xt_idxlist idxlist_src,
                                       xt_idxlist idxlist_dst) {

   xt_idxlist intersection;

   intersection = idxlist_src->vtable->get_intersection(idxlist_src, idxlist_dst);

   // if the get_intersection routine was not able to compute the intersection
   if (intersection == NULL) {

      size_t num_indices_src, num_indices_dst;

      num_indices_src = xt_idxlist_get_num_indices(idxlist_src);
      num_indices_dst = xt_idxlist_get_num_indices(idxlist_dst);

      Xt_idx *indices_src, *indices_dst;
      indices_src = malloc(num_indices_src * sizeof (indices_src[0]));
      indices_dst = malloc(num_indices_dst * sizeof (indices_dst[0]));

      xt_idxlist_get_indices(idxlist_src, indices_src);
      xt_idxlist_get_indices(idxlist_dst, indices_dst);

      xt_idxlist idxvec_src, idxvec_dst;

      idxvec_src = xt_idxvec_new(indices_src, num_indices_src);
      idxvec_dst = xt_idxvec_new(indices_dst, num_indices_dst);

      intersection = xt_idxlist_get_intersection(idxvec_src,
                                                 idxvec_dst);

      xt_idxlist_delete(idxvec_src);
      xt_idxlist_delete(idxvec_dst);
      free(indices_dst);
      free(indices_src);

      return intersection;
   }
}

xt_idxlist xt_idxlist_copy(xt_idxlist idxlist) {

   return idxlist->vtable->copy(idxlist);
}

int xt_idxlist_get_num_indices(xt_idxlist idxlist) {

   return idxlist->vtable->get_num_indices(idxlist);
}

void xt_idxlist_get_indices(xt_idxlist idxlist, Xt_idx *indices) {

   idxlist->vtable->get_indices(idxlist, indices);
}

int xt_idxlist_get_index(Xt_idxlist idxlist, int position,
                         Xt_idx * index) {

  return idxlist->vtable->get_index(idxlist, position, index);
}

int xt_idxlist_get_position_of_index(Xt_idxlist idxlist, Xt_idx index,
                                     int * position) {

  return idxlist->vtable->get_position_of_index(idxlist, index, position);
}

int xt_idxlist_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                         int * position, int offset) {

  return idxlist->vtable->get_position_of_index_off(idxlist, index,
                                                    position, offset);
}
