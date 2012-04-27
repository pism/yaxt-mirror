#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "yaxt.h"
#include "xt_idxvec.h"
#include "xt_idxlist.h"
#include "xt_mpi.h"

void xt_idxvec_delete(struct xt_idxvec *idxvec) {

   idxvec->vtable->delete(idxvec);
}

size_t xt_idxvec_get_pack_size(struct xt_idxvec *idxvec, MPI_Comm comm) {

   return idxvec->vtable->get_pack_size(idxvec, comm);
}

void xt_idxvec_pack(struct xt_idxvec *idxvec, void *buffer, int buffer_size,
                    int *position, MPI_Comm comm) {

   idxvec->vtable->pack(idxvec, buffer, buffer_size, position, comm);
}

xt_idxvec xt_idxvec_unpack(void *buffer, int buffer_size, int *position, MPI_Comm comm) {

   int type;

   // unpack the type of the index list to be unpacked
   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                   &type, 1, MPI_INT, comm), comm);

   unsigned i;
   extern struct xt_idxvec_unpack_routine_table xt_idxvec_unpack_routines[];
   extern unsigned xt_num_unpack_routines;

   for (i = 0; i < xt_num_unpack_routines; ++i) {

      if (type == xt_idxvec_unpack_routines[i].type) {
         return xt_idxvec_unpack_routines[i].xt_idxvec_unpack(buffer, buffer_size,
                                                              position, comm);
         break;
      }
   }

   fprintf(stderr, "xt_idxvec_unpack: unknown index vector type\n");
   exit(EXIT_FAILURE);
}

xt_idxvec xt_idxvec_get_intersection(xt_idxvec idxvec_src, xt_idxvec idxvec_dst) {

   xt_idxvec intersection;

   intersection = idxvec_src->vtable->get_intersection(idxvec_src, idxvec_dst);

   // if the get_intersection routine was not able to compute the intersection
   if (intersection == NULL) {

      size_t num_indices_src, num_indices_dst;

      num_indices_src = xt_idxvec_get_num_indices(idxvec_src);
      num_indices_dst = xt_idxvec_get_num_indices(idxvec_dst);

      Xt_idx *indices_src, *indices_dst;
      indices_src = malloc(num_indices_src * sizeof (indices_src[0]));
      indices_dst = malloc(num_indices_dst * sizeof (indices_dst[0]));

      xt_idxvec_get_indices(idxvec_src, indices_src);
      xt_idxvec_get_indices(idxvec_dst, indices_dst);

      xt_idxvec idxlist_src, idxlist_dst;

      idxlist_src = xt_idxlist_new(indices_src, num_indices_src);
      idxlist_dst = xt_idxlist_new(indices_dst, num_indices_dst);

      intersection = xt_idxvec_get_intersection(idxlist_src, idxlist_dst);

      xt_idxvec_delete(idxlist_src);
      xt_idxvec_delete(idxlist_dst);
      free(indices_dst);
      free(indices_src);

      return intersection;
   }
}

xt_idxvec xt_idxvec_copy(xt_idxvec idxvec) {

   return idxvec->vtable->copy(idxvec);
}

int xt_idxvec_get_num_indices(xt_idxvec idxvec) {

   return idxvec->vtable->get_num_indices(idxvec);
}

void xt_idxvec_get_indices(xt_idxvec idxvec, Xt_idx *indices) {

   idxvec->vtable->get_indices(idxvec, indices);
}

int xt_idxvec_get_index(Xt_idxvec idxvec, int position, Xt_idx * index) {

  return idxvec->vtable->get_index(idxvec, position, index);
}

int xt_idxvec_get_position_of_index(Xt_idxvec idxvec, Xt_idx index, int * position) {

  return idxvec->vtable->get_position_of_index(idxvec, index, position);
}
