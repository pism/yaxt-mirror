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

xt_idxvec xt_idxvec_get_intersection(xt_idxvec idxvec_a, xt_idxvec idxvec_b) {

   xt_idxvec intersection;

   intersection = idxvec_a->vtable->get_intersection(idxvec_a, idxvec_b);

   // if the get_intersection routine was not able to compute the intersection
   if (intersection == NULL) {

      int num_indices_a, num_indices_b;

      num_indices_a = xt_idxvec_get_num_indices(idxvec_a);
      num_indices_b = xt_idxvec_get_num_indices(idxvec_b);

      int indices_a[num_indices_a], indices_b[num_indices_b];

      xt_idxvec_get_indices(idxvec_a, indices_a);
      xt_idxvec_get_indices(idxvec_b, indices_b);

      xt_idxvec idxlist_a, idxlist_b;

      idxlist_a = xt_idxlist_new(indices_a, num_indices_a);
      idxlist_b = xt_idxlist_new(indices_b, num_indices_b);

      intersection = xt_idxvec_get_intersection(idxlist_a, idxlist_b);

      xt_idxvec_delete(idxlist_a);
      xt_idxvec_delete(idxlist_b);

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
