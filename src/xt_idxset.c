#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "yaxt.h"
#include "xt_idxset.h"
#include "xt_idxlist.h"
#include "xt_mpi.h"

void xt_idxset_delete(struct xt_idxset *idxset) {

   idxset->vtable->delete(idxset);
}

size_t xt_idxset_get_pack_size(struct xt_idxset *idxset, MPI_Comm comm) {

   return idxset->vtable->get_pack_size(idxset, comm);
}

void xt_idxset_pack(struct xt_idxset *idxset, void *buffer, int buffer_size,
                    int *position, MPI_Comm comm) {

   idxset->vtable->pack(idxset, buffer, buffer_size, position, comm);
}

xt_idxset xt_idxset_unpack(void *buffer, int buffer_size, int *position, MPI_Comm comm) {

   int type;

   // unpack the type of the index list to be unpacked
   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                   &type, 1, MPI_INT, comm), comm);

   unsigned i;
   extern struct xt_idxset_unpack_routine_table xt_idxset_unpack_routines[];
   extern unsigned xt_num_unpack_routines;

   for (i = 0; i < xt_num_unpack_routines; ++i) {

      if (type == xt_idxset_unpack_routines[i].type) {
         return xt_idxset_unpack_routines[i].xt_idxset_unpack(buffer, buffer_size,
                                                              position, comm);
         break;
      }
   }

   fprintf(stderr, "xt_idxset_unpack: unknown index set type\n");
   exit(EXIT_FAILURE);
}

xt_idxset xt_idxset_get_intersection(xt_idxset idxset_a, xt_idxset idxset_b) {

   xt_idxset intersection;

   intersection = idxset_a->vtable->get_intersection(idxset_a, idxset_b);

   // if the get_intersection routine was not able to compute the intersection
   if (intersection == NULL) {

      int num_indices_a, num_indices_b;

      num_indices_a = xt_idxset_get_num_indices(idxset_a);
      num_indices_b = xt_idxset_get_num_indices(idxset_b);

      int indices_a[num_indices_a], indices_b[num_indices_b];

      xt_idxset_get_indices(idxset_a, indices_a);
      xt_idxset_get_indices(idxset_b, indices_b);

      xt_idxset idxlist_a, idxlist_b;

      idxlist_a = xt_idxlist_new(indices_a, num_indices_a);
      idxlist_b = xt_idxlist_new(indices_b, num_indices_b);

      intersection = xt_idxset_get_intersection(idxlist_a, idxlist_b);

      xt_idxset_delete(idxlist_a);
      xt_idxset_delete(idxlist_b);

      return intersection;
   }
}

xt_idxset xt_idxset_copy(xt_idxset idxset) {

   return idxset->vtable->copy(idxset);
}

int xt_idxset_get_num_indices(xt_idxset idxset) {

   return idxset->vtable->get_num_indices(idxset);
}

void xt_idxset_get_indices(xt_idxset idxset, int *indices) {

   idxset->vtable->get_indices(idxset, indices);
}
