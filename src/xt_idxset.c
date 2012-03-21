#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "yaxt.h"
#include "xt_idxset.h"
#include "xt_idxlist.h"
#include "xt_idxset_unpack.h"
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

   char type;

   // unpack the type of the index list to be unpacked
   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                   &type, 1, MPI_CHAR, comm), comm);

   unsigned i;
   unsigned num_unpack_routines = sizeof(xt_idxset_unpack_routines) /
                                  sizeof(xt_idxset_unpack_routines[0]);

   for (i = 0; i < num_unpack_routines; ++i) {

      if (type = xt_idxset_unpack_routines[i].type) {
         return xt_idxset_unpack_routines[i].xt_idxset_unpack(buffer, buffer_size,
                                                              position, comm);
         break;
      }
   }

   fprintf(stderr, "xt_idxset_unpack: unknown index set type\n");
   exit(EXIT_FAILURE);
}

xt_idxset xt_idxset_get_intersection(xt_idxset idxset_a, xt_idxset idxset_b) {

   return idxset_a->vtable->get_intersection(idxset_a, idxset_b);
}
