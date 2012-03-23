#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

#include "yaxt.h"
#include "xt_idxset.h"
#include "xt_idxset_collection.h"
#include "xt_mpi.h"

void      idxset_collection_delete(xt_idxset data);
size_t    idxset_collection_get_pack_size(xt_idxset data, MPI_Comm comm);
void      idxset_collection_pack(xt_idxset data, void *buffer, int buffer_size,
                                 int *position, MPI_Comm comm);
xt_idxset idxset_collection_get_intersection(xt_idxset idxset_a, xt_idxset idxset_b);
xt_idxset idxset_collection_copy(xt_idxset idxset);
int       idxset_collection_get_num_indices(xt_idxset idxset);
void      idxset_collection_get_indices(xt_idxset idxset, int *indices);

static struct xt_idxset_vtable idxset_collection_vtable = 
      {idxset_collection_delete,
       idxset_collection_get_pack_size,
       idxset_collection_pack,
       idxset_collection_get_intersection,
       idxset_collection_copy,
       idxset_collection_get_num_indices,
       idxset_collection_get_indices};

struct xt_idxset_collection {

   struct xt_idxset_vtable *vtable;

   xt_idxset * idxsets;
   unsigned num_idxsets;
};

xt_idxset xt_idxset_collection_new(xt_idxset *idxsets, unsigned num_idxsets) {

   struct xt_idxset_collection *collectionset;

   collectionset = malloc(1 * sizeof(*collectionset));
   if (collectionset == NULL) {
      fprintf(stderr, "malloc failed for collectionset\n");
      exit(EXIT_FAILURE);
   }

   collectionset->vtable = &idxset_collection_vtable;
   collectionset->num_idxsets = num_idxsets;
   collectionset->idxsets = malloc(num_idxsets * sizeof(*(collectionset->idxsets)));
   if (collectionset->idxsets == NULL) {
      fprintf(stderr, "malloc failed for collectionset->idxsets\n");
      exit(EXIT_FAILURE);
   }

   unsigned i;

   for (i = 0; i < num_idxsets; ++i)
      collectionset->idxsets[i] = xt_idxset_copy(idxsets[i]);

   return (xt_idxset)collectionset;
}

void idxset_collection_delete(xt_idxset data) {

   struct xt_idxset_collection *collectionset = (struct xt_idxset_collection*)data;

   unsigned i;

   for (i = 0; i < collectionset->num_idxsets; ++i)
      xt_idxset_delete(collectionset->idxsets[i]);

   free(collectionset->idxsets);
   free(collectionset);
}

size_t idxset_collection_get_pack_size(xt_idxset data, MPI_Comm comm) {

   struct xt_idxset_collection *collectionset = (struct xt_idxset_collection*)data;

   int size_int, size_unsigned;
   size_t size_idxsets = 0;
   unsigned i;

   mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_int), comm);
   mpi_err_handler(MPI_Pack_size(1, MPI_UNSIGNED, comm, &size_unsigned), comm);

   for (i = 0; i < collectionset->num_idxsets; ++i)
      size_idxsets += xt_idxset_get_pack_size(collectionset->idxsets[i], comm);

   return size_int + size_unsigned + size_idxsets;
}

void idxset_collection_pack(xt_idxset data, void *buffer, int buffer_size, int *position,
                            MPI_Comm comm) {

   struct xt_idxset_collection *collectionset = (struct xt_idxset_collection*)data;
   int type = COLLECTION;

   mpi_err_handler(MPI_Pack(&(type), 1, MPI_INT, buffer,
                            buffer_size, position, comm), comm);
   mpi_err_handler(MPI_Pack(&(collectionset->num_idxsets), 1, MPI_UNSIGNED, buffer,
                            buffer_size, position, comm), comm);

   unsigned i;

   for (i = 0; i < collectionset->num_idxsets; ++i)
      xt_idxset_pack(collectionset->idxsets[i], buffer, buffer_size, position, comm);
}

xt_idxset xt_idxset_collection_unpack(void *buffer, int buffer_size, int *position,
                                      MPI_Comm comm) {

   struct xt_idxset_collection *collectionset;

   collectionset = malloc(1 * sizeof(*collectionset));
   if (collectionset == NULL) {
      fprintf(stderr, "malloc failed for collectionset\n");
      exit(EXIT_FAILURE);
   }

   collectionset->vtable = &idxset_collection_vtable;

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              &(collectionset->num_idxsets), 1, MPI_UNSIGNED, comm),
                   comm);
   
   collectionset->idxsets = malloc(collectionset->num_idxsets *
                                   sizeof(*(collectionset->idxsets)));
   if (collectionset->idxsets == NULL) {
      fprintf(stderr, "malloc failed for collectionset->idxsets\n");
      exit(EXIT_FAILURE);
   }

   unsigned i;

   for (i = 0; i < collectionset->num_idxsets; ++i)
      collectionset->idxsets[i] = xt_idxset_unpack(buffer, buffer_size, position, comm);

   return (xt_idxset)collectionset;
}

xt_idxset idxset_collection_get_intersection(xt_idxset idxset_a, xt_idxset idxset_b) {

   struct xt_idxset_collection *intersection,
                               *collectionset = (struct xt_idxset_collection*)idxset_a;

   intersection = malloc(1 * sizeof(*collectionset));
   if (intersection == NULL) {
      fprintf(stderr, "malloc failed for intersection\n");
      exit(EXIT_FAILURE);
   }

   intersection->vtable = &idxset_collection_vtable;
   intersection->num_idxsets = collectionset->num_idxsets;
   intersection->idxsets = malloc(collectionset->num_idxsets *
                                  sizeof(*(intersection->idxsets)));
   if (intersection->idxsets == NULL) {
      fprintf(stderr, "malloc failed for intersection->idxsets\n");
      exit(EXIT_FAILURE);
   }

   unsigned i;

   for (i = 0; i < intersection->num_idxsets; ++i)
      intersection->idxsets[i] = xt_idxset_get_intersection(collectionset->idxsets[i],
                                                            idxset_b);

   return (xt_idxset)intersection;
}

xt_idxset idxset_collection_copy(xt_idxset idxset) {

   struct xt_idxset_collection *collectionset = (struct xt_idxset_collection*)idxset;

   return xt_idxset_collection_new(collectionset->idxsets, collectionset->num_idxsets);
}

int idxset_collection_get_num_indices(xt_idxset idxset) {

   struct xt_idxset_collection *collectionset = (struct xt_idxset_collection*)idxset;

   unsigned i;
   int num_indices = 0;

   for (i = 0; i < collectionset->num_idxsets; ++i)
      num_indices += xt_idxset_get_num_indices(collectionset->idxsets[i]);

   return num_indices;
}

void idxset_collection_get_indices(xt_idxset idxset, int *indices) {

   struct xt_idxset_collection *collectionset = (struct xt_idxset_collection*)idxset;

   unsigned i;
   int offset = 0;

   for (i = 0; i < collectionset->num_idxsets; ++i) {

      xt_idxset_get_indices(collectionset->idxsets[i], indices+offset);
      offset += xt_idxset_get_num_indices(collectionset->idxsets[i]);
   }
}
