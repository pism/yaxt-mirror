#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

#include "yaxt.h"
#include "xt_idxvec.h"
#include "xt_idxvec_collection.h"
#include "xt_mpi.h"

static void
idxvec_collection_delete(xt_idxvec data);

static size_t
idxvec_collection_get_pack_size(xt_idxvec data, MPI_Comm comm);

static void
idxvec_collection_pack(xt_idxvec data, void *buffer, int buffer_size,
                       int *position, MPI_Comm comm);

static xt_idxvec
idxvec_collection_get_intersection(xt_idxvec idxvec_src, xt_idxvec idxvec_dst);

static xt_idxvec
idxvec_collection_copy(xt_idxvec idxvec);

static size_t
idxvec_collection_get_num_indices(xt_idxvec idxvec);

static void
idxvec_collection_get_indices(xt_idxvec idxvec, Xt_idx *indices);

static const struct xt_idxvec_vtable idxvec_collection_vtable =
      {idxvec_collection_delete,
       idxvec_collection_get_pack_size,
       idxvec_collection_pack,
       idxvec_collection_get_intersection,
       idxvec_collection_copy,
       idxvec_collection_get_num_indices,
       idxvec_collection_get_indices};

struct xt_idxvec_collection {

   const struct xt_idxvec_vtable *vtable;

   xt_idxvec * idxvecs;
   unsigned num_idxvecs;
};

xt_idxvec xt_idxvec_collection_new(xt_idxvec *idxvecs, unsigned num_idxvecs) {

   struct xt_idxvec_collection *collectionvector;

   collectionvector = malloc(1 * sizeof(*collectionvector));
   if (collectionvector == NULL) {
      fprintf(stderr, "malloc failed for collectionvector\n");
      exit(EXIT_FAILURE);
   }

   collectionvector->vtable = &idxvec_collection_vtable;
   collectionvector->num_idxvecs = num_idxvecs;
   collectionvector->idxvecs = malloc(num_idxvecs * sizeof(*(collectionvector->idxvecs)));
   if (collectionvector->idxvecs == NULL) {
      fprintf(stderr, "malloc failed for collectionvector->idxvecs\n");
      exit(EXIT_FAILURE);
   }

   unsigned i;

   for (i = 0; i < num_idxvecs; ++i)
      collectionvector->idxvecs[i] = xt_idxvec_copy(idxvecs[i]);

   return (xt_idxvec)collectionvector;
}

void idxvec_collection_delete(xt_idxvec data) {

   struct xt_idxvec_collection *collectionvector = (struct xt_idxvec_collection*)data;

   unsigned i;

   for (i = 0; i < collectionvector->num_idxvecs; ++i)
      xt_idxvec_delete(collectionvector->idxvecs[i]);

   free(collectionvector->idxvecs);
   free(collectionvector);
}

size_t idxvec_collection_get_pack_size(xt_idxvec data, MPI_Comm comm) {

   struct xt_idxvec_collection *collectionvector = (struct xt_idxvec_collection*)data;

   int size_int, size_unsigned;
   size_t size_idxvecs = 0;
   unsigned i;

   mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_int), comm);
   mpi_err_handler(MPI_Pack_size(1, MPI_UNSIGNED, comm, &size_unsigned), comm);

   for (i = 0; i < collectionvector->num_idxvecs; ++i)
      size_idxvecs += xt_idxvec_get_pack_size(collectionvector->idxvecs[i], comm);

   return size_int + size_unsigned + size_idxvecs;
}

void idxvec_collection_pack(xt_idxvec data, void *buffer, int buffer_size, int *position,
                            MPI_Comm comm) {

   struct xt_idxvec_collection *collectionvector = (struct xt_idxvec_collection*)data;
   int type = COLLECTION;

   mpi_err_handler(MPI_Pack(&(type), 1, MPI_INT, buffer,
                            buffer_size, position, comm), comm);
   mpi_err_handler(MPI_Pack(&(collectionvector->num_idxvecs), 1, MPI_UNSIGNED, buffer,
                            buffer_size, position, comm), comm);

   unsigned i;

   for (i = 0; i < collectionvector->num_idxvecs; ++i)
      xt_idxvec_pack(collectionvector->idxvecs[i], buffer, buffer_size, position, comm);
}

xt_idxvec xt_idxvec_collection_unpack(void *buffer, int buffer_size, int *position,
                                      MPI_Comm comm) {

   struct xt_idxvec_collection *collectionvector;

   collectionvector = malloc(1 * sizeof(*collectionvector));
   if (collectionvector == NULL) {
      fprintf(stderr, "malloc failed for collectionvector\n");
      exit(EXIT_FAILURE);
   }

   collectionvector->vtable = &idxvec_collection_vtable;

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              &(collectionvector->num_idxvecs), 1, MPI_UNSIGNED, comm),
                   comm);
   
   collectionvector->idxvecs = malloc(collectionvector->num_idxvecs *
                                   sizeof(*(collectionvector->idxvecs)));
   if (collectionvector->idxvecs == NULL) {
      fprintf(stderr, "malloc failed for collectionvector->idxvecs\n");
      exit(EXIT_FAILURE);
   }

   unsigned i;

   for (i = 0; i < collectionvector->num_idxvecs; ++i)
      collectionvector->idxvecs[i] = xt_idxvec_unpack(buffer, buffer_size, position, comm);

   return (xt_idxvec)collectionvector;
}

xt_idxvec idxvec_collection_get_intersection(xt_idxvec idxvec_src, xt_idxvec idxvec_dst) {

   return NULL;
}

xt_idxvec idxvec_collection_copy(xt_idxvec idxvec) {

   struct xt_idxvec_collection *collectionvector = (struct xt_idxvec_collection*)idxvec;

   return xt_idxvec_collection_new(collectionvector->idxvecs, collectionvector->num_idxvecs);
}

static size_t
idxvec_collection_get_num_indices(xt_idxvec idxvec) {

   struct xt_idxvec_collection *collectionvector = (struct xt_idxvec_collection*)idxvec;

   unsigned i;
   int num_indices = 0;

   for (i = 0; i < collectionvector->num_idxvecs; ++i)
      num_indices += xt_idxvec_get_num_indices(collectionvector->idxvecs[i]);

   return num_indices;
}

static void
idxvec_collection_get_indices(xt_idxvec idxvec, Xt_idx *indices) {

   struct xt_idxvec_collection *collectionvector = (struct xt_idxvec_collection*)idxvec;

   unsigned i;
   int offvector = 0;

   for (i = 0; i < collectionvector->num_idxvecs; ++i) {

      xt_idxvec_get_indices(collectionvector->idxvecs[i], indices+offvector);
      offvector += xt_idxvec_get_num_indices(collectionvector->idxvecs[i]);
   }
}
