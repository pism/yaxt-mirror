#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

#include "yaxt.h"
#include "xt_idxset.h"
#include "xt_idxset_multi.h"
#include "xt_mpi.h"

void      idxset_multi_delete(xt_idxset data);
size_t    idxset_multi_get_pack_size(xt_idxset data, MPI_Comm comm);
void      idxset_multi_pack(xt_idxset data, void *buffer, int buffer_size,
                            int *position, MPI_Comm comm);
xt_idxset idxset_multi_get_intersection(xt_idxset idxset_a, xt_idxset idxset_b);
xt_idxset idxset_multi_copy(xt_idxset idxset);
int       idxset_multi_get_num_indices(xt_idxset idxset);
void      idxset_multi_get_indices(xt_idxset idxset, int *indices);

static struct xt_idxset_vtable idxset_multi_vtable = {idxset_multi_delete,
                                                      idxset_multi_get_pack_size,
                                                      idxset_multi_pack,
                                                      idxset_multi_get_intersection,
                                                      idxset_multi_copy,
                                                      idxset_multi_get_num_indices,
                                                      idxset_multi_get_indices};

struct xt_idxset_multi {

   struct xt_idxset_vtable *vtable;

   xt_idxset * idxsets;
   unsigned num_idxsets;
};

xt_idxset xt_idxset_multi_new(xt_idxset *idxsets, unsigned num_idxsets) {

   struct xt_idxset_multi *multiset;

   multiset = malloc(1 * sizeof(*multiset));
   if (multiset == NULL) {
      fprintf(stderr, "malloc failed for multiset\n");
      exit(EXIT_FAILURE);
   }

   multiset->vtable = &idxset_multi_vtable;
   multiset->num_idxsets = num_idxsets;
   multiset->idxsets = malloc(num_idxsets * sizeof(*(multiset->idxsets)));
   if (multiset->idxsets == NULL) {
      fprintf(stderr, "malloc failed for multiset->idxsets\n");
      exit(EXIT_FAILURE);
   }

   unsigned i;

   for (i = 0; i < num_idxsets; ++i)
      multiset->idxsets[i] = xt_idxset_copy(idxsets[i]);

   return (xt_idxset)multiset;
}

void idxset_multi_delete(xt_idxset data) {

   struct xt_idxset_multi *multiset = (struct xt_idxset_multi*)data;

   unsigned i;

   for (i = 0; i < multiset->num_idxsets; ++i)
      xt_idxset_delete(multiset->idxsets[i]);

   free(multiset->idxsets);
   free(multiset);
}

size_t idxset_multi_get_pack_size(xt_idxset data, MPI_Comm comm) {

   struct xt_idxset_multi *multiset = (struct xt_idxset_multi*)data;

   int size_char, size_unsigned;
   size_t size_idxsets = 0;
   unsigned i;

   mpi_err_handler(MPI_Pack_size(1, MPI_CHAR, comm, &size_char), comm);
   mpi_err_handler(MPI_Pack_size(1, MPI_UNSIGNED, comm, &size_unsigned), comm);

   for (i = 0; i < multiset->num_idxsets; ++i)
      size_idxsets += xt_idxset_get_pack_size(multiset->idxsets[i], comm);

   return size_char + size_unsigned + size_idxsets;
}

void idxset_multi_pack(xt_idxset data, void *buffer, int buffer_size, int *position,
                       MPI_Comm comm) {

   struct xt_idxset_multi *multiset = (struct xt_idxset_multi*)data;
   char type = MULTI;

   mpi_err_handler(MPI_Pack(&(type), 1, MPI_CHAR, buffer,
                            buffer_size, position, comm), comm);
   mpi_err_handler(MPI_Pack(&(multiset->num_idxsets), 1, MPI_UNSIGNED, buffer,
                            buffer_size, position, comm), comm);

   unsigned i;

   for (i = 0; i < multiset->num_idxsets; ++i)
      xt_idxset_pack(multiset->idxsets[i], buffer, buffer_size, position, comm);
}

xt_idxset xt_idxset_multi_unpack(void *buffer, int buffer_size, int *position,
                                 MPI_Comm comm) {

   struct xt_idxset_multi *multiset;

   multiset = malloc(1 * sizeof(*multiset));
   if (multiset == NULL) {
      fprintf(stderr, "malloc failed for multiset\n");
      exit(EXIT_FAILURE);
   }

   multiset->vtable = &idxset_multi_vtable;

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              &(multiset->num_idxsets), 1, MPI_UNSIGNED, comm), comm);
   
   multiset->idxsets = malloc(multiset->num_idxsets * sizeof(*(multiset->idxsets)));
   if (multiset->idxsets == NULL) {
      fprintf(stderr, "malloc failed for multiset->idxsets\n");
      exit(EXIT_FAILURE);
   }

   unsigned i;

   for (i = 0; i < multiset->num_idxsets; ++i)
      multiset->idxsets[i] = xt_idxset_unpack(buffer, buffer_size, position, comm);

   return (xt_idxset)multiset;
}

xt_idxset idxset_multi_get_intersection(xt_idxset idxset_a, xt_idxset idxset_b) {

   struct xt_idxset_multi *intersection, *multiset = (struct xt_idxset_multi*)idxset_a;

   intersection = malloc(1 * sizeof(*multiset));
   if (intersection == NULL) {
      fprintf(stderr, "malloc failed for intersection\n");
      exit(EXIT_FAILURE);
   }

   intersection->vtable = &idxset_multi_vtable;
   intersection->num_idxsets = multiset->num_idxsets;
   intersection->idxsets = malloc(multiset->num_idxsets * sizeof(*(intersection->idxsets)));
   if (intersection->idxsets == NULL) {
      fprintf(stderr, "malloc failed for intersection->idxsets\n");
      exit(EXIT_FAILURE);
   }

   unsigned i;

   for (i = 0; i < intersection->num_idxsets; ++i)
      intersection->idxsets[i] = xt_idxset_get_intersection(multiset->idxsets[i],
                                                            idxset_b);

   return (xt_idxset)intersection;
}

xt_idxset idxset_multi_copy(xt_idxset idxset) {

   struct xt_idxset_multi *multiset = (struct xt_idxset_multi*)idxset;

   return xt_idxset_multi_new(multiset->idxsets, multiset->num_idxsets);
}

int idxset_multi_get_num_indices(xt_idxset idxset) {

   struct xt_idxset_multi *multiset = (struct xt_idxset_multi*)idxset;

   unsigned i;
   int num_indices = 0;

   for (i = 0; i < multiset->num_idxsets; ++i)
      num_indices += xt_idxset_get_num_indices(multiset->idxsets[i]);

   return num_indices;
}

void idxset_multi_get_indices(xt_idxset idxset, int *indices) {

   struct xt_idxset_multi *multiset = (struct xt_idxset_multi*)idxset;

   unsigned i;
   int offset = 0;

   for (i = 0; i < multiset->num_idxsets; ++i) {

      xt_idxset_get_indices(multiset->idxsets[i], indices+offset);
      offset += xt_idxset_get_num_indices(multiset->idxsets[i]);
   }
}
