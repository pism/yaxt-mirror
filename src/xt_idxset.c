#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xt_idxset.h"
#include "xt_mpi.h"

// index list declarations

void*     idxlist_new(int* idxlist, unsigned num_indices);
void      idxlist_delete(void* data);
size_t    idxlist_get_pack_size(void* data, MPI_Comm comm);
void      idxlist_pack(void* data, void* buffer, int buffer_size,
                       int* position, MPI_Comm comm);
void*     idxlist_unpack(void* buffer, int buffer_size, int* position,
                         MPI_Comm comm);
xt_idxset idxlist_get_intersection(void* data, xt_idxset idxset);

struct idxset_vtable {
   void* (*new)(int*, unsigned);
   void  (*delete)(void *);

   size_t (*get_pack_size)(void*, MPI_Comm);
   void (*pack)(void*, void*, int, int*, MPI_Comm);
   void* (*unpack)(void*, int, int*, MPI_Comm);
   xt_idxset (*get_intersection)(void*,xt_idxset);
};

static struct idxset_vtable idxlist_vtable = {idxlist_new,
                                              idxlist_delete,
                                              idxlist_get_pack_size,
                                              idxlist_pack,
                                              idxlist_unpack,
                                              idxlist_get_intersection};

// xt_idxset data structure

enum xt_idxset_type {
   LIST = 1,
};

struct xt_idxset {
   struct idxset_vtable vtable;
   void * data;
};

// xt_idxset routines

struct xt_idxset* xt_idxset_new(int* idxlist, unsigned num_indices) {

   struct xt_idxset* idxset;

   idxset = malloc(1 * sizeof(*idxset));
   if (idxset == NULL) {
      fprintf(stderr, "malloc failed for idxset\n");
      exit(EXIT_FAILURE);
   }

   idxset->vtable = idxlist_vtable;

   idxset->data = idxset->vtable.new(idxlist, num_indices);

   return idxset;
}

void xt_idxset_delete(struct xt_idxset* idxset) {

   idxset->vtable.delete(idxset->data);
   free(idxset);
}

size_t xt_idxset_get_pack_size(struct xt_idxset* idxset, MPI_Comm comm) {

   return idxset->vtable.get_pack_size(idxset->data, comm);
}

void xt_idxset_pack(struct xt_idxset* idxset, void* buffer, int buffer_size,
                     int* position, MPI_Comm comm) {

   idxset->vtable.pack(idxset->data, buffer, buffer_size, position, comm);
}

xt_idxset xt_idxset_unpack(void* buffer, int buffer_size, int* position, MPI_Comm comm) {

   struct xt_idxset * idxset;

   idxset = malloc(1 * sizeof(*idxset));
   if (idxset == NULL) {
      fprintf(stderr, "malloc failed for idxset\n");
      exit(EXIT_FAILURE);
   }

   char type;

   // unpack the type of the index list to be unpacked
   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                   &type, 1, MPI_CHAR, comm), comm);

   switch(type) {

      case LIST:
         idxset->vtable = idxlist_vtable;
         break;
      default:
         fprintf(stderr, "xt_idxset_unpack: unknown index list type\n");
         exit(EXIT_FAILURE);
   };

   idxset->data = idxset->vtable.unpack(buffer, buffer_size, position, comm);

   return idxset;
}

xt_idxset xt_idxset_get_intersection(xt_idxset idxset_a, xt_idxset idxset_b) {

   return idxset_a->vtable.get_intersection(idxset_a->data, idxset_b);
}

// index list data structure

struct idxlist_data {

   int* list;
   unsigned num_indices;
};

// index list routines

static int compare_int (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

void* idxlist_new(int* idxlist, unsigned num_indices) {

   struct idxlist_data* data;

   data = malloc(1 * sizeof(*data));
   if (data == NULL) {
      fprintf(stderr, "malloc failed for idxlist_data\n");
      exit(EXIT_FAILURE);
   }

   size_t list_size;

   list_size = num_indices * sizeof(*(data->list));

   data->list = malloc(list_size);
   if (data->list == NULL) {
      fprintf(stderr, "malloc failed for idxlist_data->data\n");
      exit(EXIT_FAILURE);
   }

   memcpy(data->list, idxlist, list_size);
   data->num_indices = num_indices;

   // sort the list
   qsort(data->list, data->num_indices, sizeof(int), compare_int);

   // check for double entries
   unsigned i;
   for (i = 1; i < data->num_indices; ++i) {
      if (data->list[i] == data->list[i-1]) {

         fprintf(stderr, "idxlist_new: found non-unique value in index list\n");
         exit(EXIT_FAILURE);
      }
   }

   return (void*) data;
}

void idxlist_delete(void * data) {

   free(((struct idxlist_data*)data)->list);
   free(data);
}

size_t idxlist_get_pack_size(void* data, MPI_Comm comm) {

   struct idxlist_data * idxlist = data;

   int size_unsigned, size_int, size_char;

   mpi_err_handler(MPI_Pack_size(1, MPI_CHAR, comm, &size_char), comm);
   mpi_err_handler(MPI_Pack_size(idxlist->num_indices, MPI_INTEGER, comm, &size_int), comm);
   mpi_err_handler(MPI_Pack_size(1, MPI_UNSIGNED, comm, &size_unsigned), comm);

   return size_unsigned + size_int + size_char;
}

void   idxlist_pack(void* data, void* buffer, int buffer_size,
                    int* position, MPI_Comm comm) {

   struct idxlist_data * idxlist = data;
   enum xt_idxset_type type = LIST;

   mpi_err_handler(MPI_Pack(&(type), 1, MPI_CHAR, buffer,
                            buffer_size, position, comm), comm);
   mpi_err_handler(MPI_Pack(&(idxlist->num_indices), 1, MPI_UNSIGNED, buffer,
                            buffer_size, position, comm), comm);
   mpi_err_handler(MPI_Pack(idxlist->list, idxlist->num_indices, MPI_INTEGER, buffer,
                            buffer_size, position, comm), comm);
}

void* idxlist_unpack(void* buffer, int buffer_size, int* position,
                     MPI_Comm comm) {

   struct idxlist_data * data;

   data = malloc(1 * sizeof(*data));
   if (data == NULL) {
      fprintf(stderr, "malloc failed for idxlist_data\n");
      exit(EXIT_FAILURE);
   }

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              &(data->num_indices), 1, MPI_UNSIGNED, comm), comm);

   data->list = malloc(data->num_indices * sizeof(*(data->list)));
   if (data->list == NULL) {
      fprintf(stderr, "malloc failed for idxlist_data->data\n");
      exit(EXIT_FAILURE);
   }

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              data->list, data->num_indices, MPI_INTEGER, comm), comm);

   return data;
}

xt_idxset idxlist_get_intersection(void* data, struct xt_idxset* idxset) {

   // if the index set is a index list
   if (idxset->vtable.new == idxlist_vtable.new) {

      struct idxlist_data* data_a, * data_b;

      data_a = data;
      data_b = idxset->data;

      struct idxlist_data* inter_list;

      inter_list = malloc(1 * sizeof(*inter_list));
      if (inter_list == NULL) {
         fprintf(stderr, "malloc failed for intersection list\n");
         exit(EXIT_FAILURE);
      }

      inter_list->num_indices = 0;
      inter_list->list =
         malloc(((data_a->num_indices < data_b->num_indices)?(data_a->num_indices):(data_b->num_indices))
                * sizeof(*(inter_list->list)));
      if (inter_list->list == NULL) {
         fprintf(stderr, "malloc failed for intersection list list\n");
         exit(EXIT_FAILURE);
      }

      unsigned i, j;

      for (i = 0; i < data_a->num_indices; ++i)
         for (j = 0; j < data_b->num_indices; ++j)
            if (data_a->list[i] == data_b->list[j])
               inter_list->list[inter_list->num_indices++] = data_a->list[i];

      if (inter_list->num_indices) {

         inter_list->list = realloc(inter_list->list,
                                    inter_list->num_indices * sizeof(*(inter_list->list)));
         if (inter_list->list == NULL) {
            fprintf(stderr, "malloc failed for intersection list list\n");
            exit(EXIT_FAILURE);
         }
      } else
         inter_list->list = NULL;

      struct xt_idxset* idxset;

      idxset = malloc(1 * sizeof(*idxset));
      if (idxset == NULL) {
         fprintf(stderr, "malloc failed for index set\n");
         exit(EXIT_FAILURE);
      }

      idxset->vtable = idxlist_vtable;
      idxset->data = inter_list;

      return idxset;

   } else {

      fprintf(stderr, "idxlist_get_intersection: incompatible index set\n");
      exit(EXIT_FAILURE);
   }
}
