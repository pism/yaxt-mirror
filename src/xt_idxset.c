#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xt_idxset.h"

// index list declarations

void* idxlist_new(int* idxlist, unsigned num_indices);
void  idxlist_delete(void * data);

struct idxset_vtable {
   void* (*idxset_new)(int*, unsigned);
   void  (*idxset_delete)(void *);
};

static struct idxset_vtable idxlist_vtable = {idxlist_new,
                                              idxlist_delete};

// xt_idxset data structure

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

   idxset->data = idxset->vtable.idxset_new(idxlist, num_indices);

   return idxset;
}

void xt_idxset_delete(struct xt_idxset* idxset) {

   idxset->vtable.idxset_delete(idxset->data);
   free(idxset);
}

// index list data structure

struct idxlist_data {

   int* list;
   unsigned num_indices;
};

// index list routines

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

   return (void*) data;
}

void idxlist_delete(void * data) {

   free(((struct idxlist_data*)data)->list);
   free(data);
}
