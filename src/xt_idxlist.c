#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "yaxt.h"
#include "xt_idxvec.h"
#include "xt_idxlist.h"
#include "xt_mpi.h"

static void
idxlist_delete(xt_idxvec data);

static size_t
idxlist_get_pack_size(xt_idxvec data, MPI_Comm comm);

static void
idxlist_pack(xt_idxvec data, void *buffer, int buffer_size,
             int *position, MPI_Comm comm);

static xt_idxvec
idxlist_get_intersection(xt_idxvec idxvec_a, xt_idxvec idxvec_b);

static xt_idxvec
idxlist_copy(xt_idxvec idxvec);

static size_t
idxlist_get_num_indices(xt_idxvec idxvec);

static void
idxlist_get_indices(xt_idxvec idxvec, Xt_idx *indices);

static const struct xt_idxvec_vtable idxlist_vtable = {idxlist_delete,
                                                 idxlist_get_pack_size,
                                                 idxlist_pack,
                                                 idxlist_get_intersection,
                                                 idxlist_copy,
                                                 idxlist_get_num_indices,
                                                 idxlist_get_indices};

// index list data structure

struct xt_idxlist {

   const struct xt_idxvec_vtable *vtable;

   Xt_idx *list;
   size_t num_indices;
};

static int compare_xtidx(const void * a, const void * b)
{
  return (*(Xt_idx*)a > *(Xt_idx*)b) - (*(Xt_idx*)a < *(Xt_idx*)b);
}

struct xt_idxvec *xt_idxlist_new(Xt_idx *idxlist, size_t num_indices) {

   struct xt_idxlist *idxlist_data;

   idxlist_data = malloc(1 * sizeof(*idxlist_data));
   if (idxlist_data == NULL) {
      fprintf(stderr, "malloc failed for idxlist_data\n");
      exit(EXIT_FAILURE);
   }

   idxlist_data->vtable = &idxlist_vtable;

   size_t list_size;

   list_size = num_indices * sizeof(*(idxlist_data->list));

   idxlist_data->list = malloc(list_size);
   if (idxlist_data->list == NULL) {
      fprintf(stderr, "malloc failed for idxlist_data\n");
      exit(EXIT_FAILURE);
   }

   memcpy(idxlist_data->list, idxlist, list_size);
   idxlist_data->num_indices = num_indices;

   return (void*)idxlist_data;
}

static void idxlist_delete(struct xt_idxvec *data) {

   free(((struct xt_idxlist*)data)->list);
   free(data);
}

static size_t idxlist_get_pack_size(struct xt_idxvec *data, MPI_Comm comm) {

   struct xt_idxlist *idxlist = (struct xt_idxlist*)data;

   int size_unsigned, size_int, size_int_type;

   mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);
   mpi_err_handler(MPI_Pack_size(idxlist->num_indices, Xt_idx_dt, comm,
                                 &size_int), comm);
   mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_unsigned), comm);

   return size_unsigned + size_int + size_int_type;
}

void idxlist_pack(struct xt_idxvec *data, void *buffer, int buffer_size,
                  int *position, MPI_Comm comm) {

  assert(data);
  struct xt_idxlist *idxlist = (struct xt_idxlist*)data;
  int type = LIST;
  int num_indices = (int)(idxlist->num_indices);

  assert(idxlist->num_indices <= INT_MAX);

  mpi_err_handler(MPI_Pack(&(type), 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(&num_indices, 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(idxlist->list, (int)idxlist->num_indices,
                           Xt_idx_dt, buffer,
                           buffer_size, position, comm), comm);
}

xt_idxvec xt_idxlist_unpack(void *buffer, int buffer_size, int *position,
                            MPI_Comm comm) {

   struct xt_idxlist *idxlist;
   int num_indices;

   idxlist = malloc(1 * sizeof(*idxlist));
   if (idxlist == NULL) {
      fprintf(stderr, "malloc failed for idxlist_data\n");
      exit(EXIT_FAILURE);
   }

   idxlist->vtable = &idxlist_vtable;

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              &num_indices, 1, MPI_INT, comm), comm);
   assert(num_indices >= 0);
   idxlist->num_indices = (size_t)num_indices;

   idxlist->list = malloc(idxlist->num_indices * sizeof(*(idxlist->list)));
   if (idxlist->list == NULL) {
      fprintf(stderr, "malloc failed for idxlist\n");
      exit(EXIT_FAILURE);
   }

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              idxlist->list, idxlist->num_indices,
                              Xt_idx_dt, comm),comm);

   return (xt_idxvec)idxlist;
}

static xt_idxvec
idxlist_get_intersection(xt_idxvec idxvec_src, xt_idxvec idxvec_dst) {


  if (idxvec_dst->vtable != &idxlist_vtable) return NULL;

  // both vectors are index lists:

  struct xt_idxlist *idxlist_src, *idxlist_dst;

  idxlist_src = (struct xt_idxlist*)idxvec_src;
  idxlist_dst = (struct xt_idxlist*)idxvec_dst;

  struct xt_idxlist *inter_list;

  inter_list = malloc(1 * sizeof(*inter_list));
  if (inter_list == NULL) {
    fprintf(stderr, "malloc failed for intersection list\n");
    exit(EXIT_FAILURE);
  }

  inter_list->num_indices = 0;
  inter_list->list =
    malloc(idxlist_dst->num_indices * sizeof(*(inter_list->list)));
  if (inter_list->list == NULL) {
    fprintf(stderr, "malloc failed for intersection list list\n");
    exit(EXIT_FAILURE);
  }

  size_t i, j;

  for (j = 0; j < idxlist_dst->num_indices; ++j)
    for (i = 0; i < idxlist_src->num_indices; ++i)
      if (idxlist_src->list[i] == idxlist_dst->list[j]) {
        inter_list->list[inter_list->num_indices++] = idxlist_src->list[i];
         break;
      }

  if (inter_list->num_indices) {

    // sort the list
    qsort(inter_list->list, inter_list->num_indices,
          sizeof (inter_list->list[0]), compare_xtidx);
    inter_list->list = realloc(inter_list->list,
                               inter_list->num_indices *
                               sizeof(*(inter_list->list)));
    if (inter_list->list == NULL) {
      fprintf(stderr, "malloc failed for intersection list list\n");
      exit(EXIT_FAILURE);
    }
  } else {
    free(inter_list->list);
    inter_list->list = NULL;
  }

  inter_list->vtable = &idxlist_vtable;

  return (xt_idxvec)inter_list;
}

static xt_idxvec
idxlist_copy(xt_idxvec idxvec) {

   struct xt_idxlist *idxlist_data;

   idxlist_data = (struct xt_idxlist *)idxvec;

   return xt_idxlist_new(idxlist_data->list, idxlist_data->num_indices);
}

static size_t
idxlist_get_num_indices(xt_idxvec idxvec) {

   struct xt_idxlist *idxlist_data;

   idxlist_data = (struct xt_idxlist *)idxvec;

   return idxlist_data->num_indices;
}

static void
idxlist_get_indices(xt_idxvec idxvec, Xt_idx *indices) {

   struct xt_idxlist *idxlist_data;

   idxlist_data = (struct xt_idxlist *)idxvec;

   memcpy(indices, idxlist_data->list,
          idxlist_data->num_indices * sizeof(*indices));
}
