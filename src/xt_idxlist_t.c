#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "yaxt.h"
#include "xt_idxvec.h"
#include "xt_idxlist_t.h"
#include "xt_mpi.h"

static void
idxlist_t_delete(xt_idxvec data);

static size_t
idxlist_t_get_pack_size(xt_idxvec data, MPI_Comm comm);

static void
idxlist_t_pack(xt_idxvec data, void *buffer, int buffer_size,
               int *position, MPI_Comm comm);

static xt_idxvec
idxlist_t_get_intersection(xt_idxvec idxvec_a, xt_idxvec idxvec_b);

static xt_idxvec
idxlist_t_copy(xt_idxvec idxvec);

static size_t
idxlist_t_get_num_indices(xt_idxvec idxvec);

static void
idxlist_t_get_indices(xt_idxvec idxvec, Xt_idx *indices);

static int
idxlist_t_get_index(Xt_idxvec idxvec, int position, Xt_idx * index);

static int
idxlist_t_get_position_of_index(Xt_idxvec idxvec, Xt_idx index, int * position);

static int
idxlist_t_get_position_of_index_off(Xt_idxvec idxvec, Xt_idx index,
                                    int * position, int offset);

static const struct xt_idxvec_vtable idxlist_t_vtable = {
              idxlist_t_delete,
              idxlist_t_get_pack_size,
              idxlist_t_pack,
              idxlist_t_get_intersection,
              idxlist_t_copy,
              idxlist_t_get_num_indices,
              idxlist_t_get_indices,
              idxlist_t_get_index,
              idxlist_t_get_position_of_index,
              idxlist_t_get_position_of_index_off};

// index list data structure

struct xt_idxlist_t {

   const struct xt_idxvec_vtable *vtable;

   Xt_idx *list;
   size_t num_indices;
};

static int compare_xtidx(const void * a, const void * b)
{
  return (*(Xt_idx*)a > *(Xt_idx*)b) - (*(Xt_idx*)a < *(Xt_idx*)b);
}

struct xt_idxvec *xt_idxlist_t_new(Xt_idx *idxlist_t, size_t num_indices) {

   struct xt_idxlist_t *idxlist_t_data;

   idxlist_t_data = malloc(1 * sizeof(*idxlist_t_data));
   if (idxlist_t_data == NULL) {
      fprintf(stderr, "malloc failed for idxlist_t_data\n");
      exit(EXIT_FAILURE);
   }

   idxlist_t_data->vtable = &idxlist_t_vtable;

   size_t list_size;

   list_size = num_indices * sizeof(*(idxlist_t_data->list));

   idxlist_t_data->list = malloc(list_size);
   if (idxlist_t_data->list == NULL) {
      fprintf(stderr, "malloc failed for idxlist_t_data\n");
      exit(EXIT_FAILURE);
   }

   memcpy(idxlist_t_data->list, idxlist_t, list_size);
   idxlist_t_data->num_indices = num_indices;

   return (void*)idxlist_t_data;
}

static void idxlist_t_delete(struct xt_idxvec *data) {

   free(((struct xt_idxlist_t*)data)->list);
   free(data);
}

static size_t idxlist_t_get_pack_size(struct xt_idxvec *data, MPI_Comm comm) {

   struct xt_idxlist_t *idxlist_t = (struct xt_idxlist_t*)data;

   int size_unsigned, size_int, size_int_type;

   mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);
   mpi_err_handler(MPI_Pack_size(idxlist_t->num_indices, Xt_idx_dt, comm,
                                 &size_int), comm);
   mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_unsigned), comm);

   return size_unsigned + size_int + size_int_type;
}

void idxlist_t_pack(struct xt_idxvec *data, void *buffer, int buffer_size,
                  int *position, MPI_Comm comm) {

  assert(data);
  struct xt_idxlist_t *idxlist_t = (struct xt_idxlist_t*)data;
  int type = LIST;
  int num_indices = (int)(idxlist_t->num_indices);

  assert(idxlist_t->num_indices <= INT_MAX);

  mpi_err_handler(MPI_Pack(&(type), 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(&num_indices, 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(idxlist_t->list, (int)idxlist_t->num_indices,
                           Xt_idx_dt, buffer,
                           buffer_size, position, comm), comm);
}

xt_idxvec xt_idxlist_t_unpack(void *buffer, int buffer_size, int *position,
                            MPI_Comm comm) {

   struct xt_idxlist_t *idxlist_t;
   int num_indices;

   idxlist_t = malloc(1 * sizeof(*idxlist_t));
   if (idxlist_t == NULL) {
      fprintf(stderr, "malloc failed for idxlist_t_data\n");
      exit(EXIT_FAILURE);
   }

   idxlist_t->vtable = &idxlist_t_vtable;

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              &num_indices, 1, MPI_INT, comm), comm);
   assert(num_indices >= 0);
   idxlist_t->num_indices = (size_t)num_indices;

   idxlist_t->list = malloc(idxlist_t->num_indices * sizeof(*(idxlist_t->list)));
   if (idxlist_t->list == NULL) {
      fprintf(stderr, "malloc failed for idxlist_t\n");
      exit(EXIT_FAILURE);
   }

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              idxlist_t->list, idxlist_t->num_indices,
                              Xt_idx_dt, comm),comm);

   return (xt_idxvec)idxlist_t;
}

static xt_idxvec
idxlist_t_get_intersection(xt_idxvec idxvec_src, xt_idxvec idxvec_dst) {


  if (idxvec_dst->vtable != &idxlist_t_vtable) return NULL;

  // both vectors are index lists:

  struct xt_idxlist_t *idxlist_t_src, *idxlist_t_dst;

  idxlist_t_src = (struct xt_idxlist_t*)idxvec_src;
  idxlist_t_dst = (struct xt_idxlist_t*)idxvec_dst;

  struct xt_idxlist_t *inter_list;

  inter_list = malloc(1 * sizeof(*inter_list));
  if (inter_list == NULL) {
    fprintf(stderr, "malloc failed for intersection list\n");
    exit(EXIT_FAILURE);
  }

  inter_list->num_indices = 0;
  inter_list->list =
    malloc(idxlist_t_dst->num_indices * sizeof(*(inter_list->list)));
  if (inter_list->list == NULL) {
    fprintf(stderr, "malloc failed for intersection list list\n");
    exit(EXIT_FAILURE);
  }

  size_t i, j;

  for (j = 0; j < idxlist_t_dst->num_indices; ++j)
    for (i = 0; i < idxlist_t_src->num_indices; ++i)
      if (idxlist_t_src->list[i] == idxlist_t_dst->list[j]) {
        inter_list->list[inter_list->num_indices++] = idxlist_t_src->list[i];
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

  inter_list->vtable = &idxlist_t_vtable;

  return (xt_idxvec)inter_list;
}

static xt_idxvec
idxlist_t_copy(xt_idxvec idxvec) {

   struct xt_idxlist_t *idxlist_t_data;

   idxlist_t_data = (struct xt_idxlist_t *)idxvec;

   return xt_idxlist_t_new(idxlist_t_data->list, idxlist_t_data->num_indices);
}

static size_t
idxlist_t_get_num_indices(xt_idxvec idxvec) {

   struct xt_idxlist_t *idxlist_t_data;

   idxlist_t_data = (struct xt_idxlist_t *)idxvec;

   return idxlist_t_data->num_indices;
}

static void
idxlist_t_get_indices(xt_idxvec idxvec, Xt_idx *indices) {

   struct xt_idxlist_t *idxlist_t_data;

   idxlist_t_data = (struct xt_idxlist_t *)idxvec;

   memcpy(indices, idxlist_t_data->list,
          idxlist_t_data->num_indices * sizeof(*indices));
}

static int
idxlist_t_get_index(Xt_idxvec idxvec, int position, Xt_idx * index) {

  struct xt_idxlist_t *idxlist_t_data;

  idxlist_t_data = (struct xt_idxlist_t *)idxvec;

  if (position < 0 || position >= idxlist_t_data->num_indices)
    return 1;

  *index = idxlist_t_data->list[position];

  return 0;
}

static int
idxlist_t_get_position_of_index_off(Xt_idxvec idxvec, Xt_idx index,
                                  int * position, int offset) {

  struct xt_idxlist_t *idxlist_t_data;

  idxlist_t_data = (struct xt_idxlist_t *)idxvec;

  size_t i;

  assert(offset >= 0);

  for (i = offset; i < idxlist_t_data->num_indices; ++i)
    if (idxlist_t_data->list[i] == index) {
      *position = i;
      return 0;
    }

  return 1;
}

static int
idxlist_t_get_position_of_index(Xt_idxvec idxvec, Xt_idx index, int * position) {

  return idxlist_t_get_position_of_index_off(idxvec, index, position, 0);
}
