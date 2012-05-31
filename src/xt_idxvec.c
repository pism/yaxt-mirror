#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "yaxt.h"
#include "xt_idxlist.h"
#include "xt_idxvec.h"
#include "xt_mpi.h"

static void
idxvec_delete(xt_idxlist data);

static size_t
idxvec_get_pack_size(xt_idxlist data, MPI_Comm comm);

static void
idxvec_pack(xt_idxlist data, void *buffer, int buffer_size,
            int *position, MPI_Comm comm);

static xt_idxlist
idxvec_get_intersection(xt_idxlist idxlist_a, xt_idxlist idxlist_b);

static xt_idxlist
idxvec_copy(xt_idxlist idxlist);

static size_t
idxvec_get_num_indices(xt_idxlist idxlist);

static void
idxvec_get_indices(xt_idxlist idxlist, Xt_idx *indices);

static int
idxvec_get_index(Xt_idxlist idxlist, int position, Xt_idx * index);

static int
idxvec_get_position_of_index(Xt_idxlist idxlist, Xt_idx index, int * position);

static int
idxvec_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                 int * position, int offset);

static const struct xt_idxlist_vtable idxvec_vtable = {
              idxvec_delete,
              idxvec_get_pack_size,
              idxvec_pack,
              idxvec_get_intersection,
              idxvec_copy,
              idxvec_get_num_indices,
              idxvec_get_indices,
              idxvec_get_index,
              idxvec_get_position_of_index,
              idxvec_get_position_of_index_off};

// index vector data structure
struct xt_idxvec {

   const struct xt_idxlist_vtable *vtable;

   Xt_idx *vector;
   size_t num_indices;
};

static int compare_xtidx(const void * a, const void * b)
{
  return (*(Xt_idx*)a > *(Xt_idx*)b) - (*(Xt_idx*)a < *(Xt_idx*)b);
}

struct xt_idxlist *xt_idxvec_new(Xt_idx *idxvec, size_t num_indices) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = malloc(1 * sizeof(*idxvec_data));
   if (idxvec_data == NULL) {
      fprintf(stderr, "malloc failed for idxvec_data\n");
      exit(EXIT_FAILURE);
   }

   idxvec_data->vtable = &idxvec_vtable;

   if (num_indices != 0) {

      size_t vector_size;

      vector_size = num_indices * sizeof(*(idxvec_data->vector));

      idxvec_data->vector = malloc(vector_size);
      if (idxvec_data->vector == NULL && vector_size != 0) {
         fprintf(stderr, "malloc failed for idxvec_data\n");
         exit(EXIT_FAILURE);
      }

      memcpy(idxvec_data->vector, idxvec, vector_size);

   } else {

      idxvec_data->vector = NULL;
   }

   idxvec_data->num_indices = num_indices;

   return (void*)idxvec_data;
}

static void idxvec_delete(struct xt_idxlist *data) {

   free(((struct xt_idxvec*)data)->vector);
   free(data);
}

static size_t idxvec_get_pack_size(struct xt_idxlist *data, MPI_Comm comm) {

   struct xt_idxvec *idxvec = (struct xt_idxvec*)data;

   int size_unsigned, size_int, size_int_type;

   mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_int_type), comm);
   mpi_err_handler(MPI_Pack_size(idxvec->num_indices, Xt_idx_dt, comm,
                                 &size_int), comm);
   mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &size_unsigned), comm);

   return size_unsigned + size_int + size_int_type;
}

void idxvec_pack(struct xt_idxlist *data, void *buffer, int buffer_size,
                 int *position, MPI_Comm comm) {

  assert(data);
  struct xt_idxvec *idxvec = (struct xt_idxvec*)data;
  int type = VECTOR;
  int num_indices = (int)(idxvec->num_indices);

  assert(idxvec->num_indices <= INT_MAX);

  mpi_err_handler(MPI_Pack(&(type), 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  mpi_err_handler(MPI_Pack(&num_indices, 1, MPI_INT, buffer,
                           buffer_size, position, comm), comm);
  if (num_indices != 0)
    mpi_err_handler(MPI_Pack(idxvec->vector, (int)idxvec->num_indices,
                             Xt_idx_dt, buffer,
                             buffer_size, position, comm), comm);
}

xt_idxlist xt_idxvec_unpack(void *buffer, int buffer_size, int *position,
                            MPI_Comm comm) {

   struct xt_idxvec *idxvec;
   int num_indices;

   idxvec = malloc(1 * sizeof(*idxvec));
   if (idxvec == NULL) {
      fprintf(stderr, "malloc failed for idxvec_data\n");
      exit(EXIT_FAILURE);
   }

   idxvec->vtable = &idxvec_vtable;

   mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                              &num_indices, 1, MPI_INT, comm), comm);
   assert(num_indices >= 0);
   idxvec->num_indices = (size_t)num_indices;

   if (num_indices != 0) {

      idxvec->vector = malloc(idxvec->num_indices * sizeof(*(idxvec->vector)));
      if (idxvec->vector == NULL) {
         fprintf(stderr, "malloc failed for idxvec\n");
         exit(EXIT_FAILURE);
      }
      mpi_err_handler(MPI_Unpack(buffer, buffer_size, position,
                                 idxvec->vector, idxvec->num_indices,
                                 Xt_idx_dt, comm),comm);
   } else {

      idxvec->vector = NULL;
   }

   return (xt_idxlist)idxvec;
}

static xt_idxlist
idxvec_get_intersection(xt_idxlist idxlist_src, xt_idxlist idxlist_dst) {


  if (idxlist_dst->vtable != &idxvec_vtable) return NULL;

  // both lists are index vectors:

  struct xt_idxvec *idxvec_src, *idxvec_dst;

  idxvec_src = (struct xt_idxvec*)idxlist_src;
  idxvec_dst = (struct xt_idxvec*)idxlist_dst;

  struct xt_idxvec *inter_vector;

  inter_vector = malloc(1 * sizeof(*inter_vector));
  if (inter_vector == NULL) {
    fprintf(stderr, "malloc failed for intersection vector\n");
    exit(EXIT_FAILURE);
  }

  inter_vector->num_indices = 0;
  inter_vector->vector =
    malloc(idxvec_dst->num_indices * sizeof(*(inter_vector->vector)));
  if (inter_vector->vector == NULL) {
    fprintf(stderr, "malloc failed for intersection vector\n");
    exit(EXIT_FAILURE);
  }

  size_t i, j;

  for (j = 0; j < idxvec_dst->num_indices; ++j)
    for (i = 0; i < idxvec_src->num_indices; ++i)
      if (idxvec_src->vector[i] == idxvec_dst->vector[j]) {
        inter_vector->vector[inter_vector->num_indices++] = idxvec_src->vector[i];
         break;
      }

  if (inter_vector->num_indices) {

    // sort the vector
    qsort(inter_vector->vector, inter_vector->num_indices,
          sizeof (inter_vector->vector[0]), compare_xtidx);
    inter_vector->vector = realloc(inter_vector->vector,
                                   inter_vector->num_indices *
                                   sizeof(*(inter_vector->vector)));
    if (inter_vector->vector == NULL) {
      fprintf(stderr, "malloc failed for intersection vector\n");
      exit(EXIT_FAILURE);
    }
  } else {
    free(inter_vector->vector);
    inter_vector->vector = NULL;
  }

  inter_vector->vtable = &idxvec_vtable;

  return (xt_idxlist)inter_vector;
}

static xt_idxlist
idxvec_copy(xt_idxlist idxlist) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = (struct xt_idxvec *)idxlist;

   return xt_idxvec_new(idxvec_data->vector, idxvec_data->num_indices);
}

static size_t
idxvec_get_num_indices(xt_idxlist idxlist) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = (struct xt_idxvec *)idxlist;

   return idxvec_data->num_indices;
}

static void
idxvec_get_indices(xt_idxlist idxlist, Xt_idx *indices) {

   struct xt_idxvec *idxvec_data;

   idxvec_data = (struct xt_idxvec *)idxlist;

   memcpy(indices, idxvec_data->vector,
          idxvec_data->num_indices * sizeof(*indices));
}

static int
idxvec_get_index(Xt_idxlist idxlist, int position, Xt_idx * index) {

  struct xt_idxvec *idxvec_data;

  idxvec_data = (struct xt_idxvec *)idxlist;

  if (position < 0 || position >= idxvec_data->num_indices)
    return 1;

  *index = idxvec_data->vector[position];

  return 0;
}

static int
idxvec_get_position_of_index_off(Xt_idxlist idxlist, Xt_idx index,
                                 int * position, int offset) {

  struct xt_idxvec *idxvec_data;

  idxvec_data = (struct xt_idxvec *)idxlist;

  size_t i;

  assert(offset >= 0);

  for (i = offset; i < idxvec_data->num_indices; ++i)
    if (idxvec_data->vector[i] == index) {
      *position = i;
      return 0;
    }

  return 1;
}

static int
idxvec_get_position_of_index(Xt_idxlist idxlist, Xt_idx index, int * position) {

  return idxvec_get_position_of_index_off(idxlist, index, position, 0);
}
