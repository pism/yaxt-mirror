#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"

//taken from http://beige.ucs.indiana.edu/I590/node85.html
void mpi_err_handler(int error_code, MPI_Comm comm) {


   if (error_code != MPI_SUCCESS) {

      int rank;

      MPI_Comm_rank(comm, &rank);

      char error_string[1024];
      int length_of_error_string, error_class;

      MPI_Error_class(error_code, &error_class);
      MPI_Error_string(error_class, error_string, &length_of_error_string);
      fprintf(stderr, "%3d: %s\n", rank, error_string);
      MPI_Error_string(error_code, error_string, &length_of_error_string);
      fprintf(stderr, "%3d: %s\n", rank, error_string);
      MPI_Abort(comm, error_code);
   }
}

static MPI_Datatype copy_mpi_datatype(MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  mpi_err_handler(MPI_Type_dup(old_type, &datatype), comm);

  return datatype;
}

static MPI_Datatype gen_mpi_datatype_simple(int displacement, MPI_Datatype old_type,
                                            MPI_Comm comm) {

  MPI_Datatype datatype;

  mpi_err_handler(MPI_Type_create_indexed_block(1, 1, &displacement, old_type,
                                                &datatype), comm);
  mpi_err_handler(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype gen_mpi_datatype_contiguous(int displacement, int blocklength,
                                                MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  if (displacement == 0)
    mpi_err_handler(MPI_Type_contiguous(blocklength, old_type, &datatype), comm);
  else
    mpi_err_handler(MPI_Type_create_indexed_block(1, blocklength, &displacement, old_type,
                                                  &datatype), comm);

  mpi_err_handler(MPI_Type_commit(&datatype), comm);

  return datatype;

}

static MPI_Datatype gen_mpi_datatype_vector(int stride, int blocklength, int count,
                                            int offset, MPI_Datatype old_type,
                                            MPI_Comm comm) {

  MPI_Datatype datatype;

  mpi_err_handler(MPI_Type_vector(count, blocklength, stride, old_type, &datatype), comm);
  if (offset != 0) {

    MPI_Datatype datatype_;
    int blocklength = 1;
    MPI_Aint old_type_size, old_type_lb;

    mpi_err_handler(MPI_Type_get_extent(old_type, &old_type_lb, &old_type_size), comm);

    MPI_Aint displacement = offset * old_type_size;

    mpi_err_handler(MPI_Type_create_hindexed(1, &blocklength, &displacement,
                                             datatype, &datatype_),
                    comm);
    mpi_err_handler(MPI_Type_free(&datatype), comm);
    datatype = datatype_;
  }
  mpi_err_handler(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype gen_mpi_datatype_indexed_block(int const * displacements, int blocklength,
                                                   int count, MPI_Datatype old_type,
                                                   MPI_Comm comm) {

  MPI_Datatype datatype;

  mpi_err_handler(MPI_Type_create_indexed_block(count, blocklength, (void*)displacements,
                                                old_type, &datatype), comm);
  mpi_err_handler(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype gen_mpi_datatype_indexed(int const * displacements, int * blocklengths,
                                             int count, MPI_Datatype old_type,
                                             MPI_Comm comm) {

  MPI_Datatype datatype;

  mpi_err_handler(MPI_Type_indexed(count, blocklengths, (void*)displacements, old_type,
                                   &datatype), comm);
  mpi_err_handler(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static inline int check_for_vector_type(int const * displacements, int * blocklengths,
                                        int count) {

  int blocklength = blocklengths[0];

  for (int i = 1; i < count; ++i)
    if (blocklengths[i] != blocklength)
      return 0;

  int stride = displacements[1] - displacements[0];

  for (int i = 1; i + 1 < count; ++i)
    if (displacements[i+1] - displacements[i] != stride)
      return 0;

  return 1;
}

static inline int check_for_indexed_block_type(int * blocklengths, int count) {

  int blocklength = blocklengths[0];

  for (int i = 1; i < count; ++i)
    if (blocklengths[i] != blocklength)
      return 0;

  return 1;
}

MPI_Datatype xt_mpi_generate_datatype_block(int const * displacements, int * blocklengths,
                                            int count, MPI_Datatype old_type,
                                            MPI_Comm comm) {

  MPI_Datatype datatype;

  if (count == 0)
    datatype = MPI_DATATYPE_NULL;
  else if (count == 1 && blocklengths[0] == 1 && displacements[0] == 0)
    datatype = copy_mpi_datatype(old_type, comm);
  else if (count == 1 && blocklengths[0] == 1)
    datatype = gen_mpi_datatype_simple(displacements[0], old_type, comm);
  else if (count == 1)
    datatype = gen_mpi_datatype_contiguous(displacements[0], blocklengths[0],
                                           old_type, comm);
  else if (check_for_vector_type(displacements, blocklengths, count))
    datatype = gen_mpi_datatype_vector(displacements[1]-displacements[0], blocklengths[0],
                                       count, displacements[0], old_type, comm);
  else if (check_for_indexed_block_type(blocklengths, count))
    datatype = gen_mpi_datatype_indexed_block(displacements, blocklengths[0],
                                              count, old_type, comm);
  else
    datatype = gen_mpi_datatype_indexed(displacements, blocklengths, count,
                                        old_type, comm);

  return datatype;
}

MPI_Datatype xt_mpi_generate_datatype(int const * displacements, int count,
                                      MPI_Datatype old_type, MPI_Comm comm) {

  if (count == 0)
    return MPI_DATATYPE_NULL;

  int * blocklengths;

  blocklengths = xmalloc(count * sizeof(*blocklengths));

  int new_count = 0;
  int i = 0;

  do {

    int j = 1;
  
    while (i + j < count && displacements[i] + j == displacements[i + j]) ++j;    

    blocklengths[new_count++] = j;

    i += j;
  } while (i < count);

  int * tmp_displ = NULL;
  int const * displ;

  if (new_count != count) {

    tmp_displ = xmalloc(new_count * sizeof(*tmp_displ));

    int offset = 0;

    for (i = 0; i < new_count; ++i) {

      tmp_displ[i] = displacements[offset];
      offset += blocklengths[i];
    }

    displ = tmp_displ;
  } else
    displ = displacements;

  MPI_Datatype datatype;

  datatype = xt_mpi_generate_datatype_block(displ, blocklengths, new_count,
                                            old_type, comm);

  free(blocklengths);

  free(tmp_displ);

  return datatype;
}
