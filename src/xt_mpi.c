/**
 * @file xt_mpi.c
 *
 * @copyright Copyright  (C)  2012 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://redmine.dkrz.de/doc/yaxt/html/index.html
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"

//taken from http://beige.ucs.indiana.edu/I590/node85.html
void xt_mpi_error(int error_code, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  char error_string[MPI_MAX_ERROR_STRING];
  int length_of_error_string, error_class;

  MPI_Error_class(error_code, &error_class);
  MPI_Error_string(error_class, error_string, &length_of_error_string);
  fprintf(stderr, "%3d: %s\n", rank, error_string);
  MPI_Error_string(error_code, error_string, &length_of_error_string);
  fprintf(stderr, "%3d: %s\n", rank, error_string);
  MPI_Abort(comm, error_code);
}

static MPI_Datatype copy_mpi_datatype(MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_dup(old_type, &datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_simple(int displacement, MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_create_indexed_block(1, 1, &displacement, old_type,
                                                &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_contiguous(int displacement, int blocklength,
                            MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  if (displacement == 0)
    xt_mpi_call(MPI_Type_contiguous(blocklength, old_type, &datatype),
                    comm);
  else
    xt_mpi_call(MPI_Type_create_indexed_block(1, blocklength,
                                                  &displacement, old_type,
                                                  &datatype), comm);

  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;

}

static MPI_Datatype
gen_mpi_datatype_vector(int count, int blocklength, int stride,
                        int offset, MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_vector(count, blocklength, stride, old_type,
                                  &datatype), comm);
  if (offset != 0) {

    MPI_Datatype datatype_;
    int blocklength = 1;
    MPI_Aint old_type_size, old_type_lb;

    xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                        &old_type_size), comm);

    MPI_Aint displacement = offset * old_type_size;

    xt_mpi_call(MPI_Type_create_hindexed(1, &blocklength, &displacement,
                                             datatype, &datatype_),
                    comm);
    xt_mpi_call(MPI_Type_free(&datatype), comm);
    datatype = datatype_;
  }
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_indexed_block(int const * displacements, int blocklength,
                               int count, MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_create_indexed_block(count, blocklength,
                                                (void *)displacements,
                                                old_type, &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_indexed(int const * displacements, int * blocklengths,
                         int count, MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_indexed(count, blocklengths, (void*)displacements,
                                   old_type, &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static inline int
check_for_vector_type(int const * displacements, int * blocklengths,
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

MPI_Datatype
xt_mpi_generate_datatype_block(int const * displacements, int * blocklengths,
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
    datatype = gen_mpi_datatype_vector(count, blocklengths[0],
                                       displacements[1] - displacements[0],
                                       displacements[0], old_type, comm);
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

  if (count <= 0)
    return MPI_DATATYPE_NULL;

  int * blocklengths = xmalloc((size_t)count * sizeof(*blocklengths));

  int new_count = 0;
  {
    int i = 0;
    do {
      int j = 1;
      while (i + j < count && displacements[i] + j == displacements[i + j])
        ++j;
      blocklengths[new_count++] = j;
      i += j;
    } while (i < count);
  }
  int * tmp_displ = NULL;
  const int *displ;

  if (new_count != count) {

    tmp_displ = xmalloc((size_t)new_count * sizeof(*tmp_displ));

    int offset = 0;

    for (unsigned i = 0; i < (unsigned)new_count; ++i) {

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
