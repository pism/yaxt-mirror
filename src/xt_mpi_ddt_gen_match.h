/**
 * @file xt_mpi_ddt_gen_match.h
 * Simple matching generator routines for MPI datatypes.
 *
 * @copyright Copyright  (C)  2025 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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


#define XT_TOKEN_PASTE2_(a,b) a##b
#define XT_TOKEN_PASTE2(a,b) XT_TOKEN_PASTE2_(a,b)
#define XT_TOKEN_PASTE3_(a,b,c) a##b##c
#define XT_TOKEN_PASTE3(a,b,c) XT_TOKEN_PASTE3_(a,b,c)

#define XT_GEN_MPI_DATATYPE_SIMPLE \
  XT_TOKEN_PASTE2(gen_mpi_datatype_simple,XT_MPI_DDT_GEN_MATCH_SUFFIX)
#define XT_GEN_MPI_DATATYPE_CONTIGUOUS \
  XT_TOKEN_PASTE2(gen_mpi_datatype_contiguous,XT_MPI_DDT_GEN_MATCH_SUFFIX)
#define XT_GEN_MPI_DATATYPE_VECTOR \
  XT_TOKEN_PASTE2(gen_mpi_datatype_vector,XT_MPI_DDT_GEN_MATCH_SUFFIX)
#define XT_GEN_MPI_DATATYPE_INDEXED_BLOCK \
  XT_TOKEN_PASTE2(gen_mpi_datatype_indexed_block,XT_MPI_DDT_GEN_MATCH_SUFFIX)
#define XT_GEN_MPI_DATATYPE_INDEXED \
  XT_TOKEN_PASTE2(gen_mpi_datatype_indexed,XT_MPI_DDT_GEN_MATCH_SUFFIX)
#define XT_CHECK_FOR_VECTOR_TYPE \
  XT_TOKEN_PASTE2(check_for_vector_type,XT_MPI_DDT_GEN_MATCH_SUFFIX)
#define XT_MPI_GENERATE_DATATYPE_BLOCK \
  XT_TOKEN_PASTE2(xt_mpi_generate_datatype_block_match,XT_MPI_DDT_GEN_MATCH_SUFFIX)

static MPI_Datatype
XT_GEN_MPI_DATATYPE_SIMPLE(XT_MPI_DDT_GEN_MATCH_AOFS_TYPE displacement,
                           MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;

  XT_MPI_TYPE_CREATE_INDEXED_BLOCK(1, 1, &displacement, old_type,
                                   &datatype, comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
XT_GEN_MPI_DATATYPE_CONTIGUOUS(XT_MPI_DDT_GEN_MATCH_AOFS_TYPE displacement,
                               int blocklength,
                               MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;
  if (displacement == 0)
    xt_mpi_call(MPI_Type_contiguous(blocklength, old_type, &datatype),
                comm);
  else
    XT_MPI_TYPE_CREATE_INDEXED_BLOCK(1, blocklength,
                                     &displacement, old_type,
                                     &datatype, comm);

  xt_mpi_call(MPI_Type_commit(&datatype), comm);
  return datatype;
}

static MPI_Datatype
XT_GEN_MPI_DATATYPE_VECTOR(int count, int blocklength,
                           XT_MPI_DDT_GEN_MATCH_AOFS_TYPE stride,
                           XT_MPI_DDT_GEN_MATCH_AOFS_TYPE offset,
                           MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;
  xt_mpi_call(XT_MPI_TYPE_CREATE_VECTOR(count, blocklength, stride, old_type,
                                        &datatype), comm);
  if (offset != 0) {

    MPI_Datatype datatype_;
    int hindexed_blocklength = 1;
#ifdef XT_MPI_DDT_GEN_MATCH_COMPUTE_VECTOR_DISPLACEMENT
    MPI_Aint old_type_size, old_type_lb;

    xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                    &old_type_size), comm);

    MPI_Aint displacement = offset * old_type_size;
#else
    MPI_Aint displacement = offset;
#endif
    Xt_Type_create_hindexed(1, &hindexed_blocklength,
                            &displacement, datatype, &datatype_, comm);
    xt_mpi_call(MPI_Type_free(&datatype), comm);
    datatype = datatype_;
  }
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
XT_GEN_MPI_DATATYPE_INDEXED_BLOCK(
  const XT_MPI_DDT_GEN_MATCH_AOFS_TYPE *displacements,
  int blocklength, int count, MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;
  XT_MPI_TYPE_CREATE_INDEXED_BLOCK(count, blocklength,
                                   displacements,
                                   old_type, &datatype, comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);
  return datatype;
}

static MPI_Datatype
XT_GEN_MPI_DATATYPE_INDEXED(
  const XT_MPI_DDT_GEN_MATCH_AOFS_TYPE *displacements,
  const int *blocklengths,
  int count, MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;
  XT_MPI_TYPE_INDEXED(count, blocklengths,
                      displacements,
                      old_type, &datatype, comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);
  return datatype;
}

static inline int
XT_CHECK_FOR_VECTOR_TYPE(
  const XT_MPI_DDT_GEN_MATCH_AOFS_TYPE *displacements,
  const int *blocklengths,
  int count) {

  int blocklength = blocklengths[0];
  XT_MPI_DDT_GEN_MATCH_AOFS_TYPE stride = displacements[1] - displacements[0];

  for (int i = 1; i < count; ++i)
    if ((blocklengths[i] != blocklength)
        | (displacements[i] - displacements[i-1] != stride))
      return 0;
  return 1;
}

#ifndef XT_MPI_DDT_GEN_MATCH_H_INCLUDED
static inline int
check_for_indexed_block_type(const int *blocklengths, int count) {

  int blocklength = blocklengths[0];

  for (int i = 1; i < count; ++i)
    if (blocklengths[i] != blocklength)
      return 0;

  return 1;
}

static MPI_Datatype
copy_mpi_datatype(MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_dup(old_type, &datatype), comm);

  return datatype;
}

#define XT_MPI_DDT_GEN_MATCH_H_INCLUDED
#endif

MPI_Datatype
XT_MPI_GENERATE_DATATYPE_BLOCK(
  const XT_MPI_DDT_GEN_MATCH_AOFS_TYPE *displacements,
  const int *blocklengths,
  int count,
  MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;

  if (count == 0)
    datatype = MPI_DATATYPE_NULL;
  else if (count == 1 && blocklengths[0] == 1 && displacements[0] == 0)
    datatype = copy_mpi_datatype(old_type, comm);
  else if (count == 1 && blocklengths[0] == 1)
    datatype = XT_GEN_MPI_DATATYPE_SIMPLE(displacements[0], old_type, comm);
  else if (count == 1)
    datatype = XT_GEN_MPI_DATATYPE_CONTIGUOUS(displacements[0], blocklengths[0],
                                              old_type, comm);
  else if (XT_CHECK_FOR_VECTOR_TYPE(displacements, blocklengths, count))
    datatype = XT_GEN_MPI_DATATYPE_VECTOR(count, blocklengths[0],
                                       displacements[1] - displacements[0],
                                       displacements[0], old_type, comm);
  else if (check_for_indexed_block_type(blocklengths, count))
    datatype = XT_GEN_MPI_DATATYPE_INDEXED_BLOCK(displacements, blocklengths[0],
                                                 count, old_type, comm);
  else
    datatype = XT_GEN_MPI_DATATYPE_INDEXED(displacements, blocklengths, count,
                                           old_type, comm);

  return datatype;
}

#undef XT_TOKEN_PASTE2_
#undef XT_TOKEN_PASTE2
#undef XT_TOKEN_PASTE3_
#undef XT_TOKEN_PASTE3

#undef XT_GEN_MPI_DATATYPE_SIMPLE
#undef XT_GEN_MPI_DATATYPE_CONTIGUOUS
#undef XT_GEN_MPI_DATATYPE_VECTOR
#undef XT_GEN_MPI_DATATYPE_INDEXED_BLOCK
#undef XT_GEN_MPI_DATATYPE_INDEXED
#undef XT_CHECK_FOR_VECTOR_TYPE
#undef XT_MPI_GENERATE_DATATYPE_BLOCK


/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
