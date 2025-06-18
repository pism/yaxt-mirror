/**
 * @file xt_mpi_ddt_gen.c
 *
 * @copyright Copyright  (C)  2023 Jörg Behrens <behrens@dkrz.de>
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdlib.h>

#include <mpi.h>
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_config_internal.h"
#include "xt_mpi_internal.h"
#include "xt_mpi_ddt_cache.h"

#include <xt_mpi_ddt_wrap.h>

#define XT_MPI_DDT_GEN_MATCH_SUFFIX _i
#define XT_MPI_DDT_GEN_MATCH_AOFS_TYPE int
#define XT_MPI_TYPE_CREATE_INDEXED_BLOCK Xt_Type_create_indexed_block
#define XT_MPI_TYPE_CREATE_VECTOR MPI_Type_vector
#define XT_MPI_DDT_GEN_MATCH_COMPUTE_VECTOR_DISPLACEMENT
#define XT_MPI_TYPE_INDEXED Xt_Type_indexed
#include <xt_mpi_ddt_gen_match.h>
#undef XT_MPI_DDT_GEN_MATCH_AOFS_TYPE
#undef XT_MPI_DDT_GEN_MATCH_SUFFIX
#undef XT_MPI_TYPE_CREATE_INDEXED_BLOCK
#undef XT_MPI_TYPE_CREATE_VECTOR
#undef XT_MPI_DDT_GEN_MATCH_COMPUTE_VECTOR_DISPLACEMENT
#undef XT_MPI_TYPE_INDEXED

#define XT_MPI_DDT_GEN_MATCH_SUFFIX _a
#define XT_MPI_DDT_GEN_MATCH_AOFS_TYPE MPI_Aint
#define XT_MPI_TYPE_CREATE_INDEXED_BLOCK Xt_Type_create_hindexed_block
#define XT_MPI_TYPE_CREATE_VECTOR MPI_Type_create_hvector
#define XT_MPI_TYPE_INDEXED Xt_Type_create_hindexed
#include <xt_mpi_ddt_gen_match.h>
#undef XT_MPI_DDT_GEN_MATCH_AOFS_TYPE
#undef XT_MPI_DDT_GEN_MATCH_SUFFIX
#undef XT_MPI_TYPE_CREATE_INDEXED_BLOCK
#undef XT_MPI_TYPE_CREATE_VECTOR
#undef XT_MPI_DDT_GEN_MATCH_COMPUTE_VECTOR_DISPLACEMENT
#undef XT_MPI_TYPE_INDEXED

static struct Xt_mpi_dt_ofs
xt_mpi_generate_compact_datatype_block(const int *disp, const int *blocklengths,
                                       int count, MPI_Datatype old_type,
                                       MPI_Comm comm, Xt_config config);

static struct Xt_mpi_dt_ofs
xt_mpi_generate_datatype_block_custom(const int *displacements,
                                      const int *blocklengths,
                                      int count, MPI_Datatype old_type,
                                      MPI_Comm comm, Xt_config config);

MPI_Datatype
xt_mpi_generate_datatype_block(const int *displacements,
                               const int *blocklengths,
                               int count, MPI_Datatype old_type,
                               MPI_Comm comm) {
  struct Xt_mpi_dt_ofs dt_ofs = xt_mpi_generate_datatype_block_custom(
    displacements, blocklengths, count, old_type, comm, &xt_default_config);
  if (dt_ofs.displacement != 0) {
    MPI_Datatype datatype = dt_ofs.dt;
    Xt_Type_create_hindexed(1, (int [1]){ 1 },
                            &dt_ofs.displacement, datatype, &dt_ofs.dt, comm);
    if (!xt_mpi_ddt_is_named(datatype))
      xt_mpi_call(MPI_Type_free(&datatype), comm);
    xt_mpi_call(MPI_Type_commit(&dt_ofs.dt), comm);
  } else if (dt_ofs.dt != MPI_DATATYPE_NULL && xt_mpi_ddt_is_named(dt_ofs.dt))
    xt_mpi_call(MPI_Type_dup(dt_ofs.dt, &dt_ofs.dt), comm);
  return dt_ofs.dt;
}

static struct Xt_mpi_dt_ofs
xt_mpi_generate_datatype_block_custom(const int *displacements,
                                      const int *blocklengths,
                                      int count, MPI_Datatype old_type,
                                      MPI_Comm comm, Xt_config config) {

  if (config->ddt_cache_limit >= -1)
    return xt_mpi_generate_compact_datatype_block(displacements, blocklengths,
                                                  count, old_type, comm, config);
  else
    return xt_mpi_generate_datatype_block_match_i(displacements, blocklengths,
                                                  count, old_type, comm);
}

static struct Xt_mpi_dt_ofs
xt_mpi_generate_datatype_custom(int const * displacements, int count,
                                MPI_Datatype old_type, MPI_Comm comm,
                                Xt_config config);

MPI_Datatype
xt_mpi_generate_datatype(int const * displacements, int count,
                         MPI_Datatype old_type, MPI_Comm comm)
{
  if (count <= 0)
    return MPI_DATATYPE_NULL;
  struct Xt_mpi_dt_ofs dt_ofs
    = xt_mpi_generate_datatype_custom(displacements, count, old_type,
                                      comm, &xt_default_config);
  if (dt_ofs.displacement != 0) {
    MPI_Datatype datatype = dt_ofs.dt;
    Xt_Type_create_hindexed(1, (int [1]){ 1 },
                            &dt_ofs.displacement, datatype, &dt_ofs.dt, comm);
    if (!xt_mpi_ddt_is_named(datatype))
      xt_mpi_call(MPI_Type_free(&datatype), comm);
    xt_mpi_call(MPI_Type_commit(&dt_ofs.dt), comm);
  } else if (dt_ofs.dt != MPI_DATATYPE_NULL && xt_mpi_ddt_is_named(dt_ofs.dt))
    xt_mpi_call(MPI_Type_dup(dt_ofs.dt, &dt_ofs.dt), comm);
  return dt_ofs.dt;
}

static struct Xt_mpi_dt_ofs
xt_mpi_generate_compact_datatype(int const *disp, int disp_len,
                                 MPI_Datatype old_type, MPI_Comm comm,
                                 Xt_config config);

struct Xt_mpi_dt_ofs
xt_mpi_generate_datatype_custom(int const * displacements, int count,
                                MPI_Datatype old_type, MPI_Comm comm,
                                Xt_config config)
{
  if (config->ddt_cache_limit >= -1)
    return xt_mpi_generate_compact_datatype(displacements, count, old_type,
                                            comm, config);
  size_t num_blk = xt_disp2blk_count((size_t)count, displacements);
  int *restrict blocklengths = xmalloc((size_t)num_blk * sizeof(*blocklengths)
                                       * (size_t)(1 + (num_blk != (size_t)count)));
  const int *displ;
  if (num_blk != (size_t)count) {
    xt_disp2blk((size_t)count, displacements,
                blocklengths + num_blk, blocklengths);
    displ = blocklengths + num_blk;
  } else {
    for (size_t i = 0; i < num_blk; ++i)
      blocklengths[i] = 1;
    displ = displacements;
  }
  struct Xt_mpi_dt_ofs dt_ofs = xt_mpi_generate_datatype_block_match_i(
    displ, blocklengths, (int)num_blk, old_type, comm);

  free(blocklengths);

  return dt_ofs;
}

#define XT_MPI_STRP_PRS_PREFIX
#define XT_MPI_STRP_PRS_UNITSTRIDE 1
#define XT_MPI_STRP_PRS_AOFS_TYPE int
#define XT_MPI_STRP_PRS_DISP_ADJUST(val) ((val) * params->old_type_extent)
#define XT_MPI_STRP_PRS_BLOCK_VEC_CREATE Xt_mpi_ddt_cache_acquire_vector
#define XT_MPI_STRP_PRS_INDEXED_BLOCK_CREATE \
  Xt_mpi_ddt_cache_acquire_indexed_block
#define XT_MPI_STRP_PRS_INDEXED_CREATE Xt_mpi_ddt_cache_acquire_indexed
#include "xt_mpi_stripe_parse_func.h"
#undef XT_MPI_STRP_PRS_PREFIX
#undef XT_MPI_STRP_PRS_UNITSTRIDE
#undef XT_MPI_STRP_PRS_AOFS_TYPE
#undef XT_MPI_STRP_PRS_DISP_ADJUST
#undef XT_MPI_STRP_PRS_BLOCK_VEC_CREATE
#undef XT_MPI_STRP_PRS_INDEXED_BLOCK_CREATE
#undef XT_MPI_STRP_PRS_INDEXED_CREATE


#define XT_MPI_STRP_PRS_PREFIX a
#define XT_MPI_STRP_PRS_UNITSTRIDE params->old_type_extent
#define XT_MPI_STRP_PRS_AOFS_TYPE MPI_Aint
#define XT_MPI_STRP_PRS_DISP_ADJUST(val) (val)
#define XT_MPI_STRP_PRS_BLOCK_VEC_CREATE Xt_mpi_ddt_cache_acquire_hvector
#define XT_MPI_STRP_PRS_INDEXED_BLOCK_CREATE \
  Xt_mpi_ddt_cache_acquire_hindexed_block
#define XT_MPI_STRP_PRS_INDEXED_CREATE Xt_mpi_ddt_cache_acquire_hindexed
#include "xt_mpi_stripe_parse_func.h"

static struct Xt_mpi_dt_ofs
xt_mpi_ddt_block_gen(size_t count, const int *disp, const int *blocklengths,
                     struct Xt_mpi_strp_prs_params *params);

static struct Xt_mpi_dt_ofs
xt_mpi_generate_compact_datatype_block(const int *disp, const int *blocklengths,
                                       int count, MPI_Datatype old_type,
                                       MPI_Comm comm, Xt_config config)
{
  if (count < 1)
    return (struct Xt_mpi_dt_ofs){ .dt = MPI_DATATYPE_NULL, .displacement = 0 };
  struct Xt_mpi_strp_prs_params params;
  xt_init_mpi_strp_prs_params(&params, old_type, config, comm);
  struct Xt_mpi_dt_ofs dt_ofs
    = xt_mpi_ddt_block_gen((size_t)count, disp, blocklengths, &params);
  xt_destroy_mpi_strp_prs_params(&params);
  return dt_ofs;
}

static struct Xt_mpi_dt_ofs
xt_mpi_ddt_block_gen(size_t count, const int *disp, const int *blocklengths,
                     struct Xt_mpi_strp_prs_params *params)
{
  size_t count_ = (size_t)0;
  for (size_t i=0; i < count; ++i)
    count_ += (size_t)(blocklengths[i] > 0);
  int *restrict v_disp, *restrict v_size, *restrict v_stride;
  v_stride = xmalloc(sizeof(int) * count_ * (count_ != count ? 3 : 1));
  if (count_ != count) {
    v_size = v_stride + count_;
    v_disp = v_size + count_;
    size_t j=0;
    for (size_t i=0; i < count; ++i) {
      v_disp[j] = disp[i];
      v_stride[j] = 1;
      int bl = blocklengths[i];
      v_size[j] = bl;
      j += (size_t)(bl > 0);
    }
    disp = v_disp;
    blocklengths = v_size;
  } else {
    for (size_t i = 0; i < count_; ++i)
      v_stride[i] = 1;
  }
  struct Xt_mpi_dt_ofs dt_ofs
    = parse_stripe(disp, blocklengths, v_stride, count_, params);
  free(v_stride);
  return dt_ofs;
}

struct Xt_mpi_dt_ofs
xt_mpi_generate_compact_datatype(const int *disp, int disp_len,
                                 MPI_Datatype old_type, MPI_Comm comm,
                                 Xt_config config)
{
  if (disp_len < 1)
    return (struct Xt_mpi_dt_ofs){ .dt = MPI_DATATYPE_NULL,
                                   .displacement = 0 };
  size_t vlen = xt_disp2ext_count((size_t)disp_len, disp);
  int *v_disp = xmalloc(vlen * 3U * sizeof (int)),
    *v_size = v_disp + vlen, *v_stride = v_size + vlen;
  xt_disp2ext((size_t)disp_len, disp, v_disp, v_size, v_stride);
  struct Xt_mpi_strp_prs_params params;
  xt_init_mpi_strp_prs_params(&params, old_type, config, comm);
  struct Xt_mpi_dt_ofs dt_ofs
    = parse_stripe(v_disp, v_size, v_stride, vlen, &params);
  xt_destroy_mpi_strp_prs_params(&params);
  free(v_disp);
  return dt_ofs;
}



/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
