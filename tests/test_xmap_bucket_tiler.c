/**
 * @file test_xmap_bucket_tiler.c
 *
 * @copyright Copyright  (C)  2026 Jörg Behrens <behrens@dkrz.de>
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
#include <config.h>
#endif

#include <assert.h>
#include <string.h>

#include <mpi.h>

#include <yaxt.h>

#include "test_xmap_bucket_tiler.h"

#include "core/ppm_xfuncs.h"


struct bucket_gen_1d_tile_param {
  Xt_int global_sizes[1];
};

struct bucket_gen_1d_tile_state {
  Xt_int global_sizes[1];
  int comm_size, last_rank_generated;
  Xt_idxlist last_tile;
};

static int
bucket_gen_1d_tile_init(
  void *gen_state_,
  Xt_idxlist src_idxlist,
  Xt_idxlist dst_idxlist,
  Xt_config config,
  const struct Xt_xmdd_bucket_gen_comms *comms,
  void *init_params_)
{
  (void)src_idxlist;
  (void)dst_idxlist;
  /* todo: query idxlists for minimum and maximum to assess
     range of intersecting bucket ranks */
  (void)config;
  struct bucket_gen_1d_tile_param *init_params = init_params_;
  struct bucket_gen_1d_tile_state *gen_state = gen_state_;
  memcpy(gen_state->global_sizes,
         init_params->global_sizes,
         sizeof (gen_state->global_sizes));
  MPI_Comm_size(comms->intra_comm,
                &gen_state->comm_size);
  gen_state->last_rank_generated = -1;
  gen_state->last_tile = NULL;
  return 1;
}

static void
bucket_gen_1d_tile_destroy(void *gen_state_)
{
  struct bucket_gen_1d_tile_state *gen_state = gen_state_;
  if (gen_state->last_tile)
    xt_idxlist_delete(gen_state->last_tile);
}

static int
bucket_gen_1d_tile_get_intersect_max_num(
  void *gen_state_, int type)
{
  struct bucket_gen_1d_tile_state *gen_state = gen_state_;
  /* todo: properly handle type */
  (void)type;
  return gen_state->comm_size;
}

static struct Xt_com_list
bucket_gen_1d_tile_next(void *gen_state_, int type)
{
  struct bucket_gen_1d_tile_state *gen_state = gen_state_;
  /* todo: properly handle type */
  (void)type;
  if (gen_state->last_tile) {
    xt_idxlist_delete(gen_state->last_tile);
    gen_state->last_tile = NULL;
  }
  int next_rank = gen_state->last_rank_generated + 1,
    comm_size = gen_state->comm_size;
  while (next_rank < comm_size) {
    const Xt_int *global_sizes = gen_state->global_sizes,
      next_rank_start[1] = { (Xt_int)(next_rank*global_sizes[0] / comm_size) };
    int next_rank_size[1] = {
      (int)((next_rank + 1) * global_sizes[0] / comm_size - next_rank_start[0]) };
    int next_rank_size_prod = next_rank_size[0];
    for (size_t i = 1; i < 1; ++i)
      next_rank_size_prod *= next_rank_size[i];
    if (next_rank_size_prod) {
      gen_state->last_tile
        = xt_idxsection_new(1, 1, global_sizes,
                            next_rank_size, next_rank_start);
      gen_state->last_rank_generated = next_rank;
      return (struct Xt_com_list){
        .list = gen_state->last_tile,
        .rank = next_rank
      };
    }
  }
  return (struct Xt_com_list){
    .list = NULL,
    .rank = -1,
  };
}

Xt_config
tile_bucket_config_new(int ndims,
                       const Xt_int global_dims[ndims],
                       const int num_part[ndims])
{
  Xt_xmdd_bucket_gen bucket_gen = xt_xmdd_bucket_gen_new();
  Xt_config custom_bucket_config = xt_config_new();
  /* todo: generalize for arbitrary number of dimensions */
  (void)num_part;
  assert(ndims == 1);
  struct bucket_gen_1d_tile_param *bucket_params
    = xmalloc(sizeof (*bucket_params));
  bucket_params->global_sizes[0] = global_dims[0];
  xt_xmdd_bucket_gen_define_interface(
    bucket_gen,
    bucket_gen_1d_tile_init,
    bucket_gen_1d_tile_destroy,
    bucket_gen_1d_tile_get_intersect_max_num,
    bucket_gen_1d_tile_next,
    sizeof (struct bucket_gen_1d_tile_state),
    bucket_params);
  xt_xmdd_bucket_gen_set_extended_init_params(
    bucket_gen, sizeof (struct bucket_gen_1d_tile_param),
    (Xt_xmdd_bucket_gen_iparm_destroy)0,
    (Xt_xmdd_bucket_gen_iparm_copy)0);
  xt_config_set_xmdd_bucket_gen(custom_bucket_config,
                                bucket_gen);
  xt_xmdd_bucket_gen_delete(bucket_gen);
  return custom_bucket_config;
}

void
tile_bucket_config_set_global_dims(Xt_config custom_bucket_config,
                                   const Xt_int global_dims[])
{
  Xt_xmdd_bucket_gen bucket_gen
    = xt_config_get_xmdd_bucket_gen(custom_bucket_config);
#ifndef NDEBUG
  Xt_xmdd_bucket_gen_init_state ini_func
    = xt_xmdd_bucket_gen_get_init(bucket_gen);
#endif
  assert(ini_func == bucket_gen_1d_tile_init);
  struct bucket_gen_1d_tile_param *init_params
    = xt_xmdd_bucket_gen_get_init_params(bucket_gen);
  memcpy(init_params->global_sizes,
         global_dims, sizeof (Xt_int) * 1);
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
