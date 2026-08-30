/**
 * @file xt_redist_collection.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/xt_redist_collection.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt/xt_request.h"
#include "xt_redist_internal.h"
#include "xt_exchanger.h"
#include "xt_exchanger_internal.h"
#include "xt_config_internal.h"



static const char filename[] = "xt_redist_collection.c";

enum { DEFFAULT_DATATYPE_CACHE_SIZE=16 };

static void
redist_collection_delete(Xt_redist redist);

static Xt_redist
redist_collection_copy(Xt_redist redist);

static void
redist_collection_s_exchange(Xt_redist redist, int num_src_arrays,
                             const void *const *src_data,
                             void *const *dst_data);

static void
redist_collection_a_exchange(Xt_redist redist, int num_src_arrays,
                             const void *const *src_data,
                             void *const *dst_data,
                             Xt_request *request);

static void
redist_collection_s_exchange1(Xt_redist redist,
                              const void *src_data, void *dst_data);

static void
redist_collection_a_exchange1(Xt_redist redist,
                              const void *src_data, void *dst_data,
                              Xt_request *request);

static int redist_collection_get_num_msg(Xt_redist redist,
                                         enum xt_msg_direction direction);

static MPI_Datatype
redist_collection_get_MPI_Datatype(Xt_redist redist, int rank,
                                   enum xt_msg_direction direction, bool do_dup);

static int
redist_collection_get_msg_ranks(Xt_redist redist,
                                enum xt_msg_direction direction,
                                int *restrict *ranks);

static MPI_Comm
redist_collection_get_MPI_Comm(Xt_redist redist);

static const struct xt_redist_vtable redist_collection_vtable = {
  .copy                  = redist_collection_copy,
  .delete                = redist_collection_delete,
  .s_exchange            = redist_collection_s_exchange,
  .a_exchange            = redist_collection_a_exchange,
  .s_exchange1           = redist_collection_s_exchange1,
  .a_exchange1           = redist_collection_a_exchange1,
  .get_num_msg           = redist_collection_get_num_msg,
  .get_msg_MPI_Datatype  = redist_collection_get_MPI_Datatype,
  .get_exchanger         = 0,
  .get_msg_ranks         = redist_collection_get_msg_ranks,
  .get_MPI_Comm          = redist_collection_get_MPI_Comm
};

struct exchanger_cache
{
  size_t token;
  /* array structure: [cache_size][2][num_components] */
  MPI_Aint *displacements;
  Xt_exchanger * exchangers;
};

typedef struct Xt_redist_collection_ *Xt_redist_collection;

struct Xt_redist_collection_ {

  const struct xt_redist_vtable *vtable;

  struct exchanger_cache cache;

  unsigned num_components, cache_size;

  unsigned nmsg[2];
  int *msg_ranks, *blocklengths;
  MPI_Aint *disp_ofs;

  struct Xt_config_ config;
  MPI_Comm comm;
  int tag_offset;

  MPI_Datatype all_component_dt[];
};

static void align_component_dt(
  unsigned num_redists,
  size_t num_components,
  unsigned nmsgs,
  const struct Xt_msg_param *const msg_param[num_redists],
  const unsigned num_redist_components[num_redists],
  int *restrict in_ranks[num_redists],
  const unsigned num_ranks[num_redists],
  int *restrict out_ranks,
  MPI_Datatype *restrict component_dt,
  MPI_Aint *restrict disp_ofs,
  int *restrict blocklengths)
{
  size_t rank_pos[num_redists];
  for (size_t j = 0; j < num_redists; ++j)
    rank_pos[j] = 0;
  if (nmsgs) {
    /* find ranks and corresponding component datatypes */
    for (size_t i = 0; i < nmsgs; ++i) {
      int min_rank = INT_MAX;
      for (size_t j = 0; j < num_redists; ++j)
        if (rank_pos[j] < num_ranks[j] && in_ranks[j][rank_pos[j]] < min_rank)
          min_rank = in_ranks[j][rank_pos[j]];

      for (size_t j = 0, component_offset = 0;
           j < num_redists;
           component_offset += num_redist_components[j], ++j) {
        size_t num_redist_components_j = num_redist_components[j];
        if (rank_pos[j] < num_ranks[j] && in_ranks[j][rank_pos[j]] == min_rank) {
          size_t msg_pos = rank_pos[j] * num_redist_components_j;
          for (size_t k = 0; k < num_redist_components_j; ++k) {
            component_dt[i * num_components + component_offset + k]
              = msg_param[j][msg_pos + k].datatype;
            disp_ofs[i * num_components + component_offset + k]
              = msg_param[j][msg_pos + k].displacement;
            blocklengths[i * num_components + component_offset + k]
              = msg_param[j][msg_pos + k].count;
            assert(msg_param[j][msg_pos + k].rank == min_rank);
          }
        } else {
          for (size_t k = 0; k < num_redist_components_j; ++k) {
            component_dt[i * num_components + component_offset + k]
              = MPI_DATATYPE_NULL;
            disp_ofs[i * num_components + component_offset + k] = 0;
            blocklengths[i * num_components + component_offset + k] = 0;
          }
        }
      }
      out_ranks[i] = min_rank;
      for (size_t j = 0; j < num_redists; ++j)
        rank_pos[j]
          += (rank_pos[j] < num_ranks[j] && in_ranks[j][rank_pos[j]] == min_rank);
    }
  }
  free(in_ranks[0]);
}

/* not yet used cache entries are marked with -1 as first displacement,
 * which becomes 0 later on through use */
static inline void
init_cache(struct exchanger_cache *cache, size_t cache_size,
           unsigned num_redists)
{
  cache->exchangers = xcalloc(cache_size, sizeof(*(cache->exchangers)));
  size_t num_displ = cache_size * num_redists;
  MPI_Aint *restrict q = cache->displacements
    = xmalloc(2 * num_displ * sizeof (*q));
  for (size_t i = 0; i < 2 * num_displ; i += num_redists)
    q[i] = (MPI_Aint)-1;
  cache->token = 0;
}

static inline void
destruct_cache(struct exchanger_cache *cache,
               size_t cache_size)
{
  for (size_t i = 0; i < cache_size; ++i)
    if (cache->exchangers[i] != NULL)
      xt_exchanger_delete(cache->exchangers[i]);
  free(cache->exchangers);
  free(cache->displacements);
}

Xt_redist xt_redist_collection_new(Xt_redist *redists, int num_redists,
                                   int cache_size, MPI_Comm comm)
{
  return xt_redist_collection_custom_new(redists, num_redists, cache_size,
                                         comm, (Xt_config)&xt_default_config);
}

static Xt_redist_collection
alloc_redist_coll(size_t num_components, size_t nmsg_send, size_t nmsg_recv,
                  Xt_config config, MPI_Comm comm)
{
  size_t nmsg = nmsg_recv + nmsg_send,
    header_size = sizeof (struct Xt_redist_collection_),
    all_component_dt_size = sizeof (MPI_Datatype) * num_components * nmsg,
    disp_ofs_size = sizeof (MPI_Aint) * num_components * nmsg,
    ranks_size = sizeof (int) * nmsg,
    blocklengths_size = sizeof (int) * num_components * nmsg,
    team_data_size,
    own_size = header_size + disp_ofs_size + all_component_dt_size
    + ranks_size + blocklengths_size;
  Xt_exchanger_new exchanger_new = (Xt_exchanger_new)0;
  if (!config->exchanger_team_share) {
    exchanger_new
      = xt_config_get_exchange_new_by_comm(config, comm);
    team_data_size
      = xt_exchanger_new_team_get_share_size(exchanger_new);
    /* round up own size to align beginning of team_data */
    own_size = ((own_size + team_share_align - 1)
                / team_share_align) * team_share_align;
  } else
    team_data_size = 0;
  Xt_redist_collection redist_coll
    = xmalloc(own_size + team_data_size);
  redist_coll->vtable = &redist_collection_vtable;
  redist_coll->num_components = (unsigned)num_components;
  redist_coll->nmsg[SEND] = (unsigned)nmsg_send;
  redist_coll->nmsg[RECV] = (unsigned)nmsg_recv;
  redist_coll->msg_ranks
    = (int *)(redist_coll->all_component_dt + nmsg * num_components);
  redist_coll->disp_ofs
    = (void *)(intptr_t)(redist_coll->msg_ranks + nmsg);
  redist_coll->blocklengths
    = (void *)(intptr_t)(redist_coll->disp_ofs + nmsg * num_components);
  redist_coll->config = *config;
  /* all mpi datatypes are created for the exclusive use by the
   * underlying exchanger, no need for dup'ing them */
  redist_coll->config.flags |= exch_no_dt_dup;
  if (!config->exchanger_team_share) {
    redist_coll->config.exchanger_team_share
      = (unsigned char *)redist_coll + own_size;
    xt_exchanger_new_team_share_default_init(
      exchanger_new, redist_coll->config.exchanger_team_share);
  }
  redist_coll->comm
    = xt_mpi_comm_smart_dup(comm, &redist_coll->tag_offset);
  return redist_coll;
}


#if defined __ICC && (__INTEL_COMPILER_BUILD_DATE == 20110811   \
                       || __INTEL_COMPILER == 1300)
PPM_DSO_INTERNAL
#else
static
#endif
size_t
xt_redist_collection_get_component_count(
  size_t num_redists,
  const unsigned num_ranks[num_redists],
  const unsigned num_redist_components[num_redists])
{
  size_t num_components_total = 0;
  for (size_t j = 0; j < num_redists; ++j)
    num_components_total += (size_t)num_ranks[j] * num_redist_components[j];
  return num_components_total;
}

static inline Xt_redist_collection
xrc(void *redist)
{
  return (Xt_redist_collection)redist;
}

#define XT_SORTFUNC_DECL static
#define SORT_TYPE struct Xt_msg_param
#define SORT_TYPE_SUFFIX msg_param_by_rank
#define SORT_TYPE_CMP_LT(u,v,i,j) ((u).rank < (v).rank)
#define SORT_TYPE_CMP_LE(u,v,i,j) ((u).rank <=(v).rank)
#define SORT_TYPE_CMP_EQ(u,v,i,j) ((u).rank == (v).rank)
#include "xt_mergesort_base.h"

static void
xrc_component_copy(
  struct Xt_msg_param *restrict msg_param[2],
  size_t msg_idx[2],
  const Xt_redist redist,
  size_t num_redist_components,
  enum xt_msg_direction direction,
  MPI_Comm comm
  )
{
  /* source component_redist => src_cr */
  Xt_redist_collection src_cr = xrc(redist);
  unsigned num_src_msgs = src_cr->nmsg[direction];
  struct Xt_msg_param *restrict p_msg_param = msg_param[direction];
  size_t out_idx = msg_idx[direction];
  size_t adj = direction == xt_direction_recv
    ? src_cr->nmsg[xt_direction_send] : 0;
  const int *msg_ranks = src_cr->msg_ranks + adj;
  const int *blocklengths = src_cr->blocklengths + adj * num_redist_components;
  const MPI_Aint *disp_ofs = src_cr->disp_ofs + adj * num_redist_components;
  const MPI_Datatype *dt = src_cr->all_component_dt + adj * num_redist_components;
  for (size_t i = 0; i < num_src_msgs; ++i) {
    for (size_t j = 0; j < num_redist_components; ++j, ++out_idx) {
      MPI_Datatype dt_ij = dt[i * num_redist_components + j];
      if (!xt_mpi_ddt_is_named(dt_ij))
        xt_mpi_call(MPI_Type_dup(dt_ij, &dt_ij), comm);
      p_msg_param[out_idx] = (struct Xt_msg_param){
        .rank = msg_ranks[i],
        .datatype = dt_ij,
        .displacement = disp_ofs[i * num_redist_components + j],
        .count = blocklengths[i * num_redist_components + j],
      };
    }
  }
  msg_idx[direction] = out_idx;
}

static void
setup_component_msg_params(
  size_t num_redists,
  const Xt_redist *redists,
  const unsigned num_ranks[2][num_redists],
  const unsigned num_redist_components[num_redists],
  struct Xt_msg_param *msg_param[2][num_redists],
  MPI_Comm comm, Xt_config config)
{
  size_t num_send_components_total
    = xt_redist_collection_get_component_count(
      num_redists, num_ranks[xt_direction_send], num_redist_components),
    num_recv_components_total
    = xt_redist_collection_get_component_count(
      num_redists, num_ranks[xt_direction_recv], num_redist_components);
  struct Xt_msg_param *mp_buf = xmalloc(
    (num_send_components_total + num_recv_components_total)
    * sizeof (*mp_buf));
  struct Xt_msg_param *restrict p_msg_param[2]
    = { [xt_direction_send] = mp_buf,
        [xt_direction_recv] = mp_buf + num_send_components_total };
  size_t msg_idx[2] = { 0, 0 };
  size_t old_msg_idx[2];
  for (size_t j = 0; j < num_redists; ++j) {
    for (size_t direction = 0; direction < 2; ++direction) {
      old_msg_idx[direction] = msg_idx[direction];
      msg_param[direction][j]
        = p_msg_param[direction]+msg_idx[direction];
    }
    if (redists[j]->vtable == &redist_collection_vtable) {
      xrc_component_copy(p_msg_param, msg_idx, redists[j],
                         num_redist_components[j],
                         xt_direction_send, comm);
      xrc_component_copy(p_msg_param, msg_idx, redists[j],
                         num_redist_components[j],
                         xt_direction_recv, comm);
    } else {
      Xt_exchanger exchanger = xt_redist_get_exchanger(redists[j]);
      Xt_exchanger_iter exch_iter
        = xt_exchanger_get_iterator(exchanger, config);
      while (exch_iter) {
        int direction = xt_exchanger_iter_get_direction(exch_iter);
        p_msg_param[direction][msg_idx[direction]++] = (struct Xt_msg_param){
          .datatype = xt_exchanger_iter_get_MPI_Datatype(exch_iter),
          .displacement = xt_exchanger_iter_get_displacement(exch_iter),
          .rank = xt_exchanger_iter_get_rank(exch_iter),
          .count = xt_exchanger_iter_get_count(exch_iter),
        };
        if (!xt_exchanger_iter_next(exch_iter)) {
          xt_exchanger_iter_delete(exch_iter);
          exch_iter = NULL;
        }
      }
    }
    assert(msg_idx[xt_direction_send] - old_msg_idx[xt_direction_send]
           == num_ranks[xt_direction_send][j] * num_redist_components[j]
           && msg_idx[xt_direction_recv] - old_msg_idx[xt_direction_recv]
           == num_ranks[xt_direction_recv][j] * num_redist_components[j]);
    for (size_t direction = 0; direction < 2; ++direction)
      xt_mergesort_msg_param_by_rank(p_msg_param[direction]
                                     + old_msg_idx[direction],
                                     msg_idx[direction]
                                     - old_msg_idx[direction]);
  }
}

Xt_redist xt_redist_collection_custom_new(Xt_redist *redists, int num_redists,
                                          int cache_size, MPI_Comm comm,
                                          Xt_config config)
{
  // ensure that yaxt is initialized
  assert(xt_initialized());
#ifndef NDEBUG
  if (num_redists > 0) ; else
    Xt_abort(comm, "ERROR: invalid number of redists passed",
             "xt_redist_collection_custom_new", __LINE__);
  if (cache_size >= -1) ; else
    Xt_abort(comm, "ERROR: invalid cache size in xt_redist_collection_new",
             filename, __LINE__);
#endif
  unsigned num_redists_ = (unsigned)num_redists;
  unsigned num_ranks[2][num_redists_];
  int *restrict ranks[2][num_redists_];
  unsigned num_redist_components[num_redists_];
  unsigned nmsg_send = xt_redist_agg_msg_count(
      num_redists_, xt_direction_send, redists,
      num_ranks[SEND], ranks[SEND], config),
    nmsg_recv = xt_redist_agg_msg_count(
      num_redists_, xt_direction_recv, redists,
      num_ranks[RECV], ranks[RECV], config);
  size_t num_components = 0;
  for (size_t i = 0; i < num_redists_; ++i)
    num_components +=
      num_redist_components[i] =
      (redists[i]->vtable == &redist_collection_vtable)
      ? ((Xt_redist_collection)redists[i])->num_components : 1;
  struct Xt_msg_param *msg_param[2][num_redists_];
  setup_component_msg_params(num_redists_, redists, (void *)num_ranks,
                             num_redist_components, msg_param, comm, config);
  Xt_redist_collection redist_coll
    = alloc_redist_coll(num_components, nmsg_send, nmsg_recv, config, comm);
  redist_coll->cache_size
    = (cache_size == -1)
    ?((unsigned)DEFFAULT_DATATYPE_CACHE_SIZE):(unsigned)cache_size;

  xt_redist_check_comms(redists, num_redists, comm);

  MPI_Datatype *all_component_dt = redist_coll->all_component_dt;
  align_component_dt(num_redists_, num_components, nmsg_send,
                     (void *)msg_param[xt_direction_send], num_redist_components,
                     ranks[SEND], num_ranks[SEND], redist_coll->msg_ranks,
                     all_component_dt, redist_coll->disp_ofs,
                     redist_coll->blocklengths);
  align_component_dt(num_redists_, num_components, nmsg_recv,
                     (void *)msg_param[xt_direction_recv], num_redist_components,
                     ranks[xt_direction_recv], num_ranks[xt_direction_recv],
                     redist_coll->msg_ranks+nmsg_send,
                     all_component_dt + nmsg_send * num_components,
                     redist_coll->disp_ofs + nmsg_send * num_components,
                     redist_coll->blocklengths + nmsg_send * num_components);
  free(msg_param[0][0]);
  init_cache(&redist_coll->cache, redist_coll->cache_size, num_components);

  return (Xt_redist)redist_coll;
}


static void
create_all_dt_for_dir(
  unsigned num_messages, size_t num_redists,
  const int ranks[num_messages],
  const MPI_Datatype component_dt[num_messages][num_redists],
  MPI_Aint *restrict displacements,
  const MPI_Aint disp_ofs[num_messages][num_redists],
  const int blocklengths[num_messages][num_redists],
  struct Xt_msg_param msgs[num_messages],
  struct Xt_mpiddt_list *ddt_list,
  MPI_Comm comm)
{
  for (size_t i = 0; i < num_messages; ++i) {
    for (size_t j = 0; j < num_redists; ++j)
      displacements[j] += disp_ofs[i][j];
    struct Xt_mpi_dt_ofs dt_ofs
      = xt_create_compound_datatype(num_redists, displacements,
                                    component_dt[i],
                                    blocklengths[i], ddt_list, comm);
    msgs[i].datatype = dt_ofs.dt;
    for (size_t j = 0; j < num_redists; ++j)
      displacements[j] -= disp_ofs[i][j];
    msgs[i].rank = ranks[i];
    msgs[i].displacement = dt_ofs.displacement;
    msgs[i].count = 1;
  }
}

static void
compute_displ(const void *const *data, unsigned num_redists,
              MPI_Aint displacements[num_redists])
{
  if (num_redists) {
    MPI_Aint base_addr, offset;
    base_addr = (MPI_Aint)(intptr_t)(const void *)data[0];
    displacements[0] = 0;
    for (size_t i = 1; i < num_redists; ++i) {
      offset = (MPI_Aint)(intptr_t)(const void *)data[i];
      displacements[i] = offset - base_addr;
    }
  }
}

static size_t
lookup_cache_index(unsigned num_redists,
                   const MPI_Aint displacements[2][num_redists],
                   const MPI_Aint (*cached_displacements)[2][num_redists],
                   size_t cache_size)
{
  for (size_t i = 0; i < cache_size &&
       cached_displacements[i][SEND][0] == (MPI_Aint)0 &&
       cached_displacements[i][RECV][0] == (MPI_Aint)0; ++i) {
    bool mismatch = false;
    for (size_t d = 0; d < 2; ++d)
      for (size_t j = 0; j < num_redists; ++j)
        mismatch |= (displacements[d][j] != cached_displacements[i][d][j]);
    if (!mismatch) return i;
  }
  return cache_size;
}

enum { redist_msg_stack_alloc_lim = 16 };

static Xt_exchanger
create_exchanger(struct Xt_redist_collection_ *redist_coll,
                 size_t num_redists,
                 MPI_Aint displacements[2][num_redists])
{
  unsigned nmsg_send = redist_coll->nmsg[SEND],
    nmsg_recv = redist_coll->nmsg[RECV];
  size_t nmsg = (size_t)nmsg_send + nmsg_recv;
  struct Xt_msg_param *restrict msgs = xmalloc(nmsg * sizeof (*msgs));

  const MPI_Datatype (*all_component_dt)[num_redists]
    = (void *)redist_coll->all_component_dt;
  const MPI_Aint (*disp_ofs)[num_redists]
    = (void *)redist_coll->disp_ofs;
  const int (*blocklengths)[num_redists]
    = (void *)redist_coll->blocklengths;
  struct Xt_mpiddt_list ddt_list;
  Xt_mpi_ddt_cache_init(&ddt_list, &redist_coll->config);
  MPI_Comm comm = redist_coll->comm;
  create_all_dt_for_dir(nmsg_send, num_redists,
                        redist_coll->msg_ranks, all_component_dt,
                        displacements[SEND], disp_ofs, blocklengths,
                        msgs, &ddt_list, comm);
  /* the table view needs to be adjusted by this to switch to the
     receive direction */
  /* unfortunately, older Intel compilers flunk the offset computation */
#if defined __ICC && ( __INTEL_COMPILER_BUILD_DATE == 20110811   \
                       || __INTEL_COMPILER == 1300) \
  || defined __PGI && __PGIC__ <= 15
  all_component_dt = (void *)(redist_coll->all_component_dt
                              + num_redists * nmsg_send);
  disp_ofs = (void *)(redist_coll->disp_ofs + num_redists * nmsg_send);
  blocklengths = (void *)(redist_coll->blocklengths
                          + num_redists * nmsg_send);
#else
  all_component_dt += nmsg_send;
  disp_ofs += nmsg_send;
  blocklengths += nmsg_send;
#endif
  create_all_dt_for_dir(nmsg_recv, num_redists,
                        redist_coll->msg_ranks+nmsg_send,
                        all_component_dt, displacements[RECV], disp_ofs,
                        blocklengths, msgs + nmsg_send, &ddt_list, comm);
  Xt_mpi_ddt_cache_free(&ddt_list, comm);
  Xt_exchanger_new exchanger_new
    = xt_config_get_exchange_new_by_comm(&redist_coll->config, comm);
  Xt_exchanger exch
    = exchanger_new((int)nmsg_send, (int)nmsg_recv, msgs, msgs + nmsg_send,
                    comm, redist_coll->tag_offset, &redist_coll->config);
  free(msgs);
  return exch;
}


static Xt_exchanger
get_exchanger(struct Xt_redist_collection_ *redist_coll,
              const void *const * src_data, const void *const * dst_data,
              unsigned num_redists)
{
  MPI_Aint displacements[2][num_redists];
  compute_displ(src_data, num_redists, displacements[SEND]);
  compute_displ(dst_data, num_redists, displacements[RECV]);

  Xt_exchanger exchanger;

  struct exchanger_cache *restrict cache = &redist_coll->cache;
  size_t cache_size = redist_coll->cache_size,
    cache_index = 0;
  if (cache_size > 0)
  {
    cache_index
      = lookup_cache_index(num_redists,
                           (const MPI_Aint (*)[num_redists])displacements,
                           (const MPI_Aint (*)[2][num_redists])
                           cache->displacements,
                           cache_size);

    if (cache_index != cache_size)
      return cache->exchangers[cache_index];

    cache_index = cache->token;
    cache->token = (cache_index + 1) % cache_size;
    memcpy(cache->displacements + cache_index * 2 * num_redists,
           displacements, sizeof (displacements));

    if (cache->exchangers[cache_index] != NULL)
      xt_exchanger_delete(cache->exchangers[cache_index]);
  }
  exchanger =
    create_exchanger(redist_coll, num_redists, displacements);
  if (cache_size > 0)
    cache->exchangers[cache_index] = exchanger;
  return exchanger;
}

static void
redist_collection_s_exchange(Xt_redist redist, int num_arrays,
                             const void *const *src_data,
                             void *const *dst_data)
{
  Xt_redist_collection redist_coll = xrc(redist);

  if (num_arrays != (int)redist_coll->num_components)
    Xt_abort(redist_coll->comm, "ERROR: wrong number of arrays in "
             "redist_collection_s_exchange", filename, __LINE__);


  Xt_exchanger exchanger = get_exchanger(redist_coll,
                                         src_data,
                                         (const void *const (*))dst_data,
                                         redist_coll->num_components);

  xt_exchanger_s_exchange(exchanger, src_data[0], dst_data[0]);

  if (redist_coll->cache_size == 0)
    xt_exchanger_delete(exchanger);
}

static void
redist_collection_a_exchange(Xt_redist redist, int num_arrays,
                             const void *const *src_data, void *const *dst_data,
                             Xt_request *request) {

  Xt_redist_collection redist_coll = xrc(redist);

  if (num_arrays != (int)redist_coll->num_components)
    Xt_abort(redist_coll->comm, "ERROR: wrong number of arrays in "
             "redist_collection_a_exchange", filename, __LINE__);


  Xt_exchanger exchanger = get_exchanger(redist_coll,
                                         src_data,
                                         (const void *const (*))dst_data,
                                         redist_coll->num_components);

  xt_exchanger_a_exchange(exchanger, src_data[0], dst_data[0], request);

  if (redist_coll->cache_size == 0)
    xt_exchanger_delete(exchanger);

}

static void
copy_component_dt(size_t num_component_dt,
                  const MPI_Datatype *component_dt_orig,
                  MPI_Datatype *component_dt_copy,
                  MPI_Comm comm)
{
  for (size_t i = 0; i < num_component_dt; ++i)
  {
    MPI_Datatype orig_dt = component_dt_orig[i];
    if (orig_dt != MPI_DATATYPE_NULL && !xt_mpi_ddt_is_named(orig_dt))
      xt_mpi_call(MPI_Type_dup(orig_dt, component_dt_copy + i), comm);
    else
      component_dt_copy[i] = orig_dt;
  }
}

static Xt_redist
redist_collection_custom_copy(Xt_redist redist, Xt_config config);

static Xt_redist
redist_collection_copy(Xt_redist redist)
{
  Xt_redist_collection redist_coll = xrc(redist);
  unsigned nmsg_send = redist_coll->nmsg[SEND],
    nmsg_recv = redist_coll->nmsg[RECV];
  size_t nmsg = (size_t)nmsg_recv + nmsg_send;
  struct Xt_config_ config = redist_coll->config;
  /* if redist comes with an integrated team share, cause creation of
   * a new one */
  void *team_share = config.exchanger_team_share;
  if ((unsigned char *)team_share > (unsigned char *)redist_coll
      && (unsigned char *)team_share <
      (unsigned char *)(redist_coll->blocklengths
                        + nmsg * redist_coll->num_components) + team_share_align)
    config.exchanger_team_share = NULL;
  return redist_collection_custom_copy(redist, &config);
}

static Xt_redist
redist_collection_custom_copy(Xt_redist redist, Xt_config config)
{
  Xt_redist_collection redist_coll = xrc(redist);
  unsigned num_components = redist_coll->num_components,
    nmsg_send = redist_coll->nmsg[SEND],
    nmsg_recv = redist_coll->nmsg[RECV];
  size_t nmsg = (size_t)nmsg_recv + nmsg_send;
  Xt_redist_collection redist_copy
    = alloc_redist_coll(num_components, nmsg_send, nmsg_recv, config,
                        redist_coll->comm);
  memcpy(redist_copy->msg_ranks, redist_coll->msg_ranks,
         sizeof (*redist_copy->msg_ranks) * nmsg);
  memcpy(redist_copy->blocklengths, redist_coll->blocklengths,
         sizeof (*redist_copy->blocklengths) * nmsg * num_components);
  memcpy(redist_copy->disp_ofs, redist_coll->disp_ofs,
         sizeof (*redist_copy->disp_ofs) * nmsg * num_components);
  copy_component_dt(num_components * nmsg,
                    redist_coll->all_component_dt,
                    redist_copy->all_component_dt, redist_copy->comm);
  unsigned cache_size = redist_coll->cache_size;
  redist_copy->cache_size = cache_size;
  init_cache(&redist_copy->cache, cache_size, num_components);
  return (Xt_redist)redist_copy;
}

static void
free_component_dt(size_t num_dt, MPI_Datatype *all_component_dt, MPI_Comm comm)
{
  for (size_t i = 0; i < num_dt; ++i)
    if (all_component_dt[i] != MPI_DATATYPE_NULL
        && !xt_mpi_ddt_is_named(all_component_dt[i]))
      xt_mpi_call(MPI_Type_free(all_component_dt + i), comm);
}

static void
redist_collection_delete(Xt_redist redist) {

  Xt_redist_collection redist_coll = xrc(redist);

  unsigned num_components = redist_coll->num_components;
  size_t nmsg = (size_t)redist_coll->nmsg[RECV] + redist_coll->nmsg[SEND];
  free_component_dt(nmsg * num_components, redist_coll->all_component_dt,
                    redist_coll->comm);

  destruct_cache(&redist_coll->cache, redist_coll->cache_size);

  void *team_share = redist_coll->config.exchanger_team_share;
  if ((unsigned char *)team_share > (unsigned char *)redist_coll
      && (unsigned char *)team_share
      < ((unsigned char *)(redist_coll->blocklengths + nmsg * num_components)
         + team_share_align)) {
    Xt_exchanger_new exchanger_new
      = xt_config_get_exchange_new_by_comm(&redist_coll->config,
                                           redist_coll->comm);
    xt_exchanger_new_team_share_destroy(exchanger_new, team_share);
  }
  xt_mpi_comm_smart_dedup(&(redist_coll->comm), redist_coll->tag_offset);

  free(redist_coll);
}

static int redist_collection_get_num_msg(Xt_redist redist,
                                         enum xt_msg_direction direction)
{
  return (int)(xrc(redist)->nmsg[direction]);
}

static MPI_Datatype
redist_collection_get_MPI_Datatype(Xt_redist redist, int rank,
                                   enum xt_msg_direction direction,
                                   bool do_dup)
{
  Xt_redist_collection redist_coll = xrc(redist);

  MPI_Datatype dt = MPI_DATATYPE_NULL;
  if (redist_coll->num_components > 1)
    Xt_abort(redist_coll->comm, "ERROR: datatype retrieval is not"
             " supported for this xt_redist type (Xt_redist_collection)",
             filename, __LINE__);
  else {
    size_t adj = direction == xt_direction_recv
      ? redist_coll->nmsg[xt_direction_send] : 0;
    const int *msg_ranks = redist_coll->msg_ranks + adj;
    const MPI_Datatype *all_dt = redist_coll->all_component_dt + adj;
    for (size_t i = 0, num_msg = redist_coll->nmsg[direction]; i < num_msg; ++i)
    {
      if (msg_ranks[i] < rank)
        goto ret;
      else if (msg_ranks[i] == rank) {
        dt = all_dt[i];
        goto check_dup;
      }
    }
  check_dup:
    if (do_dup)
      xt_mpi_call(MPI_Type_dup(dt, &dt), redist_coll->comm);
  }
ret:
  return dt;
}

static void
redist_collection_s_exchange1(Xt_redist redist,
                              const void *src_data, void *dst_data)
{

  Xt_redist_collection redist_coll = xrc(redist);
  if (redist_coll->num_components == 1)
    redist_collection_s_exchange(redist, 1, &src_data, &dst_data);
  else
    Xt_abort(redist_coll->comm, "ERROR: s_exchange1 is not implemented for"
             " this xt_redist type (Xt_redist_collection)", filename, __LINE__);
}

static void
redist_collection_a_exchange1(Xt_redist redist,
                              const void *src_data, void *dst_data,
                              Xt_request *request)
{

  Xt_redist_collection redist_coll = xrc(redist);
  if (redist_coll->num_components == 1)
    redist_collection_a_exchange(redist, 1, &src_data, &dst_data, request);
  else
    Xt_abort(redist_coll->comm, "ERROR: a_exchange1 is not implemented for"
             " this xt_redist type (Xt_redist_collection)", filename, __LINE__);
}

static int
redist_collection_get_msg_ranks(Xt_redist redist,
                                enum xt_msg_direction direction,
                                int *restrict *ranks)
{
  Xt_redist_collection redist_coll = xrc(redist);
  unsigned nmsg_direction = redist_coll->nmsg[direction],
    nmsg_send = redist_coll->nmsg[SEND];
  size_t nmsg = (size_t)nmsg_direction + redist_coll->nmsg[!direction];
  int *ranks_orig
    = (int *)(redist_coll->all_component_dt + nmsg * redist_coll->num_components)
    + (((unsigned)direction-1) & nmsg_send);
  if (nmsg_direction > 0) {
    int *ranks_ = *ranks;
    if (!ranks_)
      ranks_ = *ranks = xmalloc(nmsg_direction * sizeof (*ranks_));
    memcpy(ranks_, ranks_orig, nmsg_direction * sizeof (*ranks_));
  }
  return (int)nmsg_direction;
}


static MPI_Comm
redist_collection_get_MPI_Comm(Xt_redist redist) {

  Xt_redist_collection redist_coll = xrc(redist);

  return redist_coll->comm;
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
