/**
 * @file xt_redist_collection.c
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

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt/xt_redist_collection.h"
#include "ensure_array_size.h"
#include "xt/xt_redist.h"
#include "xt_redist_internal.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

/**
 * @brief round size to next multiple of mult
 * @param size size to round up from
 * @param mult
 * @return q with q >= size and q mod mult == 0
 */
static inline size_t sizeRoundUp(size_t size, size_t mult)
{
  return ((size + mult - 1)/mult) * mult;
}

#define DEFFAULT_DATATYPE_CACHE_SIZE (16)

static void
redist_collection_delete(Xt_redist redist);

static void
redist_collection_s_exchange(Xt_redist redist, unsigned num_src_arrays,
                             void **src_data, void **dst_data);

static void
redist_collection_s_exchange1(Xt_redist redist, void *src_data, void *dst_data);

static MPI_Datatype
redist_collection_get_send_MPI_Datatype(Xt_redist redist, int rank);

static MPI_Datatype
redist_collection_get_recv_MPI_Datatype(Xt_redist redist, int rank);

static MPI_Comm
redist_collection_get_MPI_Comm(Xt_redist redist);

static const struct xt_redist_vtable redist_collection_vtable = {
  .delete                = redist_collection_delete,
  .s_exchange            = redist_collection_s_exchange,
  .s_exchange1           = redist_collection_s_exchange1,
  .get_send_MPI_Datatype = redist_collection_get_send_MPI_Datatype,
  .get_recv_MPI_Datatype = redist_collection_get_recv_MPI_Datatype,
  .get_MPI_Comm          = redist_collection_get_MPI_Comm
};

struct redist_collection_msg {

  int rank;
  MPI_Datatype *component_dt; // datatypes of the redists (size == num_redists)
};

struct dt_cache
{
  size_t token;
  MPI_Aint *displacements;
  MPI_Datatype *dt;
};

typedef struct Xt_redist_collection_ *Xt_redist_collection;

struct Xt_redist_collection_ {

  const struct xt_redist_vtable *vtable;

  unsigned num_redists;

  struct dt_cache src_cache, dst_cache;

  int ndst, nsrc;
  struct redist_collection_msg * send_msgs;
  struct redist_collection_msg * recv_msgs;

  size_t cache_size;

  MPI_Comm comm;
};

static void copy_component_dt(struct redist_collection_msg ** msgs, int * nmsgs,
                              Xt_redist * redists, unsigned num_redists,
                              MPI_Comm comm,
                              MPI_Datatype (*get_MPI_datatype)(Xt_redist,int))
{
  size_t msgs_array_size = 0;

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  MPI_Datatype datatypes[num_redists];
  assert(*nmsgs >= 0);
  size_t num_messages = (size_t)*nmsgs;
  struct redist_collection_msg *p = *msgs;

  for (int i = 0; i < comm_size; ++i) {

    int flag = 0;

    for (unsigned j = 0; j < num_redists; ++j)
      flag |= ((datatypes[j] = get_MPI_datatype(redists[j], i))
               != MPI_DATATYPE_NULL);

    if (flag) {

        ENSURE_ARRAY_SIZE(p, msgs_array_size, num_messages+1);

        p[num_messages].rank = i;
        p[num_messages].component_dt = xmalloc(num_redists *
          sizeof(*(p[num_messages].component_dt)));
        memcpy(p[num_messages].component_dt, datatypes,
               num_redists * sizeof(*datatypes));

        ++num_messages;
    }
  }

  if (num_messages > 0)
    p = xrealloc(p, num_messages * sizeof(**msgs));
  *msgs = p;
  *nmsgs = (int)num_messages;
}

/* not yet used cache entries are marked with -1 as first displacement,
 * which becomes 0 later on through use */
static inline void
init_cache(struct dt_cache *cache,
           size_t cache_size, size_t ntx, unsigned num_redists)
{
  size_t num_cache_dt = ntx * cache_size;
  MPI_Datatype *p
    = cache->dt = xmalloc(num_cache_dt * sizeof (*p));
  for (size_t i = 0; i < num_cache_dt; ++i)
    p[i] = MPI_DATATYPE_NULL;
  size_t num_displ = cache_size * num_redists;
  MPI_Aint *q
    = cache->displacements = xmalloc(num_displ * sizeof (*q));
  for (size_t i = 0; i < num_displ; i += num_redists)
    q[i] = (MPI_Aint)-1;
  cache->token = 0;
}

static inline void
destruct_cache(struct dt_cache *cache,
               size_t cache_size, size_t ntx, MPI_Comm comm)
{
  size_t num_cache_dt = ntx * cache_size;
  for (size_t i = 0; i < num_cache_dt; ++i)
    if (cache->dt[i] != MPI_DATATYPE_NULL)
      xt_mpi_call(MPI_Type_free(cache->dt + i), comm);
  free(cache->dt);
  free(cache->displacements);
}

static void check_comms(Xt_redist * redists, unsigned num_redists,
                        MPI_Comm comm) {

  int result;

  for (int i = 0; i < num_redists; ++i) {

    xt_mpi_call(MPI_Comm_compare(xt_redist_get_MPI_Comm(redists[i]),
                                 comm, &result), comm);

    if ((result != MPI_IDENT) && (result != MPI_CONGRUENT))
      Xt_abort(comm, "ERROR: MPI communicators do not match; cannot build "
               "redist_collection_static\n", __FILE__, __LINE__);
  }
}

Xt_redist xt_redist_collection_new(Xt_redist * redists, unsigned num_redists,
                                   int cache_size, MPI_Comm comm) {

  Xt_redist_collection redist_coll = xmalloc(sizeof (*redist_coll));

  redist_coll->vtable = &redist_collection_vtable;
  redist_coll->num_redists = num_redists;
  redist_coll->ndst = 0;
  redist_coll->nsrc = 0;
  redist_coll->send_msgs = NULL;
  redist_coll->recv_msgs = NULL;
  if (cache_size < -1)
    Xt_abort(comm, "ERROR: invalid cache size in xt_redist_collection_new",
             __FILE__, __LINE__);
  redist_coll->cache_size
    = (cache_size == -1)?(DEFFAULT_DATATYPE_CACHE_SIZE):(size_t)cache_size;

  xt_mpi_call(MPI_Comm_dup(comm, &(redist_coll->comm)), comm);

  check_comms(redists, num_redists, comm);

  copy_component_dt(&(redist_coll->send_msgs), &(redist_coll->nsrc), redists,
                    num_redists, redist_coll->comm,
                    xt_redist_get_send_MPI_Datatype);
  init_cache(&redist_coll->src_cache,
             redist_coll->cache_size, (size_t)redist_coll->nsrc, num_redists);

  copy_component_dt(&(redist_coll->recv_msgs), &(redist_coll->ndst), redists,
                    num_redists, redist_coll->comm,
                    xt_redist_get_recv_MPI_Datatype);
  init_cache(&redist_coll->dst_cache,
             redist_coll->cache_size, (size_t)redist_coll->ndst, num_redists);


  return (Xt_redist)redist_coll;
}

static MPI_Datatype
create_compound_dt(MPI_Aint *displacements, int *block_lengths,
                   struct redist_collection_msg *msg,
                   unsigned num_redists, MPI_Comm comm)
{
  MPI_Datatype datatype;

  unsigned num_datatypes = 0;

  for (unsigned i = 0; i < num_redists; ++i)
    num_datatypes += (msg->component_dt[i] != MPI_DATATYPE_NULL);

  MPI_Datatype * datatypes;
  MPI_Aint * displacements_;

  if (num_datatypes != num_redists) {

    datatypes = xmalloc(num_datatypes * sizeof(*datatypes));
    displacements_ = xmalloc(num_datatypes * sizeof(*displacements));

    num_datatypes = 0;

    for (unsigned i = 0; i < num_redists; ++i) {
      if (msg->component_dt[i] != MPI_DATATYPE_NULL) {
        datatypes[num_datatypes] = msg->component_dt[i];
        displacements_[num_datatypes] = displacements[i];
        ++num_datatypes;
      }
    }
  } else {
    datatypes = msg->component_dt;
    displacements_ = displacements;
  }

  assert(num_datatypes <= INT_MAX);
  xt_mpi_call(MPI_Type_create_struct((int)num_datatypes, block_lengths,
                                     displacements_, datatypes, &datatype),
              comm);

  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  if (num_datatypes != num_redists) {
    free(datatypes);
    free(displacements_);
  }
  return datatype;
}

static void
create_all_dt_for_dir(struct redist_collection_msg *msgs,
                      int num_messages, unsigned num_redists,
                      MPI_Aint displacements[num_redists],
                      MPI_Datatype dt[num_messages], MPI_Comm comm)
{
  int block_lengths[num_redists];

  for (unsigned i = 0; i < num_redists; ++i)
    block_lengths[i] = 1;
  for (int i = 0; i < num_messages; ++i)
  {
    if (dt[i] != MPI_DATATYPE_NULL)
      xt_mpi_call(MPI_Type_free(dt + i), comm);
    dt[i] = create_compound_dt(displacements, block_lengths,
                               msgs+i, num_redists, comm);
  }
}

static void
compute_displ(void **data, unsigned num_redists,
              MPI_Aint displacements[num_redists],
              MPI_Comm comm)
{
  if (num_redists)
  {
    MPI_Aint base_addr, offset;
    xt_mpi_call(MPI_Get_address(data[0], &base_addr), comm);
    displacements[0] = 0;
    for (unsigned i = 1; i < num_redists; ++i)
    {
      xt_mpi_call(MPI_Get_address(data[i], &offset), comm);
      displacements[i] = offset - base_addr;
    }
  }
}

static size_t
lookup_cache_index(unsigned num_redists,
                   MPI_Aint displacements[num_redists],
                   MPI_Aint (*cached_displacements)[num_redists],
                   size_t cache_size)
{
  for (size_t i = 0;
       i < cache_size && cached_displacements[i][0] == (MPI_Aint)0;
       ++i) {
    int mismatch = 0;
    for (unsigned j = 0; j < num_redists; ++j)
      mismatch |= (displacements[j] != cached_displacements[i][j]);
    if (!mismatch) return i;
  }
  return cache_size;
}

static MPI_Datatype *
get_compound_datatype(void ** data, struct redist_collection_msg * msgs,
                      int num_messages, unsigned num_redists,
                      struct dt_cache *cache, size_t cache_size,
                      MPI_Datatype temp_dt[num_messages],
                      MPI_Comm comm)
{
  MPI_Aint displacements[num_redists];
  compute_displ(data, num_redists, displacements, comm);
  if (cache_size > 0)
  {
    size_t cache_index
      = lookup_cache_index(num_redists, displacements,
                           (MPI_Aint (*)[num_redists])cache->displacements,
                           cache_size);
    MPI_Datatype *dt;
    if (cache_index == cache_size)
    {
      cache_index = cache->token;
      dt = cache->dt + (size_t)num_messages * cache_index;
      create_all_dt_for_dir(msgs, num_messages, num_redists,
                            displacements, dt, comm);
      memcpy(cache->displacements + cache_index * num_redists, displacements,
             sizeof (displacements));
      cache->token = (cache->token + 1) % cache_size;
    }
    else
      dt = cache->dt + (size_t)num_messages * cache_index;
    return dt;
  }
  else
  {
    create_all_dt_for_dir(msgs, num_messages, num_redists,
                          displacements, temp_dt, comm);
    return temp_dt;
  }
}

static void clear_temp_dt(int num_messages,
                          MPI_Datatype dt[num_messages], MPI_Comm comm) {

  for (int i = 0; i < num_messages; ++i)
    xt_mpi_call(MPI_Type_free(dt + i), comm);
}

static inline Xt_redist_collection
xrc(void *redist)
{
  return (Xt_redist_collection)redist;
}

static void
redist_collection_s_exchange(Xt_redist redist, unsigned num_arrays,
                             void **src_data, void **dst_data) {

  Xt_redist_collection redist_coll = xrc(redist);

  if (num_arrays != redist_coll->num_redists)
    Xt_abort(redist_coll->comm, "ERROR: wrong number of arrays in "
             "redist_collection_s_exchange", __FILE__, __LINE__);

  MPI_Request * recv_requests;
  MPI_Datatype *temp_dt_src, *temp_dt_dst;
  size_t num_dt = (size_t)redist_coll->ndst + (size_t)redist_coll->nsrc;
  {
    size_t rq_size = (size_t)redist_coll->ndst * sizeof (*recv_requests);
    size_t temp_dt_size = num_dt * sizeof (MPI_Datatype);
    size_t block_ofs = sizeRoundUp(rq_size, sizeof (MPI_Datatype));
    size_t allocSize = block_ofs + temp_dt_size;
    recv_requests = xmalloc(allocSize);
    temp_dt_dst = (void *)((unsigned char *)recv_requests + block_ofs);
    temp_dt_src = temp_dt_dst + (size_t)redist_coll->ndst;
  }


  size_t cache_size = redist_coll->cache_size;
  if (cache_size == 0)
    for (size_t i = 0; i < num_dt; ++i)
      temp_dt_dst[i] = MPI_DATATYPE_NULL;

  MPI_Datatype *dst_dt
    = get_compound_datatype(dst_data, redist_coll->recv_msgs,
                            redist_coll->ndst, num_arrays,
                            &redist_coll->dst_cache, cache_size,
                            temp_dt_dst, redist_coll->comm);

  for (int i = 0; i < redist_coll->ndst; ++i)
    xt_mpi_call(MPI_Irecv(
                  dst_data[0], 1,
                  dst_dt[i],
                  redist_coll->recv_msgs[i].rank, 0, redist_coll->comm,
                  recv_requests+i), redist_coll->comm);

  MPI_Datatype *src_dt
    = get_compound_datatype(src_data, redist_coll->send_msgs,
                            redist_coll->nsrc, num_arrays,
                            &redist_coll->src_cache, cache_size,
                            temp_dt_src, redist_coll->comm);

  for (int i = 0; i < redist_coll->nsrc; ++i)
    xt_mpi_call(MPI_Send(
                  src_data[0], 1,
                  src_dt[i],
                  redist_coll->send_msgs[i].rank, 0, redist_coll->comm),
                redist_coll->comm);

  xt_mpi_call(MPI_Waitall(redist_coll->ndst, recv_requests,
                          MPI_STATUSES_IGNORE), redist_coll->comm);

  if (cache_size == 0) {
    clear_temp_dt(redist_coll->ndst + redist_coll->nsrc,
                  temp_dt_dst, redist_coll->comm);
  }

  free(recv_requests);
}

static void
free_redist_collection_msgs(struct redist_collection_msg * msgs,
                            int nmsgs, unsigned num_redists,
                            MPI_Comm comm) {

  for (int i = 0; i < nmsgs; ++i) {

    for (unsigned j = 0; j < num_redists; ++j)
      if (msgs[i].component_dt[j] != MPI_DATATYPE_NULL)
        xt_mpi_call(MPI_Type_free(msgs[i].component_dt+j), comm);
    free(msgs[i].component_dt);
  }
  free(msgs);
}

static void
redist_collection_delete(Xt_redist redist) {

  Xt_redist_collection redist_coll = xrc(redist);

  free_redist_collection_msgs(redist_coll->send_msgs, redist_coll->nsrc,
                              redist_coll->num_redists,
                              redist_coll->comm);

  free_redist_collection_msgs(redist_coll->recv_msgs, redist_coll->ndst,
                              redist_coll->num_redists,
                              redist_coll->comm);

  destruct_cache(&redist_coll->src_cache, redist_coll->cache_size,
                 (size_t)redist_coll->nsrc, redist_coll->comm);
  destruct_cache(&redist_coll->dst_cache, redist_coll->cache_size,
                 (size_t)redist_coll->ndst, redist_coll->comm);

  xt_mpi_call(MPI_Comm_free(&(redist_coll->comm)), MPI_COMM_WORLD);

  free(redist_coll);
}

static MPI_Datatype
redist_collection_get_send_MPI_Datatype(Xt_redist redist, int XT_UNUSED(rank))
{
  Xt_redist_collection redist_coll = xrc(redist);

  Xt_abort(redist_coll->comm, "ERROR: get_send_MPI_Datatype is not"
           " supported for this xt_redist type (Xt_redist_collection)",
           __FILE__, __LINE__);

  return MPI_DATATYPE_NULL;
}

static MPI_Datatype
redist_collection_get_recv_MPI_Datatype(Xt_redist redist, int XT_UNUSED(rank)) {

  Xt_redist_collection redist_coll = xrc(redist);

  Xt_abort(redist_coll->comm, "ERROR: get_recv_MPI_Datatype is not"
           " supported for this xt_redist type (Xt_redist_collection)",
           __FILE__, __LINE__);

  return MPI_DATATYPE_NULL;
}

static void
redist_collection_s_exchange1(Xt_redist redist, void *src_data, void *dst_data)
{

  Xt_redist_collection redist_coll = xrc(redist);
  if (redist_coll->num_redists == 1U)
    redist_collection_s_exchange(redist, 1U, &src_data, &dst_data);
  else
    Xt_abort(redist_coll->comm, "ERROR: s_exchange1 is not implemented for"
             " this xt_redist type (Xt_redist_collection)", __FILE__, __LINE__);
}

static MPI_Comm
redist_collection_get_MPI_Comm(Xt_redist redist) {

  Xt_redist_collection redist_coll = xrc(redist);

  return redist_coll->comm;
}
