/**
 * @file xt_xmap_dist_dir.c
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
 * URL: https://doc.redmine.dkrz.de/yaxt/html/
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

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <mpi.h>

#include "xt/xt_idxlist.h"
#include "xt/xt_idxlist_collection.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_xmap.h"
#include "xt/xt_xmap_dist_dir.h"
#include "xt/xt_mpi.h"
#include "xt_arithmetic_util.h"
#include "xt_idxstripes_internal.h"
#include "xt_mpi_internal.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "ensure_array_size.h"
#include "xt/xt_xmap_intersection.h"
#include "instr.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

enum {
  SEND_SIZE_SRC = 0,
  SEND_SIZE_DST = 1,
  SEND_NUM_SRC = 2,
  SEND_NUM_DST = 3,
  SEND_SIZE_ASIZE,
};

static int
xt_com_list_rank_cmp(const void *a_, const void *b_)
{
  const struct Xt_com_list *a = a_, *b = b_;
  return a->rank - b->rank;
}

struct dist_dir {
  int num_entries;
  struct Xt_com_list entries[];
};


static inline Xt_int
get_dist_dir_global_interval_size(Xt_idxlist idxlist, MPI_Comm comm,
                                  int comm_size) {

  unsigned long long num_indices, global_num_indices;

  num_indices = (unsigned)xt_idxlist_get_num_indices(idxlist);

  xt_mpi_call(MPI_Allreduce(&num_indices, &global_num_indices, 1,
                            MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm), comm);

  return (Xt_int)(((global_num_indices + (unsigned)comm_size - 1)
                   / (unsigned)comm_size) * (unsigned)comm_size);
}

static inline Xt_int get_min_idxlist_index(Xt_idxlist a, Xt_idxlist b) {

  int num_a = xt_idxlist_get_num_indices(a),
    num_b = xt_idxlist_get_num_indices(b);
  Xt_int min_index_a = num_a ? xt_idxlist_get_min_index(a) : XT_INT_MAX,
    min_index_b = num_b ? xt_idxlist_get_min_index(b) : XT_INT_MAX,
    min_index = (Xt_int)MIN(min_index_a, min_index_b);
  return min_index;
}

static inline Xt_int get_max_idxlist_index(Xt_idxlist a, Xt_idxlist b) {

  int num_a = xt_idxlist_get_num_indices(a),
    num_b = xt_idxlist_get_num_indices(b);
  Xt_int max_index_a = num_a ? xt_idxlist_get_max_index(a) : XT_INT_MIN,
    max_index_b = num_b ? xt_idxlist_get_max_index(b) : XT_INT_MIN,
    max_index = MAX(max_index_a, max_index_b);
  return max_index;
}

struct bucket_params {
  Xt_int global_interval, local_interval,
    local_index_range_lbound, local_index_range_ubound;
  int num_src, num_dst;
};

static struct bucket_params
get_bucket_params(Xt_idxlist src_idxlist,
                  Xt_idxlist dst_idxlist, MPI_Comm comm, int comm_size)
{
  Xt_int global_interval
    = get_dist_dir_global_interval_size(src_idxlist, comm, comm_size);
  /* guard vs. comm_size being larger than number of indices */
  Xt_int local_interval = MAX((Xt_int)1, (Xt_int)(global_interval / comm_size));
  int num_src = xt_idxlist_get_num_indices(src_idxlist),
    num_dst = xt_idxlist_get_num_indices(dst_idxlist);
  Xt_int local_index_range_lbound
    = get_min_idxlist_index(src_idxlist, dst_idxlist);
  Xt_int local_index_range_ubound
    = get_max_idxlist_index(src_idxlist, dst_idxlist);
  return (struct bucket_params){
    .global_interval = global_interval,
    .local_interval = local_interval,
    .num_src = num_src,
    .num_dst = num_dst,
    .local_index_range_lbound = local_index_range_lbound,
    .local_index_range_ubound = local_index_range_ubound,
  };
}

static Xt_idxlist
get_bucket(const struct bucket_params *bucket_params,
           struct Xt_stripe **stripes_, size_t *stripes_array_size,
           int dist_dir_rank)
{
  Xt_int global_interval = bucket_params->global_interval;
  Xt_int local_interval = bucket_params->local_interval;
  Xt_int local_index_range_lbound = bucket_params->local_index_range_lbound;
  Xt_int local_index_range_ubound = bucket_params->local_index_range_ubound;
  int num_stripes = 0;

  /* find first index in bucket of dist_dir_rank
     <= local_index_range_lbound */
  Xt_int start = (Xt_int)(0 + dist_dir_rank * local_interval);
  {
    long long start_correction
      = (long long)local_index_range_lbound - (long long)start;
    Xt_int corr_steps
      = (Xt_int)((start_correction
                  - (llsign_mask(start_correction)
                     & (long long)(global_interval - 1)))
                 / (long long)global_interval);
    start = (Xt_int)(start + corr_steps * global_interval);
  }
  /* next find last stripe in bucket of dist_dir_rank
   * <= local_index_range_ubound */
  Xt_int end
    = (Xt_int)(start
               + (((long long)local_index_range_ubound - (long long)start)
                  / global_interval) * global_interval);
  Xt_int use_start_stripe
    = (Xt_int)(start + local_interval > local_index_range_lbound);
  num_stripes = (int)(((long long)end - (long long)start)/global_interval)
    + (int)use_start_stripe;
  start = (Xt_int)(start
                   + ((Xt_int)((Xt_uint)use_start_stripe - (Xt_uint)1)
                      & global_interval));
  if (!num_stripes)
    return NULL;

  struct Xt_stripe *restrict stripes = *stripes_;
  ENSURE_ARRAY_SIZE(stripes, *stripes_array_size, (size_t)num_stripes);

  for (int j = 0; j < num_stripes; ++j) {

    stripes[j].start = (Xt_int)(start + j * global_interval);
    stripes[j].stride = 1;
    stripes[j].nstrides = (int)local_interval;
  }

  *stripes_ = stripes;
  return xt_idxstripes_prealloc_new(stripes, num_stripes);
}

struct idxlist_dst_src_pair {
  Xt_idxlist src, dst;
};

static void
rank_no_send(size_t rank, int (*restrict send_size)[SEND_SIZE_ASIZE],
             struct idxlist_dst_src_pair *restrict send_list)
{
  send_size[rank][SEND_SIZE_SRC] = 0;
  send_size[rank][SEND_NUM_SRC] = 0;
  send_size[rank][SEND_SIZE_DST] = 0;
  send_size[rank][SEND_NUM_DST] = 0;
  send_list[rank].src = NULL;
  send_list[rank].dst = NULL;
}

static int
compute_and_pack_bucket_intersections(struct bucket_params *bucket_params,
                                      Xt_idxlist src_idxlist,
                                      Xt_idxlist dst_idxlist,
                                      int (*send_size)[SEND_SIZE_ASIZE],
                                      void ** send_buffer,
                                      MPI_Comm comm, int comm_size) {


  int num_msg = 0;
  if (bucket_params->local_index_range_lbound
      <= bucket_params->local_index_range_ubound) {
    int send_buffer_size = 0;
    struct idxlist_dst_src_pair *send_list
      = xmalloc((size_t)comm_size * sizeof(*send_list));
    struct Xt_stripe *stripes = NULL;
    size_t stripes_array_size = 0;
    size_t i = 0;
    /* is it impossible for early buckets to overlap our lists? */
    if (bucket_params->local_index_range_lbound >= 0
        && (bucket_params->local_index_range_ubound
            < bucket_params->global_interval)) {
      size_t non_overlapping_bucket_prefix
        = bucket_params->local_index_range_lbound
        / bucket_params->local_interva;
      for (; i < non_overlapping_bucket_prefix; ++i)
        rank_no_send(i, send_size, send_list);
    }
    /* is it impossible for later ranks to overlap our lists? */
    size_t start_of_non_overlapping_bucket_suffix
      = (bucket_params->local_index_range_ubound
         + bucket_params->local_interval - 1) / bucket_params->local_interval;
    if (bucket_params->local_index_range_lbound < 0
        || start_of_non_overlapping_bucket_suffix > comm_size)
      start_of_non_overlapping_bucket_suffix = comm_size;

    for (; i < (size_t)start_of_non_overlapping_bucket_suffix; ++i) {
      Xt_idxlist bucket
        = get_bucket(bucket_params, &stripes, &stripes_array_size, (int)i);
      if (bucket) {
        send_list[i].src = xt_idxlist_get_intersection(src_idxlist, bucket);
        if (xt_idxlist_get_num_indices(send_list[i].src) > 0)
        {
          send_size[i][SEND_SIZE_SRC]
            = (int)xt_idxlist_get_pack_size(send_list[i].src, comm);
          send_size[i][SEND_NUM_SRC] = 1;
          ++num_msg;
        }
        else
        {
          send_size[i][SEND_SIZE_SRC] = 0;
          send_size[i][SEND_NUM_SRC] = 0;
          xt_idxlist_delete(send_list[i].src);
          send_list[i].src = NULL;
        }

        send_list[i].dst = xt_idxlist_get_intersection(bucket, dst_idxlist);
        xt_idxlist_delete(bucket);

        if (xt_idxlist_get_num_indices(send_list[i].dst) > 0)
        {
          send_size[i][SEND_SIZE_DST]
            = (int)xt_idxlist_get_pack_size(send_list[i].dst, comm);
          send_size[i][SEND_NUM_DST] = 1;
          ++num_msg;
        }
        else
        {
          send_size[i][SEND_SIZE_DST] = 0;
          send_size[i][SEND_NUM_DST] = 0;
          xt_idxlist_delete(send_list[i].dst);
          send_list[i].dst = NULL;
        }

      } else
        rank_no_send(i, send_size, send_list);

      send_buffer_size
        += send_size[i][SEND_SIZE_SRC] + send_size[i][SEND_SIZE_DST];
    }
    for (; i < (size_t)comm_size; ++i)
      rank_no_send(i, send_size, send_list);

    *send_buffer = xrealloc(stripes, (size_t)send_buffer_size);

    int position = 0;

    for (int i = 0; i < comm_size; ++i) {

      if (send_size[i][SEND_SIZE_SRC] > 0)
      {
        xt_idxlist_pack(send_list[i].src, *send_buffer, send_buffer_size,
                        &position, comm);
        xt_idxlist_delete(send_list[i].src);
      }
    }

    for (int i = 0; i < comm_size; ++i) {

      if (send_size[i][SEND_SIZE_DST] > 0)
      {
        xt_idxlist_pack(send_list[i].dst, *send_buffer, send_buffer_size,
                        &position, comm);
        xt_idxlist_delete(send_list[i].dst);
      }
    }

    free(send_list);
  } else
    memset(send_size, 0, (size_t)comm_size * sizeof (*send_size));
  return num_msg;
}

static void
get_remote_packed_intersection_size(int recv_size[SEND_SIZE_ASIZE],
                                    int (*send_size)[SEND_SIZE_ASIZE],
                                    MPI_Comm comm) {

#if MPI_VERSION > 2 || ( MPI_VERSION == 2 && MPI_SUBVERSION >= 2)
  xt_mpi_call(MPI_Reduce_scatter_block((int *)send_size, (int *)recv_size,
                                       SEND_SIZE_ASIZE, MPI_INT, MPI_SUM,
                                       comm), comm);
#else
  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int *recv_count = xmalloc((size_t)comm_size * sizeof(*recv_count));
  for (int i = 0; i < comm_size; ++i) recv_count[i] = SEND_SIZE_ASIZE;

  xt_mpi_call(MPI_Reduce_scatter(send_size, recv_size, recv_count, MPI_INT,
                                 MPI_SUM, comm), comm);

  free(recv_count);
#endif
}

static void send_intersections(void * send_buffer,
                               int (*send_size)[SEND_SIZE_ASIZE],
                               MPI_Request *dir_init_send_requests,
                               int tag_offset, MPI_Comm comm, int comm_size) {
  int offset = 0;
  size_t reqOfs = 0;

  // pack the intersections into the send buffer
  for (int i = 0; i < comm_size; ++i)
    if (send_size[i][SEND_SIZE_SRC] > 0) {
      xt_mpi_call(MPI_Isend((char *)send_buffer + offset,
                            send_size[i][SEND_SIZE_SRC],
                            MPI_PACKED, i, tag_offset
                            + xt_mpi_tag_xmap_dist_dir_src_send,
                            comm, dir_init_send_requests + reqOfs),
                  comm);
      ++reqOfs;
      offset += send_size[i][SEND_SIZE_SRC];
    }

  for (int i = 0; i < comm_size; ++i)
    if (send_size[i][SEND_SIZE_DST] > 0) {
      xt_mpi_call(MPI_Isend((char *)send_buffer + offset,
                            send_size[i][SEND_SIZE_DST],
                            MPI_PACKED, i, tag_offset
                            + xt_mpi_tag_xmap_dist_dir_dst_send,
                            comm, dir_init_send_requests + reqOfs),
                  comm);
      ++reqOfs;
      offset += send_size[i][SEND_SIZE_DST];
    }
}

static void
recv_and_unpack_intersection(struct dist_dir *dist_dir, int recv_size,
                             int recv_count, void * recv_buffer, int tag,
                             MPI_Comm comm) {

  // initialize distributed directories
  int total_recv_size = 0;

  for (int i = 0; i < recv_count; ++i)
  {
    MPI_Status status;

    xt_mpi_call(MPI_Recv(recv_buffer, recv_size, MPI_PACKED, MPI_ANY_SOURCE,
                         tag, comm, &status), comm);

    int received_count;
    xt_mpi_call(MPI_Get_count(&status, MPI_PACKED, &received_count), comm);

    int position = 0;

    dist_dir->entries[i].rank = status.MPI_SOURCE;
    dist_dir->entries[i].list =
      xt_idxlist_unpack(recv_buffer, received_count, &position, comm);

    total_recv_size += received_count;
  }

  if (total_recv_size != recv_size)
    Xt_abort(comm, "ERROR: recv_intersections received wrong number of bytes",
             __FILE__, __LINE__);
  dist_dir->num_entries = recv_count;
}

static void
recv_and_unpack_intersections(int recv_size[SEND_SIZE_ASIZE],
                              struct dist_dir **src_dist_dir,
                              struct dist_dir **dst_dist_dir,
                              int tag_offset, MPI_Comm comm) {

  *src_dist_dir = xmalloc(sizeof (struct dist_dir)
                          + (sizeof (struct Xt_com_list)
                             * (size_t)recv_size[SEND_NUM_SRC]));
  *dst_dist_dir = xmalloc(sizeof (struct dist_dir)
                          + (sizeof (struct Xt_com_list)
                             * (size_t)recv_size[SEND_NUM_DST]));


  void * recv_buffer = xmalloc((size_t)MAX(recv_size[SEND_SIZE_SRC],
                                           recv_size[SEND_SIZE_DST]));

  recv_and_unpack_intersection(*src_dist_dir, recv_size[SEND_SIZE_SRC],
                               recv_size[SEND_NUM_SRC], recv_buffer,
                               tag_offset + xt_mpi_tag_xmap_dist_dir_src_send,
                               comm);
  recv_and_unpack_intersection(*dst_dist_dir, recv_size[SEND_SIZE_DST],
                               recv_size[SEND_NUM_DST], recv_buffer,
                               tag_offset + xt_mpi_tag_xmap_dist_dir_dst_send,
                               comm);

  free(recv_buffer);
}


struct isect
{
  int rank;
  int dst_list, src_list;
  Xt_idxlist idxlist;
};

static size_t
match_src_dst_dist_dirs(struct dist_dir * src_dist_dir,
                        struct dist_dir * dst_dist_dir,
                        struct isect **src_dst_intersections)
{
  struct isect (*src_dst_intersections_)
    = (*src_dst_intersections)
    = xmalloc((size_t)src_dist_dir->num_entries
              * (size_t)dst_dist_dir->num_entries
              * sizeof(**src_dst_intersections));
  size_t isect_fill = 0;

  for (int i = 0; i < src_dist_dir->num_entries; ++i)
    for (int j = 0; j < dst_dist_dir->num_entries; ++j)
    {
      Xt_idxlist intersection = xt_idxlist_get_intersection(
        src_dist_dir->entries[i].list, dst_dist_dir->entries[j].list);
      if (xt_idxlist_get_num_indices(intersection) > 0)
      {
        src_dst_intersections_[isect_fill]
          = (struct isect){ .dst_list = j, .src_list = i,
                            .idxlist = intersection };
        ++isect_fill;
      }
      else
      {
        xt_idxlist_delete(intersection);
      }
    }
  *src_dst_intersections
    = xrealloc(src_dst_intersections_,
               isect_fill * sizeof (*src_dst_intersections_));
  return isect_fill;
}

static int
cmp_isect_rank(const void *a_, const void *b_)
{
  const struct isect *a = a_, *b = b_;
  /* this is safe vs. overflow because ranks are in [0..MAX_INT) */
  return a->rank - b->rank;
}

static int
match_and_pack_src_dst_dist_dirs(struct dist_dir * src_dist_dir,
                                 struct dist_dir * dst_dist_dir,
                                 int (*send_size)[SEND_SIZE_ASIZE],
                                 void ** send_buffer,
                                 MPI_Comm comm) {

  struct isect *src_dst_intersections;

  size_t num_intersections
    = match_src_dst_dist_dirs(src_dist_dir, dst_dist_dir,
                              &src_dst_intersections);

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int total_send_size = 0;

  int rank_pack_size;

  xt_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &rank_pack_size), comm);

  for (int i = 0; i < comm_size; ++i)
    send_size[i][SEND_SIZE_SRC] = 0, send_size[i][SEND_SIZE_DST] = 0,
      send_size[i][SEND_NUM_SRC] = 0, send_size[i][SEND_NUM_DST] = 0;

  for (size_t i = 0; i < num_intersections; ++i)
  {
    int msg_size = rank_pack_size
      + (int)xt_idxlist_get_pack_size(src_dst_intersections[i].idxlist,
                                      comm);
    int src = src_dst_intersections[i].src_list,
        dst = src_dst_intersections[i].dst_list;
    send_size[src_dist_dir->entries[src].rank][SEND_SIZE_SRC] += msg_size;
    ++(send_size[src_dist_dir->entries[src].rank][SEND_NUM_SRC]);
    send_size[dst_dist_dir->entries[dst].rank][SEND_SIZE_DST] += msg_size;
    ++(send_size[dst_dist_dir->entries[dst].rank][SEND_NUM_DST]);
    total_send_size += 2*msg_size;
    src_dst_intersections[i].rank = src_dist_dir->entries[src].rank;
  }

  /* sort intersections by src rank */
  qsort(src_dst_intersections, num_intersections,
        sizeof (src_dst_intersections[0]), cmp_isect_rank);

  (*send_buffer) = xmalloc((size_t)total_send_size);

  int position = 0;
  int prev_send_rank = -1;
  int num_send_indices_requests = 0;
  for (size_t i = 0; i < num_intersections; ++i)
  {
    int dst = src_dst_intersections[i].dst_list;
    int dst_rank = dst_dist_dir->entries[dst].rank;

    /* see if this generates a new request? */
    num_send_indices_requests
      += src_dst_intersections[i].rank != prev_send_rank;
    prev_send_rank = src_dst_intersections[i].rank;

    // pack rank
    xt_mpi_call(MPI_Pack(&dst_rank, 1, MPI_INT,
                         (*send_buffer), total_send_size, &position,
                         comm), comm);
    // pack intersection
    xt_idxlist_pack(src_dst_intersections[i].idxlist, (*send_buffer),
                    total_send_size, &position, comm);

    src_dst_intersections[i].rank = dst_rank;
  }

  /* sort intersections by dst rank */
  qsort(src_dst_intersections, num_intersections,
        sizeof (src_dst_intersections[0]), cmp_isect_rank);

  prev_send_rank = -1;
  for (size_t i = 0; i < num_intersections; ++i)
  {
    int src = src_dst_intersections[i].src_list;
    int src_rank = src_dist_dir->entries[src].rank;

    /* see if this generates a new request? */
    num_send_indices_requests
      += src_dst_intersections[i].rank != prev_send_rank;
    prev_send_rank = src_dst_intersections[i].rank;

    // pack rank
    xt_mpi_call(MPI_Pack(&src_rank, 1, MPI_INT,
                         (*send_buffer), total_send_size, &position,
                         comm), comm);
    // pack intersection
    xt_idxlist_pack(src_dst_intersections[i].idxlist, (*send_buffer),
                    total_send_size, &position, comm);
    /* after this step, the idxlist is no longer needed */
    xt_idxlist_delete(src_dst_intersections[i].idxlist);
  }

  free(src_dst_intersections);
  return num_send_indices_requests;
}

static void generate_distributed_directories(struct dist_dir **src_dist_dir,
                                             struct dist_dir **dst_dist_dir,
                                             Xt_idxlist src_idxlist,
                                             Xt_idxlist dst_idxlist,
                                             int tag_offset,
                                             MPI_Comm comm) {

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  struct bucket_params bucket_params
    = get_bucket_params(src_idxlist, dst_idxlist, comm, comm_size);

  void *send_buffer = NULL;

  int (*send_size)[SEND_SIZE_ASIZE]
    = xmalloc((size_t)comm_size * sizeof(*send_size));

  int num_msg
    = compute_and_pack_bucket_intersections(&bucket_params,
                                            src_idxlist, dst_idxlist,
                                            send_size, &send_buffer,
                                            comm, comm_size);

  int recv_size[SEND_SIZE_ASIZE]; // for src and dst

  get_remote_packed_intersection_size(recv_size, send_size, comm);

  MPI_Request *dir_init_send_requests
    = xmalloc((size_t)num_msg * sizeof(*dir_init_send_requests));
  send_intersections(send_buffer, send_size, dir_init_send_requests,
                     tag_offset, comm, comm_size);

  free(send_size);

  recv_and_unpack_intersections(recv_size, src_dist_dir, dst_dist_dir,
                                tag_offset, comm);

  // wait for the sends to be completed
  xt_mpi_call(MPI_Waitall(num_msg, dir_init_send_requests,
                          MPI_STATUSES_IGNORE), comm);
  free(dir_init_send_requests);
  free(send_buffer);
}

static void
recv_and_unpack_dist_dir_result(struct dist_dir *dist_dir, int recv_size,
                                void *restrict recv_buffer, int tag,
                                MPI_Comm comm)
{

  // initiate distributed directories
  int num_entries = 0;
  while (recv_size > 0) {

    MPI_Status status;

    xt_mpi_call(MPI_Recv(recv_buffer, recv_size, MPI_PACKED,
                         MPI_ANY_SOURCE, tag, comm, &status), comm);

    int received_count;
    xt_mpi_call(MPI_Get_count(&status, MPI_PACKED, &received_count), comm);

    recv_size -= received_count;

    int position = 0;

    while (received_count > position) {

      xt_mpi_call(MPI_Unpack(recv_buffer, received_count, &position,
                             &dist_dir->entries[num_entries].rank,
                             1, MPI_INT, comm), comm);

      dist_dir->entries[num_entries].list =
        xt_idxlist_unpack(recv_buffer, received_count, &position, comm);

      ++num_entries;
    }
  }
  qsort(dist_dir->entries, (size_t)num_entries, sizeof(*dist_dir->entries),
        xt_com_list_rank_cmp);

  if (0 != recv_size)
    Xt_abort(comm, "ERROR: recv_and_unpack_dist_dir_result"
             " received wrong number of bytes", __FILE__, __LINE__);

  dist_dir->num_entries = num_entries;

}

static int
stripe_cmp(const void *a, const void *b)
{
  typedef const struct Xt_stripe *csx;
  return (((csx)a)->start > ((csx)b)->start)
    - (((csx)b)->start > ((csx)a)->start);
}

/*
 * @param dist_dir_results contains the intersections of this ranks
 * dst or src idxlist with other ranks in bucket-sized chunks, the
 * chunks belonging to the same communication partner are merged
 * in-place here, i.e. on return (*dist_dir_results)->num_entries is
 * less than or equal to the previous count
 */
static void
merge_bucket_intersections(struct dist_dir **dist_dir_results) {

  struct Xt_com_list *restrict entries = (*dist_dir_results)->entries;
  size_t num_isect_agg = 0;

  size_t i = 0, num_shards = (size_t)(*dist_dir_results)->num_entries;
  while (i < num_shards) {
    int rank = entries[i].rank;
    size_t j = i + 1;
    /* find all entries matching the currently considered rank */
    while (j < num_shards && entries[j].rank == rank)
      ++j;

    size_t num_intersections_per_rank = j - i;

    if (num_intersections_per_rank > 0) {
      struct Xt_stripe *restrict stripes = NULL;
      size_t num_stripes = 0;
      for (; i < j; ++i) {
        struct Xt_stripe *stripes_of_intersection;
        int num_stripes_of_intersection;
        xt_idxlist_get_index_stripes(entries[i].list,
                                     &stripes_of_intersection,
                                     &num_stripes_of_intersection);
        xt_idxlist_delete(entries[i].list);
        if (stripes) {
          stripes = xrealloc(stripes,
                             (num_stripes + (size_t)num_stripes_of_intersection)
                             * sizeof (*stripes));
          memcpy(stripes + num_stripes, stripes_of_intersection,
                 (size_t)num_stripes_of_intersection * sizeof (*stripes));
          free(stripes_of_intersection);
        } else
          stripes = stripes_of_intersection;
        num_stripes += (size_t)num_stripes_of_intersection;
      }
      qsort(stripes, num_stripes, sizeof (*stripes), stripe_cmp);
      entries[num_isect_agg].list
        = xt_idxstripes_new(stripes, (int)num_stripes);
      free(stripes);
      entries[num_isect_agg].rank = rank;
      ++num_isect_agg;
    }
  }
  (*dist_dir_results)->num_entries = (int)num_isect_agg;
  *dist_dir_results = xrealloc(*dist_dir_results, sizeof (struct dist_dir)
                               + (size_t)num_isect_agg
                               * sizeof(struct Xt_com_list));
}

static void free_dist_dir(struct dist_dir * dist_dir) {

  for (int i = 0; i < dist_dir->num_entries; ++i)
    xt_idxlist_delete(dist_dir->entries[i].list);
  free(dist_dir);
}

static bool
test_mpi_req(int *restrict num_req_,
             MPI_Request *restrict req,
             int *restrict ops_completed, MPI_Comm comm)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int done_count;
  size_t num_req = (size_t)*num_req_;
  xt_mpi_call(MPI_Testsome(*num_req_, req, &done_count, ops_completed,
                           MPI_STATUSES_IGNORE), comm);
  if (done_count != MPI_UNDEFINED) {
    if (num_req > (size_t)done_count) {
      size_t j = num_req;
      for (size_t i = 0;
           i < (size_t)done_count && j >= num_req - (size_t)done_count;
           ++i)
        if (ops_completed[i] < (int)num_req - done_count) {
          while (req[--j] == MPI_REQUEST_NULL);
          req[ops_completed[i]] = req[j];
        }
      num_req -= (size_t)done_count;
    }
    else
      num_req = 0;
  }
  *num_req_ = (int)num_req;
  return num_req == 0;
}


static void
recv_and_unpack_dist_dir_results(int recv_size[SEND_SIZE_ASIZE],
                                 struct dist_dir **src_intersections,
                                 struct dist_dir **dst_intersections,
                                 int *num_send_indices_requests,
                                 MPI_Request *send_indices_requests,
                                 int tag_offset,
                                 MPI_Comm comm) {

  struct dist_dir *src_dist_dir_results
    = xmalloc(sizeof (struct dist_dir)
              + (sizeof (struct Xt_com_list)
                 * (size_t)recv_size[SEND_NUM_SRC])),
    *dst_dist_dir_results
    = xmalloc(sizeof (struct dist_dir)
              + (sizeof (struct Xt_com_list)
                 * (size_t)recv_size[SEND_NUM_DST]));

  void *recv_buffer = xmalloc((size_t)MAX(recv_size[SEND_SIZE_SRC],
                                          recv_size[SEND_SIZE_DST]));

  recv_and_unpack_dist_dir_result(src_dist_dir_results,
                                  recv_size[SEND_SIZE_SRC],
                                  recv_buffer, tag_offset
                                  + xt_mpi_tag_xmap_dist_dir_src_send, comm);
  assert(src_dist_dir_results->num_entries
         == recv_size[SEND_NUM_SRC]);

  enum { ops_completed_auto_size = 16 };
  int ops_completed_auto[ops_completed_auto_size];
  int *ops_completed
    = *num_send_indices_requests > ops_completed_auto_size
    ? xmalloc((size_t)*num_send_indices_requests * sizeof (*ops_completed))
    : ops_completed_auto;
  bool all_sends_done
    = test_mpi_req(num_send_indices_requests, send_indices_requests,
                   ops_completed, comm);

  recv_and_unpack_dist_dir_result(dst_dist_dir_results,
                                  recv_size[SEND_SIZE_DST],
                                  recv_buffer, tag_offset
                                  + xt_mpi_tag_xmap_dist_dir_dst_send, comm);
  assert(dst_dist_dir_results->num_entries
         == recv_size[SEND_NUM_DST]);

  if (!all_sends_done)
    all_sends_done
      = test_mpi_req(num_send_indices_requests, send_indices_requests,
                     ops_completed, comm);
  free(recv_buffer);

  merge_bucket_intersections(&src_dist_dir_results);
  *src_intersections = src_dist_dir_results;

  if (!all_sends_done)
    all_sends_done
      = test_mpi_req(num_send_indices_requests, send_indices_requests,
                     ops_completed, comm);

  merge_bucket_intersections(&dst_dist_dir_results);
  *dst_intersections = dst_dist_dir_results;
  if (ops_completed != ops_completed_auto) free(ops_completed);
}

static void exchange_idxlists(struct dist_dir **src_intersections,
                              struct dist_dir **dst_intersections,
                              Xt_idxlist src_idxlist,
                              Xt_idxlist dst_idxlist,
                              int tag_offset,
                              MPI_Comm comm) {

  int comm_size;

  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  struct dist_dir *src_dist_dir, *dst_dist_dir;

  generate_distributed_directories(&src_dist_dir, &dst_dist_dir,
                                   src_idxlist, dst_idxlist,
                                   tag_offset, comm);

  void * send_buffer;

  int recv_size[SEND_SIZE_ASIZE], (*send_size)[SEND_SIZE_ASIZE]
    = xmalloc((size_t)comm_size * sizeof(*send_size));

  // match the source and destination entries in the local distributed
  // directories and pack the results into a sendable format
  int num_send_indices_requests
    = match_and_pack_src_dst_dist_dirs(src_dist_dir, dst_dist_dir, send_size,
                                       &send_buffer, comm);
  free_dist_dir(src_dist_dir);
  free_dist_dir(dst_dist_dir);

  // get the data size the local process will receive from other processes
  get_remote_packed_intersection_size(recv_size, send_size, comm);

  MPI_Request *send_indices_requests
    = xmalloc((size_t)num_send_indices_requests
              * sizeof(*send_indices_requests));

  send_intersections(send_buffer, send_size, send_indices_requests,
                     tag_offset, comm, comm_size);

  recv_and_unpack_dist_dir_results(recv_size,
                                   src_intersections, dst_intersections,
                                   &num_send_indices_requests,
                                   send_indices_requests,
                                   tag_offset, comm);

  xt_mpi_call(MPI_Waitall(num_send_indices_requests, send_indices_requests,
                          MPI_STATUSES_IGNORE), comm);

  free(send_buffer);
  free(send_size);
  free(send_indices_requests);
}

Xt_xmap xt_xmap_dist_dir_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                             MPI_Comm comm) {

  INSTR_DEF(this_instr,"xt_xmap_all2all_new")
  INSTR_START(this_instr);

  int tag_offset;
  MPI_Comm newcomm = xt_mpi_comm_smart_dup(comm, &tag_offset);

  struct dist_dir *src_intersections, *dst_intersections;

  exchange_idxlists(&src_intersections, &dst_intersections,
                    src_idxlist, dst_idxlist, tag_offset, newcomm);

  Xt_xmap xmap
    = xt_xmap_intersection_new(src_intersections->num_entries,
                               src_intersections->entries,
                               dst_intersections->num_entries,
                               dst_intersections->entries,
                               src_idxlist, dst_idxlist, newcomm);

  xt_mpi_comm_smart_dedup(&newcomm, tag_offset);

  free_dist_dir(src_intersections);
  free_dist_dir(dst_intersections);
  INSTR_STOP(this_instr);
  return xmap;
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
