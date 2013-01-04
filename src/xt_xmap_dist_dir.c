/**
 * @file xt_xmap_dist_dir.c
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
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <mpi.h>

#include "xt/xt_idxlist.h"
#include "xt/xt_idxlist_collection.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_xmap.h"
#include "xt/xt_xmap_dist_dir.h"
#include "xt/xt_mpi.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "ensure_array_size.h"
#include "xt/xt_xmap_intersection.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

enum {
  SEND_SIZE_SRC = 0,
  SEND_SIZE_DST = 1,
  SEND_NUM_SRC = 2,
  SEND_NUM_DST = 3,
  SEND_SIZE_ASIZE,
};

struct dist_dir_entry {
  Xt_idxlist idxlist;
  int rank;
};

struct dist_dir {
  unsigned num_entries;
  struct dist_dir_entry entries[];
};

static inline Xt_int
get_dist_dir_global_interval_size(Xt_idxlist idxlist, MPI_Comm comm) {

  Xt_int num_indices, global_num_indices;

  num_indices = xt_idxlist_get_num_indices(idxlist);

  xt_mpi_call(MPI_Allreduce(&num_indices, &global_num_indices, 1, Xt_int_dt,
                                MPI_SUM, comm), comm);

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  return ((global_num_indices + comm_size - 1) / comm_size) * comm_size;
}

static inline Xt_int get_min_idxlist_index(Xt_idxlist a, Xt_idxlist b) {

  Xt_int a_min_index, b_min_index;

  a_min_index = xt_idxlist_get_min_index(a);
  b_min_index = xt_idxlist_get_min_index(b);

  return MIN(a_min_index, b_min_index);
}

static inline Xt_int get_max_idxlist_index(Xt_idxlist a, Xt_idxlist b) {

  Xt_int a_max_index, b_max_index;

  a_max_index = xt_idxlist_get_max_index(a);
  b_max_index = xt_idxlist_get_max_index(b);

  return MAX(a_max_index, b_max_index);
}

static void generate_buckets(Xt_idxlist * buckets, Xt_idxlist src_idxlist,
                             Xt_idxlist dst_idxlist, MPI_Comm comm) {

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  Xt_int global_interval, local_interval;
  global_interval = get_dist_dir_global_interval_size(src_idxlist, comm);
  local_interval  = global_interval / comm_size;

  Xt_int local_index_range_lbound
    = get_min_idxlist_index(src_idxlist, dst_idxlist);
  Xt_int local_index_range_ubound
    = get_max_idxlist_index(src_idxlist, dst_idxlist);

  struct Xt_stripe * stripes = NULL;
  size_t stripes_array_size = 0;

  // generate buckets for each process
  for (int i = 0; i < comm_size; ++i) {

    Xt_int start = 0 + i * local_interval;
    Xt_int num_stripes = 0;

    while (start > local_index_range_lbound)
      start -= global_interval, num_stripes++;

    if (local_index_range_ubound > 0)
      num_stripes += ((local_index_range_ubound - start + global_interval)
                      / global_interval);
    else
      num_stripes++;

    ENSURE_ARRAY_SIZE(stripes, stripes_array_size, (size_t)num_stripes);

    for (Xt_int j = 0; j < num_stripes; ++j) {

      stripes[j].start = start + j * global_interval;
      stripes[j].stride = 1;
      stripes[j].nstrides = local_interval;
    }

    buckets[i] = xt_idxstripes_new(stripes, num_stripes);
  }

  free(stripes);
}

static void
compute_and_pack_bucket_intersections(Xt_idxlist * buckets,
                                      Xt_idxlist src_idxlist,
                                      Xt_idxlist dst_idxlist,
                                      int (*send_size)[SEND_SIZE_ASIZE],
                                      void ** send_buffer, MPI_Comm comm) {

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int send_buffer_size = 0;

  struct {
    Xt_idxlist src, dst;
  } *send_list = xmalloc((size_t)comm_size * sizeof(*send_list));

  for (int i = 0; i < comm_size; ++i) {

    send_list[i].src = xt_idxlist_get_intersection(src_idxlist, buckets[i]);
    send_list[i].dst = xt_idxlist_get_intersection(buckets[i], dst_idxlist);

    if (xt_idxlist_get_num_indices(send_list[i].src) > 0)
    {
      send_size[i][SEND_SIZE_SRC]
        = (int)xt_idxlist_get_pack_size(send_list[i].src, comm);
      send_size[i][SEND_NUM_SRC] = 1;
    }
    else
    {
      send_size[i][SEND_SIZE_SRC] = 0;
      send_size[i][SEND_NUM_SRC] = 0;
      xt_idxlist_delete(send_list[i].src);
      send_list[i].src = NULL;
    }

    if (xt_idxlist_get_num_indices(send_list[i].dst) > 0)
    {
      send_size[i][SEND_SIZE_DST]
        = (int)xt_idxlist_get_pack_size(send_list[i].dst, comm);
      send_size[i][SEND_NUM_DST] = 1;
    }
    else
    {
      send_size[i][SEND_SIZE_DST] = 0;
      send_size[i][SEND_NUM_DST] = 0;
      xt_idxlist_delete(send_list[i].dst);
      send_list[i].dst = NULL;
    }

    send_buffer_size += send_size[i][SEND_SIZE_SRC];
    send_buffer_size += send_size[i][SEND_SIZE_DST];
  }

  *send_buffer = xmalloc((size_t)send_buffer_size);

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
}

static inline void free_buckets(Xt_idxlist * buckets, int num_buckets) {

  for (int i = 0; i < num_buckets; ++i)
    xt_idxlist_delete(buckets[i]);

  free(buckets);
}

static void
get_remote_packed_intersection_size(int recv_size[SEND_SIZE_ASIZE],
                                    int (*send_size)[SEND_SIZE_ASIZE],
                                    MPI_Comm comm) {

#if MPI_VERSION > 2 || ( MPI_VERSION == 2 && MPI_SUBVERSION >= 2)
  xt_mpi_call(MPI_Reduce_scatter_block(send_size, recv_size, SEND_SIZE_ASIZE,
                                       MPI_INT, MPI_SUM, comm), comm);
#else
  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int * recv_count;
  recv_count = xmalloc((size_t)comm_size * sizeof(*recv_count));
  for (int i = 0; i < comm_size; ++i) recv_count[i] = SEND_SIZE_ASIZE;

  xt_mpi_call(MPI_Reduce_scatter(send_size, recv_size, recv_count, MPI_INT,
                                 MPI_SUM, comm), comm);

  free(recv_count);
#endif
}

static void send_intersections(void * send_buffer,
                               int (*send_size)[SEND_SIZE_ASIZE],
                               MPI_Request * dir_init_send_requests,
                               MPI_Comm comm) {

  int const src_send_tag = 0;
  int const dst_send_tag = 1;

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  for (int i = 0; i < 2 * comm_size; ++i)
    dir_init_send_requests[i] = MPI_REQUEST_NULL;

  int offset = 0;

  // pack the intersections into the send buffer
  for (int i = 0; i < comm_size; ++i)
    if (send_size[i][SEND_SIZE_SRC] > 0) {
      xt_mpi_call(MPI_Isend((char *)send_buffer + offset,
                            send_size[i][SEND_SIZE_SRC],
                            MPI_PACKED, i, src_send_tag, comm,
                            dir_init_send_requests+2*i+0),
                  comm);
      offset += send_size[i][SEND_SIZE_SRC];
    }

  for (int i = 0; i < comm_size; ++i)
    if (send_size[i][SEND_SIZE_DST] > 0) {
      xt_mpi_call(MPI_Isend((char *)send_buffer + offset,
                            send_size[i][SEND_SIZE_DST],
                            MPI_PACKED, i, dst_send_tag, comm,
                            dir_init_send_requests+2*i+1),
                  comm);
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
    dist_dir->entries[i].idxlist =
      xt_idxlist_unpack(recv_buffer, received_count, &position, comm);

    total_recv_size += received_count;
  }

  if (total_recv_size != recv_size)
    Xt_abort(comm, "ERROR: recv_intersections received wrong number of bytes",
             __FILE__, __LINE__);
  dist_dir->num_entries = (unsigned)recv_count;
}

static void
recv_and_unpack_intersections(int recv_size[SEND_SIZE_ASIZE],
                              struct dist_dir **src_dist_dir,
                              struct dist_dir **dst_dist_dir, MPI_Comm comm) {

  int const src_recv_tag = 0;
  int const dst_recv_tag = 1;

  *src_dist_dir = xmalloc(sizeof (struct dist_dir)
                          + (sizeof (struct dist_dir_entry)
                             * (size_t)recv_size[SEND_NUM_SRC]));
  *dst_dist_dir = xmalloc(sizeof (struct dist_dir)
                          + (sizeof (struct dist_dir_entry)
                             * (size_t)recv_size[SEND_NUM_DST]));


  void * recv_buffer = xmalloc((size_t)MAX(recv_size[SEND_SIZE_SRC],
                                           recv_size[SEND_SIZE_DST]));

  recv_and_unpack_intersection(*src_dist_dir, recv_size[SEND_SIZE_SRC],
                               recv_size[SEND_NUM_SRC], recv_buffer,
                               src_recv_tag, comm);
  recv_and_unpack_intersection(*dst_dist_dir, recv_size[SEND_SIZE_DST],
                               recv_size[SEND_NUM_DST], recv_buffer,
                               dst_recv_tag, comm);

  free(recv_buffer);
}

static void
match_src_dst_dist_dirs(struct dist_dir * src_dist_dir,
                        struct dist_dir * dst_dist_dir,
                        Xt_idxlist
                        (**src_dst_intersections)[dst_dist_dir->num_entries])
{
  Xt_idxlist (*src_dst_intersections_)[dst_dist_dir->num_entries]
    = (*src_dst_intersections)
    = xmalloc(src_dist_dir->num_entries * sizeof(**src_dst_intersections));

  for (unsigned i = 0; i < src_dist_dir->num_entries; ++i)
    for (unsigned j = 0; j < dst_dist_dir->num_entries; ++j)
    {
      src_dst_intersections_[i][j] = xt_idxlist_get_intersection(
        src_dist_dir->entries[i].idxlist, dst_dist_dir->entries[j].idxlist);
      if (!(xt_idxlist_get_num_indices(src_dst_intersections_[i][j]) > 0))
      {
        xt_idxlist_delete(src_dst_intersections_[i][j]);
        src_dst_intersections_[i][j] = NULL;
      }
    }
}

static void
match_and_pack_src_dst_dist_dirs(struct dist_dir * src_dist_dir,
                                 struct dist_dir * dst_dist_dir,
                                 int (*send_size)[SEND_SIZE_ASIZE],
                                 void ** send_buffer,
                                 MPI_Comm comm) {

  Xt_idxlist (*src_dst_intersections)[dst_dist_dir->num_entries];

  match_src_dst_dist_dirs(src_dist_dir, dst_dist_dir, &src_dst_intersections);

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int total_send_size = 0;

  int rank_pack_size;

  xt_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &rank_pack_size), comm);

  for (int i = 0; i < comm_size; ++i)
    send_size[i][SEND_SIZE_SRC] = 0, send_size[i][SEND_SIZE_DST] = 0,
      send_size[i][SEND_NUM_SRC] = 0, send_size[i][SEND_NUM_DST] = 0;

  for (unsigned i = 0; i < src_dist_dir->num_entries; ++i)
    for (unsigned j = 0; j < dst_dist_dir->num_entries; ++j)
      if (src_dst_intersections[i][j] != NULL) {

        int msg_size;

        msg_size = rank_pack_size
          + (int)xt_idxlist_get_pack_size(src_dst_intersections[i][j], comm);

        send_size[src_dist_dir->entries[i].rank][SEND_SIZE_SRC] += msg_size;
        ++(send_size[src_dist_dir->entries[i].rank][SEND_NUM_SRC]);
        send_size[dst_dist_dir->entries[j].rank][SEND_SIZE_DST] += msg_size;
        ++(send_size[dst_dist_dir->entries[j].rank][SEND_NUM_DST]);

        total_send_size += 2*msg_size;
      }

  (*send_buffer) = xmalloc((size_t)total_send_size);

  int position = 0;

  for (unsigned i = 0; i < (unsigned)comm_size; ++i) {

    if (send_size[i][SEND_SIZE_SRC] > 0) {

      // collect all src->dst data for the current rank

      for (unsigned j = 0; j < src_dist_dir->num_entries; ++j) {

        if ((unsigned)src_dist_dir->entries[j].rank != i) continue;

        for (unsigned k = 0; k < dst_dist_dir->num_entries; ++k) {

          if (src_dst_intersections[j][k] != NULL) {

            // pack rank
            xt_mpi_call(MPI_Pack(&(dst_dist_dir->entries[k].rank), 1, MPI_INT,
                                 (*send_buffer), total_send_size, &position,
                                 comm), comm);

            // pack intersection
            xt_idxlist_pack(src_dst_intersections[j][k], (*send_buffer),
                            total_send_size, &position, comm);
          }
        }
      }
    }
  }

  for (int i = 0; i < comm_size; ++i) {

    if (send_size[i][SEND_SIZE_DST] > 0) {

      // collect all dst->src data for the current rank

      for (unsigned j = 0; j < dst_dist_dir->num_entries; ++j) {

        if (dst_dist_dir->entries[j].rank != i) continue;

        for (unsigned k = 0; k < src_dist_dir->num_entries; ++k) {

          if (src_dst_intersections[k][j] != NULL) {

            // pack rank
            xt_mpi_call(MPI_Pack(&(src_dist_dir->entries[k].rank), 1, MPI_INT,
                                 (*send_buffer), total_send_size, &position,
                                 comm), comm);

            // pack intersection
            xt_idxlist_pack(src_dst_intersections[k][j], (*send_buffer),
                            total_send_size, &position, comm);
          }
        }
      }
    }
  }

  for (unsigned i = 0; i < src_dist_dir->num_entries; ++i)
    for (unsigned j = 0; j < dst_dist_dir->num_entries; ++j)
      if (src_dst_intersections[i][j] != NULL)
        xt_idxlist_delete(src_dst_intersections[i][j]);
  free(src_dst_intersections);
}

static void generate_distributed_directories(struct dist_dir **src_dist_dir,
                                             struct dist_dir **dst_dist_dir,
                                             Xt_idxlist src_idxlist,
                                             Xt_idxlist dst_idxlist,
                                             MPI_Comm comm) {

  int comm_size;
  Xt_idxlist * buckets;

  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  buckets = xmalloc((size_t)comm_size * sizeof(*buckets));
  generate_buckets(buckets, src_idxlist, dst_idxlist, comm);

  int (*send_size)[SEND_SIZE_ASIZE];
  void * send_buffer;

  send_size = xmalloc((size_t)comm_size * sizeof(*send_size));
  compute_and_pack_bucket_intersections(buckets, src_idxlist, dst_idxlist,
                                        send_size, &send_buffer, comm);


  free_buckets(buckets, comm_size);

  int recv_size[SEND_SIZE_ASIZE]; // for src and dst

  get_remote_packed_intersection_size(recv_size, send_size, comm);

  MPI_Request * dir_init_send_requests;

  dir_init_send_requests
    = xmalloc(2 * (size_t)comm_size * sizeof(*dir_init_send_requests));
  send_intersections(send_buffer, send_size, dir_init_send_requests, comm);

  free(send_size);

  recv_and_unpack_intersections(recv_size, src_dist_dir, dst_dist_dir, comm);

  // wait for the sends to be completed
  xt_mpi_call(MPI_Waitall(2 * comm_size, dir_init_send_requests,
                             MPI_STATUSES_IGNORE), comm);
  free(dir_init_send_requests);
  free(send_buffer);
}

static void
recv_and_unpack_dist_dir_result(struct dist_dir * dist_dir, int recv_size,
                                void * recv_buffer, int tag, MPI_Comm comm) {

  // initiate distributed directories
  dist_dir->num_entries = 0;

  while (recv_size > 0) {

    MPI_Status status;

    xt_mpi_call(MPI_Recv(recv_buffer, recv_size, MPI_PACKED,
                             MPI_ANY_SOURCE, tag, comm, &status), comm);

    int received_count;
    xt_mpi_call(MPI_Get_count(&status, MPI_PACKED, &received_count), comm);

    recv_size -= received_count;

    int position = 0;

    while (received_count > position) {

      int rank;

      xt_mpi_call(MPI_Unpack(recv_buffer, received_count, &position,
                                 &rank, 1, MPI_INT, comm), comm);

      dist_dir->entries[dist_dir->num_entries].rank = rank;
      dist_dir->entries[dist_dir->num_entries].idxlist =
        xt_idxlist_unpack(recv_buffer, received_count, &position, comm);

      dist_dir->num_entries++;
    }
  }

  if (0 != recv_size)
    Xt_abort(comm, "ERROR: recv_and_unpack_dist_dir_result"
             " received wrong number of bytes", __FILE__, __LINE__);
}

static void
generate_intersection_from_dist_dir_results(struct dist_dir *dist_dir_results,
                                            struct Xt_com_list **src_com,
                                            Xt_int *num_src_intersections,
                                            MPI_Comm comm) {

  int comm_size, comm_rank;

  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

  Xt_idxlist * intersections = NULL;
  size_t intersections_array_size = 0;

  struct Xt_com_list *p
    = xmalloc(dist_dir_results->num_entries * sizeof(*p));
  size_t num_src_ix = 0;

  for (int i = 0; i < comm_size; ++i) {

    unsigned num_intersections_per_rank = 0;
    for (unsigned j = 0; j < dist_dir_results->num_entries; ++j)
      num_intersections_per_rank += (dist_dir_results->entries[j].rank == i);

    if (num_intersections_per_rank > 0) {

      ENSURE_ARRAY_SIZE(intersections, intersections_array_size,
                        num_intersections_per_rank);

      unsigned offset = 0;
      int j = 0;

      while (num_intersections_per_rank != offset) {

        if (dist_dir_results->entries[j].rank == i) {
          intersections[offset] = dist_dir_results->entries[j].idxlist;
          offset++;
        }
        ++j;
      }

      p[num_src_ix].list =
        xt_idxlist_collection_new(intersections, num_intersections_per_rank);
      p[num_src_ix].rank = i;
      ++num_src_ix;
    }
  }

  p = xrealloc(p, num_src_ix * sizeof(*p));
  *src_com = p;
  *num_src_intersections = (Xt_int)num_src_ix;
  free(intersections);
}

static void free_dist_dir(struct dist_dir * dist_dir) {

  unsigned i;

  for (i = 0; i < dist_dir->num_entries; ++i)
    xt_idxlist_delete(dist_dir->entries[i].idxlist);
  free(dist_dir);
}

static void recv_and_unpack_dist_dir_results(int recv_size[SEND_SIZE_ASIZE],
                                             struct Xt_com_list **src_com,
                                             Xt_int *num_src_intersections,
                                             struct Xt_com_list **dst_com,
                                             Xt_int *num_dst_intersections,
                                             MPI_Comm comm) {

  int const src_recv_tag = 0;
  int const dst_recv_tag = 1;

  void * recv_buffer;

  recv_buffer = xmalloc((size_t)MAX(recv_size[SEND_SIZE_SRC],
                                    recv_size[SEND_SIZE_DST]));

  struct dist_dir *src_dist_dir_results, *dst_dist_dir_results;

  src_dist_dir_results = xmalloc(sizeof (struct dist_dir)
                                 + (sizeof (struct dist_dir_entry)
                                    * (size_t)recv_size[SEND_NUM_SRC]));
  dst_dist_dir_results = xmalloc(sizeof (struct dist_dir)
                                 + (sizeof (struct dist_dir_entry)
                                    * (size_t)recv_size[SEND_NUM_DST]));

  recv_and_unpack_dist_dir_result(src_dist_dir_results,
                                  recv_size[SEND_SIZE_SRC],
                                  recv_buffer, src_recv_tag, comm);
  assert(src_dist_dir_results->num_entries
         == (unsigned)recv_size[SEND_NUM_SRC]);

  recv_and_unpack_dist_dir_result(dst_dist_dir_results,
                                  recv_size[SEND_SIZE_DST],
                                  recv_buffer, dst_recv_tag, comm);
  assert(dst_dist_dir_results->num_entries
         == (unsigned)recv_size[SEND_NUM_DST]);

  free(recv_buffer);

  generate_intersection_from_dist_dir_results(src_dist_dir_results,
                                              src_com, num_src_intersections,
                                              comm);

  generate_intersection_from_dist_dir_results(dst_dist_dir_results,
                                              dst_com, num_dst_intersections,
                                              comm);

  free_dist_dir(src_dist_dir_results);
  free_dist_dir(dst_dist_dir_results);
}

static void exchange_idxlists(struct Xt_com_list **src_com,
                              Xt_int *num_src_intersections,
                              struct Xt_com_list **dst_com,
                              Xt_int *num_dst_intersections,
                              Xt_idxlist src_idxlist,
                              Xt_idxlist dst_idxlist,
                              MPI_Comm comm) {

  int comm_size;

  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  struct dist_dir *src_dist_dir, *dst_dist_dir;

  generate_distributed_directories(&src_dist_dir, &dst_dist_dir,
                                   src_idxlist, dst_idxlist, comm);

  void * send_buffer;
  int (*send_size)[SEND_SIZE_ASIZE];
  int recv_size[SEND_SIZE_ASIZE];

  send_size = xmalloc((size_t)comm_size * sizeof(*send_size));

  // match the source and destination entries in the local distributed
  // directories and pack the results into a sendable format
  match_and_pack_src_dst_dist_dirs(src_dist_dir, dst_dist_dir, send_size,
                                   &send_buffer, comm);
  free_dist_dir(src_dist_dir);
  free_dist_dir(dst_dist_dir);

  // get the data size the local process will receive from other processes
  get_remote_packed_intersection_size(recv_size, send_size, comm);

  MPI_Request * send_indices_requests;
  send_indices_requests
    = xmalloc(2 * (size_t)comm_size * sizeof(*send_indices_requests));

  send_intersections(send_buffer, send_size, send_indices_requests, comm);

  recv_and_unpack_dist_dir_results(recv_size, src_com, num_src_intersections,
                                   dst_com, num_dst_intersections, comm);

  xt_mpi_call(MPI_Waitall(2 * comm_size, send_indices_requests,
                              MPI_STATUSES_IGNORE), comm);

  free(send_buffer);
  free(send_size);
  free(send_indices_requests);
}

Xt_xmap xt_xmap_dist_dir_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                             MPI_Comm comm) {

  MPI_Comm exchange_comm;
  xt_mpi_call(MPI_Comm_dup(comm, &exchange_comm), comm);

  struct Xt_com_list * src_intersections, * dst_intersections;
  Xt_int num_src_intersections, num_dst_intersections;

  exchange_idxlists(&src_intersections, &num_src_intersections,
                    &dst_intersections, &num_dst_intersections,
                    src_idxlist, dst_idxlist, exchange_comm);

  MPI_Comm_free(&exchange_comm);

  Xt_xmap xmap
    = xt_xmap_intersection_new(num_src_intersections, src_intersections,
                               num_dst_intersections, dst_intersections,
                               src_idxlist, dst_idxlist, comm);

  for (Xt_int i = 0; i < num_src_intersections; ++i)
    if (src_intersections[i].list != NULL)
      xt_idxlist_delete(src_intersections[i].list);
  for (Xt_int i = 0; i < num_dst_intersections; ++i)
    if (dst_intersections[i].list != NULL)
      xt_idxlist_delete(dst_intersections[i].list);
  free(src_intersections), free(dst_intersections);

  return xmap;
}
