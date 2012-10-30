#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <mpi.h>

#include "xt/xt_idxlist.h"
#include "xt/xt_idxlist_collection.h"
#include "xt_idxvec.h"
#include "xt_idxstripes.h"
#include "xt_xmap.h"
#include "xt_mpi.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "ensure_array_size.h"
#include "xt_xmap_intersection.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


struct dist_dir_entry {
  Xt_idxlist idxlist;
  int rank;
};

struct dist_dir {

  struct dist_dir_entry * entries;
  unsigned num_entries;
};

static inline Xt_count get_dist_dir_global_interval_size(Xt_idxlist idxlist, MPI_Comm comm) {

  Xt_count num_indices, global_num_indices;

  num_indices = xt_idxlist_get_num_indices(idxlist);

  mpi_err_handler(MPI_Allreduce(&num_indices, &global_num_indices, 1, Xt_count_dt,
                                MPI_SUM, comm), comm);

  int comm_size;
  mpi_err_handler(MPI_Comm_size(comm, &comm_size), comm);

  return ((global_num_indices + comm_size - 1) / comm_size) * comm_size;
}

static inline Xt_idx get_min_idxlist_index(Xt_idxlist a, Xt_idxlist b) {

  Xt_idx a_min_index, b_min_index;

  a_min_index = xt_idxlist_get_min_index(a);
  b_min_index = xt_idxlist_get_min_index(b);

  return MIN(a_min_index, b_min_index);
}

static inline Xt_idx get_max_idxlist_index(Xt_idxlist a, Xt_idxlist b) {

  Xt_idx a_max_index, b_max_index;

  a_max_index = xt_idxlist_get_max_index(a);
  b_max_index = xt_idxlist_get_max_index(b);

  return MAX(a_max_index, b_max_index);
}

static void generate_buckets(Xt_idxlist * buckets, Xt_idxlist src_idxlist,
                             Xt_idxlist dst_idxlist, MPI_Comm comm) {

  int comm_size;
  mpi_err_handler(MPI_Comm_size(comm, &comm_size), comm);

  Xt_count global_interval, local_interval;
  global_interval = get_dist_dir_global_interval_size(src_idxlist, comm);
  local_interval  = global_interval / comm_size;

  Xt_idx local_index_range[2];
  local_index_range[0] = get_min_idxlist_index(src_idxlist, dst_idxlist);
  local_index_range[1] = get_max_idxlist_index(src_idxlist, dst_idxlist);

  struct Xt_stripe * stripes = NULL;
  Xt_count stripes_array_size = 0;

  // generate buckets for each process
  for (int i = 0; i < comm_size; ++i) {

    Xt_idx start = 0 + i * local_interval;
    Xt_count num_stripes = 0;

    while (start > local_index_range[0])
      start -= global_interval, num_stripes++;

    if (local_index_range[1] > 0)
      num_stripes += (local_index_range[1] - start + global_interval) / global_interval;
    else
      num_stripes++;

    ENSURE_ARRAY_SIZE(stripes, stripes_array_size, num_stripes);

    for (Xt_count j = 0; j < num_stripes; ++j) {

      stripes[j].start = start + j * global_interval;
      stripes[j].stride = 1;
      stripes[j].nstrides = local_interval;
    }

    buckets[i] = xt_idxstripes_new(stripes, num_stripes);
  }

  free(stripes);
}

static void
compute_and_pack_bucket_intersections(Xt_idxlist * buckets, Xt_idxlist src_idxlist,
                                      Xt_idxlist dst_idxlist, int * send_size,
                                      void ** send_buffer, MPI_Comm comm) {

  int comm_size;
  mpi_err_handler(MPI_Comm_size(comm, &comm_size), comm);

  int send_buffer_size = 0;

  Xt_idxlist * src_send_list, * dst_send_list;

  src_send_list = xmalloc(comm_size * sizeof(*src_send_list));
  dst_send_list = xmalloc(comm_size * sizeof(*dst_send_list));

  for (int i = 0; i < comm_size; ++i) {

    src_send_list[i] = xt_idxlist_get_intersection(src_idxlist, buckets[i]);
    dst_send_list[i] = xt_idxlist_get_intersection(buckets[i], dst_idxlist);

    if (xt_idxlist_get_num_indices(src_send_list[i]) > 0)
      send_size[2*i+0] = xt_idxlist_get_pack_size(src_send_list[i], comm);
    else
      send_size[2*i+0] = 0;

    if (xt_idxlist_get_num_indices(dst_send_list[i]) > 0)
      send_size[2*i+1] = xt_idxlist_get_pack_size(dst_send_list[i], comm);
    else
      send_size[2*i+1] = 0;

    send_buffer_size += send_size[2*i+0];
    send_buffer_size += send_size[2*i+1];
  }

  *send_buffer = xmalloc(send_buffer_size);

  int position = 0;

  for (int i = 0; i < comm_size; ++i) {

    if (xt_idxlist_get_num_indices(src_send_list[i]) > 0)
      xt_idxlist_pack(src_send_list[i], *send_buffer, send_buffer_size,
                      &position, comm);

    xt_idxlist_delete(src_send_list[i]);
  }

  for (int i = 0; i < comm_size; ++i) {

    if (xt_idxlist_get_num_indices(dst_send_list[i]) > 0)
      xt_idxlist_pack(dst_send_list[i], *send_buffer, send_buffer_size,
                      &position, comm);

    xt_idxlist_delete(dst_send_list[i]);
  }

  free(src_send_list);
  free(dst_send_list);
}

static inline void free_buckets(Xt_idxlist * buckets, int num_buckets) {

  for (int i = 0; i < num_buckets; ++i)
    xt_idxlist_delete(buckets[i]);

  free(buckets);
}

static void
get_remote_packed_intersection_size(int * recv_size, int * send_size, MPI_Comm comm) {

  int comm_size;
  mpi_err_handler(MPI_Comm_size(comm, &comm_size), comm);

  int * recv_count;
  recv_count = xmalloc(comm_size * sizeof(*recv_count));
  for (int i = 0; i < comm_size; ++i) recv_count[i] = 2;

  mpi_err_handler(MPI_Reduce_scatter(send_size, recv_size, recv_count, MPI_INT,
                                     MPI_SUM, comm), comm);

  free(recv_count);
}

static void send_intersections(void * send_buffer, int * send_size,
                               MPI_Request * dir_init_send_requests,
                               MPI_Comm comm) {

  int const src_send_tag = 0;
  int const dst_send_tag = 1;

  int comm_size;
  mpi_err_handler(MPI_Comm_size(comm, &comm_size), comm);

  for (int i = 0; i < 2 * comm_size; ++i)
    dir_init_send_requests[i] = MPI_REQUEST_NULL;

  int offset = 0;

  // pack the intersections into the send buffer
  for (int i = 0; i < comm_size; ++i)
    if (send_size[i*2+0] > 0) {
      mpi_err_handler(MPI_Isend((char *)send_buffer + offset, send_size[i*2+0], MPI_PACKED,
                                i, src_send_tag, comm, dir_init_send_requests+2*i+0),
                      comm);
      offset += send_size[i*2+0];
    }

  for (int i = 0; i < comm_size; ++i)
    if (send_size[i*2+1] > 0) {
      mpi_err_handler(MPI_Isend((char *)send_buffer + offset, send_size[i*2+1], MPI_PACKED,
                                i, dst_send_tag, comm, dir_init_send_requests+2*i+1),
                      comm);
      offset += send_size[i*2+1];
    }
}

static void
recv_and_unpack_intersection(struct dist_dir * dist_dir, int recv_size,
                             void * recv_buffer, int tag, MPI_Comm comm) {

  // initiate distributed directories
  dist_dir->num_entries = 0;
  dist_dir->entries = NULL;

  int dist_dir_entries_array_size = 0;

  while (recv_size > 0) {

    MPI_Status status;

    mpi_err_handler(MPI_Recv(recv_buffer, recv_size, MPI_PACKED, MPI_ANY_SOURCE,
                             tag, comm, &status), comm);

    int received_count;
    mpi_err_handler(MPI_Get_count(&status, MPI_PACKED, &received_count), comm);

    ENSURE_ARRAY_SIZE(dist_dir->entries, dist_dir_entries_array_size,
                      dist_dir->num_entries + 1);

    int position = 0;

    dist_dir->entries[dist_dir->num_entries].rank = status.MPI_SOURCE;
    dist_dir->entries[dist_dir->num_entries].idxlist =
      xt_idxlist_unpack(recv_buffer, received_count, &position, comm);

    recv_size -= received_count;

    dist_dir->num_entries++;
  }

  if (0 != recv_size)
    Xt_abort(comm, "ERROR: recv_intersections received wrong number of bytes",
             __FILE__, __LINE__);
}

static void
recv_and_unpack_intersections(int * recv_size, struct dist_dir * src_dist_dir,
                              struct dist_dir * dst_dist_dir, MPI_Comm comm) {

  int const src_recv_tag = 0;
  int const dst_recv_tag = 1;

  void * recv_buffer;

  recv_buffer = xmalloc(MAX(recv_size[0],recv_size[1]));

  recv_and_unpack_intersection(src_dist_dir, recv_size[0], recv_buffer,
                               src_recv_tag, comm);
  recv_and_unpack_intersection(dst_dist_dir, recv_size[1], recv_buffer,
                               dst_recv_tag, comm);

  free(recv_buffer);
}

static void
match_src_dst_dist_dirs(struct dist_dir * src_dist_dir, struct dist_dir * dst_dist_dir,
                        Xt_idxlist *** src_dst_intersections) {

  Xt_idxlist * src_dst_intersections_;

  src_dst_intersections_ = xmalloc(src_dist_dir->num_entries * dst_dist_dir->num_entries *
                                   sizeof(*src_dst_intersections_));
  (*src_dst_intersections) = xmalloc(src_dist_dir->num_entries *
                                     sizeof(**src_dst_intersections));
  for (unsigned i = 0; i < src_dist_dir->num_entries; ++i)
    (*src_dst_intersections)[i] = src_dst_intersections_ + i * dst_dist_dir->num_entries;

  for (unsigned i = 0; i < src_dist_dir->num_entries; ++i)
    for (unsigned j = 0; j < dst_dist_dir->num_entries; ++j)
      (*src_dst_intersections)[i][j] = xt_idxlist_get_intersection(
        src_dist_dir->entries[i].idxlist, dst_dist_dir->entries[j].idxlist);
}

static void
match_and_pack_src_dst_dist_dirs(struct dist_dir * src_dist_dir,
                                 struct dist_dir * dst_dist_dir,
                                 int * send_size, void ** send_buffer,
                                 MPI_Comm comm) {

  Xt_idxlist ** src_dst_intersections;

  match_src_dst_dist_dirs(src_dist_dir, dst_dist_dir, &src_dst_intersections);

  int comm_size;
  mpi_err_handler(MPI_Comm_size(comm, &comm_size), comm);

  int total_send_size = 0;

  int rank_pack_size;

  mpi_err_handler(MPI_Pack_size(1, MPI_INT, comm, &rank_pack_size), comm);

  for (int i = 0; i < 2 * comm_size; ++i)
    send_size[i] = 0;

  for (int i = 0; i < src_dist_dir->num_entries; ++i)
    for (int j = 0; j < dst_dist_dir->num_entries; ++j)
      if (xt_idxlist_get_num_indices(src_dst_intersections[i][j]) > 0) {

        int msg_size;

        msg_size = rank_pack_size + xt_idxlist_get_pack_size(src_dst_intersections[i][j],
                                                             comm);
      
        send_size[2*src_dist_dir->entries[i].rank+0] += msg_size;
        send_size[2*dst_dist_dir->entries[j].rank+1] += msg_size;

        total_send_size += 2*msg_size;
      }

  (*send_buffer) = xmalloc(total_send_size);

  int position = 0;

  for (int i = 0; i < comm_size; ++i) {

    if (send_size[2*i+0] > 0) {

      // collect all src->dst data for the current rank

      for (int j = 0; j < src_dist_dir->num_entries; ++j) {

        if (src_dist_dir->entries[j].rank != i) continue;

        for (int k = 0; k < dst_dist_dir->num_entries; ++k) {

          if (xt_idxlist_get_num_indices(src_dst_intersections[j][k]) > 0) {

            // pack rank
            mpi_err_handler(MPI_Pack(&(dst_dist_dir->entries[k].rank), 1, MPI_INT,
                                     (*send_buffer), total_send_size, &position, comm),
                            comm);

            // pack intersection
            xt_idxlist_pack(src_dst_intersections[j][k], (*send_buffer),
                            total_send_size, &position, comm);
          }
        }
      }
    }
  }

  for (int i = 0; i < comm_size; ++i) {

    if (send_size[2*i+1] > 0) {

      // collect all dst->src data for the current rank

      for (int j = 0; j < dst_dist_dir->num_entries; ++j) {

        if (dst_dist_dir->entries[j].rank != i) continue;

        for (int k = 0; k < src_dist_dir->num_entries; ++k) {

          if (xt_idxlist_get_num_indices(src_dst_intersections[k][j]) > 0) {

            // pack rank
            mpi_err_handler(MPI_Pack(&(src_dist_dir->entries[k].rank), 1, MPI_INT,
                                     (*send_buffer), total_send_size, &position, comm),
                            comm);

            // pack intersection
            xt_idxlist_pack(src_dst_intersections[k][j], (*send_buffer),
                            total_send_size, &position, comm);
          }
        }
      }
    }
  }

  for (int i = 0; i < src_dist_dir->num_entries; ++i)
    for (int j = 0; j < dst_dist_dir->num_entries; ++j)
      if (src_dst_intersections[i][j] != NULL)
        xt_idxlist_delete(src_dst_intersections[i][j]);
  if (*src_dst_intersections != NULL) free(*src_dst_intersections);
  free(src_dst_intersections);
}

static void generate_distributed_directories(struct dist_dir * src_dist_dir,
                                             struct dist_dir * dst_dist_dir,
                                             Xt_idxlist src_idxlist,
                                             Xt_idxlist dst_idxlist,
                                             MPI_Comm comm) {

  int comm_size;
  Xt_idxlist * buckets;

  mpi_err_handler(MPI_Comm_size(comm, &comm_size), comm);
  buckets = xmalloc(comm_size * sizeof(*buckets));
  generate_buckets(buckets, src_idxlist, dst_idxlist, comm);

  int * send_size;
  void * send_buffer;

  send_size = xmalloc(2 * comm_size * sizeof(*send_size));
  compute_and_pack_bucket_intersections(buckets, src_idxlist, dst_idxlist,
                                        send_size, &send_buffer, comm);


  free_buckets(buckets, comm_size);

  int recv_size[2]; // for src and dst

  get_remote_packed_intersection_size(recv_size, send_size, comm);

  MPI_Request * dir_init_send_requests;

  dir_init_send_requests = xmalloc(2 * comm_size * sizeof(*dir_init_send_requests));
  send_intersections(send_buffer, send_size, dir_init_send_requests, comm);

  free(send_size);

  recv_and_unpack_intersections(recv_size, src_dist_dir, dst_dist_dir, comm);

  // wait for the sends to be completed
  mpi_err_handler(MPI_Waitall(2 * comm_size, dir_init_send_requests,
                             MPI_STATUSES_IGNORE), comm);
  free(dir_init_send_requests);
  free(send_buffer);
}

static void
recv_and_unpack_dist_dir_result(struct dist_dir * dist_dir, int recv_size,
                                void * recv_buffer, int tag, MPI_Comm comm) {

  // initiate distributed directories
  dist_dir->num_entries = 0;
  dist_dir->entries = NULL;

  int dist_dir_entries_array_size = 0;

  while (recv_size > 0) {

    MPI_Status status;

    mpi_err_handler(MPI_Recv(recv_buffer, recv_size, MPI_PACKED,
                             MPI_ANY_SOURCE, tag, comm, &status), comm);

    int received_count;
    mpi_err_handler(MPI_Get_count(&status, MPI_PACKED, &received_count), comm);

    recv_size -= received_count;

    int position = 0;

    while (received_count > position) {

      ENSURE_ARRAY_SIZE(dist_dir->entries, dist_dir_entries_array_size,
                        dist_dir->num_entries + 1);

      int rank;

      mpi_err_handler(MPI_Unpack(recv_buffer, received_count, &position,
                                 &rank, 1, MPI_INT, comm), comm);

      dist_dir->entries[dist_dir->num_entries].rank = rank;
      dist_dir->entries[dist_dir->num_entries].idxlist =
        xt_idxlist_unpack(recv_buffer, received_count, &position, comm);

      dist_dir->num_entries++;
    }
  }

  if (0 != recv_size)
    Xt_abort(comm, "ERROR: recv_and_unpack_dist_dir_result received wrong number of bytes",
             __FILE__, __LINE__);
}

static void
generate_intersection_from_dist_dir_results(struct dist_dir dist_dir_results,
                                            Xt_idxlist ** src_intersections,
                                            int ** src_ranks,
                                            Xt_count * num_src_intersections,
                                            MPI_Comm comm) {

  int comm_size, comm_rank;

  mpi_err_handler(MPI_Comm_size(comm, &comm_size), comm);
  mpi_err_handler(MPI_Comm_rank(comm, &comm_rank), comm);

  Xt_idxlist * intersections = NULL;
  unsigned intersections_array_size = 0;

  (*src_intersections) = xmalloc(dist_dir_results.num_entries * sizeof(**src_intersections));
  (*src_ranks) = xmalloc(dist_dir_results.num_entries * sizeof(**src_ranks));
  (*num_src_intersections) = 0;

  for (int i = 0; i < comm_size; ++i) {

    unsigned num_intersections_per_rank = 0;
    for (unsigned j = 0; j < dist_dir_results.num_entries; ++j)
      if (dist_dir_results.entries[j].rank == i) num_intersections_per_rank++;

    if (num_intersections_per_rank > 0) {

      ENSURE_ARRAY_SIZE(intersections, intersections_array_size,
                        num_intersections_per_rank);

      unsigned offset = 0;
      int j = 0;

      while (num_intersections_per_rank != offset) {

        if (dist_dir_results.entries[j].rank == i) {
          intersections[offset] = dist_dir_results.entries[j].idxlist;
          offset++;
        }
        ++j;
      }

      (*src_intersections)[(*num_src_intersections)] =
        xt_idxlist_collection_new(intersections, num_intersections_per_rank);
      (*src_ranks)[(*num_src_intersections)] = i;
      (*num_src_intersections)++;
    }
  }

  (*src_intersections) = xrealloc((*src_intersections), (*num_src_intersections) * sizeof(**src_intersections));
  (*src_ranks) = xrealloc((*src_ranks), (*num_src_intersections) * sizeof(**src_ranks));

  free(intersections);
}

void free_dist_dir(struct dist_dir * dist_dir) {

  unsigned i;

  for (i = 0; i < dist_dir->num_entries; ++i)
    xt_idxlist_delete(dist_dir->entries[i].idxlist);
  free(dist_dir->entries);

  dist_dir->entries = NULL;
  dist_dir->num_entries = 0;
}

static void recv_and_unpack_dist_dir_results(int * recv_size,
                                             Xt_idxlist ** src_intersections,
                                             int ** src_ranks,
                                             Xt_count * num_src_intersections,
                                             Xt_idxlist ** dst_intersections,
                                             int ** dst_ranks,
                                             Xt_count * num_dst_intersections,
                                             MPI_Comm comm) {

  int const src_recv_tag = 0;
  int const dst_recv_tag = 1;

  void * recv_buffer;

  recv_buffer = xmalloc(MAX(recv_size[0],recv_size[1]));

  struct dist_dir src_dist_dir_results, dst_dist_dir_results;

  recv_and_unpack_dist_dir_result(&src_dist_dir_results, recv_size[0], recv_buffer,
                                  src_recv_tag, comm);
  recv_and_unpack_dist_dir_result(&dst_dist_dir_results, recv_size[1], recv_buffer,
                                  dst_recv_tag, comm);

  free(recv_buffer);

  generate_intersection_from_dist_dir_results(src_dist_dir_results, src_intersections,
                                              src_ranks, num_src_intersections, comm);

  generate_intersection_from_dist_dir_results(dst_dist_dir_results, dst_intersections,
                                              dst_ranks, num_dst_intersections, comm);

  free_dist_dir(&src_dist_dir_results);
  free_dist_dir(&dst_dist_dir_results);
}

static void exchange_idxlists(Xt_idxlist ** src_intersections,
                              int ** src_ranks,
                              Xt_count * num_src_intersections,
                              Xt_idxlist ** dst_intersections,
                              int ** dst_ranks,
                              Xt_count * num_dst_intersections,
                              Xt_idxlist src_idxlist,
                              Xt_idxlist dst_idxlist,
                              MPI_Comm comm) {

  int comm_size;

  mpi_err_handler(MPI_Comm_size(comm, &comm_size), comm);

  struct dist_dir src_dist_dir, dst_dist_dir;

  generate_distributed_directories(&src_dist_dir, &dst_dist_dir,
                                   src_idxlist, dst_idxlist, comm);

  void * send_buffer;
  int * send_size;
  int recv_size[2];

  send_size = xmalloc(2 * comm_size * sizeof(*send_size));

  // match the source and destination entries in the local distributed directories
  // and pack the results into a sendable format
  match_and_pack_src_dst_dist_dirs(&src_dist_dir, &dst_dist_dir, send_size, &send_buffer,
                                   comm);

  // get the data size the local process will receive from other processes
  get_remote_packed_intersection_size(recv_size, send_size, comm);

  MPI_Request * send_indices_requests;
  send_indices_requests = xmalloc(2 * comm_size * sizeof(send_indices_requests));

  send_intersections(send_buffer, send_size, send_indices_requests, comm);

  recv_and_unpack_dist_dir_results(recv_size, src_intersections, src_ranks,
                                   num_src_intersections, dst_intersections,
                                   dst_ranks, num_dst_intersections, comm);

  mpi_err_handler(MPI_Waitall(2 * comm_size, send_indices_requests,
                              MPI_STATUSES_IGNORE), comm);

  free(send_buffer);
  free_dist_dir(&src_dist_dir);
  free_dist_dir(&dst_dist_dir);
  free(send_size);
  free(send_indices_requests);
}

Xt_xmap xt_xmap_dist_dir_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist, MPI_Comm comm) {

  MPI_Comm exchange_comm;
  mpi_err_handler(MPI_Comm_dup(comm, &exchange_comm), comm);

  Xt_idxlist * src_intersections, * dst_intersections;
  int * src_ranks, * dst_ranks;
  Xt_count num_src_intersections, num_dst_intersections;

  exchange_idxlists(&src_intersections, &src_ranks, &num_src_intersections,
                    &dst_intersections, &dst_ranks, &num_dst_intersections,
                    src_idxlist, dst_idxlist, exchange_comm);

  MPI_Comm_free(&exchange_comm);

  Xt_xmap xmap;

  xmap = xt_xmap_intersection_new(src_intersections, src_ranks, num_src_intersections,
                                  dst_intersections, dst_ranks, num_dst_intersections,
                                  src_idxlist, dst_idxlist, comm);

  for (Xt_count i = 0; i < num_src_intersections; ++i)
    if (src_intersections[i] != NULL)
      xt_idxlist_delete(src_intersections[i]);
  for (Xt_count i = 0; i < num_dst_intersections; ++i)
    if (dst_intersections[i] != NULL)
      xt_idxlist_delete(dst_intersections[i]);
  free(src_intersections), free(dst_intersections);
  free(src_ranks), free(dst_ranks);

  return xmap;
}
