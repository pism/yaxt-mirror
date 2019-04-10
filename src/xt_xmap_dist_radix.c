/**
 * @file xt_xmap_dist_radix.c
 *
 * @copyright Copyright  (C)  2019 Jörg Behrens <behrens@dkrz.de>
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <mpi.h>

#include "xt/xt_xmap.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_xmap_dist_radix.h"
#include "xt/xt_xmap_intersection.h"
#include "xt_memcmp.h"
#include "xt_id_buckets.h"
#include "instr.h"

#define MAX(a,b) ((a) >= (b) ? (a) : (b))

struct exchange_data {
  int * pos;
  int * counts;
  int * ranks;
  int num_msgs;
  int buffer[]; // buffer for pos, counts, ranks (this ensures that no alignment
                // issue occur when allocating everything in a single call
};

static int compare_3int(const void * a, const void * b) {
  int const * a_int = (int const *)a;
  int const * b_int = (int const *)b;

  int ret = a_int[0] - b_int[0];
  if (ret) return ret;
  ret = a_int[1] - b_int[1];
  if (ret) return ret;
  return a_int[2] - b_int[2];
}

// Allocates an instance of struct exchange data and sets num_msgs, counts, and
// ranks. The array pos needs to be set elsewhere.
static struct exchange_data * generate_basic_exchange_data (
  int num_ids, int comm_size, int *restrict num_id_per_rank) {

  int num_msgs = 0;
  for (int i = 0; i < comm_size; ++i) if (num_id_per_rank[i]) ++num_msgs;

  struct exchange_data * exchange_data =
    xmalloc(sizeof(*exchange_data) +
            (2 * (size_t)num_msgs + (size_t)num_ids) * sizeof(int));

  int *restrict pos = &(exchange_data->buffer[0]);
  int *restrict counts = &(exchange_data->buffer[num_ids]);
  int *restrict ranks =
    &(exchange_data->buffer[(size_t)num_ids + (size_t)num_msgs]);
  exchange_data->pos = pos;
  exchange_data->counts = counts;
  exchange_data->ranks = ranks;
  exchange_data->num_msgs = num_msgs;

  for (int i = 0, j = 0; i < comm_size; ++i) {
    if (num_id_per_rank[i]) {
      counts[j] = num_id_per_rank[i];
      ranks[j] = i;
      ++j;
    }
  }

  return exchange_data;
}

static struct exchange_data * generate_basic_send_exchange_data (
  int * basic_exchange_data, int basic_exchange_data_size, int comm_size,
  int *restrict num_id_per_rank, int * remote_pos) {

  // basic exchange data contains blocks of three integers
  // each of these block is associated with an id and contains
  // the following data: [remote rank, local position, remote position]

  // sort the data first by rank, then by local location, and
  // followed by remote position
  int num_ids = basic_exchange_data_size / 3;
  qsort(basic_exchange_data, num_ids, 3 * sizeof(int), compare_3int);

  for (int i = 0; i < comm_size; ++i) num_id_per_rank[i] = 0;

  for (int i = 0; i < num_ids; ++i)
    num_id_per_rank[basic_exchange_data[3*i+0]]++;

  struct exchange_data * exchange_data =
    generate_basic_exchange_data(num_ids, comm_size, num_id_per_rank);

  int *restrict local_pos = &(exchange_data->buffer[0]);
  for (int i = 0; i < num_ids; ++i) {
    local_pos[i] = basic_exchange_data[3*i+1];
    remote_pos[i] = basic_exchange_data[3*i+2];
  }

  return exchange_data;
}

// computes recv_count, send_displ, recv_displ, and total_recv_count based on
// the contents of send_count
// these arrays are used for a call to MPI_Alltoallv
static void compute_alltoallv_parameters(
  int *restrict send_count, int *restrict send_displ,
  int *restrict recv_count, int *restrict recv_displ,
  int * total_recv_count, MPI_Comm comm) {

  int comm_size, is_inter;
  xt_mpi_call(MPI_Comm_test_inter(comm, &is_inter), comm);
  int (*get_comm_size)(MPI_Comm, int *)
    = is_inter ? MPI_Comm_remote_size : MPI_Comm_size;
  xt_mpi_call(get_comm_size(comm, &comm_size), comm);

  xt_mpi_call(
    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm), comm);

  int total_recv_count_ = 0;

  send_displ[0] = 0; // we hack send_displ slightly to make the packing of the
                     // sendbuffer easier
  for (int i = 0, send_accu = 0, recv_accu = 0; i < comm_size; ++i) {
    send_displ[i+1] = send_accu;
    recv_displ[i] = recv_accu;
    send_accu += send_count[i];
    recv_accu += recv_count[i];
    total_recv_count_ += recv_count[i];
  }
  *total_recv_count = total_recv_count_;
}

static void generate_ubounds(void const * src_ids, size_t num_src_ids,
                             void const * dst_ids, size_t num_dst_ids,
                             size_t id_size, MPI_Comm comm, size_t num_buckets,
                             void * ubounds) {

  // merge source and destination ids into single array
  size_t num_all_ids = num_src_ids + num_dst_ids;
  size_t src_ids_size = num_src_ids * id_size;
  size_t dst_ids_size = num_dst_ids * id_size;
  void * all_ids = xmalloc(src_ids_size + dst_ids_size);
  memcpy(all_ids, src_ids, src_ids_size);
  memcpy((char*)all_ids + src_ids_size, dst_ids, dst_ids_size);

  // sort all ids
  qsort(all_ids, num_all_ids, id_size, xt_memcmp_get_func(id_size));

  // generate upper bounds for the buckets
  Xt_id_buckets buckets =
    xt_id_buckets_new(all_ids, id_size, num_all_ids, num_buckets, comm);
  xt_id_buckets_get_ubounts(buckets, ubounds);
  xt_id_buckets_delete(buckets);

  free(all_ids);
}

static void redistribute_ids(
  void const * ids, int num_ids, size_t id_size,
  void ** redist_ids, int ** ranks, int ** orig_pos,
  int * num_redist_ids, void * ubounds, MPI_Datatype ID_TYPE, MPI_Comm comm) {

  int comm_size, is_inter;
  xt_mpi_call(MPI_Comm_test_inter(comm, &is_inter), comm);
  int (*get_comm_size)(MPI_Comm, int *)
    = is_inter ? MPI_Comm_remote_size : MPI_Comm_size;
  xt_mpi_call(get_comm_size(comm, &comm_size), comm);

  Xt_memcmp_func compare_id = xt_memcmp_get_func(id_size);

  // determine buckets for each id
  int * int_buffer = xmalloc(((size_t)num_ids + 4 * (size_t)comm_size + 1) *
                             sizeof(*int_buffer));
  int * send_count = int_buffer + 0 * (size_t)comm_size;
  int * recv_count = int_buffer + 1 * (size_t)comm_size;
  int * send_displ = int_buffer + 2 * (size_t)comm_size;
  int * recv_displ = int_buffer + 3 * (size_t)comm_size + 1;
  int * bucket_ids = int_buffer + 4 * (size_t)comm_size + 1;
  for (int i = 0; i < comm_size; ++i) send_count[i] = 0;
  for (int i = 0, j; i < num_ids; ++i) {
    void const * curr_id = (char*)ids + (size_t)i * id_size;
    for (j = 0; j < comm_size; ++j)
      if (compare_id(curr_id, (char*)ubounds + (size_t)j * id_size) <= 0) break;
    bucket_ids[i] = j;
    send_count[j]++;
  }

  int total_recv_count;
  compute_alltoallv_parameters(
    send_count, send_displ, recv_count, recv_displ, &total_recv_count, comm);

  // id, rank, orig_pos
  size_t id_meta_size = id_size + 2 * sizeof(int);
  size_t send_buffer_size = (size_t)num_ids * id_size;
  size_t recv_buffer_size = (size_t)total_recv_count * id_size;
  size_t sorted_id_meta_size = (size_t)total_recv_count * id_meta_size;

  // pack the send buffer, first source ids then destination ids
  void * id_buffer =
    xmalloc(send_buffer_size + recv_buffer_size + sorted_id_meta_size);
  void * id_send_buffer = id_buffer;
  void * id_recv_buffer = (char*)id_buffer + send_buffer_size;
  void * sorted_id_meta = (char*)id_buffer + send_buffer_size + recv_buffer_size;
  int * orig_pos_send_buffer =
    xmalloc(((size_t)num_ids + (size_t)total_recv_count) * sizeof(int));
  int * orig_pos_recv_buffer = orig_pos_send_buffer + (size_t)num_ids;
  for (int i = 0; i < num_ids; ++i) {
    // get current displacement and adjust send_displ
    size_t curr_displ = (size_t)(send_displ[bucket_ids[i]+1]++);
    memcpy((char*)id_send_buffer + curr_displ * id_size,
           (char*)ids + (size_t)i * id_size, id_size);
    orig_pos_send_buffer[curr_displ] = i;
  }

  // redistribute ids and orig_pos
  xt_mpi_call(MPI_Alltoallv(id_send_buffer, send_count, send_displ, ID_TYPE,
                            id_recv_buffer, recv_count, recv_displ, ID_TYPE,
                            comm), comm);
  xt_mpi_call(MPI_Alltoallv(orig_pos_send_buffer, send_count, send_displ, MPI_INT,
                            orig_pos_recv_buffer, recv_count, recv_displ, MPI_INT,
                            comm), comm);


  // add meta data (rank and orig_pos) to each id and sort them all
  // (first by id, then by rank, and by orig_pos)
  size_t id_meta_offset = 0, id_offset = 0;
  for (int rank = 0, j = 0; rank < comm_size; ++rank) {

    int curr_num_ids = recv_count[rank];
    for (int i = 0; i < curr_num_ids;
         ++i, ++j) {
      void * curr_id_meta = (char*)sorted_id_meta + (size_t)j * id_meta_size;
      memcpy(curr_id_meta, (char*)id_recv_buffer + (size_t)j * id_size, id_size);
      memcpy((char*)curr_id_meta + id_size, &rank, sizeof(rank));
      memcpy((char*)curr_id_meta + id_size + sizeof(rank),
             orig_pos_recv_buffer + j, sizeof(int));
    }
  }

  // sort id_meta data
  qsort(sorted_id_meta, total_recv_count, id_meta_size,
        xt_memcmp_get_func(id_meta_size));

  void *restrict redist_ids_ = xmalloc(total_recv_count * id_size);
  int *restrict ranks_ = xmalloc(total_recv_count * sizeof(ranks_));
  int *restrict orig_pos_ = xmalloc(total_recv_count * sizeof(*orig_pos_));

  *redist_ids = redist_ids_;
  *ranks = ranks_;
  *orig_pos = orig_pos_;
  *num_redist_ids = total_recv_count;

  for (int i = 0; i < total_recv_count; ++i) {
    void *restrict curr_id_meta =
      (char*)sorted_id_meta + (size_t)i * id_meta_size;
    void *restrict curr_id = (char*)redist_ids_ + (size_t)i * id_size;
    memcpy(curr_id, curr_id_meta, id_size);
    memcpy(ranks_ + i, (char*)curr_id_meta + id_size, sizeof(int));
    memcpy(
      orig_pos_ + i, (char*)curr_id_meta + id_size + sizeof(int), sizeof(int));
  }

  free(int_buffer);
  free(orig_pos_send_buffer);
  free(id_buffer);
}

// generates group ids for an intercomm
// for all processes within one group it will return 0, while for all processes
// of the other group 1 is returned
// this makes coordination between the groups easier
static int get_intercomm_group_id(MPI_Comm intercomm) {

  MPI_Group local_group, remote_group, world_group;
  xt_mpi_call(MPI_Comm_group(intercomm, &local_group), intercomm);
  xt_mpi_call(MPI_Comm_remote_group(intercomm, &remote_group), intercomm);
  xt_mpi_call(MPI_Comm_group(MPI_COMM_WORLD, &world_group), MPI_COMM_WORLD);
  int local_root_world_rank, remote_root_world_rank, root = 0;
  xt_mpi_call(
    MPI_Group_translate_ranks(local_group, 1, &root, world_group,
                              &local_root_world_rank), MPI_COMM_WORLD);
  xt_mpi_call(
    MPI_Group_translate_ranks(remote_group, 1, &root, world_group,
                              &remote_root_world_rank), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Group_free(&world_group), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Group_free(&remote_group), intercomm);
  xt_mpi_call(MPI_Group_free(&local_group), intercomm);

  return local_root_world_rank > remote_root_world_rank;
}

static void match_ids(void * redist_src_ids, int num_redist_src_ids,
                      void * redist_dst_ids, int num_redist_dst_ids,
                      int ** matching_src_id_idx_, size_t id_size,
                      MPI_Comm comm) {

  Xt_memcmp_func compare_id = xt_memcmp_get_func(id_size);

  // contains for each redistributed destination id the index of the matching
  // redistributed source id
  int * matching_src_id_idx =
    xmalloc((size_t)num_redist_dst_ids * sizeof(*matching_src_id_idx));
  *matching_src_id_idx_ = matching_src_id_idx;

  // match source and destination ids
  for (int i = 0, src_idx = 0; i < num_redist_dst_ids; ++i) {

    void * curr_dst_id = (char*)redist_dst_ids + (size_t)i * id_size;

    // find for the current destination id a matching source id
    int ret = -1;
    while ((src_idx < num_redist_src_ids) &&
           ((ret = compare_id((char*)redist_src_ids +
                              (size_t)src_idx * id_size, curr_dst_id))))
      ++src_idx;

    // if the current destination id has no matching source id
    if (ret) Xt_abort(comm, "ERROR(generate_exchange_data): unmatched "
                      "destination id", __FILE__, __LINE__);

    matching_src_id_idx[i] = src_idx;
  }
}

static void generate_exchange_data(
  void const * src_ids, int num_src_ids,
  void const * dst_ids, int num_dst_ids, size_t id_size,
  struct exchange_data ** send_data, struct exchange_data ** recv_data,
  MPI_Comm comm) {

  //----------------------------------------------------------------------------
  // in case comm is an intercommunicator we need some additional communication
  // data variables, which are generated here
  //----------------------------------------------------------------------------

  int is_inter, intercomm_group_id;
  xt_mpi_call(MPI_Comm_test_inter(comm, &is_inter), comm);

  // intracomm: local_comm_size = remote_comm_size = comm_size
  int local_comm_size, remote_comm_size;
  MPI_Comm ubound_comm; // comm used to generate ubounds
  MPI_Comm src_comm; // comm used to redistribute source ids

  // if comm is an intercommunicator
  if (is_inter) {

    // determine intercomm group id
    intercomm_group_id = get_intercomm_group_id(comm);

    xt_mpi_call(
      MPI_Intercomm_merge(comm, intercomm_group_id, &ubound_comm), comm);
    xt_mpi_call(
      MPI_Comm_split(ubound_comm, intercomm_group_id, 0, &src_comm),
      ubound_comm);

    xt_mpi_call(MPI_Comm_size(comm, &local_comm_size), comm);
    xt_mpi_call(MPI_Comm_remote_size(comm, &remote_comm_size), comm);

  } else {
    intercomm_group_id = 0;
    xt_mpi_call(MPI_Comm_size(comm, &local_comm_size), comm);
    remote_comm_size = local_comm_size;
    ubound_comm = comm;
    src_comm = comm;
  }

  //----------------------------------------------------------------------------
  // Redistribute source and destination ids such that identical ids are on
  // the same processes. This is done by assigning buckets to each process.
  // All processes locally have the upper bounds for each bucket.
  //----------------------------------------------------------------------------

  // generate upper bounds for all buckets
  void * src_ubounds, * dst_ubounds;
  if (is_inter) {

    // source ids of one group are matched with the destination ids of the other
    // group. Therefore we need two ubounds in case comm is an intercomm.
    // The src_ubounds of the local group are the dst_ubounds of the remote one.
    src_ubounds = xmalloc(((size_t)local_comm_size + (size_t)remote_comm_size) *
                          id_size);
    dst_ubounds = (char*)src_ubounds + (size_t)local_comm_size * id_size;
    if (intercomm_group_id) {
      generate_ubounds(src_ids, (size_t)num_src_ids, NULL, 0, id_size,
                       ubound_comm, (size_t)local_comm_size, src_ubounds);
      generate_ubounds(NULL, 0, dst_ids, (size_t)num_dst_ids, id_size,
                       ubound_comm, (size_t)remote_comm_size, dst_ubounds);
    } else {
      generate_ubounds(NULL, 0, dst_ids, (size_t)num_dst_ids, id_size,
                       ubound_comm, (size_t)remote_comm_size, dst_ubounds);
      generate_ubounds(src_ids, (size_t)num_src_ids, NULL, 0, id_size,
                       ubound_comm, (size_t)local_comm_size, src_ubounds);
    }
    xt_mpi_call(MPI_Comm_free(&ubound_comm), comm);

  } else {

      void * ubounds = xmalloc((size_t)local_comm_size * id_size);
      generate_ubounds(src_ids, (size_t)num_src_ids,
                       dst_ids, (size_t)num_dst_ids,
                       id_size, ubound_comm, (size_t)local_comm_size, ubounds);
      src_ubounds = dst_ubounds = ubounds;
  }

    // generate data type for id's
  MPI_Datatype ID_TYPE;
  xt_mpi_call(MPI_Type_contiguous((int)id_size, MPI_BYTE, &ID_TYPE), comm);
  xt_mpi_call(MPI_Type_commit(&ID_TYPE), comm);

  // redistribute source and destination ids according to previously computed
  // buckets
  // intercomm: source ids are redistributed within the local group, while the
  //            destination ids are redistributed within the remote group
  void * redist_src_ids, * redist_dst_ids;
  int * redist_src_id_ranks, * redist_dst_id_ranks;
  int * redist_src_id_orig_pos, * redist_dst_id_orig_pos;
  int num_redist_src_ids, num_redist_dst_ids;
  redistribute_ids(
    src_ids, num_src_ids, id_size, &redist_src_ids, &redist_src_id_ranks,
    &redist_src_id_orig_pos, &num_redist_src_ids, src_ubounds, ID_TYPE,
    src_comm);
  redistribute_ids(
    dst_ids, num_dst_ids, id_size, &redist_dst_ids, &redist_dst_id_ranks,
    &redist_dst_id_orig_pos, &num_redist_dst_ids, dst_ubounds, ID_TYPE,
    comm);

  xt_mpi_call(MPI_Type_free(&ID_TYPE), comm);
  free(src_ubounds);

  // determine for each destination id the matching source id
  int * matching_src_id_idx;
  match_ids(
    redist_src_ids, num_redist_src_ids, redist_dst_ids, num_redist_dst_ids,
    &matching_src_id_idx, id_size, comm);

  free(redist_src_ids);
  free(redist_dst_ids);

  //----------------------------------------------------------------------------
  // inform all owner of source ids that match with destination ids about the
  // matches
  // intercomm: all source id owners are in the local group
  //----------------------------------------------------------------------------

  size_t max_num_procs = (size_t)(MAX(local_comm_size, remote_comm_size));
  int * send_count = xmalloc((4 * max_num_procs + 1) * sizeof(*send_count));
  int * recv_count = send_count + 1 * max_num_procs;
  int * send_displ = send_count + 2 * max_num_procs;
  int * recv_displ = send_count + 3 * max_num_procs + 1;
  int total_recv_count;

  // compute the number of source elements each process has to send in an
  // exchange (per element we send dst_rank, orig_src_pos, and orig_dst_pos
  for (int i = 0; i < local_comm_size; ++i) send_count[i] = 0;
  for (int i = 0; i < num_redist_dst_ids; ++i)
    send_count[redist_src_id_ranks[matching_src_id_idx[i]]] += 3;

  compute_alltoallv_parameters(
    send_count, send_displ, recv_count, recv_displ, &total_recv_count,
    src_comm);

  int * send_buffer =
    xmalloc((3 * (size_t)num_redist_dst_ids + total_recv_count) * sizeof(int));
  int * recv_buffer = send_buffer + 3 * (size_t)num_redist_dst_ids;

  // pack the send buffer
  for (int i = 0; i < num_redist_dst_ids; ++i) {

    int src_id_idx = matching_src_id_idx[i];

    int src_rank = redist_src_id_ranks[src_id_idx];
    int dst_rank = redist_dst_id_ranks[i];
    int src_orig_pos = redist_src_id_orig_pos[src_id_idx];
    int dst_orig_pos = redist_dst_id_orig_pos[i];

    int * curr_buffer = send_buffer + send_displ[src_rank+1];
    send_displ[src_rank+1] += 3;
    curr_buffer[0] = dst_rank;
    curr_buffer[1] = src_orig_pos;
    curr_buffer[2] = dst_orig_pos;
  }

  free(matching_src_id_idx);
  free(redist_src_id_ranks);
  free(redist_dst_id_ranks);
  free(redist_src_id_orig_pos);
  free(redist_dst_id_orig_pos);

  // exchange results
  xt_mpi_call(MPI_Alltoallv(send_buffer, send_count, send_displ, MPI_INT,
                            recv_buffer, recv_count, recv_displ, MPI_INT,
                            src_comm), src_comm);

  if (is_inter) xt_mpi_call(MPI_Comm_free(&src_comm), comm);

  //----------------------------------------------------------------------------
  // generate send data
  // in addition, the new send count and send_buffer (containing the positions
  // of the destination ids) are filled
  //----------------------------------------------------------------------------

  *send_data = generate_basic_send_exchange_data(
    recv_buffer, total_recv_count, remote_comm_size, send_count, send_buffer);

  //----------------------------------------------------------------------------
  // owners of source ids that match with destination ids, will now inform the
  // respective processes from which the recv data is generated
  // intercomm: the owners of the destination ids are in the remote group
  //----------------------------------------------------------------------------

  compute_alltoallv_parameters(
    send_count, send_displ, recv_count, recv_displ, &total_recv_count, comm);

  // we need to generate recv_data first, because we receive the id positions
  // directly into the (*recv_data)->pos
  *recv_data =
    generate_basic_exchange_data(num_dst_ids, remote_comm_size, recv_count);

  // inform processes about data they will receive from the local process in an
  // exchange operation
  xt_mpi_call(MPI_Alltoallv(send_buffer, send_count, send_displ+1, MPI_INT,
                            &((*recv_data)->pos[0]), recv_count, recv_displ,
                            MPI_INT, comm), comm);

  free(send_buffer);
  free(send_count);
}

Xt_xmap xt_xmap_dist_radix_new(void const * src_ids, int num_src_ids,
                               void const * dst_ids, int num_dst_ids,
                               size_t id_size, MPI_Comm comm) {

  INSTR_DEF(t_xt_xmap_dist_radix_new,"xt_xmap_dist_radix_new")
  INSTR_START(t_xt_xmap_dist_radix_new);

  // ensure that yaxt is initialized
  assert(xt_initialized());

  int tag_offset;
  MPI_Comm newcomm = xt_mpi_comm_smart_dup(comm, &tag_offset);

  struct exchange_data * send_data, * recv_data;
  generate_exchange_data(src_ids, num_src_ids, dst_ids, num_dst_ids,
                         id_size, &send_data, &recv_data, newcomm);

  Xt_xmap xmap =
    xt_xmap_base_new(
      recv_data->pos, recv_data->num_msgs, recv_data->counts, recv_data->ranks,
      send_data->pos, send_data->num_msgs, send_data->counts, send_data->ranks,
      newcomm);

  free(send_data);
  free(recv_data);

  xt_mpi_comm_smart_dedup(&newcomm, tag_offset);
  INSTR_STOP(t_xt_xmap_dist_radix_new);
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
