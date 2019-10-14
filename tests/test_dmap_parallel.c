/**
 * @file test_dmap_parallel.c
 *
 * @copyright Copyright  (C)  2018 Jörg Behrens <behrens@dkrz.de>
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

#include <string.h>
#include <mpi.h>

#include <yaxt.h>

#include "../src/xt/quicksort.h"

#include "tests.h"

enum {
  NUM_RANKS = 4,
  NUM_INDICES_PER_RANK = 4,
};

static void test_dmap(Xt_xmap xmap, enum xt_reorder_type reorder_type) {

  int rank;
  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);

  const int src_offsets[NUM_RANKS][NUM_INDICES_PER_RANK] =
    {{0,1,2,3}, {3,2,1,0}, {0,2,3,1}, {3,1,2,0}};
  const int dst_offsets[NUM_RANKS][NUM_RANKS*NUM_INDICES_PER_RANK] =
    {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15},
     {15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0},
     {10,15,13,4,3,1,8,14,12,9,0,2,6,7,5,11},
     {7,12,6,9,5,4,11,14,1,8,10,2,13,3,15,0}};

  Xt_dmap dmap =
    xt_dmap_offset_new(xmap, src_offsets[rank], dst_offsets[rank]);

  xt_dmap_reorder(dmap, reorder_type);

  Xt_redist redist = xt_dmap_generate_redist(dmap, MPI_INT);

  // test mpi datatypes in redist
  for (int other_rank = 0; other_rank < NUM_RANKS; ++other_rank) {

    { // test send datatype (used to extract requested data from the local
      // source data)
      int local_send_offsets[NUM_INDICES_PER_RANK];  // offsets of the local points that need to be sent
      int remote_recv_offsets[NUM_INDICES_PER_RANK]; // offsets on the remote process for receiving the data
      int send_buffer[NUM_INDICES_PER_RANK];         // local send data
      int recv_buffer[NUM_INDICES_PER_RANK];         // contains data as it would arrive at the remote process
      int ref_recv_buffer[NUM_INDICES_PER_RANK];     // expected contents of recv_buffer

      memcpy(local_send_offsets, src_offsets[rank],
             NUM_INDICES_PER_RANK * sizeof(local_send_offsets[0]));
      for (int i = 0; i < NUM_INDICES_PER_RANK; ++i)
        send_buffer[src_offsets[rank][i]] = rank * NUM_INDICES_PER_RANK + i;
      for (int i = 0; i < NUM_INDICES_PER_RANK; ++i)
        remote_recv_offsets[i] =
          dst_offsets[other_rank][rank * NUM_INDICES_PER_RANK + i];
      for (int i = 0; i < NUM_INDICES_PER_RANK; ++i) recv_buffer[i] = -1;
      for (int i = 0; i < NUM_INDICES_PER_RANK; ++i)
        ref_recv_buffer[i] = rank * NUM_INDICES_PER_RANK + i;

      // by default yaxt sends data in an order that is based on the sorted
      // indices
      // apply reordering to ref_recv_buffer (based on send or receive offsets)
      if (reorder_type != XT_REORDER_NONE) {
        int *offsets_to_sort
          = reorder_type == XT_REORDER_SEND_UP ? local_send_offsets
          : remote_recv_offsets;
        xt_quicksort_int_permutation(
          offsets_to_sort, NUM_INDICES_PER_RANK, ref_recv_buffer);
      }

      MPI_Datatype send_datatype =
        xt_redist_get_send_MPI_Datatype(redist, other_rank);

      xt_mpi_call(
        MPI_Sendrecv(
          send_buffer, 1, send_datatype, 0, 0,
          recv_buffer, 4, MPI_INT, 0, 0,
          MPI_COMM_SELF, MPI_STATUS_IGNORE), MPI_COMM_WORLD);

      for (int i = 0; i < NUM_INDICES_PER_RANK; ++i)
        if (ref_recv_buffer[i] != recv_buffer[i])
          PUT_ERR("ERROR: wrong send datatype");

      xt_mpi_call(MPI_Type_free(&send_datatype), MPI_COMM_WORLD);
    }

    { // test recv datatype
      int local_recv_offsets[NUM_INDICES_PER_RANK];  // offsets on the local process for receiving the data
      int remote_send_offsets[NUM_INDICES_PER_RANK]; // offsets of the remote points that need to be sent
      int send_buffer[NUM_INDICES_PER_RANK];         // contains data as it would arrive from the remote process
      int recv_buffer[NUM_RANKS*NUM_INDICES_PER_RANK];     // contains the destination array
      int ref_recv_buffer[NUM_RANKS*NUM_INDICES_PER_RANK]; // expected contents of recv_buffer

      for (int i = 0; i < NUM_RANKS*NUM_INDICES_PER_RANK; ++i)
        ref_recv_buffer[i] = -1;
      for (int i = 0; i < NUM_INDICES_PER_RANK; ++i) {
        local_recv_offsets[i] =
          dst_offsets[rank][other_rank * NUM_INDICES_PER_RANK + i];
        ref_recv_buffer[local_recv_offsets[i]] =
          other_rank * NUM_INDICES_PER_RANK + i;
      }
      memcpy(remote_send_offsets, src_offsets[other_rank],
             NUM_INDICES_PER_RANK * sizeof(remote_send_offsets[0]));
      for (int i = 0; i < NUM_INDICES_PER_RANK; ++i)
        send_buffer[i] = other_rank * NUM_INDICES_PER_RANK + i;
      for (int i = 0; i < NUM_RANKS*NUM_INDICES_PER_RANK; ++i)
        recv_buffer[i] = -1;


      // by default yaxt sends data in an order that is based on the sorted
      // indices
      // apply reordering to send_buffer (based on send or receive offsets)
      if (reorder_type != XT_REORDER_NONE) {
        int *offsets_to_sort
          = reorder_type == XT_REORDER_SEND_UP ? remote_send_offsets
          : local_recv_offsets;
        xt_quicksort_int_permutation(
          offsets_to_sort, NUM_INDICES_PER_RANK, send_buffer);
      }

      MPI_Datatype recv_datatype =
        xt_redist_get_recv_MPI_Datatype(redist, other_rank);

      xt_mpi_call(
        MPI_Sendrecv(
          send_buffer, 4, MPI_INT, 0, 0,
          recv_buffer, 1, recv_datatype, 0, 0,
          MPI_COMM_SELF, MPI_STATUS_IGNORE), MPI_COMM_WORLD);

      for (int i = 0; i < NUM_RANKS*NUM_INDICES_PER_RANK; ++i)
        if (ref_recv_buffer[i] != recv_buffer[i])
          PUT_ERR("ERROR: wrong receive datatype");

      xt_mpi_call(MPI_Type_free(&recv_datatype), MPI_COMM_WORLD);
    }
  }

  xt_redist_delete(redist);
  xt_dmap_delete(dmap);
}

// test setup:
// The source array of each process contains the indices in the range
//   [rank*NUM_INDICES_PER_RANK;(rank+1)*NUM_INDICES_PER_RANK[
// The destination array contains all available indices, so that test generates
// a allgather communication pattern
//
// The indices in the source and destination arrays are sorted in ascending
// order. However, the offset arrays differ for all processes.
static void test1(MPI_Comm comm) {

  int rank;

  xt_mpi_call(MPI_Comm_rank(comm, &rank), comm);

  static Xt_int src_index_vector[NUM_INDICES_PER_RANK];
  static Xt_int dst_index_vector[NUM_RANKS*NUM_INDICES_PER_RANK];

  for (int i = 0; i < NUM_INDICES_PER_RANK; ++i)
    src_index_vector[i] = (Xt_int)(NUM_INDICES_PER_RANK * rank + i);
  for (int i = 0; i < NUM_RANKS*NUM_INDICES_PER_RANK; ++i)
    dst_index_vector[i] = (Xt_int)i;

  Xt_idxlist src_index_list =
    xt_idxvec_new(src_index_vector, NUM_INDICES_PER_RANK);
  Xt_idxlist dst_index_list =
    xt_idxvec_new(dst_index_vector, NUM_RANKS*NUM_INDICES_PER_RANK);

  Xt_xmap xmap =
    xt_xmap_all2all_new(src_index_list, dst_index_list, comm);

  enum xt_reorder_type reorder_types[] =
    {XT_REORDER_NONE, XT_REORDER_SEND_UP, XT_REORDER_RECV_UP};

  for (size_t i = 0; i < sizeof(reorder_types)/sizeof(reorder_types[0]); ++i)
    test_dmap(xmap, reorder_types[i]);

  xt_xmap_delete(xmap);
  xt_idxlist_delete(dst_index_list);
  xt_idxlist_delete(src_index_list);
}

// test setup:
// The source array of each process contains the indices in the range
//   [rank*NPROMA*NBLK;(rank+1)*NPROMA*NBLK[
// The destination array contains all available indices, so that test generates
// a allgather communication pattern.
// To simulate the data layout of ICON, a level dimension is introduced.
// Processes always have the whole level column for each 2D point.
// The exchange map is generated for the 2D domain, the third dimension is
// through the use xt_dmap_repeat.
static void test2(MPI_Comm comm) {

  int rank;
  xt_mpi_call(MPI_Comm_rank(comm, &rank), comm);

  enum {
    NPROMA = 4,
    NLEV = 20,
    NBLK = 6,
  };

  static Xt_int src_index_vector[NPROMA*NBLK];
  static Xt_int dst_index_vector[NPROMA*NBLK*NUM_RANKS];

  for (int i = 0; i < NPROMA*NBLK; ++i)
    src_index_vector[i] = (Xt_int)(NPROMA*NBLK * rank + i);
  for (int i = 0; i < NPROMA*NBLK*NUM_RANKS; ++i)
    dst_index_vector[i] = (Xt_int)i;

  Xt_idxlist src_index_list =
    xt_idxvec_new(src_index_vector, NPROMA*NBLK);
  Xt_idxlist dst_index_list =
    xt_idxvec_new(dst_index_vector, NPROMA*NBLK*NUM_RANKS);

  Xt_xmap xmap =
    xt_xmap_all2all_new(src_index_list, dst_index_list, comm);

  int src_offsets[NPROMA*NBLK];
  int dst_offsets[NPROMA*NBLK*NUM_RANKS];
  for (int i = 0, k = 0; i < NBLK; ++i)
    for (int j = 0; j < NPROMA; ++j, ++k)
      src_offsets[k] = j + i * NPROMA * NLEV;
  for (int i = 0, k = 0; i < NBLK*NUM_RANKS; ++i)
    for (int j = 0; j < NPROMA; ++j, ++k)
      dst_offsets[k] = j + i * NPROMA * NLEV;

  int lev_displacements[NLEV];
  for (int i = 0; i < NLEV; ++i) lev_displacements[i] = i * NPROMA;

  int src_data[NBLK][NLEV][NPROMA];
  int ref_dst_data[NBLK*NUM_RANKS][NLEV][NPROMA];
  for (int i = 0; i < NLEV; ++i)
    for (int j = 0; j < NBLK; ++j)
      for (int k = 0; k < NPROMA; ++k)
        src_data[j][i][k] = rank * NPROMA * NBLK + k + j * NPROMA +
                            i * NUM_RANKS * NPROMA * NBLK;
  for (int i = 0; i < NLEV; ++i)
    for (int j = 0; j < NBLK * NUM_RANKS; ++j)
      for (int k = 0; k < NPROMA; ++k)
        ref_dst_data[j][i][k] = k + j * NPROMA + i * NUM_RANKS * NPROMA * NBLK;

  enum xt_reorder_type reorder_types[] =
    {XT_REORDER_NONE, XT_REORDER_SEND_UP, XT_REORDER_RECV_UP};

  for (size_t i = 0; i < sizeof(reorder_types)/sizeof(reorder_types[0]); ++i) {

    Xt_dmap dmap = xt_dmap_offset_new(xmap, src_offsets, dst_offsets);
    xt_dmap_repeat(dmap, NLEV, lev_displacements, lev_displacements);

    xt_dmap_reorder(dmap, reorder_types[i]);

    Xt_redist redist = xt_dmap_generate_redist(dmap, MPI_INT);

    int dst_data[NBLK*NUM_RANKS][NLEV][NPROMA];
    for (int j = 0; j < NBLK*NUM_RANKS; ++j)
      for (int k = 0; k < NLEV; ++k)
        for (int l = 0; l < NPROMA; ++l)
          dst_data[j][k][l] = -1;

    xt_redist_s_exchange1(redist, &src_data[0][0][0], &dst_data[0][0][0]);

    if (memcmp(&dst_data[0][0][0], &ref_dst_data[0][0][0],
               NBLK*NUM_RANKS*NLEV*NPROMA * sizeof(dst_data[0][0][0])))
      PUT_ERR("ERROR: xt_dmap_repeat failed");

    xt_redist_delete(redist);
    xt_dmap_delete(dmap);
  }

  xt_xmap_delete(xmap);
  xt_idxlist_delete(dst_index_list);
  xt_idxlist_delete(src_index_list);
}

int main(void) {

  // init mpi
  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);
  xt_initialize(MPI_COMM_WORLD);

  int size;
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &size), MPI_COMM_WORLD);

  if (size == 4) {

    test1(MPI_COMM_WORLD);
    test2(MPI_COMM_WORLD);

  } else PUT_ERR("Wrong number of processes");

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
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
