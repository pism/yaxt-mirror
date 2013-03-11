
#include <mpi.h>
#include <yaxt.h>

#include "tests.h"
#include "test_xmap_common.h"

static void
test_xmap_allgather_analog(xmap_constructor xmap_new);
static void
test_ring_1d(xmap_constructor xmap_new);
static void
test_pair(xmap_constructor xmap_new);
static void
test_ping_pong(xmap_constructor xmap_new);

static int comm_rank, comm_size;

int
xt_xmap_parallel_test_main(xmap_constructor xmap_new)
{
  // init mpi
  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  test_xmap_allgather_analog(xmap_new);

  if (comm_size > 2) // skip test if there are not enough processes
    test_ring_1d(xmap_new);

  if (comm_size == 2)
    test_pair(xmap_new);

  if (comm_size > 1)
    test_ping_pong(xmap_new);

  xt_finalize();
  MPI_Finalize();
  return TEST_EXIT_CODE;
}

static void
test_xmap_allgather_analog(xmap_constructor xmap_new)
{
  // test in which every process requests data from all processes
  // source index list
  Xt_int src_index_list[1] = {comm_rank};

  Xt_idxlist src_idxlist = xt_idxvec_new(src_index_list, 1);

  // destination index list
  Xt_int dst_index_list[comm_size];

  for (int i = 0; i < comm_size; ++i)
    dst_index_list[i] = i;

  Xt_idxlist dst_idxlist;

  dst_idxlist = xt_idxvec_new(dst_index_list, comm_size);

  // test of exchange map
  Xt_xmap xmap;

  xmap = xmap_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

  // test results

  if (xt_xmap_get_num_destinations(xmap) != comm_size)
    PUT_ERR("error in xt_xmap_get_num_destinations\n");

  if (xt_xmap_get_num_sources(xmap) != comm_size)
    PUT_ERR("error in xt_xmap_get_num_sources\n");

  int ranks[comm_size];

  xt_xmap_get_destination_ranks(xmap, ranks);
  for (int i = 0; i < comm_size; ++i)
    if (ranks[i] != i)
      PUT_ERR("error in xt_xmap_get_destination_ranks\n");

  xt_xmap_get_source_ranks(xmap, ranks);
  for (int i = 0; i < comm_size; ++i)
    if (ranks[i] != i)
      PUT_ERR("error in xt_xmap_get_source_ranks\n");

  // clean up
  xt_xmap_delete(xmap);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);
}

static void
test_ring_1d(xmap_constructor xmap_new)
{
  // test in which each process talks with two other
  // processes

  Xt_int src_index_list[1] = {comm_rank};

  Xt_idxlist src_idxlist;

  src_idxlist = xt_idxvec_new(src_index_list, 1);

  // destination index list
  Xt_int dst_index_list[2] = {(comm_rank + comm_size - 1)%comm_size,
                              (comm_rank             + 1)%comm_size};

  if (dst_index_list[0] > dst_index_list[1]) {
    dst_index_list[0] ^= dst_index_list[1];
    dst_index_list[1] ^= dst_index_list[0];
    dst_index_list[0] ^= dst_index_list[1];
  }

  Xt_idxlist dst_idxlist;

  dst_idxlist = xt_idxvec_new(dst_index_list, 2);

  // test of exchange map
  Xt_xmap xmap;

  xmap = xmap_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

  // test results

  if (xt_xmap_get_num_destinations(xmap) != 2)
    PUT_ERR("error in xt_xmap_get_num_destinations\n");

  if (xt_xmap_get_num_sources(xmap) != 2)
    PUT_ERR("error in xt_xmap_get_num_sources\n");

  int ranks[2];

  xt_xmap_get_destination_ranks(xmap, ranks);
  if (ranks[0] != dst_index_list[0] ||
      ranks[1] != dst_index_list[1])
    PUT_ERR("error in xt_xmap_get_destination_ranks\n");

  xt_xmap_get_source_ranks(xmap, ranks);
  if (ranks[0] != dst_index_list[0] ||
      ranks[1] != dst_index_list[1])
    PUT_ERR("error in xt_xmap_get_source_ranks\n");

  // clean up
  xt_xmap_delete(xmap);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);
}

static void
test_pair(xmap_constructor xmap_new)
{
  Xt_int src_index_list[2][20] = { //src_index_list[rank][index];
    {1,2,3,4,5,9,10,11,12,13,17,18,19,20,21,25,26,27,28,29},
    {4,5,6,7,8,12,13,14,15,16,20,21,22,23,24,28,29,30,31,32} };

  Xt_idxlist src_idxlist;

  src_idxlist = xt_idxvec_new(src_index_list[comm_rank], 20);

  // destination index list
  Xt_int dst_index_list[2][20] = { //dst_index_list[rank][index];
    {10,15,14,13,12,15,10,11,12,13,23,18,19,20,21,31,26,27,28,29},
    {13,12,11,10,15,12,13,14,15,10,20,21,22,23,18,28,29,30,31,26}};

  Xt_idxlist dst_idxlist;

  dst_idxlist = xt_idxvec_new(dst_index_list[comm_rank], 20);

  // test of exchange map
  Xt_xmap xmap = xmap_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);

  // test results

  if (xt_xmap_get_num_destinations(xmap) != 2)
    PUT_ERR("error in xt_xmap_get_num_destinations\n");

  if (xt_xmap_get_num_sources(xmap) != 2)
    PUT_ERR("error in xt_xmap_get_num_sources\n");

  int ranks[2];

  xt_xmap_get_destination_ranks(xmap, ranks);
  if (ranks[0] != 0 || ranks[1] != 1)
    PUT_ERR("error in xt_xmap_get_destination_ranks\n");

  xt_xmap_get_source_ranks(xmap, ranks);
  if (ranks[0] != 0 || ranks[1] != 1)
    PUT_ERR("error in xt_xmap_get_source_ranks\n");

  // clean up
  xt_xmap_delete(xmap);
}

static void
test_ping_pong(xmap_constructor xmap_new)
{
  Xt_int index_list[5] = {0,1,2,3,4};

  Xt_idxlist src_idxlist = (comm_rank == 0)?
    xt_idxvec_new(index_list, 5) : xt_idxempty_new();

  Xt_idxlist dst_idxlist = (comm_rank == comm_size-1)?
    xt_idxvec_new(index_list, 5) : xt_idxempty_new();

  // test of exchange map
  Xt_xmap xmap = xmap_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);

  // test results

  if (comm_rank == 0) {
    if (xt_xmap_get_num_destinations(xmap) != 1)
      PUT_ERR("error in xt_xmap_get_num_destinations (rank == 0)\n");
  } else {
    if (xt_xmap_get_num_destinations(xmap) != 0)
      PUT_ERR("error in xt_xmap_get_num_destinations (rank != 0)\n");
  }

  if (comm_rank == comm_size-1) {
    if (xt_xmap_get_num_sources(xmap) != 1)
      PUT_ERR("error in xt_xmap_get_num_sources (rank == 0)\n");
  } else {
    if (xt_xmap_get_num_sources(xmap) != 0)
      PUT_ERR("error in xt_xmap_get_num_sources (rank != 0)\n");
  }

  if (comm_rank == 0) {
    int dst_rank;
    xt_xmap_get_destination_ranks(xmap, &dst_rank);
    if (dst_rank != comm_size-1)
      PUT_ERR("error in xt_xmap_get_destination_ranks\n");
  }

  if (comm_rank == comm_size) {
    int src_rank;
    xt_xmap_get_source_ranks(xmap, &src_rank);
    if (src_rank != 0)
      PUT_ERR("error in xt_xmap_get_source_ranks\n");
  }

  // clean up
  xt_xmap_delete(xmap);
}
