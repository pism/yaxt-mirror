#include <mpi.h>

#include <yaxt.h>

#define VERBOSE
#include "tests.h"
#include "test_xmap_common.h"

static void
test_xmap1(xmap_constructor new_xmap);
static void
test_xmap2(xmap_constructor new_xmap);

static int my_rank;

int
xt_xmap_self_test_main(xmap_constructor new_xmap)
{
  // init mpi
  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  test_xmap1(new_xmap);
  test_xmap2(new_xmap);

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
test_xmap(const Xt_int *src_index_list, Xt_int num_src_indices,
          const Xt_int *dst_index_list, Xt_int num_dst_indices,
          xmap_constructor new_xmap)
{
  Xt_idxlist src_idxlist = xt_idxvec_new(src_index_list, num_src_indices);
  Xt_idxlist dst_idxlist = xt_idxvec_new(dst_index_list, num_dst_indices);

  // test of exchange map
  Xt_xmap xmap = new_xmap(src_idxlist, dst_idxlist, MPI_COMM_WORLD);
  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);

  // test results
  if (xt_xmap_get_num_destinations(xmap) != 1)
    PUT_ERR("error in xt_xmap_get_num_destinations\n");

  if (xt_xmap_get_num_sources(xmap) != 1)
    PUT_ERR("error in xt_xmap_get_num_sources\n");

  int rank;

  xt_xmap_get_destination_ranks(xmap, &rank);
  if (rank != my_rank)
    PUT_ERR("error in xt_xmap_get_destination_ranks\n");

  xt_xmap_get_source_ranks(xmap, &rank);
  if (rank != my_rank)
    PUT_ERR("error in xt_xmap_get_source_ranks\n");
  // clean up
  xt_xmap_delete(xmap);
}

static inline void shift_idx(Xt_int idx[], int num, int offset)  {
  for (int i=0; i<num; i++) {
    idx[i] = (Xt_int)(idx[i] + my_rank * offset);
  }
}

static void
test_xmap1(xmap_constructor new_xmap)
{
  // source index list
  Xt_int src_index_list[] = {1,2,3,4,5,6,7};
  int num_src_indices
    = sizeof(src_index_list) / sizeof(src_index_list[0]);
  shift_idx(src_index_list, num_src_indices, 7);

  // destination index list
  Xt_int dst_index_list[] = {7,6,5,4,3,2,1};
  int num_dst_indices
    = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
  shift_idx(dst_index_list, num_dst_indices, 7);

  test_xmap(src_index_list, (Xt_int)num_src_indices, dst_index_list,
            (Xt_int)num_dst_indices, new_xmap);
}

static void
test_xmap2(xmap_constructor new_xmap)
{
  // source index list
  Xt_int src_index_list[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
  int num_src_indices
    = sizeof(src_index_list) / sizeof(src_index_list[0]);
  shift_idx(src_index_list, num_src_indices, 100);

  // destination index list
  Xt_int dst_index_list[] = {5,4,3,96,1,5,4,5,4,3,13,2,1};
  int num_dst_indices
    = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
  shift_idx(dst_index_list, num_dst_indices, 100);

  test_xmap(src_index_list, (Xt_int)num_src_indices,
            dst_index_list, (Xt_int)num_dst_indices,
            new_xmap);
}
