#include <mpi.h>

#include "tests.h"

#include "yaxt.h"
#include "xt_mpi.h"
#include "xt_idxvec.h"
#include "xt_xmap.h"
#include "xt_xmap_dist_dir.h"

static void
check_idxdata(Xt_idxlist idxlist, Xt_idx * ref_indices,
              Xt_count ref_num_indices);

int main(void) {

  // init mpi
  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  unsigned i;

  { // test in which every process requests data from all processes
    // source index list

    Xt_idx src_index_list[1] = {comm_rank};

    Xt_idxlist src_idxlist;

    src_idxlist = xt_idxvec_new(src_index_list, 1);

    // destination index list
    Xt_idx dst_index_list[comm_size];

    for (i = 0; i < comm_size; ++i)
      dst_index_list[i] = i;

    Xt_idxlist dst_idxlist;

    dst_idxlist = xt_idxvec_new(dst_index_list, comm_size);

    // test of exchange map
    Xt_xmap xmap;

    xmap = xt_xmap_dist_dir_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // test results

    if (xt_xmap_get_num_destinations(xmap) != comm_size)
      PUT_ERR("error in xt_xmap_get_num_destinations\n");

    if (xt_xmap_get_num_sources(xmap) != comm_size)
      PUT_ERR("error in xt_xmap_get_num_destinations\n");

    int ranks[comm_size];

    xt_xmap_get_destination_ranks(xmap, ranks);
    for (i = 0; i < comm_size; ++i)
      if (ranks[i] != i)
        PUT_ERR("error in xt_xmap_get_destination_ranks\n");

    xt_xmap_get_source_ranks(xmap, ranks);
    for (i = 0; i < comm_size; ++i)
      if (ranks[i] != i)
        PUT_ERR("error in xt_xmap_get_destination_ranks\n");

    // clean up
    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);
  }

  if (comm_size > 2) // skip test if there are not enough processes
  { // test in which each process talks with two other
    // processes

    Xt_idx src_index_list[1] = {comm_rank};

    Xt_idxlist src_idxlist;

    src_idxlist = xt_idxvec_new(src_index_list, 1);

    // destination index list
    Xt_idx dst_index_list[2] = {(comm_rank + comm_size - 1)%comm_size,
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

    xmap = xt_xmap_dist_dir_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // test results

    if (xt_xmap_get_num_destinations(xmap) != 2)
      PUT_ERR("error in xt_xmap_get_num_destinations\n");

    if (xt_xmap_get_num_sources(xmap) != 2)
      PUT_ERR("error in xt_xmap_get_num_destinations\n");

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

  if (comm_size == 2) {

    Xt_idx src_index_list[2][20] = { //src_index_list[rank][index];
      {1,2,3,4,5,9,10,11,12,13,17,18,19,20,21,25,26,27,28,29},
      {4,5,6,7,8,12,13,14,15,16,20,21,22,23,24,28,29,30,31,32} };

    Xt_idxlist src_idxlist;

    src_idxlist = xt_idxvec_new(src_index_list[comm_rank], 20);

    // destination index list
    Xt_idx dst_index_list[2][20] = { //dst_index_list[rank][index];
      {10,15,14,13,12,15,10,11,12,13,23,18,19,20,21,31,26,27,28,29},
      {13,12,11,10,15,12,13,14,15,10,20,21,22,23,18,28,29,30,31,26}};

    Xt_idxlist dst_idxlist;

    dst_idxlist = xt_idxvec_new(dst_index_list[comm_rank], 20);

    // test of exchange map
    Xt_xmap xmap;

    xmap = xt_xmap_dist_dir_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

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
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);

  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
check_idxdata(Xt_idxlist idxlist, Xt_idx * ref_indices,
              Xt_count ref_num_indices) {

  unsigned i;

  if (xt_idxlist_get_num_indices(idxlist) != ref_num_indices)
    PUT_ERR("wrong number of indices in idxlist\n");

  Xt_idx indices[ref_num_indices];

  xt_idxlist_get_indices(idxlist, indices);

  for (i = 0; i < ref_num_indices; ++i)
    if (indices[i] != ref_indices[i])
      PUT_ERR("index list contains wrong indices\n");
}
