#include "mpi.h"

#define VERBOSE

#include "tests.h"

#include "yaxt.h"
#include "xt_mpi.h"
#include "xt_idxvec.h"
#include "xt_xmap.h"

void check_idxdata(xt_idxlist idxlist, Xt_idx * ref_indices, unsigned ref_num_indices);

int main(void) {

  // init mpi
  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  Xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  unsigned i;

  { // test in which every process requests data from all processes
    // source index list

    Xt_idx src_index_list[1] = {comm_rank};

    xt_idxlist src_idxlist;

    src_idxlist = xt_idxvec_new(src_index_list, 1);

    // destination index list
    Xt_idx dst_index_list[comm_size];

    for (i = 0; i < comm_size; ++i)
      dst_index_list[i] = i;

    xt_idxlist dst_idxlist;

    dst_idxlist = xt_idxvec_new(dst_index_list, comm_size);

    // test of exchange map
    xt_xmap xmap;

    xmap = xt_xmap_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

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

    for (i = 0; i < comm_size; ++i) {

      Xt_idxlist src_xlist, dst_xlist;

      dst_xlist = xt_xmap_get_dst_xlist(xmap, i);
      src_xlist = xt_xmap_get_src_xlist(xmap, i);

      Xt_idx ref_src_xlist[1] = {0};
      Xt_idx ref_dst_xlist[1] = {i};

      check_idxdata(src_xlist, ref_src_xlist, 1);
      check_idxdata(dst_xlist, ref_dst_xlist, 1);

      xt_idxlist_delete(dst_xlist);
      xt_idxlist_delete(src_xlist);
    }

    // clean up
    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);
  }
  
  if (comm_size > 2) // skip test if there are not enough processes
  { // test in which each process talks with two other
    // processes

    Xt_idx src_index_list[1] = {comm_rank};

    xt_idxlist src_idxlist;

    src_idxlist = xt_idxvec_new(src_index_list, 1);

    // destination index list
    Xt_idx dst_index_list[2] = {(comm_rank + comm_size - 1)%comm_size,
                                (comm_rank             + 1)%comm_size};

    if (dst_index_list[0] > dst_index_list[1]) {
      dst_index_list[0] ^= dst_index_list[1];
      dst_index_list[1] ^= dst_index_list[0];
      dst_index_list[0] ^= dst_index_list[1];
    }

    xt_idxlist dst_idxlist;

    dst_idxlist = xt_idxvec_new(dst_index_list, 2);

    // test of exchange map
    xt_xmap xmap;

    xmap = xt_xmap_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

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

    for (i = 0; i < 2; ++i) {

      Xt_idxlist src_xlist, dst_xlist;

      dst_xlist = xt_xmap_get_dst_xlist(xmap, dst_index_list[i]);
      src_xlist = xt_xmap_get_src_xlist(xmap, dst_index_list[i]);

      Xt_idx ref_src_xlist[1] = {0};
      Xt_idx ref_dst_xlist[1] = {i};

      check_idxdata(src_xlist, ref_src_xlist, 1);
      check_idxdata(dst_xlist, ref_dst_xlist, 1);

      xt_idxlist_delete(dst_xlist);
      xt_idxlist_delete(src_xlist);
    }

    // clean up
    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

void check_idxdata(xt_idxlist idxlist, Xt_idx * ref_indices, unsigned ref_num_indices) {

  unsigned i;

  if (xt_idxlist_get_num_indices(idxlist) != ref_num_indices)
    PUT_ERR("wrong number of indices in idxlist\n");

  Xt_idx indices[ref_num_indices];

  xt_idxlist_get_indices(idxlist, indices);

  for (i = 0; i < ref_num_indices; ++i)
    if (indices[i] != ref_indices[i])
      PUT_ERR("index list contains wrong indices\n");
}
