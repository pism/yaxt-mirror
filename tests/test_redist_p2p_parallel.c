#include <stdlib.h>

#include <mpi.h>

#include "tests.h"

#include "yaxt.h"
#include "xt_idxvec.h"
#include "xt_redist_p2p.h"
#include "xt_idxlist.h"
#include "xt_mpi.h"

int main(void) {

  // init mpi

  Xt_count i;
  int rank, size;

  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  Xt_initialize(MPI_COMM_WORLD);

  mpi_err_handler(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);
  mpi_err_handler(MPI_Comm_size(MPI_COMM_WORLD, &size), MPI_COMM_WORLD);

  {
    // source index list
    Xt_idx src_index_list[10];
    Xt_count src_num_indices = 10;
    for (i = 0; i < src_num_indices; ++i)
      src_index_list[i] = rank * 10 + i;

    xt_idxlist src_idxlist;

    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    // destination index list
    Xt_idx dst_index_list[10];
    Xt_count dst_num_indices = 10;
    for (i = 0; i < dst_num_indices; ++i)
      dst_index_list[i] = (rank * 10 + i + 2) % (size * 10);

    xt_idxlist dst_idxlist;

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

    // xmap
    xt_xmap xmap;

    xmap = xt_xmap_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // redist_p2p
    Xt_redist redist;

    redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    // test exchange

    double src_data[10];
    double dst_data[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

    for (i = 0; i < src_num_indices; ++i)
      src_data[i] = rank * 10 + i;

    double * src_data_p = &src_data[0];
    double * dst_data_p = &dst_data[0];

    xt_redist_s_exchange(redist, (void**)(&src_data_p), 1, (void**)(&dst_data_p), 1);

    unsigned i;

    for (i = 0; i < sizeof(dst_data) / sizeof(dst_data[0]); ++i)
      if (dst_index_list[i] != dst_data[i])
        PUT_ERR("error in xt_redist_s_exchange\n");

    // clean up

    xt_redist_delete(redist);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);
  }

  MPI_Finalize();

 return TEST_EXIT_CODE;
}
