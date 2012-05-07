#include <stdlib.h>

#include <mpi.h>

#include "tests.h"

#include "yaxt.h"
#include "xt_idxlist.h"
#include "xt_redist_p2p.h"

int main(void) {

  // init mpi

  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  {
    // source index vector
    Xt_idx src_index_list[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
    unsigned src_num_indices = sizeof(src_index_list) / sizeof(src_index_list[0]);

    xt_idxvec src_idxvec;

    src_idxvec = xt_idxlist_new(src_index_list, src_num_indices);

    // destination index vector
    Xt_idx dst_index_list[] = {41,5,4,7,3,96,1,5,48,54,4,5,4,6,3,13,7,2,1};
    unsigned dst_num_indices = sizeof(dst_index_list) / sizeof(dst_index_list[0]);

    xt_idxvec dst_idxvec;

    dst_idxvec = xt_idxlist_new(dst_index_list, dst_num_indices);

    // xmap
    xt_xmap xmap;

    xmap = xt_xmap_new(src_idxvec, dst_idxvec, MPI_COMM_WORLD);

    // redist_p2p
    Xt_redist redist;

    redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    // test exchange

    double src_data[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
    double dst_data[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

    double * src_data_p = &src_data[0];
    double * dst_data_p = &dst_data[0];

    xt_redist_s_exchange(redist, (void**)(&src_data_p), 1, (void**)(&dst_data_p), 1);

    // clean up

    xt_redist_delete(redist);
    xt_xmap_delete(xmap);
    xt_idxvec_delete(src_idxvec);
    xt_idxvec_delete(dst_idxvec);
  }

  MPI_Finalize();

 return TEST_EXIT_CODE;
}