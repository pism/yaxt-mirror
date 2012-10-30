#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

#include "tests.h"

#include "yaxt.h"
#include "xt/xt_idxvec.h"
#include "xt_redist_p2p.h"
#include "xt/xt_idxlist.h"
#include "xt_mpi.h"
#include "xt_xmap_all2all.h"

int main(void) {

  // init mpi

  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  // offset-free test:
  {
    // source index list
    Xt_idx src_index_list[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
    Xt_count src_num_indices
      = sizeof(src_index_list) / sizeof(src_index_list[0]);

    Xt_idxlist src_idxlist;

    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    // destination index list
    Xt_idx dst_index_list[] = {41,5,4,7,3,96,1,5,48,54,4,5,4,6,3,13,7,2,1};
    Xt_count dst_num_indices
      = sizeof(dst_index_list) / sizeof(dst_index_list[0]);

    Xt_idxlist dst_idxlist;

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

    // xmap
    Xt_xmap xmap;

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // redist_p2p
    Xt_redist redist;

    redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    // test exchange

    double src_data[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
    double dst_data[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                         -1,-1,-1,-1,-1,-1,-1,-1,-1};

    double * src_data_p = &src_data[0];
    double * dst_data_p = &dst_data[0];

    xt_redist_s_exchange(redist, (void**)(&src_data_p), 1,
                         (void**)(&dst_data_p), 1);

    double ref_dst_data[] = {-1,0,2,-1,13,9,7,0,-1,-1,2,0,2,-1,13,4,-1,6,7};

    unsigned i;

    for (i = 0; i < sizeof(dst_data) / sizeof(dst_data[0]); ++i)
      if (ref_dst_data[i] != dst_data[i])
        PUT_ERR("error in xt_redist_s_exchange\n");

    // clean up

    xt_redist_delete(redist);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);
  }

  // check offsets
  {
    // source index list
    Xt_idx src_index_list[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
    unsigned src_num = sizeof(src_index_list) / sizeof(src_index_list[0]);

    Xt_idxlist src_idxlist;

    src_idxlist = xt_idxvec_new(src_index_list, src_num);

    // destination index list
    Xt_idx dst_index_list[] = {41,5,4,7,3,96,1,5,48,54,4,5,4,6,3,13,7,2,1};
    unsigned dst_num = sizeof(dst_index_list) / sizeof(dst_index_list[0]);

    Xt_idxlist dst_idxlist;

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num);

    // xmap
    Xt_xmap xmap;

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // redist_p2p with offsets

    unsigned i;

    int src_pos[src_num];
    int dst_pos[dst_num];

    for (i = 0; i < src_num; i++)
      src_pos[i] = i;
    for (i = 0; i < dst_num; i++)
      dst_pos[i] = dst_num-1-i;

    Xt_redist redist;
    redist = xt_redist_p2p_off_new(xmap, src_pos, dst_pos, MPI_DOUBLE);

    // test exchange

    double src_data[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
    double dst_data[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                         -1,-1,-1,-1,-1,-1,-1,-1,-1};
    assert(sizeof(src_data)/sizeof(src_data[0]) == src_num);
    assert(sizeof(dst_data)/sizeof(dst_data[0]) == dst_num);

    double * src_data_p = &src_data[0];
    double * dst_data_p = &dst_data[0];

    xt_redist_s_exchange(redist, (void**)(&src_data_p), 1,
                         (void**)(&dst_data_p), 1);

    double ref_dst_data[] = {-1,0,2,-1,13,9,7,0,-1,-1,
                             2,0,2,-1,13,4,-1,6,7};

    for (i = 0; i < dst_num; ++i)
      if (ref_dst_data[dst_pos[i]] != dst_data[i])
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
