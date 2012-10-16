#include <string.h>

#include <mpi.h>

#include "tests.h"

#include "yaxt.h"
#include "xt_idxlist.h"
#include "xt_idxvec.h"
#include "xt_redist_p2p.h"
#include "xt_redist_collection_static.h"
#include "xt_mpi.h"
#include "xt_xmap_all2all.h"

int main(void) {

  // init mpi

  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // general test with one redist
    // set up data
    Xt_idx src_index_list[] = {1,2,3,4,5};
    unsigned src_num_indices = sizeof(src_index_list) / sizeof(src_index_list[0]);
    Xt_idxlist src_idxlist;
    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    Xt_idx dst_index_list[] = {1,3,5};
    unsigned dst_num_indices = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
    Xt_idxlist dst_idxlist;
    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

    Xt_xmap xmap;
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);

    Xt_redist redist;
    redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    xt_xmap_delete(xmap);

    // generate redist_collection

    Xt_redist redist_coll;
    MPI_Aint src_displacements[1] = {0};
    MPI_Aint dst_displacements[1] = {0};

    redist_coll = xt_redist_collection_static_new(&redist, 1, src_displacements,
                                                  dst_displacements, MPI_COMM_WORLD);

    xt_redist_delete(redist);

    // test exchange

    double src_data[] = {1,2,3,4,5};
    double dst_data[] = {-1,-1,-1};

    xt_redist_s_exchange1(redist_coll, (void*)src_data, (void*)dst_data);

    double ref_dst_data[] = {1,3,5};

    unsigned i;

    for (i = 0; i < sizeof(dst_data) / sizeof(dst_data[0]); ++i)
      if (ref_dst_data[i] != dst_data[i])
        PUT_ERR("error in xt_redist_s_exchange\n");

    // clean up

    xt_redist_delete(redist_coll);
  }

  { // test with one redist used three times (two exchanges)
    // set up data
    Xt_idx src_index_list[] = {1,2,3,4,5};
    unsigned src_num_indices = sizeof(src_index_list) / sizeof(src_index_list[0]);
    Xt_idxlist src_idxlist;
    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    Xt_idx dst_index_list[] = {1,3,5};
    unsigned dst_num_indices = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
    Xt_idxlist dst_idxlist;
    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

    Xt_xmap xmap;
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);

    Xt_redist redist;
    redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    xt_xmap_delete(xmap);

    // generate redist_collection

    Xt_redist redist_coll;
    Xt_redist redists[3] = {redist, redist, redist};
    double src_data[3][5];
    double dst_data[3][3];
    MPI_Aint src_displacements[3] = {0,(src_data[1]-src_data[0])*sizeof(double),
                                     (src_data[2]-src_data[0])*sizeof(double)};
    MPI_Aint dst_displacements[3] = {0,(dst_data[1]-dst_data[0])*sizeof(double),
                                     (dst_data[2]-dst_data[0])*sizeof(double)};

    redist_coll = xt_redist_collection_static_new(redists, 3, src_displacements,
                                                  dst_displacements, MPI_COMM_WORLD);

    xt_redist_delete(redist);

    // test exchange
    {
      double src_data[3][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
      double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

      xt_redist_s_exchange1(redist_coll, (void*)src_data, (void*)dst_data);

      double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

      unsigned i, j;

      for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
          if (ref_dst_data[i][j] != dst_data[i][j])
            PUT_ERR("error in xt_redist_s_exchange\n");
    }

    {
      double src_data[3][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
      double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

      xt_redist_s_exchange1(redist_coll, (void*)src_data, (void*)dst_data);

      double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

      unsigned i, j;

      for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
          if (ref_dst_data[i][j] != dst_data[i][j])
            PUT_ERR("error in xt_redist_s_exchange\n");
    }

    // clean up

    xt_redist_delete(redist_coll);
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
