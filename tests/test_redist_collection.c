#include <string.h>

#include <mpi.h>

#include "tests.h"

#include "yaxt.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_redist_p2p.h"
#include "xt/xt_redist_collection.h"
#include "xt/xt_mpi.h"
#include "xt_xmap_all2all.h"

int main(void) {

  // init mpi

  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // general test with one redist
    // set up data
    Xt_idx src_index_list[] = {1,2,3,4,5};
    unsigned src_num_indices
      = sizeof(src_index_list) / sizeof(src_index_list[0]);
    Xt_idxlist src_idxlist;
    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    Xt_idx dst_index_list[] = {1,3,5};
    unsigned dst_num_indices
      = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
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

    redist_coll = xt_redist_collection_new(&redist, 1, -1, MPI_COMM_WORLD);

    xt_redist_delete(redist);

    // test exchange

    double src_data[] = {1,2,3,4,5};
    double dst_data[] = {-1,-1,-1};

    double * src_data_p = &src_data[0];
    double * dst_data_p = &dst_data[0];

    xt_redist_s_exchange(redist_coll, (void**)(&src_data_p), 1,
                         (void**)(&dst_data_p), 1);

    double ref_dst_data[] = {1,3,5};

    unsigned i;

    for (i = 0; i < sizeof(dst_data) / sizeof(dst_data[0]); ++i)
      if (ref_dst_data[i] != dst_data[i])
        PUT_ERR("error in xt_redist_s_exchange\n");

    // clean up

    xt_redist_delete(redist_coll);
  }

  { // test with one redist used three times (with two different input data
    // displacements -> test of cache) (with default cache size)
    // set up data
    Xt_idx src_index_list[] = {1,2,3,4,5};
    unsigned src_num_indices
      = sizeof(src_index_list) / sizeof(src_index_list[0]);
    Xt_idxlist src_idxlist;
    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    Xt_idx dst_index_list[] = {1,3,5};
    unsigned dst_num_indices
      = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
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

    redist_coll = xt_redist_collection_new(redists, 3, -1, MPI_COMM_WORLD);

    xt_redist_delete(redist);

    // test exchange
    {
      double src_data[3][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
      double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

      double *src_data_p[3] = {src_data[0],src_data[1],src_data[2]};
      double *dst_data_p[3] = {dst_data[0],dst_data[1],dst_data[2]};

      xt_redist_s_exchange(redist_coll, (void**)src_data_p, 3,
                           (void**)dst_data_p, 3);

      double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

      unsigned i, j;

      for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
          if (ref_dst_data[i][j] != dst_data[i][j])
            PUT_ERR("error in xt_redist_s_exchange\n");
    }

    // test exchange with changed displacements
    {
      double src_data[3][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
      double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

      double *src_data_p[3] = {src_data[1],src_data[0],src_data[2]};
      double *dst_data_p[3] = {dst_data[1],dst_data[0],dst_data[2]};

      xt_redist_s_exchange(redist_coll, (void**)src_data_p, 3,
                           (void**)dst_data_p, 3);

      double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

      unsigned i, j;

      for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
          if (ref_dst_data[i][j] != dst_data[i][j])
            PUT_ERR("error in xt_redist_s_exchange\n");
    }

    // test exchange with original displacements
    {
      double src_data[3][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
      double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

      double *src_data_p[3] = {src_data[0],src_data[1],src_data[2]};
      double *dst_data_p[3] = {dst_data[0],dst_data[1],dst_data[2]};

      xt_redist_s_exchange(redist_coll, (void**)src_data_p, 3,
                           (void**)dst_data_p, 3);

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

  { // test with one redist used three times (with two different input data
    // displacements -> test of cache) (with cache size == 0)
    // set up data
    Xt_idx src_index_list[] = {1,2,3,4,5};
    unsigned src_num_indices
      = sizeof(src_index_list) / sizeof(src_index_list[0]);
    Xt_idxlist src_idxlist;
    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    Xt_idx dst_index_list[] = {1,3,5};
    unsigned dst_num_indices
      = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
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

    redist_coll = xt_redist_collection_new(redists, 3, 0, MPI_COMM_WORLD);

    xt_redist_delete(redist);

    // test exchange
    {
      double src_data[3][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
      double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

      double *src_data_p[3] = {src_data[0],src_data[1],src_data[2]};
      double *dst_data_p[3] = {dst_data[0],dst_data[1],dst_data[2]};

      xt_redist_s_exchange(redist_coll, (void**)src_data_p, 3,
                           (void**)dst_data_p, 3);

      double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

      unsigned i, j;

      for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
          if (ref_dst_data[i][j] != dst_data[i][j])
            PUT_ERR("error in xt_redist_s_exchange\n");
    }

    // test exchange with changed displacements
    {
      double src_data[3][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
      double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

      double *src_data_p[3] = {src_data[1],src_data[0],src_data[2]};
      double *dst_data_p[3] = {dst_data[1],dst_data[0],dst_data[2]};

      xt_redist_s_exchange(redist_coll, (void**)src_data_p, 3,
                           (void**)dst_data_p, 3);

      double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

      unsigned i, j;

      for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
          if (ref_dst_data[i][j] != dst_data[i][j])
            PUT_ERR("error in xt_redist_s_exchange\n");
    }

    // test exchange with original displacements
    {
      double src_data[3][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
      double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

      double *src_data_p[3] = {src_data[0],src_data[1],src_data[2]};
      double *dst_data_p[3] = {dst_data[0],dst_data[1],dst_data[2]};

      xt_redist_s_exchange(redist_coll, (void**)src_data_p, 3,
                           (void**)dst_data_p, 3);

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

  { // test with one redist used three times (with different input
    // data displacements until the cache is full)
    // set up data
    Xt_idx src_index_list[] = {1,2,3,4,5};
    unsigned src_num_indices
      = sizeof(src_index_list) / sizeof(src_index_list[0]);
    Xt_idxlist src_idxlist;
    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    Xt_idx dst_index_list[] = {1,3,5};
    unsigned dst_num_indices
      = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
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

    redist_coll = xt_redist_collection_new(redists, 3, -1, MPI_COMM_WORLD);

    xt_redist_delete(redist);

    double src_data[3][5] = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

    const int cache_size = 16;

    double * src_data_, * dst_data_;

    src_data_ = malloc((5 + cache_size + 2) * sizeof(*src_data_));
    dst_data_ = malloc((3 + cache_size + 2) * sizeof(*dst_data_));

    double *src_data_p[3] = {src_data[0],src_data[1],NULL};
    double *dst_data_p[3] = {dst_data[0],dst_data[1],NULL};

    // test exchange
    for (int k = 0; k < cache_size + 2; ++k) {

      memcpy(src_data_+k, src_data[2], 5 * sizeof(*src_data_));
      memcpy(dst_data_+k, dst_data[2], 3 * sizeof(*dst_data_));

      src_data_p[2] = src_data_+k;
      dst_data_p[2] = dst_data_+k;

      xt_redist_s_exchange(redist_coll, (void**)src_data_p, 3,
                           (void**)dst_data_p, 3);

      double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

      unsigned i, j;

      for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
          if (ref_dst_data[i][j] != dst_data_p[i][j])
            PUT_ERR("error in xt_redist_s_exchange\n");
    }

    // clean up

    free(src_data_);
    free(dst_data_);
    xt_redist_delete(redist_coll);
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
