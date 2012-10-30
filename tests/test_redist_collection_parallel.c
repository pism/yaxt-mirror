#include <stdlib.h>

#include <mpi.h>
#include "tests.h"

#include "yaxt.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxsection.h"
#include "xt/xt_redist_p2p.h"
#include "xt/xt_redist_collection.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_idxlist_collection.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_mpi.h"
#include "xt/xt_xmap_all2all.h"

static inline int imin(int a, int b) {
  return  (((a)<(b))?(a):(b));
}

static inline int imax(int a, int b) {
  return  (((a)>(b))?(a):(b));
}

int main(void) {

  // init mpi

  int rank, size;

  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  mpi_err_handler(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);
  mpi_err_handler(MPI_Comm_size(MPI_COMM_WORLD, &size), MPI_COMM_WORLD);

  if (size > 1) {

    { // redist test with four different redists

      Xt_idxlist indices_a, indices_b, indices_all;

      {
        Xt_idxlist indices_a_[2];

        Xt_idx start = 0;
        unsigned global_size[2] = {2*size, size*size};
        unsigned local_size[2] = {size,size};
        unsigned local_start[2][2]
          = {{0,rank*size},{size,size*size-(rank+1)*size}};

        indices_a_[0] = xt_idxsection_new(start, 2, global_size, local_size,
                                            local_start[0]);
        indices_a_[1] = xt_idxsection_new(start, 2, global_size, local_size,
                                            local_start[1]);

        indices_a = xt_idxlist_collection_new(indices_a_, 2);

        xt_idxlist_delete(indices_a_[0]);
        xt_idxlist_delete(indices_a_[1]);
      }

      {
        struct Xt_stripe stripe = {rank*2*size*size, 2*size*size, 1};

        indices_b = xt_idxstripes_new(&stripe, 1);
      }

      {
        struct Xt_stripe stripe = {0, 2*size*size*size, 1};

        indices_all = xt_idxstripes_new(&stripe, 1);
      }

      Xt_idx index_vector_a[2*size*size], index_vector_b[2*size*size],
             index_vector_all[2*size*size*size];

      xt_idxlist_get_indices(indices_a, index_vector_a);
      xt_idxlist_get_indices(indices_b, index_vector_b);
      xt_idxlist_get_indices(indices_all, index_vector_all);

      Xt_xmap xmaps[4] = {xt_xmap_all2all_new(indices_a, indices_b,
                                              MPI_COMM_WORLD),
                          xt_xmap_all2all_new(indices_b, indices_a,
                                              MPI_COMM_WORLD),
                          xt_xmap_all2all_new(indices_a, indices_all,
                                              MPI_COMM_WORLD),
                          xt_xmap_all2all_new(indices_b, indices_all,
                                              MPI_COMM_WORLD)};

      xt_idxlist_delete(indices_a);
      xt_idxlist_delete(indices_b);
      xt_idxlist_delete(indices_all);

      Xt_redist redists[4] = {xt_redist_p2p_new(xmaps[0], Xt_idx_dt),
                              xt_redist_p2p_new(xmaps[1], Xt_idx_dt),
                              xt_redist_p2p_new(xmaps[2], Xt_idx_dt),
                              xt_redist_p2p_new(xmaps[3], Xt_idx_dt)};

      xt_xmap_delete(xmaps[0]), xt_xmap_delete(xmaps[1]);
      xt_xmap_delete(xmaps[2]), xt_xmap_delete(xmaps[3]);

      Xt_redist redist = xt_redist_collection_new(redists, 4, -1,
                                                  MPI_COMM_WORLD);

      xt_redist_delete(redists[0]), xt_redist_delete(redists[1]);
      xt_redist_delete(redists[2]), xt_redist_delete(redists[3]);

      Xt_idx results_1[2*size*size], results_2[2*size*size];
      Xt_idx results_3[2*size*size*size], results_4[2*size*size*size];

      Xt_idx * results[4] = {results_1, results_2, results_3, results_4};

      Xt_idx * input[4] = {index_vector_a, index_vector_b,
                           index_vector_a, index_vector_b};

      xt_redist_s_exchange(redist, (void**)input, 4, (void**)results, 4);

      // check results

      for (int i = 0; i < 2*size*size; ++i)
        if (results[0][i] != index_vector_b[i])
          PUT_ERR("error on xt_redist_s_exchange\n");

      for (int i = 0; i < 2*size*size; ++i)
        if (results[1][i] != index_vector_a[i])
          PUT_ERR("error on xt_redist_s_exchange\n");

      for (int i = 0; i < 2*size*size*size; ++i)
        if (results[2][i] != i)
          PUT_ERR("error on xt_redist_s_exchange\n");

      for (int i = 0; i < 2*size*size*size; ++i)
        if (results[3][i] != i)
          PUT_ERR("error on xt_redist_s_exchange\n");

      // clean up

      xt_redist_delete(redist);
    }

    { // redist test with two redists that do a round robin exchange in
      // different directions

      Xt_idxlist src_indices, dst_indices[2];

      Xt_idx src_indices_[5], dst_indices_[2][5];

      for (int i = 0; i < 5; ++i) {
        src_indices_[i] = rank * 5 + i;
        dst_indices_[0][i] = (src_indices_[i] + 1) % (size * 5);
        int temp = src_indices_[i] - 1;
        dst_indices_[1][i] = (temp < 0)?(size * 5 - 1):temp;
      }

      src_indices = xt_idxvec_new(src_indices_, 5);
      dst_indices[0] = xt_idxvec_new(dst_indices_[0], 5);
      dst_indices[1] = xt_idxvec_new(dst_indices_[1], 5);

      Xt_xmap xmaps[2] = {xt_xmap_all2all_new(src_indices, dst_indices[0],
                                              MPI_COMM_WORLD),
                          xt_xmap_all2all_new(src_indices, dst_indices[1],
                                              MPI_COMM_WORLD)};

      xt_idxlist_delete(src_indices);
      xt_idxlist_delete(dst_indices[0]);
      xt_idxlist_delete(dst_indices[1]);

      Xt_redist redists[2] = {xt_redist_p2p_new(xmaps[0], Xt_idx_dt),
                              xt_redist_p2p_new(xmaps[1], Xt_idx_dt)};

      xt_xmap_delete(xmaps[0]), xt_xmap_delete(xmaps[1]);

      Xt_redist redist = xt_redist_collection_new(redists, 2, -1,
                                                  MPI_COMM_WORLD);

      xt_redist_delete(redists[0]), xt_redist_delete(redists[1]);

      Xt_idx results_1[5] = {-1,-1,-1,-1,-1}, results_2[5] = {-1,-1,-1,-1,-1};

      Xt_idx * results[2] = {results_1, results_2};

      Xt_idx * input[2] = {src_indices_, src_indices_};

      xt_redist_s_exchange(redist, (void**)input, 2, (void**)results, 2);

      // check results

      for (int i = 0; i < 5; ++i) {
        if (results[0][i] != dst_indices_[0][i])
          PUT_ERR("error on xt_redist_s_exchange\n");
        if (results[1][i] != dst_indices_[1][i])
          PUT_ERR("error on xt_redist_s_exchange\n");
      }

      // clean up

      xt_redist_delete(redist);
    }
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
