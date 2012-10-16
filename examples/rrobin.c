#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "mpi.h"

#include "yaxt.h"
#include "xt_mpi.h"
#include "xt_xmap_all2all.h"
#include "xt_idxlist.h"
#include "xt_redist_p2p.h"
#include "xt_idxstripes.h"

int main(void) {

  //init mpi

  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize (MPI_COMM_WORLD);

  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  {
    // init source and destination data array, local data 5 elements length
    int leng = 5;
    int src_array[leng];
    int dst_array[leng];

    for (int i = 0; i < leng; i++) {
      // first position = 1, second position = rank, third position element_index; for rank < 10, leng < 10
      src_array[i] = 100 + rank*10 + i;
      dst_array[i] = -1;

      //print source
      printf("SOURCEvalue: %d, element_index: %d, rank: %d \n", src_array[i], i, rank);
    }

    // output barrier
    MPI_Barrier (MPI_COMM_WORLD);

    // source index list by stripes
    Xt_idxlist src_idxlist;

    struct Xt_stripe src_stripes = {rank*leng, leng, 1};
    src_idxlist = xt_idxstripes_new(&src_stripes, 1);

    // destination index list by stripes
    Xt_idxlist dst_idxlist;

    struct Xt_stripe dst_stripes = {((rank+1)*leng)%(size*leng), leng, 1};
    dst_idxlist = xt_idxstripes_new(&dst_stripes, 1);

    // xmap
    Xt_xmap xmap;

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // redist_p2p
    Xt_redist redist;

    redist = xt_redist_p2p_new(xmap, MPI_INTEGER);

    //array poiter, especially necessary for data array number > 1
    int* src_array_p = &src_array[0];
    int* dst_array_p = &dst_array[0];

    //Exchange
    xt_redist_s_exchange(redist, (void**)(&src_array_p), 1, (void**)(&dst_array_p), 1);

    // print result
    // for emelent_index < 10, leng < 10:
    // first and last positions: should stay; second position: should be +1, except last rank should have 0
    for (int p = 0; p < leng; p++)
      printf("DESTvalue: %d, element_index: %d, rank: %d \n", dst_array[p], p, rank);


    xt_redist_delete(redist);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(dst_idxlist);
    xt_idxlist_delete(src_idxlist);
  }

  MPI_Finalize();

  return 0;
}
