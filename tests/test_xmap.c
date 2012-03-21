#include "mpi.h"

#include "tests.h"

#include "yaxt.h"
#include "xt_mpi.h"
#include "xt_idxlist.h"

int main(void) {

   // init mpi

   mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

   {
      // source index set
      int src_index_list[] = {1,2,3,4,5,6,7};
      unsigned src_num_indices = sizeof(src_index_list) / sizeof(src_index_list[0]);

      xt_idxset src_idxset;

      src_idxset = xt_idxlist_new(src_index_list, src_num_indices);

      // destination index set
      int dst_index_list[] = {7,6,5,4,3,2,1};
      unsigned dst_num_indices = sizeof(dst_index_list) / sizeof(dst_index_list[0]);

      xt_idxset dst_idxset;

      dst_idxset = xt_idxlist_new(dst_index_list, dst_num_indices);

      // test of exchange map
      xt_xmap xmap;

      xmap = xt_xmap_new_single(src_idxset, dst_idxset, MPI_COMM_WORLD);

      // clean up
      xt_xmap_delete(xmap);
      xt_idxset_delete(src_idxset);
      xt_idxset_delete(dst_idxset);
   }

   MPI_Finalize();

   return TEST_EXIT_CODE;
}
