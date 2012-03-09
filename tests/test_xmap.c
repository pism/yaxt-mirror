#include "mpi.h"

#include "tests.h"

#include "yaxt.h"

int main(void) {

   // source index set
   int src_index_list[] = {1,2,3,4,5,6,7};
   unsigned src_num_indices = sizeof(src_index_list) / sizeof(src_index_list[0]);

   xt_idxset src_idxset;

   src_idxset = xt_idxset_new(src_index_list, src_num_indices);

   // destination index set
   int dst_index_list[] = {7,6,5,4,3,2,1};
   unsigned dst_num_indices = sizeof(dst_index_list) / sizeof(dst_index_list[0]);

   xt_idxset dst_idxset;

   dst_idxset = xt_idxset_new(dst_index_list, dst_num_indices);

   // test of exchange map
   xt_xmap xmap;

   xmap = xt_xmap_new(src_idxset, dst_idxset, MPI_COMM_WORLD);

   // clean up
   xt_xmap_delete(xmap);
   xt_idxset_delete(src_idxset);
   xt_idxset_delete(dst_idxset);

   return TEST_EXIT_CODE;
}
