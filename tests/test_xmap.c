#include "mpi.h"

#include "tests.h"

#include "yaxt.h"
#include "xt_mpi.h"
#include "xt_idxlist.h"
#include "xt_xmap.h"

void check_idxdata(xt_idxvec idxvec, Xt_idx * ref_indices, unsigned ref_num_indices);

int main(void) {

   // init mpi

   mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

   {
      // source index vector
      Xt_idx src_index_list[] = {1,2,3,4,5,6,7};
      unsigned src_num_indices = sizeof(src_index_list) / sizeof(src_index_list[0]);

      xt_idxvec src_idxvec;

      src_idxvec = xt_idxlist_new(src_index_list, src_num_indices);

      // destination index vector
      Xt_idx dst_index_list[] = {7,6,5,4,3,2,1};
      unsigned dst_num_indices = sizeof(dst_index_list) / sizeof(dst_index_list[0]);

      xt_idxvec dst_idxvec;

      dst_idxvec = xt_idxlist_new(dst_index_list, dst_num_indices);

      // test of exchange map
      xt_xmap xmap;

      xmap = xt_xmap_new(src_idxvec, dst_idxvec, MPI_COMM_WORLD);

      // test results

      Xt_idxvec src_xlist, dst_xlist;

      dst_xlist = xt_xmap_get_dst_xlist(xmap, 0);
      src_xlist = xt_xmap_get_src_xlist(xmap, 0);

      Xt_idx ref_src_xlist[7] = {0,1,2,3,4,5,6};
      Xt_idx ref_dst_xlist[7] = {6,5,4,3,2,1,0};

      check_idxdata(src_xlist, ref_src_xlist, 7);
      check_idxdata(dst_xlist, ref_dst_xlist, 7);

      // clean up
      xt_idxvec_delete(dst_xlist);
      xt_idxvec_delete(src_xlist);
      xt_xmap_delete(xmap);
      xt_idxvec_delete(src_idxvec);
      xt_idxvec_delete(dst_idxvec);
   }

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

      xt_xmap xmap;

      xmap = xt_xmap_new(src_idxvec, dst_idxvec, MPI_COMM_WORLD);

      // test results

      // the global id intersection should look as follows:
      // 1,1,2,3,3,4,4,4,5,5,5,13,96

      Xt_idxvec src_xlist, dst_xlist;

      dst_xlist = xt_xmap_get_dst_xlist(xmap, 0);
      src_xlist = xt_xmap_get_src_xlist(xmap, 0);

      Xt_idx ref_src_xlist[] = {7,7,6,13,13,2,2,2,0,0,0,4,9};
      Xt_idx ref_dst_xlist[] = {6,18,17,4,14,2,10,12,1,7,11,15,5};

      check_idxdata(src_xlist, ref_src_xlist,
                    sizeof(ref_src_xlist) / sizeof(ref_src_xlist[0]));
      check_idxdata(dst_xlist, ref_dst_xlist,
                    sizeof(ref_dst_xlist) / sizeof(ref_dst_xlist[0]));

      // clean up
      xt_idxvec_delete(dst_xlist);
      xt_idxvec_delete(src_xlist);
      xt_xmap_delete(xmap);
      xt_idxvec_delete(src_idxvec);
      xt_idxvec_delete(dst_idxvec);
   }

   MPI_Finalize();

   return TEST_EXIT_CODE;
}

void check_idxdata(xt_idxvec idxvec, Xt_idx * ref_indices, unsigned ref_num_indices) {

   unsigned i;

   if (xt_idxvec_get_num_indices(idxvec) != ref_num_indices)
      PUT_ERR("wrong number of indices in idxvec\n");

   Xt_idx indices[ref_num_indices];

   xt_idxvec_get_indices(idxvec, indices);

   for (i = 0; i < ref_num_indices; ++i)
      if (indices[i] != ref_indices[i])
         PUT_ERR("index vector contains wrong indices\n");
}
