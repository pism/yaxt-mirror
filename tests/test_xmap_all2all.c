#include "mpi.h"

#define VERBOSE

#include "tests.h"

#include "yaxt.h"
#include "xt_mpi.h"
#include "xt_idxvec.h"
#include "xt_xmap.h"
#include "xt_xmap_all2all.h"

static void
check_idxdata(Xt_idxlist idxlist, Xt_idx * ref_indices,
              Xt_count ref_num_indices);

static void
shift_idx(Xt_idx idx[], int num, int offset)  {
  int rank, i;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for (i=0; i<num; i++) {
    idx[i] = idx[i] + rank*offset;
  }
}

int main(void) {

   // init mpi
   mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

   xt_initialize(MPI_COMM_WORLD);

   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   {
      // source index list
      Xt_idx src_index_list[] = {1,2,3,4,5,6,7};
      Xt_count src_num_indices
        = sizeof(src_index_list) / sizeof(src_index_list[0]);
      shift_idx(src_index_list, src_num_indices, 7);

      Xt_idxlist src_idxlist;

      src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

      // destination index list
      Xt_idx dst_index_list[] = {7,6,5,4,3,2,1};
      Xt_count dst_num_indices
        = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
      shift_idx(dst_index_list, dst_num_indices, 7);

      Xt_idxlist dst_idxlist;

      dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

      // test of exchange map
      Xt_xmap xmap;

      xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

      // test results

      if (xt_xmap_get_num_destinations(xmap) != 1)
        PUT_ERR("error in xt_xmap_get_num_destinations\n");

      if (xt_xmap_get_num_sources(xmap) != 1)
        PUT_ERR("error in xt_xmap_get_num_destinations\n");

      int rank;

      xt_xmap_get_destination_ranks(xmap, &rank);
      if (rank != my_rank)
        PUT_ERR("error in xt_xmap_get_destination_ranks\n");

      xt_xmap_get_source_ranks(xmap, &rank);
      if (rank != my_rank)
        PUT_ERR("error in xt_xmap_get_source_ranks\n");

      // clean up
      xt_xmap_delete(xmap);
      xt_idxlist_delete(src_idxlist);
      xt_idxlist_delete(dst_idxlist);
   }

   {
      // source index list
      Xt_idx src_index_list[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
      Xt_count src_num_indices
        = sizeof(src_index_list) / sizeof(src_index_list[0]);
      shift_idx(src_index_list, src_num_indices, 100);

      Xt_idxlist src_idxlist;

      src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

      // destination index list
      Xt_idx dst_index_list[] = {41,5,4,7,3,96,1,5,48,54,4,5,4,6,3,13,7,2,1};
      Xt_count dst_num_indices
        = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
      shift_idx(dst_index_list, dst_num_indices, 100);

      Xt_idxlist dst_idxlist;

      dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

      Xt_xmap xmap;

      xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

      // test results

      // the global id intersection should look as follows:
      // 1,1,2,3,3,4,4,4,5,5,5,13,96

      if (xt_xmap_get_num_destinations(xmap) != 1)
        PUT_ERR("error in xt_xmap_get_num_destinations\n");

      if (xt_xmap_get_num_sources(xmap) != 1)
        PUT_ERR("error in xt_xmap_get_num_destinations\n");

      int rank;

      xt_xmap_get_destination_ranks(xmap, &rank);
      if (rank != my_rank)
        PUT_ERR("error in xt_xmap_get_destination_ranks\n");

      xt_xmap_get_source_ranks(xmap, &rank);
      if (rank != my_rank)
        PUT_ERR("error in xt_xmap_get_source_ranks\n");

      // clean up
      xt_xmap_delete(xmap);
      xt_idxlist_delete(src_idxlist);
      xt_idxlist_delete(dst_idxlist);
   }

   MPI_Finalize();

   return TEST_EXIT_CODE;
}

static void
check_idxdata(Xt_idxlist idxlist, Xt_idx * ref_indices,
              Xt_count ref_num_indices) {

   Xt_count i;

   if (xt_idxlist_get_num_indices(idxlist) != ref_num_indices)
      PUT_ERR("wrong number of indices in idxlist\n");

   Xt_idx indices[ref_num_indices];

   xt_idxlist_get_indices(idxlist, indices);

   for (i = 0; i < ref_num_indices; ++i)
     if (indices[i] != ref_indices[i])
       PUT_ERR("index list contains wrong indices\n");
}
