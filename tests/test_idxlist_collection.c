#include <stdlib.h>

#include "mpi.h"

#define VERBOSE

#include "tests.h"

#include "yaxt.h"
#include "xt_idxvec.h"
#include "xt_idxlist_collection.h"
#include "test_idxlist_utils.h"
#include "xt_mpi.h"
#include "core/ppm_xfuncs.h"

void check_stripes(struct Xt_stripe const * stripes, Xt_count num_stripes,
                   struct Xt_stripe const * ref_stripes, Xt_count ref_num_stripes);

int main(void) {

   // init mpi

   mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

   { // test packing and unpacking of collection index lists

      Xt_idx index_list[2][7] = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1}};
      Xt_count num_indices = 7;

      Xt_idxlist idxlists[2];

      // create two index lists

      idxlists[0] = xt_idxvec_new(index_list[0], num_indices);
      idxlists[1] = xt_idxvec_new(index_list[1], num_indices);

      // generate a collection index list

      Xt_idxlist collectionlist;

      collectionlist = xt_idxlist_collection_new(idxlists, 2);

      xt_idxlist_delete(idxlists[0]);
      xt_idxlist_delete(idxlists[1]);

      // test generated collection list

      check_idxlist(collectionlist, &(index_list[0][0]), 2*num_indices);

      // get the size of the index list

      size_t buffer_size;

      buffer_size = xt_idxlist_get_pack_size(collectionlist, MPI_COMM_WORLD);

      // allocate send buffers

      void* send_buffer, * recv_buffer;

      send_buffer = xmalloc(buffer_size);
      recv_buffer = xmalloc(buffer_size);

      // pack the index list

      int position = 0;

      xt_idxlist_pack(collectionlist, send_buffer, buffer_size, &position, MPI_COMM_WORLD);

      // send the buffer

      int rank;

      mpi_err_handler(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);

      MPI_Request request;

      mpi_err_handler(MPI_Isend(send_buffer, buffer_size, MPI_PACKED, rank, 0,
                                MPI_COMM_WORLD, &request), MPI_COMM_WORLD);
      mpi_err_handler(MPI_Request_free(&request), MPI_COMM_WORLD);

      // receive the buffer

      mpi_err_handler(MPI_Recv(recv_buffer, buffer_size, MPI_PACKED, rank, 0,
                               MPI_COMM_WORLD, MPI_STATUS_IGNORE), MPI_COMM_WORLD);

      // unpack the buffer

      Xt_idxlist collectionlist_copy;

      position = 0;
      collectionlist_copy = xt_idxlist_unpack(recv_buffer, buffer_size, &position,
                                              MPI_COMM_WORLD);

      // check received collection list

      check_idxlist(collectionlist_copy, &(index_list[0][0]), 2*num_indices);

      // compute intersection between the two index lists

      Xt_idxlist intersection;

      intersection = xt_idxlist_get_intersection(collectionlist, collectionlist_copy);

      // check intersection data
      Xt_idx ref_index_data[14] = {1,1,2,2,3,3,4,4,5,5,6,6,7,7};
      check_idxlist(intersection, ref_index_data, 2*num_indices);

      // check the conversion to stripes

      struct Xt_stripe * stripes;
      Xt_count num_stripes;

      xt_idxlist_get_index_stripes(collectionlist, &stripes, &num_stripes);

      struct Xt_stripe ref_stripes[8] = {{1,7,1}, {7,1,1}, {6,1,1}, {5,1,1},
                                         {4,1,1}, {3,1,1}, {2,1,1}, {1,1,1}};
      Xt_count ref_num_stripes = 8;

      check_stripes(stripes, num_stripes, ref_stripes, ref_num_stripes);

      // clean up

      free(stripes);

      free(send_buffer);
      free(recv_buffer);

      xt_idxlist_delete(collectionlist);
      xt_idxlist_delete(collectionlist_copy);
      xt_idxlist_delete(intersection);
   }

   { // test copying of collection index lists

      Xt_idx index_list[2][7] = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1}};
      Xt_count num_indices = 7;

      Xt_idxlist idxlists[2];

      // create two index lists

      idxlists[0] = xt_idxvec_new(index_list[0], num_indices);
      idxlists[1] = xt_idxvec_new(index_list[1], num_indices);

      // generate a collection index list

      Xt_idxlist collectionlist, collectionlist_copy;

      collectionlist = xt_idxlist_collection_new(idxlists, 2);

      xt_idxlist_delete(idxlists[0]);
      xt_idxlist_delete(idxlists[1]);

      // test generated collection list

      check_idxlist(collectionlist, &(index_list[0][0]), 2*num_indices);

      // copy the index list

      collectionlist_copy = xt_idxlist_copy(collectionlist);

      // check copy of collection list

      check_idxlist(collectionlist_copy, &(index_list[0][0]), 2*num_indices);

      // check the conversion to stripes

      struct Xt_stripe * stripes;
      Xt_count num_stripes;

      xt_idxlist_get_index_stripes(collectionlist, &stripes, &num_stripes);

      struct Xt_stripe ref_stripes[8] = {{1,7,1}, {7,1,1}, {6,1,1}, {5,1,1},
                                         {4,1,1}, {3,1,1}, {2,1,1}, {1,1,1}};
      Xt_count ref_num_stripes = 8;

      check_stripes(stripes, num_stripes, ref_stripes, ref_num_stripes);

      // clean up

      free(stripes);

      xt_idxlist_delete(collectionlist);
      xt_idxlist_delete(collectionlist_copy);
   }

   MPI_Finalize();

   return TEST_EXIT_CODE;
}

void check_stripes(struct Xt_stripe const * stripes, Xt_count num_stripes,
                   struct Xt_stripe const * ref_stripes, Xt_count ref_num_stripes) {

  Xt_count i;

  if (num_stripes != ref_num_stripes)
    PUT_ERR("wrong number of stripes\n");

  for (i = 0; i < ref_num_stripes; ++i)
    if (stripes[i].start != ref_stripes[i].start ||
        stripes[i].nstrides != ref_stripes[i].nstrides ||
        stripes[i].stride != ref_stripes[i].stride)
      PUT_ERR("error in stripe (start, nstrides and/or stride)\n");
}
