#include <stdlib.h>

#include "mpi.h"

#define VERBOSE

#include "tests.h"

#include "yaxt.h"
#include "xt_idxlist.h"
#include "xt_idxset_collection.h"
#include "xt_mpi.h"

int main(void) {

   // init mpi

   mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

   { // test packing and unpacking of collection index sets

      int index_list[2][7] = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1}};
      unsigned num_indices = 7;

      xt_idxset idxsets[2];

      // create two index lists

      idxsets[0] = xt_idxlist_new(index_list[0], num_indices);
      idxsets[1] = xt_idxlist_new(index_list[1], num_indices);

      // generate a collection index set

      xt_idxset collectionset;

      collectionset = xt_idxset_collection_new(idxsets, 2);

      // get the size of the index list

      size_t buffer_size;

      buffer_size = xt_idxset_get_pack_size(collectionset, MPI_COMM_WORLD);

      // allocate send buffers

      void* send_buffer, * recv_buffer;

      send_buffer = malloc(buffer_size);
      recv_buffer = malloc(buffer_size);

      // pack the index set

      int position = 0;

      xt_idxset_pack(collectionset, send_buffer, buffer_size, &position, MPI_COMM_WORLD);

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

      xt_idxset collectionset_copy;

      position = 0;
      collectionset_copy = xt_idxset_unpack(recv_buffer, buffer_size, &position,
                                            MPI_COMM_WORLD);

      // compute intersection between the two index sets

      xt_idxset intersection;
      
      intersection = xt_idxset_get_intersection(collectionset, collectionset_copy);

      // clean up

      free(send_buffer);
      free(recv_buffer);

      xt_idxset_delete(collectionset);
      xt_idxset_delete(collectionset_copy);
      xt_idxset_delete(intersection);
   }

   { // test copying of collection index sets

      int index_list[2][7] = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1}};
      unsigned num_indices = 7;

      xt_idxset idxsets[2];

      // create two index lists

      idxsets[0] = xt_idxlist_new(index_list[0], num_indices);
      idxsets[1] = xt_idxlist_new(index_list[1], num_indices);

      // generate a collection index set

      xt_idxset collectionset, collectionset_copy;

      collectionset = xt_idxset_collection_new(idxsets, 2);

      // copy the index list

      collectionset_copy = xt_idxset_copy(collectionset);

      // clean up

      xt_idxset_delete(collectionset);
      xt_idxset_delete(collectionset_copy);
   }

   { // test get_num_indices and get_indices

      int index_list[2][4] = {{1,2,3,4},{5,6,7,8}};
      unsigned num_indices = 4;

      xt_idxset idxsets[2];

      // create two index lists

      idxsets[0] = xt_idxlist_new(index_list[0], num_indices);
      idxsets[1] = xt_idxlist_new(index_list[1], num_indices);

      // generate a collection index set

      xt_idxset collectionset, collectionset_copy;

      collectionset = xt_idxset_collection_new(idxsets, 2);

      // testing

      if (xt_idxset_get_num_indices(collectionset) != 2 * num_indices)
         PUT_ERR("error in xt_idxset_get_num_indices\n");

      unsigned i;
      int indices[2*num_indices];

      xt_idxset_get_indices(collectionset, indices);

      for (i = 0; i < num_indices; ++i)
         if (indices[i] != index_list[i/4][i&3]) {
            PUT_ERR("error in xt_idxset_get_indices\n");
            break;
         }
   }

   MPI_Finalize();

   return TEST_EXIT_CODE;
}
