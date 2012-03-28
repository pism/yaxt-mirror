#include <stdlib.h>

#include "mpi.h"

#define VERBOSE

#include "tests.h"

#include "yaxt.h"
#include "xt_idxlist.h"
#include "xt_idxvec_collection.h"
#include "xt_mpi.h"

int main(void) {

   // init mpi

   mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

   { // test packing and unpacking of collection index vectors

      Xt_idx index_list[2][7] = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1}};
      size_t num_indices = 7;

      xt_idxvec idxvecs[2];

      // create two index lists

      idxvecs[0] = xt_idxlist_new(index_list[0], num_indices);
      idxvecs[1] = xt_idxlist_new(index_list[1], num_indices);

      // generate a collection index vector

      xt_idxvec collectionvector;

      collectionvector = xt_idxvec_collection_new(idxvecs, 2);

      // get the size of the index list

      size_t buffer_size;

      buffer_size = xt_idxvec_get_pack_size(collectionvector, MPI_COMM_WORLD);

      // allocate send buffers

      void* send_buffer, * recv_buffer;

      send_buffer = malloc(buffer_size);
      recv_buffer = malloc(buffer_size);

      // pack the index vector

      int position = 0;

      xt_idxvec_pack(collectionvector, send_buffer, buffer_size, &position, MPI_COMM_WORLD);

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

      xt_idxvec collectionvector_copy;

      position = 0;
      collectionvector_copy = xt_idxvec_unpack(recv_buffer, buffer_size, &position,
                                            MPI_COMM_WORLD);

      // compute intersection between the two index vectors

      xt_idxvec intersection;
      
      intersection = xt_idxvec_get_intersection(collectionvector, collectionvector_copy);

      // clean up

      free(send_buffer);
      free(recv_buffer);

      xt_idxvec_delete(collectionvector);
      xt_idxvec_delete(collectionvector_copy);
      xt_idxvec_delete(intersection);
   }

   { // test copying of collection index vectors

      Xt_idx index_list[2][7] = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1}};
      size_t num_indices = 7;

      xt_idxvec idxvecs[2];

      // create two index lists

      idxvecs[0] = xt_idxlist_new(index_list[0], num_indices);
      idxvecs[1] = xt_idxlist_new(index_list[1], num_indices);

      // generate a collection index vector

      xt_idxvec collectionvector, collectionvector_copy;

      collectionvector = xt_idxvec_collection_new(idxvecs, 2);

      // copy the index list

      collectionvector_copy = xt_idxvec_copy(collectionvector);

      // clean up

      xt_idxvec_delete(collectionvector);
      xt_idxvec_delete(collectionvector_copy);
   }

   { // test get_num_indices and get_indices

      Xt_idx index_list[2][4] = {{1,2,3,4},{5,6,7,8}};
      size_t num_indices = 4;

      xt_idxvec idxvecs[2];

      // create two index lists

      idxvecs[0] = xt_idxlist_new(index_list[0], num_indices);
      idxvecs[1] = xt_idxlist_new(index_list[1], num_indices);

      // generate a collection index vector

      xt_idxvec collectionvector, collectionvector_copy;

      collectionvector = xt_idxvec_collection_new(idxvecs, 2);

      // testing

      if (xt_idxvec_get_num_indices(collectionvector) != 2 * num_indices)
         PUT_ERR("error in xt_idxvec_get_num_indices\n");

      size_t i;
      Xt_idx indices[2*num_indices];

      xt_idxvec_get_indices(collectionvector, indices);

      for (i = 0; i < num_indices; ++i)
         if (indices[i] != index_list[i/4][i&3]) {
            PUT_ERR("error in xt_idxvec_get_indices\n");
            break;
         }
   }

   MPI_Finalize();

   return TEST_EXIT_CODE;
}
