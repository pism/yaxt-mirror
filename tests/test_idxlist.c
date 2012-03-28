#include <stdlib.h>

#include "mpi.h"

#define VERBOSE

#include "tests.h"

#include "yaxt.h"
#include "xt_idxlist.h"
#include "xt_mpi.h"

int main(void) {

   // init mpi

   mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

   { // test packing and unpacking of index lists

      int index_list[] = {1,2,3,4,5,6,7};
      unsigned num_indices = sizeof(index_list) / sizeof(index_list[0]);

      xt_idxvec idxvec;

      // create an index list

      idxvec = xt_idxlist_new(index_list, num_indices);

      // get the size of the index list

      size_t buffer_size;

      buffer_size = xt_idxvec_get_pack_size(idxvec, MPI_COMM_WORLD);

      // allocate send buffers

      void* send_buffer, * recv_buffer;

      send_buffer = malloc(buffer_size);
      recv_buffer = malloc(buffer_size);

      // pack the index vector

      int position = 0;

      xt_idxvec_pack(idxvec, send_buffer, buffer_size, &position, MPI_COMM_WORLD);

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

      xt_idxvec idxvec_copy;

      position = 0;
      idxvec_copy = xt_idxvec_unpack(recv_buffer, buffer_size, &position, MPI_COMM_WORLD);

      // compute intersection between the two index vectors

      xt_idxvec intersection;
      
      intersection = xt_idxvec_get_intersection(idxvec, idxvec_copy);

      // clean up

      free(send_buffer);
      free(recv_buffer);

      xt_idxvec_delete(idxvec);
      xt_idxvec_delete(idxvec_copy);
      xt_idxvec_delete(intersection);
   }

   { // test copying of index list

      int index_list[] = {1,2,3,4,5,6,7};
      unsigned num_indices = sizeof(index_list) / sizeof(index_list[0]);

      xt_idxvec idxvec, idxvec_copy;

      // create an index list

      idxvec = xt_idxlist_new(index_list, num_indices);

      // copy the index list

      idxvec_copy = xt_idxvec_copy(idxvec);

      // clean up

      xt_idxvec_delete(idxvec);
      xt_idxvec_delete(idxvec_copy);
   }

   { // test get_num_indices and get_indices

      int index_list[] = {1,2,3,4,5,6,7};
      unsigned num_indices = sizeof(index_list) / sizeof(index_list[0]);

      xt_idxvec idxvec;

      // create an index list

      idxvec = xt_idxlist_new(index_list, num_indices);

      // testing

      if (xt_idxvec_get_num_indices(idxvec) != num_indices)
         PUT_ERR("error in xt_idxvec_get_num_indices\n");

      unsigned i;
      int indices[num_indices];

      xt_idxvec_get_indices(idxvec, indices);

      for (i = 0; i < num_indices; ++i)
         if (indices[i] != index_list[i]) {
            PUT_ERR("error in xt_idxvec_get_indices\n");
            break;
         }
   }

   MPI_Finalize();

   return TEST_EXIT_CODE;
}
