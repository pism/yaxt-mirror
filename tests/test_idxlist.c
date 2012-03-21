#include <stdlib.h>

#include "mpi.h"

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

      xt_idxset idxset;

      // create an index list

      idxset = xt_idxlist_new(index_list, num_indices);

      // get the size of the index list

      size_t buffer_size;

      buffer_size = xt_idxset_get_pack_size(idxset, MPI_COMM_WORLD);

      // allocate send buffers

      void* send_buffer, * recv_buffer;

      send_buffer = malloc(buffer_size);
      recv_buffer = malloc(buffer_size);

      // pack the index set

      int position = 0;

      xt_idxset_pack(idxset, send_buffer, buffer_size, &position, MPI_COMM_WORLD);

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

      xt_idxset idxset_copy;

      position = 0;
      idxset_copy = xt_idxset_unpack(recv_buffer, buffer_size, &position, MPI_COMM_WORLD);

      // compute intersection between the two index sets

      xt_idxset intersection;
      
      intersection = xt_idxset_get_intersection(idxset, idxset_copy);

      // clean up

      free(send_buffer);
      free(recv_buffer);

      xt_idxset_delete(idxset);
      xt_idxset_delete(idxset_copy);
      xt_idxset_delete(intersection);
   }

   { // test copying of index list

      int index_list[] = {1,2,3,4,5,6,7};
      unsigned num_indices = sizeof(index_list) / sizeof(index_list[0]);

      xt_idxset idxset, idxset_copy;

      // create an index list

      idxset = xt_idxlist_new(index_list, num_indices);

      // copy the index list

      idxset_copy = xt_idxset_copy(idxset);

      // clean up

      xt_idxset_delete(idxset);
      xt_idxset_delete(idxset_copy);
   }

   MPI_Finalize();

   return TEST_EXIT_CODE;
}
