#include <stdlib.h>

#include "mpi.h"

#define VERBOSE

#include "tests.h"

#include "yaxt.h"
#include "xt_idxlist_t.h"
#include "xt_mpi.h"

void check_idxlist_t(xt_idxlist idxlist, Xt_idx * ref_indices, unsigned ref_num_indices);

int main(void) {

   // init mpi

   mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

   { // test packing and unpacking of index lists

      Xt_idx index_list[] = {1,2,3,4,5,6,7};
      size_t num_indices = sizeof(index_list) / sizeof(index_list[0]);

      xt_idxlist idxlist;

      // create an index list

      idxlist = xt_idxlist_t_new(index_list, num_indices);

      // test the generated index vector

      check_idxlist_t(idxlist, index_list, num_indices);

      // get the size of the index list

      size_t buffer_size;

      buffer_size = xt_idxlist_get_pack_size(idxlist, MPI_COMM_WORLD);

      // allocate send buffers

      void* send_buffer, * recv_buffer;

      send_buffer = malloc(buffer_size);
      recv_buffer = malloc(buffer_size);

      // pack the index vector

      int position = 0;

      xt_idxlist_pack(idxlist, send_buffer, buffer_size, &position, MPI_COMM_WORLD);

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

      xt_idxlist idxlist_copy;

      position = 0;
      idxlist_copy = xt_idxlist_unpack(recv_buffer, buffer_size, &position, MPI_COMM_WORLD);

      // check the received index vector

      check_idxlist_t(idxlist_copy, index_list, num_indices);

      // compute intersection between the two index vectors

      xt_idxlist intersection;
      
      intersection = xt_idxlist_get_intersection(idxlist, idxlist_copy);

      // check the computed intersection (should be identically to the original

      check_idxlist_t(intersection, index_list, num_indices);

      // index vector)

      // clean up

      free(send_buffer);
      free(recv_buffer);

      xt_idxlist_delete(idxlist);
      xt_idxlist_delete(idxlist_copy);
      xt_idxlist_delete(intersection);
   }

   { // test copying of index list

      Xt_idx index_list[] = {1,2,3,4,5,6,7};
      size_t num_indices = sizeof(index_list) / sizeof(index_list[0]);

      xt_idxlist idxlist, idxlist_copy;

      // create an index list

      idxlist = xt_idxlist_t_new(index_list, num_indices);

      // test the generated index vector

      check_idxlist_t(idxlist, index_list, num_indices);

      // copy the index list

      idxlist_copy = xt_idxlist_copy(idxlist);

      // check the generated copy

      check_idxlist_t(idxlist_copy, index_list, num_indices);

      // clean up

      xt_idxlist_delete(idxlist);
      xt_idxlist_delete(idxlist_copy);
   }

   MPI_Finalize();

   return TEST_EXIT_CODE;
}

void check_idxlist_t(xt_idxlist idxlist, Xt_idx * ref_indices,
                     unsigned ref_num_indices) {

   unsigned i;

   if (xt_idxlist_get_num_indices(idxlist) != ref_num_indices)
      PUT_ERR("wrong number of indices in idxlist\n");

   Xt_idx indices[ref_num_indices];

   xt_idxlist_get_indices(idxlist, indices);

   for (i = 0; i < ref_num_indices; ++i)
      if (indices[i] != ref_indices[i])
         PUT_ERR("index vector contains wrong indices\n");
}
