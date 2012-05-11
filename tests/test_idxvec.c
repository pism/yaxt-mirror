#include <stdlib.h>

#include "mpi.h"

#define VERBOSE

#include "tests.h"

#include "yaxt.h"
#include "xt_idxvec.h"
#include "xt_mpi.h"

void check_idxvec(xt_idxlist idxlist, Xt_idx * ref_indices, unsigned ref_num_indices);

int main(void) {

   // init mpi

   mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

   { // test packing and unpacking of index vectors

      Xt_idx index_vector[] = {1,2,3,4,5,6,7};
      size_t num_indices = sizeof(index_vector) / sizeof(index_vector[0]);

      xt_idxlist idxvector;

      // create an index vector

      idxvector = xt_idxvec_new(index_vector, num_indices);

      // test the generated index vector

      check_idxvec(idxvector, index_vector, num_indices);

      // get the size of the index vector

      size_t buffer_size;

      buffer_size = xt_idxlist_get_pack_size(idxvector, MPI_COMM_WORLD);

      // allocate send buffers

      void* send_buffer, * recv_buffer;

      send_buffer = malloc(buffer_size);
      recv_buffer = malloc(buffer_size);

      // pack the index vector

      int position = 0;

      xt_idxlist_pack(idxvector, send_buffer, buffer_size, &position, MPI_COMM_WORLD);

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

      xt_idxlist idxvector_copy;

      position = 0;
      idxvector_copy = xt_idxlist_unpack(recv_buffer, buffer_size, &position, MPI_COMM_WORLD);

      // check the received index vector

      check_idxvec(idxvector_copy, index_vector, num_indices);

      // compute intersection between the two index vectors

      xt_idxlist intersection;
      
      intersection = xt_idxlist_get_intersection(idxvector, idxvector_copy);

      // check the computed intersection (should be identically to the original
      // index vector)

      check_idxvec(intersection, index_vector, num_indices);

      // clean up

      free(send_buffer);
      free(recv_buffer);

      xt_idxlist_delete(idxvector);
      xt_idxlist_delete(idxvector_copy);
      xt_idxlist_delete(intersection);
   }

   { // test copying of index vector

      Xt_idx index_vector[] = {1,2,3,4,5,6,7};
      size_t num_indices = sizeof(index_vector) / sizeof(index_vector[0]);

      xt_idxlist idxvector, idxvector_copy;

      // create an index vector

      idxvector = xt_idxvec_new(index_vector, num_indices);

      // test the generated index vector

      check_idxvec(idxvector, index_vector, num_indices);

      // copy the index vector

      idxvector_copy = xt_idxlist_copy(idxvector);

      // check the generated copy

      check_idxvec(idxvector_copy, index_vector, num_indices);

      // clean up

      xt_idxlist_delete(idxvector);
      xt_idxlist_delete(idxvector_copy);
   }

   MPI_Finalize();

   return TEST_EXIT_CODE;
}

void check_idxvec(xt_idxlist idxlist, Xt_idx * ref_indices,
                  unsigned ref_num_indices) {

   unsigned i;

   if (xt_idxlist_get_num_indices(idxlist) != ref_num_indices)
      PUT_ERR("wrong number of indices in idxlist\n");

   Xt_idx indices[ref_num_indices];

   xt_idxlist_get_indices(idxlist, indices);

   for (i = 0; i < ref_num_indices; ++i)
      if (indices[i] != ref_indices[i])
         PUT_ERR("index list contains wrong indices\n");
}
