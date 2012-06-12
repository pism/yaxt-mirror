#include <stdlib.h>

#include "mpi.h"

#define VERBOSE

#include "tests.h"
#include "test_idxlist_utils.h"

#include "yaxt.h"
#include "xt_idxvec.h"
#include "xt_mpi.h"
#include "core/ppm_xfuncs.h"

void check_stripes(struct Xt_stripe const * stripes, unsigned num_stripes,
                   struct Xt_stripe const * ref_stripes, unsigned ref_num_stripes);

int main(void) {

   // init mpi

   mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  Xt_initialize(MPI_COMM_WORLD);

   { // test packing and unpacking of index vectors

      Xt_idx index_vector[] = {1,2,3,4,5,6,7};
      size_t num_indices = sizeof(index_vector) / sizeof(index_vector[0]);

      xt_idxlist idxvector;

      // create an index vector

      idxvector = xt_idxvec_new(index_vector, num_indices);

      // test the generated index vector

      check_idxlist(idxvector, index_vector, num_indices);

      // get the size of the index vector

      size_t buffer_size;

      buffer_size = xt_idxlist_get_pack_size(idxvector, MPI_COMM_WORLD);

      // allocate send buffers

      void* send_buffer, * recv_buffer;

      send_buffer = xmalloc(buffer_size);
      recv_buffer = xmalloc(buffer_size);

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

      check_idxlist(idxvector_copy, index_vector, num_indices);

      // compute intersection between the two index vectors

      xt_idxlist intersection;
      
      intersection = xt_idxlist_get_intersection(idxvector, idxvector_copy);

      // check the computed intersection (should be identically to the original
      // index vector)

      check_idxlist(intersection, index_vector, num_indices);

      // check the conversion to stripes

      struct Xt_stripe * stripes;
      unsigned num_stripes;

      xt_idxlist_get_index_stripes(idxvector, &stripes, &num_stripes);

      struct Xt_stripe ref_stripes[1] = {{1,7,1}};
      unsigned ref_num_stripes = 1;

      check_stripes(stripes, num_stripes, ref_stripes, ref_num_stripes);

      // clean up

      free(stripes);

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

      check_idxlist(idxvector, index_vector, num_indices);

      // copy the index vector

      idxvector_copy = xt_idxlist_copy(idxvector);

      // check the generated copy

      check_idxlist(idxvector_copy, index_vector, num_indices);

      // clean up

      xt_idxlist_delete(idxvector);
      xt_idxlist_delete(idxvector_copy);
   }

   MPI_Finalize();

   return TEST_EXIT_CODE;
}

void check_stripes(struct Xt_stripe const * stripes, unsigned num_stripes,
                   struct Xt_stripe const * ref_stripes, unsigned ref_num_stripes) {

  unsigned i;

  if (num_stripes != ref_num_stripes)
    PUT_ERR("wrong number of stripes\n");

  for (i = 0; i < ref_num_stripes; ++i)
    if (stripes[i].start != ref_stripes[i].start ||
        stripes[i].extent != ref_stripes[i].extent ||
        stripes[i].stride != ref_stripes[i].stride)
      PUT_ERR("error in stripe (start, extent and/or stride)\n");
}
