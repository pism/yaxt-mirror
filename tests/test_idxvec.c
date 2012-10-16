#include <stdlib.h>
#include <assert.h>

#include "mpi.h"

#define VERBOSE

#include "tests.h"
#include "test_idxlist_utils.h"

#include "yaxt.h"
#include "xt_idxvec.h"
#include "xt_mpi.h"
#include "core/ppm_xfuncs.h"

void check_stripes(struct Xt_stripe const * stripes, Xt_count num_stripes,
                   struct Xt_stripe const * ref_stripes, Xt_count ref_num_stripes);

int main(void) {

   // init mpi

   mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

   { // test packing and unpacking of index vectors

      Xt_idx index_vector[] = {1,2,3,4,5,6,7};
      Xt_count num_indices = sizeof(index_vector) / sizeof(index_vector[0]);

      Xt_idxlist idxvector;

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

      Xt_idxlist idxvector_copy;

      position = 0;
      idxvector_copy = xt_idxlist_unpack(recv_buffer, buffer_size, &position, MPI_COMM_WORLD);

      // check the received index vector

      check_idxlist(idxvector_copy, index_vector, num_indices);

      // compute intersection between the two index vectors

      Xt_idxlist intersection;

      intersection = xt_idxlist_get_intersection(idxvector, idxvector_copy);

      // check the computed intersection (should be identically to the original
      // index vector)

      check_idxlist(intersection, index_vector, num_indices);

      // check the conversion to stripes

      struct Xt_stripe * stripes;
      Xt_count num_stripes;

      xt_idxlist_get_index_stripes(idxvector, &stripes, &num_stripes);

      struct Xt_stripe ref_stripes[1] = {{1,7,1}};
      Xt_count ref_num_stripes = 1;

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
      Xt_count num_indices = sizeof(index_vector) / sizeof(index_vector[0]);

      Xt_idxlist idxvector, idxvector_copy;

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

   { // test repeated equal indices

      Xt_idx index_vector[] = {1,2,3,7,5,6,7,7};
      Xt_count num_indices = sizeof(index_vector) / sizeof(index_vector[0]);

      Xt_idxlist idxvector;

      // create an index vector

      idxvector = xt_idxvec_new(index_vector, num_indices);

      // test the generated index vector

      check_idxlist(idxvector, index_vector, num_indices);

      // clean up
      xt_idxlist_delete(idxvector);
   }

   { // subtest from test_ut

     const int single_match_only = 1;
     Xt_idx index_vector[] = {10,15,14,13,12,15,10,11,12,13,23,18,19,20,21,31,26,27,28,29};
     size_t index_num = sizeof(index_vector) / sizeof(index_vector[0]);

     Xt_idx intersection_vector[] = {12,12,13,13,14,15,15,20,21,23,28,29,31};
     size_t intersection_num = sizeof(intersection_vector) / sizeof(intersection_vector[0]);

     Xt_idxlist idxvector;

     // create an index vector

     idxvector = xt_idxvec_new(index_vector, index_num);

     // get positions:
     int retval;
     int intersection_pos[intersection_num];
     int ref_intersection_pos[] = {4,8,3,9,2,1,5,13,14,10,18,19,15};
     size_t ref_intersection_num = sizeof(ref_intersection_pos) / sizeof(ref_intersection_pos[0]);

     assert(intersection_num == ref_intersection_num);

     retval = xt_idxlist_get_position_of_indices(idxvector, intersection_vector,
                                                 intersection_num, intersection_pos, single_match_only);
     assert(retval != 1);

     // check positions:
     for(int i=0; i<intersection_num; i++)
       if (intersection_pos[i] != ref_intersection_pos[i])
         PUT_ERR("unexpected results from xt_idxlist_get_position_of_indices\n");

     // clean up
     xt_idxlist_delete(idxvector);
   }

   { // check intersection
      Xt_idx index_vector[2][3] = {{1,2,3},{1,2,3}};
      Xt_count num_indices[2] = {sizeof(index_vector[0]) / sizeof(index_vector[0][0]),
                                 sizeof(index_vector[0]) / sizeof(index_vector[0][0])};

      Xt_idxlist idxvector[2];

      // create index vectors

      idxvector[0] = xt_idxvec_new(index_vector[0], num_indices[0]);
      idxvector[1] = xt_idxvec_new(index_vector[1], num_indices[1]);

      // compute intersection between the two index vectors

      Xt_idxlist intersection;

      intersection = xt_idxlist_get_intersection(idxvector[0], idxvector[1]);

      // check the computed intersection

      Xt_idx ref_intersection_indices[] = {1,2,3};
      Xt_count ref_intersection_num_indices =
         sizeof(ref_intersection_indices) / sizeof(ref_intersection_indices[0]);

      check_idxlist(intersection, ref_intersection_indices,
                    ref_intersection_num_indices);

     // clean up

     xt_idxlist_delete(idxvector[0]);
     xt_idxlist_delete(idxvector[1]);
     xt_idxlist_delete(intersection);
   }

   { // check intersection
      Xt_idx index_vector[2][3] = {{1,2,3},{2,3,4}};
      Xt_count num_indices[2] = {sizeof(index_vector[0]) / sizeof(index_vector[0][0]),
                                 sizeof(index_vector[0]) / sizeof(index_vector[0][0])};

      Xt_idxlist idxvector[2];

      // create index vectors

      idxvector[0] = xt_idxvec_new(index_vector[0], num_indices[0]);
      idxvector[1] = xt_idxvec_new(index_vector[1], num_indices[1]);

      // compute intersection between the two index vectors

      Xt_idxlist intersection;

      intersection = xt_idxlist_get_intersection(idxvector[0], idxvector[1]);

      // check the computed intersection

      Xt_idx ref_intersection_indices[] = {2,3};
      Xt_count ref_intersection_num_indices =
         sizeof(ref_intersection_indices) / sizeof(ref_intersection_indices[0]);

      check_idxlist(intersection, ref_intersection_indices,
                    ref_intersection_num_indices);

     // clean up

     xt_idxlist_delete(idxvector[0]);
     xt_idxlist_delete(idxvector[1]);
     xt_idxlist_delete(intersection);
   }

   { // check intersection
      Xt_idx index_vector[2][3] = {{2,3,4},{1,2,3}};
      Xt_count num_indices[2] = {sizeof(index_vector[0]) / sizeof(index_vector[0][0]),
                                 sizeof(index_vector[0]) / sizeof(index_vector[0][0])};

      Xt_idxlist idxvector[2];

      // create index vectors

      idxvector[0] = xt_idxvec_new(index_vector[0], num_indices[0]);
      idxvector[1] = xt_idxvec_new(index_vector[1], num_indices[1]);

      // compute intersection between the two index vectors

      Xt_idxlist intersection;

      intersection = xt_idxlist_get_intersection(idxvector[0], idxvector[1]);

      // check the computed intersection

      Xt_idx ref_intersection_indices[] = {2,3};
      Xt_count ref_intersection_num_indices =
         sizeof(ref_intersection_indices) / sizeof(ref_intersection_indices[0]);

      check_idxlist(intersection, ref_intersection_indices,
                    ref_intersection_num_indices);

     // clean up

     xt_idxlist_delete(idxvector[0]);
     xt_idxlist_delete(idxvector[1]);
     xt_idxlist_delete(intersection);
   }

   { // check intersection
      Xt_idx index_vector[2][3] = {{4,2,3},{3,1,2}};
      Xt_count num_indices[2] = {sizeof(index_vector[0]) / sizeof(index_vector[0][0]),
                                 sizeof(index_vector[0]) / sizeof(index_vector[0][0])};

      Xt_idxlist idxvector[2];

      // create index vectors

      idxvector[0] = xt_idxvec_new(index_vector[0], num_indices[0]);
      idxvector[1] = xt_idxvec_new(index_vector[1], num_indices[1]);

      // compute intersection between the two index vectors

      Xt_idxlist intersection;

      intersection = xt_idxlist_get_intersection(idxvector[0], idxvector[1]);

      // check the computed intersection

      Xt_idx ref_intersection_indices[] = {2,3};
      Xt_count ref_intersection_num_indices =
         sizeof(ref_intersection_indices) / sizeof(ref_intersection_indices[0]);

      check_idxlist(intersection, ref_intersection_indices,
                    ref_intersection_num_indices);

     // clean up

     xt_idxlist_delete(idxvector[0]);
     xt_idxlist_delete(idxvector[1]);
     xt_idxlist_delete(intersection);
   }

   { // check intersection
      Xt_idx index_vector[2][3] = {{3,1,2},{4,2,3}};
      Xt_count num_indices[2] = {sizeof(index_vector[0]) / sizeof(index_vector[0][0]),
                                 sizeof(index_vector[0]) / sizeof(index_vector[0][0])};

      Xt_idxlist idxvector[2];

      // create index vectors

      idxvector[0] = xt_idxvec_new(index_vector[0], num_indices[0]);
      idxvector[1] = xt_idxvec_new(index_vector[1], num_indices[1]);

      // compute intersection between the two index vectors

      Xt_idxlist intersection;

      intersection = xt_idxlist_get_intersection(idxvector[0], idxvector[1]);

      // check the computed intersection

      Xt_idx ref_intersection_indices[] = {2,3};
      Xt_count ref_intersection_num_indices =
         sizeof(ref_intersection_indices) / sizeof(ref_intersection_indices[0]);

      check_idxlist(intersection, ref_intersection_indices,
                    ref_intersection_num_indices);

     // clean up

     xt_idxlist_delete(idxvector[0]);
     xt_idxlist_delete(idxvector[1]);
     xt_idxlist_delete(intersection);
   }

   { // test constructor xt_idxvec_from_stripes_new

      struct Xt_stripe stripes[2] = {{5,5,1}, {4,5,-1}};

      // construct

      Xt_idxlist idxvec;

      idxvec = xt_idxvec_from_stripes_new(stripes, 2);

      // test results

      Xt_idx ref_indices[] = {5,6,7,8,9,4,3,2,1,0};
      check_idxlist(idxvec, ref_indices, 10);

      // clean up

      xt_idxlist_delete(idxvec);
   }

   { // test constructor xt_idxvec_from_stripes_new

      struct Xt_stripe stripes[3] = {{0,5,1}, {2,5,1},{4,5,1}};

      // construct

      Xt_idxlist idxvec;

      idxvec = xt_idxvec_from_stripes_new(stripes, 3);

      // test results

      Xt_idx ref_indices[] = {0,1,2,3,4,2,3,4,5,6,4,5,6,7,8};
      check_idxlist(idxvec, ref_indices, 15);

      // clean up

      xt_idxlist_delete(idxvec);
   }

   { // test constructor xt_idxvec_from_stripes_new

      struct Xt_stripe stripes[3] = {{2,5,1},{0,5,1},{4,5,1}};

      // construct

      Xt_idxlist idxvec;

      idxvec = xt_idxvec_from_stripes_new(stripes, 3);

      // test results

      Xt_idx ref_indices[] = {2,3,4,5,6,0,1,2,3,4,4,5,6,7,8};
      check_idxlist(idxvec, ref_indices, 15);

      // clean up

      xt_idxlist_delete(idxvec);
   }

   { // test constructor xt_idxvec_from_stripes_new

      struct Xt_stripe stripes[3] = {{2,5,1},{4,5,-1},{4,5,1}};

      // construct

      Xt_idxlist idxvec;

      idxvec = xt_idxvec_from_stripes_new(stripes, 3);

      // test results

      Xt_idx ref_indices[] = {2,3,4,5,6,4,3,2,1,0,4,5,6,7,8};
      check_idxlist(idxvec, ref_indices, 15);

      // clean up

      xt_idxlist_delete(idxvec);
   }

   { // test constructor xt_idxvec_from_stripes_new

      struct Xt_stripe stripes[3] = {{0,5,3},{1,5,3},{14,5,-3}};

      // construct

      Xt_idxlist idxvec;

      idxvec = xt_idxvec_from_stripes_new(stripes, 3);

      // test results

      Xt_idx ref_indices[] = {0,3,6,9,12,1,4,7,10,13,14,11,8,5,2};
      check_idxlist(idxvec, ref_indices, 15);

      // clean up

      xt_idxlist_delete(idxvec);
   }

   { // test constructor xt_idxvec_from_stripes_new

      struct Xt_stripe stripes[3] = {{0,5,3},{2,5,3},{14,5,-3}};

      // construct

      Xt_idxlist idxvec;

      idxvec = xt_idxvec_from_stripes_new(stripes, 3);

      // test results

      Xt_idx ref_indices[] = {0,3,6,9,12,2,5,8,11,14,14,11,8,5,2};
      check_idxlist(idxvec, ref_indices, 15);

      // clean up

      xt_idxlist_delete(idxvec);
   }

   { // test constructor xt_idxvec_from_stripes_new

      struct Xt_stripe stripes[4] = {{0,5,-1},{1,5,1},{-5,5,-1},{6,5,1}};

      // construct

      Xt_idxlist idxvec;

      idxvec = xt_idxvec_from_stripes_new(stripes, 4);

      // test results

      Xt_idx ref_indices[] = {0,-1,-2,-3,-4,1,2,3,4,5,-5,-6,-7,-8,-9,6,7,8,9,10};
      check_idxlist(idxvec, ref_indices, 20);

      // clean up

      xt_idxlist_delete(idxvec);
   }

   {
     // check idxvec_get_indices_at_positions
     // case: mixed valid and invalid positions
     const Xt_idx undef_idx = -1234567;
     Xt_idx indices[] = {0,3,6,9,12,1,4,7,10,13,14,11,8,5,2,1};
     int num_indices = sizeof(indices) / sizeof(indices[0]);
     Xt_idxlist idxvec = xt_idxvec_new(indices, num_indices);
     int pos[] = {0,2,7,9,11,100,11,200,9,300,7,400,5};
     int num_pos = sizeof(pos) / sizeof(pos[0]);
     Xt_idx ref_sel_idx[num_pos];
     int ref_undef_count = 0;
     for (int ip=0; ip<num_pos; ip++) {
       int p = pos[ip];
       if (xt_idxlist_get_index_at_position(idxvec, p, &ref_sel_idx[ip]) != 0) {
         ref_sel_idx[ip] = undef_idx;
         ref_undef_count++;
       }
     }
     Xt_idx sel_idx[num_pos];
     int undef_count = xt_idxlist_get_indices_at_positions(idxvec, pos, num_pos, sel_idx, undef_idx);
     if (undef_count != ref_undef_count) PUT_ERR("test_idxvec.c: (undef_count != ref_undef_count)\n");
     for (int ip=0; ip<num_pos; ip++) {
       if (sel_idx[ip] != ref_sel_idx[ip]) PUT_ERR("test_idxvec.c: (sel_idx[ip] != ref_sel_idx[ip])\n");
     }
     xt_idxlist_delete(idxvec);
   }

   {
     // check idxvec_get_indices_at_positions
     // case: only valid positions
     const Xt_idx undef_idx = -1234567;
     Xt_idx indices[] = {0,3,6,9,12,1,4,7,10,13,14,11,8,5,2,1};
     int num_indices = sizeof(indices) / sizeof(indices[0]);
     Xt_idxlist idxvec = xt_idxvec_new(indices, num_indices);
     int pos[] = {0,2,7,9,11,11,9,7,5};
     int num_pos = sizeof(pos) / sizeof(pos[0]);
     Xt_idx ref_sel_idx[num_pos];
     int ref_undef_count = 0;
     for (int ip=0; ip<num_pos; ip++) {
       int p = pos[ip];
       if (xt_idxlist_get_index_at_position(idxvec, p, &ref_sel_idx[ip]) != 0) {
         ref_sel_idx[ip] = undef_idx;
         ref_undef_count++;
       }
     }
     Xt_idx sel_idx[num_pos];
     int undef_count = xt_idxlist_get_indices_at_positions(idxvec, pos, num_pos, sel_idx, undef_idx);
     if (undef_count != ref_undef_count) PUT_ERR("test_idxvec.c: (undef_count != ref_undef_count)\n");
     for (int ip=0; ip<num_pos; ip++) {
       if (sel_idx[ip] != ref_sel_idx[ip]) PUT_ERR("test_idxvec.c: (sel_idx[ip] != ref_sel_idx[ip])\n");
     }
     xt_idxlist_delete(idxvec);
   }

   {
     // check idxvec_get_indices_at_positions
     // case: only invalid positions
     const Xt_idx undef_idx = -1234567;
     Xt_idx indices[] = {0,3,6,9,12,1,4,7,10,13,14,11,8,5,2,1};
     int num_indices = sizeof(indices) / sizeof(indices[0]);
     Xt_idxlist idxvec = xt_idxvec_new(indices, num_indices);
     int pos[] = {100,102,107,109,1011,1011,109,107,105};
     int num_pos = sizeof(pos) / sizeof(pos[0]);
     Xt_idx ref_sel_idx[num_pos];
     int ref_undef_count = 0;
     for (int ip=0; ip<num_pos; ip++) {
       int p = pos[ip];
       if (xt_idxlist_get_index_at_position(idxvec, p, &ref_sel_idx[ip]) != 0) {
         ref_sel_idx[ip] = undef_idx;
         ref_undef_count++;
       }
     }
     Xt_idx sel_idx[num_pos];
     int undef_count = xt_idxlist_get_indices_at_positions(idxvec, pos, num_pos, sel_idx, undef_idx);
     if (undef_count != ref_undef_count) PUT_ERR("test_idxvec.c: (undef_count != ref_undef_count)\n");
     for (int ip=0; ip<num_pos; ip++) {
       if (sel_idx[ip] != ref_sel_idx[ip]) PUT_ERR("test_idxvec.c: (sel_idx[ip] != ref_sel_idx[ip])\n");
     }
     xt_idxlist_delete(idxvec);
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
