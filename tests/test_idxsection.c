#include <stdlib.h>

#include "mpi.h"

#include "tests.h"
#include "test_idxlist_utils.h"

#include "yaxt.h"
#include "xt_idxsection.h"
#include "xt_mpi.h"
#include "core/ppm_xfuncs.h"

void do_tests(xt_idxlist idxlist, Xt_idx * ref_indices, Xt_count num_indices,
              struct Xt_stripe * ref_stripes, Xt_count ref_num_stripes);
void check_stripes(struct Xt_stripe const * stripes, Xt_count num_stripes,
                   struct Xt_stripe const * ref_stripes, Xt_count ref_num_stripes);

int main(void) {

  // init mpi

  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // test 1D section

    Xt_idxlist idxsection;

    Xt_idx start = 0;
    unsigned num_dimensions = 1;
    unsigned global_size[1] = {10};
    unsigned local_size [1] = {5};
    unsigned local_start[1] = {3};

    // create index section

    idxsection = xt_idxsection_new(start, num_dimensions, global_size,
                                   local_size, local_start);

    // testing

    Xt_idx ref_indices[5] = {3,4,5,6,7};
    struct Xt_stripe ref_stripes[1] = {{3,5,1}};
    do_tests(idxsection, ref_indices, 5, ref_stripes, 1);

    // clean up

    xt_idxlist_delete(idxsection);
  }

  { // test 2D section

    Xt_idxlist idxsection;

    Xt_idx start = 0;
    unsigned num_dimensions = 2;
    unsigned global_size[2] = {5,6};
    unsigned local_size [2] = {3,2};
    unsigned local_start[2] = {1,2};

    // create index section

    idxsection = xt_idxsection_new(start, num_dimensions, global_size,
                                   local_size, local_start);

    // testing

    Xt_idx ref_indices[6] = {8,9,14,15,20,21};
    struct Xt_stripe ref_stripes[3] = {{8,2,1}, {14,2,1}, {20,2,1}};
    do_tests(idxsection, ref_indices, 6, ref_stripes, 3);

    // clean up

    xt_idxlist_delete(idxsection);
  }

  { // test 3D section

    Xt_idxlist idxsection;

    Xt_idx start = 0;
    unsigned num_dimensions = 3;
    unsigned global_size[3] = {4,4,4};
    unsigned local_size[3] = {4,2,2};
    unsigned local_start[3] = {0,1,1};

    // create index section

    idxsection = xt_idxsection_new(start, num_dimensions, global_size,
                                   local_size, local_start);

    // testing

    Xt_idx ref_indices[16] = {5,6,9,10, 21,22,25,26, 37,38,41,42, 53,54,57,58};
    struct Xt_stripe ref_stripes[8] = {{5,2,1}, {9,2,1}, {21,2,1},
                                       {25,2,1}, {37,2,1}, {41,2,1},
                                       {53,2,1}, {57,2,1}};
    do_tests(idxsection, ref_indices, 16, ref_stripes, 8);

    // clean up

    xt_idxlist_delete(idxsection);
  }

  { // 1D intersection test

    Xt_idxlist idxsection_a, idxsection_b;

    Xt_idx start = 0;
    unsigned num_dimensions = 1;
    unsigned global_size_a[1] = {10};
    unsigned local_size_a [1] = {5};
    unsigned local_start_a[1] = {4};

    unsigned global_size_b[1] = {15};
    unsigned local_size_b [1] = {6};
    unsigned local_start_b[1] = {7};

    // create index sections

    idxsection_a = xt_idxsection_new(start, num_dimensions, global_size_a,
                                     local_size_a, local_start_a);
    idxsection_b = xt_idxsection_new(start, num_dimensions, global_size_b,
                                     local_size_b, local_start_b);

    // compute intersection

    Xt_idxlist intersection;

    intersection = xt_idxlist_get_intersection(idxsection_a, idxsection_b);

    // check intersection

    Xt_idx ref_indices[2] = {7,8};
    struct Xt_stripe ref_stripes[1] = {{7,2,1}};
    
    do_tests(intersection, ref_indices, 2, ref_stripes, 1);

    xt_idxlist_delete(intersection);
    xt_idxlist_delete(idxsection_a);
    xt_idxlist_delete(idxsection_b);
  }

  { // 1D intersection test (with empty intersection)

    Xt_idxlist idxsection_a, idxsection_b;

    Xt_idx start = 0;
    unsigned num_dimensions = 1;
    unsigned global_size_a[1] = {10};
    unsigned local_size_a [1] = {1};
    unsigned local_start_a[1] = {3};

    unsigned global_size_b[1] = {10};
    unsigned local_size_b [1] = {5};
    unsigned local_start_b[1] = {4};

    // create index sections

    idxsection_a = xt_idxsection_new(start, num_dimensions, global_size_a,
                                     local_size_a, local_start_a);
    idxsection_b = xt_idxsection_new(start, num_dimensions, global_size_b,
                                     local_size_b, local_start_b);

    // compute intersection

    Xt_idxlist intersection;

    intersection = xt_idxlist_get_intersection(idxsection_a, idxsection_b);

    // check intersection

    Xt_idx ref_indices[1] = {0};
    struct Xt_stripe ref_stripes[1] = {{0,1,1}};
    
    do_tests(intersection, ref_indices, 0, ref_stripes, 0);

    xt_idxlist_delete(intersection);
    xt_idxlist_delete(idxsection_a);
    xt_idxlist_delete(idxsection_b);
  }

  { // 2D intersection test

    Xt_idxlist idxsection_a, idxsection_b;

    Xt_idx start = 0;
    unsigned num_dimensions = 2;
    unsigned global_size_a[2] = {6,6};
    unsigned local_size_a [2] = {4,2};
    unsigned local_start_a[2] = {1,1};

    unsigned global_size_b[2] = {6,6};
    unsigned local_size_b [2] = {3,3};
    unsigned local_start_b[2] = {3,2};

    // create index sections

    idxsection_a = xt_idxsection_new(start, num_dimensions, global_size_a,
                                     local_size_a, local_start_a);
    idxsection_b = xt_idxsection_new(start, num_dimensions, global_size_b,
                                     local_size_b, local_start_b);

    // compute intersection

    Xt_idxlist intersection;

    intersection = xt_idxlist_get_intersection(idxsection_a, idxsection_b);

    // check intersection

    Xt_idx ref_indices[2] = {20,26};
    struct Xt_stripe ref_stripes[2] = {{20,1,1}, {26,1,1}};
    
    do_tests(intersection, ref_indices, 2, ref_stripes, 2);

    xt_idxlist_delete(intersection);
    xt_idxlist_delete(idxsection_a);
    xt_idxlist_delete(idxsection_b);
  }

  { // 2D test

    Xt_idxlist idxsection;

    Xt_idx start = 0;
    unsigned num_dimensions = 2;
    unsigned global_size[2] = {4,4};
    unsigned local_size [2] = {2,2};
    unsigned local_start[2] = {0,2};

    idxsection = xt_idxsection_new(start, num_dimensions, global_size,
                                   local_size, local_start);

    Xt_idx indices[4], ref_indices[4] = {2,3,6,7};

    xt_idxlist_get_indices(idxsection, indices);

    for (int i = 0; i < 4; ++i)
      if (indices[i] != ref_indices[i])
        PUT_ERR("error in xt_idxlist_get_indices\n");
  }

  { // 2D test

    Xt_idxlist idxsection;

    Xt_idx start = 1;
    unsigned num_dimensions = 2;
    unsigned global_size[2] = {4,4};
    unsigned local_size [2] = {2,2};
    unsigned local_start[2] = {0,2};

    idxsection = xt_idxsection_new(start, num_dimensions, global_size,
                                   local_size, local_start);

    Xt_idx indices[4], ref_indices[4] = {3,4,7,8};

    xt_idxlist_get_indices(idxsection, indices);

    for (int i = 0; i < 4; ++i)
      if (indices[i] != ref_indices[i])
        PUT_ERR("error in xt_idxlist_get_indices\n");
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

void do_tests(xt_idxlist idxlist, Xt_idx * ref_indices, Xt_count num_indices,
              struct Xt_stripe * ref_stripes, Xt_count ref_num_stripes) {

  check_idxlist(idxlist, ref_indices, num_indices);

  struct Xt_stripe * stripes;
  Xt_count num_stripes;

  xt_idxlist_get_index_stripes(idxlist, &stripes, &num_stripes);

  check_stripes(stripes, num_stripes, ref_stripes, ref_num_stripes);

  free(stripes);

  { // test packing and unpacking

    size_t buffer_size;

    buffer_size = xt_idxlist_get_pack_size(idxlist, MPI_COMM_WORLD);

    // allocate send buffers

    void* send_buffer, * recv_buffer;

    send_buffer = xmalloc(buffer_size);
    recv_buffer = xmalloc(buffer_size);

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

    // check copy

    check_idxlist(idxlist_copy, ref_indices, num_indices);

    // clean up

    free(send_buffer);
    free(recv_buffer);
    xt_idxlist_delete(idxlist_copy);
  }

  { // test copying

    xt_idxlist idxlist_copy;

    idxlist_copy = xt_idxlist_copy(idxlist);

    // check copy

    check_idxlist(idxlist_copy, ref_indices, num_indices);

    // clean up

    xt_idxlist_delete(idxlist_copy);
  }
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
