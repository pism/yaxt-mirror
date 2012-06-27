#include <stdlib.h>

#include "mpi.h"

#include "tests.h"
#include "test_idxlist_utils.h"

#include "yaxt.h"
#include "xt_idxstripes.h"
#include "xt_mpi.h"
#include "core/ppm_xfuncs.h"

void do_tests(xt_idxlist idxlist, Xt_idx * ref_indices, unsigned num_indices,
              struct Xt_stripe * ref_stripes, unsigned ref_num_stripes);
void check_stripes(struct Xt_stripe const * stripes, unsigned num_stripes,
                   struct Xt_stripe const * ref_stripes, unsigned ref_num_stripes);

int main(void) {

  // init mpi

  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  Xt_initialize(MPI_COMM_WORLD);

  { // general tests
    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[3] = {{0,5,1},{10,5,1},{20,5,1}};
    unsigned num_stripes = 3;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    Xt_idx ref_indices[15] = {0,1,2,3,4, 10,11,12,13,14, 20,21,22,23,24};
    unsigned ref_num_indices = 15;

    do_tests(idxstripes, ref_indices, ref_num_indices, stripes, num_stripes);

    xt_idxlist_delete(idxstripes);
  }

  { // general tests
    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[3] = {{0,5,1},{10,5,2},{20,5,3}};
    unsigned num_stripes = 3;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    Xt_idx ref_indices[15] = {0,1,2,3,4, 10,12,14,16,18, 20,23,26,29,32};
    unsigned ref_num_indices = 15;
    struct Xt_stripe ref_stripes[11] = {{0,5,1}, {10,1,1}, {12,1,1}, {14,1,1},
                                        {16,1,1}, {18,1,1}, {20,1,1}, {23,1,1},
                                        {26,1,1}, {29,1,1}, {32,1,1}};

    do_tests(idxstripes, ref_indices, ref_num_indices, ref_stripes, 11);

    xt_idxlist_delete(idxstripes);
  }

  { // general tests
    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[2] = {{0,5,6},{1,5,3}};
    unsigned num_stripes = 2;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    Xt_idx ref_indices[10] = {0,6,12,18,24, 1,4,7,10,13};
    unsigned ref_num_indices = 10;
    struct Xt_stripe ref_stripes[10] = {{0,1,1}, {6,1,1}, {12,1,1}, {18,1,1},
                                        {24,1,1}, {1,1,1}, {4,1,1}, {7,1,1},
                                        {10,1,1}, {13,1,1}};

    do_tests(idxstripes, ref_indices, ref_num_indices, ref_stripes, 10);

    xt_idxlist_delete(idxstripes);
  }

  { // general tests
    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[2] = {{0,5,6},{1,5,3}};
    unsigned num_stripes = 2;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    Xt_idx ref_indices[10] = {0,6,12,18,24, 1,4,7,10,13};
    unsigned ref_num_indices = 10;
    struct Xt_stripe ref_stripes[10] = {{0,1,1}, {6,1,1}, {12,1,1}, {18,1,1},
                                        {24,1,1}, {1,1,1}, {4,1,1}, {7,1,1},
                                        {10,1,1}, {13,1,1}};

    do_tests(idxstripes, ref_indices, ref_num_indices, ref_stripes, 10);

    xt_idxlist_delete(idxstripes);
  }

  { // general tests
    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[2] = {{0,5,-1},{1,5,1}};
    unsigned num_stripes = 2;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    Xt_idx ref_indices[10] = {0,-1,-2,-3,-4, 1,2,3,4,5};
    unsigned ref_num_indices = 10;
    struct Xt_stripe ref_stripes[6] = {{0,1,1}, {-1,1,1}, {-2,1,1}, {-3,1,1},
                                       {-4,1,1}, {1,5,1}};

    do_tests(idxstripes, ref_indices, ref_num_indices, ref_stripes, 6);

    xt_idxlist_delete(idxstripes);
  }

  { // general tests
    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[2] = {{9,5,-2},{0,5,2}};
    unsigned num_stripes = 2;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    Xt_idx ref_indices[10] = {9,7,5,3,1, 0,2,4,6,8};
    unsigned ref_num_indices = 10;
    struct Xt_stripe ref_stripes[10] = {{9,1,1}, {7,1,1}, {5,1,1}, {3,1,1},
                                        {1,1,1}, {0,1,1}, {2,1,1}, {4,1,1},
                                        {6,1,1}, {8,1,1}};

    do_tests(idxstripes, ref_indices, ref_num_indices, ref_stripes, 10);

    xt_idxlist_delete(idxstripes);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[2] = {{0,4,1},{6,4,1}};
    struct Xt_stripe stripes_b[1] = {{1,8,1}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 2);
    idxstripes_b = xt_idxstripes_new(stripes_b, 1);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_idx ref_indices[6] = {1,2,3, 6,7,8};
    unsigned ref_num_indices = 6;
    struct Xt_stripe ref_stripes[2] = {{1,3,1}, {6,3,1}};

    do_tests(intersection[0], ref_indices, ref_num_indices, ref_stripes, 2);
    do_tests(intersection[1], ref_indices, ref_num_indices, ref_stripes, 2);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[3] = {{0,4,1},{6,4,1}, {11,4,1}};
    struct Xt_stripe stripes_b[2] = {{1,7,1}, {9,5,1}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 3);
    idxstripes_b = xt_idxstripes_new(stripes_b, 2);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_idx ref_indices[9] = {1,2,3, 6,7, 9, 11,12,13};
    unsigned ref_num_indices = 9;
    struct Xt_stripe ref_stripes[4] = {{1,3,1}, {6,2,1}, {9,1,1}, {11,3,1}};

    do_tests(intersection[0], ref_indices, ref_num_indices,ref_stripes, 4);
    do_tests(intersection[1], ref_indices, ref_num_indices,ref_stripes, 4);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[2] = {{0,3,1},{8,3,1}};
    struct Xt_stripe stripes_b[2] = {{3,5,1}, {11,3,1}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 2);
    idxstripes_b = xt_idxstripes_new(stripes_b, 2);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_idx *ref_indices = NULL;
    unsigned ref_num_indices = 0;
    struct Xt_stripe ref_stripes[1] = {{0,1,1}};

    do_tests(intersection[0], ref_indices, ref_num_indices, ref_stripes, 0);
    do_tests(intersection[1], ref_indices, ref_num_indices, ref_stripes, 0);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[1] = {{0,10,1}};
    struct Xt_stripe stripes_b[2] = {{0,5,2}, {9,5,-2}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 1);
    idxstripes_b = xt_idxstripes_new(stripes_b, 2);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_idx ref_indices[10] = {0,1,2,3,4,5,6,7,8,9};
    unsigned ref_num_indices = 10;
    struct Xt_stripe ref_stripes[1] = {{0,10,1}};

    do_tests(intersection[0], ref_indices, ref_num_indices, ref_stripes, 1);
    do_tests(intersection[1], ref_indices, ref_num_indices, ref_stripes, 1);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[2] = {{0,5,3},{1,5,7}};
    struct Xt_stripe stripes_b[2] = {{0,7,2}, {24,10,-1}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 2);
    idxstripes_b = xt_idxstripes_new(stripes_b, 2);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_idx ref_indices[6] = {0,6,8,12,15,22};
    unsigned ref_num_indices = 6;
    struct Xt_stripe ref_stripes[6] = {{0,1,1}, {6,1,1}, {8,1,1}, {12,1,1}, {15,1,1}, {22,1,1}};

    do_tests(intersection[0], ref_indices, ref_num_indices, ref_stripes, 6);
    do_tests(intersection[1], ref_indices, ref_num_indices, ref_stripes, 6);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

void do_tests(xt_idxlist idxlist, Xt_idx * ref_indices, unsigned num_indices,
              struct Xt_stripe * ref_stripes, unsigned ref_num_stripes) {

  check_idxlist(idxlist, ref_indices, num_indices);

  struct Xt_stripe * stripes;
  unsigned num_stripes;

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

void check_stripes(struct Xt_stripe const * stripes, unsigned num_stripes,
                   struct Xt_stripe const * ref_stripes, unsigned ref_num_stripes) {

  unsigned i;

  if (num_stripes != ref_num_stripes)
    PUT_ERR("wrong number of stripes\n");

  for (i = 0; i < ref_num_stripes; ++i)
    if (stripes[i].start != ref_stripes[i].start ||
        stripes[i].nstrides != ref_stripes[i].nstrides ||
        stripes[i].stride != ref_stripes[i].stride)
      PUT_ERR("error in stripe (start, nstrides and/or stride)\n");
}
