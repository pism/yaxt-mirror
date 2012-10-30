#include <stdlib.h>

#include <mpi.h>

#define VERBOSE

#include "tests.h"

#include "yaxt.h"
#include "xt_idxlist.h"
#include "xt_idxvec.h"
#include "xt_idxstripes.h"
#include "xt_idxsection.h"
#include "xt_idxlist_collection.h"
#include "test_idxlist_utils.h"
#include "xt_mpi.h"
#include "core/ppm_xfuncs.h"

static void
check_stripes(struct Xt_stripe const * stripes, Xt_count num_stripes,
              struct Xt_stripe const * ref_stripes, Xt_count ref_num_stripes);

int main(void) {

  // init mpi

  mpi_err_handler(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // test packing and unpacking of collection index lists

    Xt_idx index_list[2][7] = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1}};
    Xt_count num_indices = 7;

    Xt_idxlist idxlists[2];

    // create two index lists

    idxlists[0] = xt_idxvec_new(index_list[0], num_indices);
    idxlists[1] = xt_idxvec_new(index_list[1], num_indices);

    // generate a collection index list

    Xt_idxlist collectionlist;

    collectionlist = xt_idxlist_collection_new(idxlists, 2);

    xt_idxlist_delete(idxlists[0]);
    xt_idxlist_delete(idxlists[1]);

    // test generated collection list

    check_idxlist(collectionlist, &(index_list[0][0]), 2*num_indices);

    // get the size of the index list

    size_t buffer_size;

    buffer_size = xt_idxlist_get_pack_size(collectionlist, MPI_COMM_WORLD);

    // allocate send buffers

    void* send_buffer, * recv_buffer;

    send_buffer = xmalloc(buffer_size);
    recv_buffer = xmalloc(buffer_size);

    // pack the index list

    int position = 0;

    xt_idxlist_pack(collectionlist, send_buffer, buffer_size, &position, MPI_COMM_WORLD);

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

    Xt_idxlist collectionlist_copy;

    position = 0;
    collectionlist_copy = xt_idxlist_unpack(recv_buffer, buffer_size, &position,
                                            MPI_COMM_WORLD);

    // check received collection list

    check_idxlist(collectionlist_copy, &(index_list[0][0]), 2*num_indices);

    // compute intersection between the two index lists

    Xt_idxlist intersection;

    intersection = xt_idxlist_get_intersection(collectionlist, collectionlist_copy);

    // check intersection data
    Xt_idx ref_index_data[14] = {1,1,2,2,3,3,4,4,5,5,6,6,7,7};
    check_idxlist(intersection, ref_index_data, 2*num_indices);

    // check the conversion to stripes

    struct Xt_stripe * stripes;
    Xt_count num_stripes;

    xt_idxlist_get_index_stripes(collectionlist, &stripes, &num_stripes);

    struct Xt_stripe ref_stripes[8] = {{1,7,1}, {7,1,1}, {6,1,1}, {5,1,1},
                                       {4,1,1}, {3,1,1}, {2,1,1}, {1,1,1}};
    Xt_count ref_num_stripes = 8;

    check_stripes(stripes, num_stripes, ref_stripes, ref_num_stripes);

    // clean up

    free(stripes);

    free(send_buffer);
    free(recv_buffer);

    xt_idxlist_delete(collectionlist);
    xt_idxlist_delete(collectionlist_copy);
    xt_idxlist_delete(intersection);
  }

  { // test copying of collection index lists

    Xt_idx index_list[2][7] = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1}};
    Xt_count num_indices = 7;

    Xt_idxlist idxlists[2];

    // create two index lists

    idxlists[0] = xt_idxvec_new(index_list[0], num_indices);
    idxlists[1] = xt_idxvec_new(index_list[1], num_indices);

    // generate a collection index list

    Xt_idxlist collectionlist, collectionlist_copy;

    collectionlist = xt_idxlist_collection_new(idxlists, 2);

    xt_idxlist_delete(idxlists[0]);
    xt_idxlist_delete(idxlists[1]);

    // test generated collection list

    check_idxlist(collectionlist, &(index_list[0][0]), 2*num_indices);

    // copy the index list

    collectionlist_copy = xt_idxlist_copy(collectionlist);

    // check copy of collection list

    check_idxlist(collectionlist_copy, &(index_list[0][0]), 2*num_indices);

    // check the conversion to stripes

    struct Xt_stripe * stripes;
    Xt_count num_stripes;

    xt_idxlist_get_index_stripes(collectionlist, &stripes, &num_stripes);

    struct Xt_stripe ref_stripes[8] = {{1,7,1}, {7,1,1}, {6,1,1}, {5,1,1},
                                       {4,1,1}, {3,1,1}, {2,1,1}, {1,1,1}};
    Xt_count ref_num_stripes = 8;

    check_stripes(stripes, num_stripes, ref_stripes, ref_num_stripes);

    // clean up

    free(stripes);

    xt_idxlist_delete(collectionlist);
    xt_idxlist_delete(collectionlist_copy);
  }

  { // test computing intersection of index list collections

    Xt_idx index_list[3][7] = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1},{2,6,1,4,7,3,0}};
    Xt_idx sorted_index_list[21] = {0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,6,6,6,7,7,7};
    Xt_count num_indices = 7;

    Xt_idxlist idxlists[3];

    // create two index lists

    idxlists[0] = xt_idxvec_new(index_list[0], num_indices);
    idxlists[1] = xt_idxvec_new(index_list[1], num_indices);
    idxlists[2] = xt_idxvec_new(index_list[2], num_indices);

    // generate a collection index list

    Xt_idxlist collectionlist;

    collectionlist = xt_idxlist_collection_new(idxlists, 3);

    xt_idxlist_delete(idxlists[0]);
    xt_idxlist_delete(idxlists[1]);
    xt_idxlist_delete(idxlists[2]);

    // test generated collection list

    check_idxlist(collectionlist, &(index_list[0][0]), 3*num_indices);

    // test computing of intersection

    Xt_idxlist intersection, ref_idxvec;

    ref_idxvec = xt_idxvec_new(&(index_list[0][0]), 3 * num_indices);

    intersection = xt_idxlist_get_intersection(ref_idxvec, collectionlist);

    check_idxlist(intersection, sorted_index_list, 3 * num_indices);

    xt_idxlist_delete(intersection);

    intersection = xt_idxlist_get_intersection(collectionlist, ref_idxvec);

    check_idxlist(intersection, sorted_index_list, 3 * num_indices);

    xt_idxlist_delete(intersection);
    xt_idxlist_delete(collectionlist);
  }

  { // test using various tpyes of index lists

    Xt_idx index_list[] = {1,3,5,7,9,11};
    struct Xt_stripe stripes[] = {{0,5,2},{1,5,2}};
    unsigned global_size[2] = {10,10};
    unsigned local_size[2] = {5,5};
    unsigned local_start[2] = {2,2};

    Xt_idx ref_index_list[] = {1,3,5,7,9,11,
                               0,2,4,6,8,1,3,5,7,9,
                               22,23,24,25,26,
                               32,33,34,35,36,
                               42,43,44,45,46,
                               52,53,54,55,56,
                               62,63,64,65,66};

    Xt_idxlist idxlists[3] = {xt_idxvec_new(index_list, sizeof(index_list)/sizeof(index_list[0])),
                              xt_idxstripes_new(stripes, sizeof(stripes)/sizeof(stripes[0])),
                              xt_idxsection_new(0, 2, global_size, local_size, local_start)};

    // generate a collection index list

    Xt_idxlist collectionlist;

    collectionlist = xt_idxlist_collection_new(idxlists, 3);

    xt_idxlist_delete(idxlists[0]);
    xt_idxlist_delete(idxlists[1]);
    xt_idxlist_delete(idxlists[2]);

    // test generated collection list

    check_idxlist(collectionlist, ref_index_list,
                  sizeof(ref_index_list)/sizeof(ref_index_list[0]));

    xt_idxlist_delete(collectionlist);
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
check_stripes(struct Xt_stripe const * stripes, Xt_count num_stripes,
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
