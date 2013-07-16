/**
 * @file test_idxempty.c
 *
 * @copyright Copyright  (C)  2012 Moritz Hanke <hanke@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 * URL: https://redmine.dkrz.de/doc/yaxt/html/index.html
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdlib.h>
#include <limits.h>
#include <assert.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "test_idxlist_utils.h"
#include "core/ppm_xfuncs.h"

int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // test packing and unpacking

    Xt_idxlist idxempty;

    idxempty = xt_idxempty_new();

    check_idxlist(idxempty, NULL, 0);

    size_t buffer_size = xt_idxlist_get_pack_size(idxempty, MPI_COMM_WORLD);

    void* send_buffer, * recv_buffer;

    send_buffer = xmalloc(buffer_size);
    recv_buffer = xmalloc(buffer_size);

    int position = 0;

    assert(buffer_size <= INT_MAX);
    xt_idxlist_pack(idxempty, send_buffer, (int)buffer_size, &position,
                    MPI_COMM_WORLD);

    // send the buffer

    int rank;

    xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);

    MPI_Request request;

    xt_mpi_call(MPI_Isend(send_buffer, (int)buffer_size, MPI_PACKED, rank, 0,
                          MPI_COMM_WORLD, &request), MPI_COMM_WORLD);
    xt_mpi_call(MPI_Request_free(&request), MPI_COMM_WORLD);

    // receive the buffer

    xt_mpi_call(MPI_Recv(recv_buffer, (int)buffer_size, MPI_PACKED, rank, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE),
                MPI_COMM_WORLD);

    // unpack the buffer

    Xt_idxlist idxempty_copy;

    position = 0;
    idxempty_copy = xt_idxlist_unpack(recv_buffer, (int)buffer_size, &position,
                                      MPI_COMM_WORLD);

    // check the received index vector

    check_idxlist(idxempty_copy, NULL, 0);

    // compute intersection between the two index vectors

    Xt_idxlist intersection;

    intersection = xt_idxlist_get_intersection(idxempty, idxempty_copy);

    // check the computed intersection (should be identically to the original
    // index vector)

    check_idxlist(intersection, NULL, 0);

    // check the conversion to stripes

    struct Xt_stripe * stripes;
    int num_stripes;

    xt_idxlist_get_index_stripes(idxempty, &stripes, &num_stripes);

    if (stripes != NULL || num_stripes != 0)
      PUT_ERR("ERROR: xt_idxlist_get_index_stripes\n");

    // check the computation of a bounding box

    unsigned ndim = 3;
    Xt_int global_size[ndim];
    Xt_int global_start_index = 0;
    struct Xt_bounds bounds[ndim];

    for (unsigned i = 0; i < ndim; ++i)
      global_size[i] = 10;

    xt_idxlist_get_bounding_box(idxempty, ndim, global_size,
                                global_start_index, bounds);

    for (unsigned i = 0; i < ndim; ++i)
      if (bounds[i].size != 0)
        PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

    // clean up

    free(stripes);

    free(send_buffer);
    free(recv_buffer);

    xt_idxlist_delete(idxempty);
    xt_idxlist_delete(idxempty_copy);
    xt_idxlist_delete(intersection);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}
