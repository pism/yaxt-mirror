/**
 * @file test_redist_single_array_base.c
 *
 * @copyright Copyright  (C)  2017 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://doc.redmine.dkrz.de/yaxt/html/
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "test_redist_common.h"

int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  // single double
  {
    // send and receive messages
    struct Xt_redist_msg send_msgs[] = { {.rank = 0, .datatype = MPI_DOUBLE} };
    struct Xt_redist_msg recv_msgs[] = { {.rank = 0, .datatype = MPI_DOUBLE} };
    enum { nsend = sizeof(send_msgs) / sizeof(send_msgs[0]) };
    enum { nrecv = sizeof(recv_msgs) / sizeof(recv_msgs[0]) };

    // redist_single_array_base
    Xt_redist redist =
      xt_redist_single_array_base_new(
        nsend, nrecv, send_msgs, recv_msgs, MPI_COMM_WORLD);

    // test communicator of redist
    if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist),
                                     MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    // test exchange
    static const double src_data[] = {0};
    static const double ref_dst_data[] = {0};
    enum { num_ref_values = sizeof(ref_dst_data) / sizeof(ref_dst_data[0]) };
    double dst_data[num_ref_values];
    for (size_t i = 0; i < num_ref_values; ++i) dst_data[i] = -1.0;
    check_redist(redist, src_data, sizeof (src_data[0]) * num_ref_values,
                 dst_data, ref_dst_data);
    Xt_redist redist_copy = xt_redist_copy(redist);
    xt_redist_delete(redist);
    for (size_t i = 0; i < num_ref_values; ++i) dst_data[i] = -1.0;
    check_redist(redist_copy, src_data, sizeof (src_data[0]) * num_ref_values,
                 dst_data, ref_dst_data);

    // clean up
    xt_redist_delete(redist_copy);
  }

  // reverse order of some doubles
  {
    // generate datatypes
    MPI_Datatype send_type, recv_type;
    int recv_displs[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
    enum { nelem = sizeof(recv_displs) / sizeof(recv_displs[0]) };
    MPI_Type_contiguous(nelem, MPI_DOUBLE, &send_type);
    MPI_Type_create_indexed_block(
      nelem, 1, recv_displs, MPI_DOUBLE, &recv_type);
    MPI_Type_commit(&send_type);
    MPI_Type_commit(&recv_type);

    // send and receive messages
    struct Xt_redist_msg send_msgs[] = { {.rank = 0, .datatype = send_type} };
    struct Xt_redist_msg recv_msgs[] = { {.rank = 0, .datatype = recv_type} };
    enum { nsend = sizeof(send_msgs) / sizeof(send_msgs[0]) };
    enum { nrecv = sizeof(recv_msgs) / sizeof(recv_msgs[0]) };

    // redist_single_array_base
    Xt_redist redist =
      xt_redist_single_array_base_new(
        nsend, nrecv, send_msgs, recv_msgs, MPI_COMM_WORLD);
    MPI_Type_free(&recv_type);
    MPI_Type_free(&send_type);

    // test communicator of redist
    if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist),
                                     MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    // test exchange
    static const double src_data[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    static const double ref_dst_data[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
    enum { num_ref_values = sizeof(ref_dst_data) / sizeof(ref_dst_data[0]) };
    double dst_data[num_ref_values];
    for (size_t i = 0; i < num_ref_values; ++i) dst_data[i] = -1.0;
    check_redist(redist, src_data, sizeof (src_data[0]) * num_ref_values,
                 dst_data, ref_dst_data);
    Xt_redist redist_copy = xt_redist_copy(redist);
    xt_redist_delete(redist);
    for (size_t i = 0; i < num_ref_values; ++i) dst_data[i] = -1.0;
    check_redist(redist_copy, src_data, sizeof (src_data[0]) * num_ref_values,
                 dst_data, ref_dst_data);

    // clean up
    xt_redist_delete(redist_copy);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
