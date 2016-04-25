/**
 * @file test_request_parallel.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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
#include <unistd.h>

#include <yaxt.h>

#include "tests.h"

int main(int argc, char **argv) {

  // init mpi
  xt_mpi_call(MPI_Init(&argc, &argv), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  int rank, size;

  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &size), MPI_COMM_WORLD);

  { // trivial test

    MPI_Request requests[] = {MPI_REQUEST_NULL};
    int n = sizeof(requests) / sizeof(requests[0]);
    MPI_Comm comm = MPI_COMM_WORLD;

    Xt_request request = xt_request_msgs_new(n, requests, comm);

    if (request == XT_REQUEST_NULL)
      PUT_ERR("request ==  XT_REQUEST_NULL after xt_request_msgs_new\n");

    xt_request_wait(&request);

    if (request != XT_REQUEST_NULL)
      PUT_ERR("request !=  XT_REQUEST_NULL after xt_request_wait\n");
  }

  { // each process receives and sends data from/to every other process

    int msg_size = 1 * 1024 * 1024; // 4 * 1 MiB
    int * send_buf = malloc(msg_size * sizeof(*send_buf));
    int * recv_buf = malloc(size * msg_size * sizeof(*recv_buf));

    if ((send_buf == NULL) || (recv_buf == NULL)) {

      free(send_buf), free(recv_buf);
      PUT_ERR("malloc failed for \"big\" test\n");

    } else {

      for (int i = 0; i < msg_size; ++i) send_buf[i] = rank;
      for (int i = 0; i < size * msg_size; ++i) recv_buf[i] = -1;

      MPI_Request requests[2 * size];
      int n = sizeof(requests) / sizeof(requests[0]);
      MPI_Comm comm = MPI_COMM_WORLD;

      for (unsigned i = 0; i < size; ++i)
        xt_mpi_call(MPI_Isend(send_buf, msg_size, MPI_INT, i, 0, MPI_COMM_WORLD,
                              &(requests[i + size])), MPI_COMM_WORLD);
      for (unsigned i = 0; i < size; ++i)
        xt_mpi_call(MPI_Irecv(&(recv_buf[i * msg_size]), msg_size, MPI_INT, i, 0,
                              MPI_COMM_WORLD, &(requests[i])), MPI_COMM_WORLD);

      Xt_request request = xt_request_msgs_new(n, requests, comm);

      if (request == XT_REQUEST_NULL)
        PUT_ERR("request ==  XT_REQUEST_NULL after xt_request_msgs_new\n");

      xt_request_wait(&request);

      if (request != XT_REQUEST_NULL)
        PUT_ERR("request !=  XT_REQUEST_NULL after xt_request_wait\n");

      for (int i = 0; i < size; ++i)
        for (int j = 0; j < msg_size; ++j)
          if (recv_buf[i * msg_size + j] != i)
            PUT_ERR("recv_buf[i * msg_size + j] != i\n");
    }
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}
