/**
 * @file test_xmap_intersection_parallel.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
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
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

#include <yaxt.h>

#include "tests.h"

struct test_message {

  int rank;       // rank of communication partner
  const int *pos;   // positions to be sent/received
  int num_pos; // number of positions
};

struct Xt_redist_msg {

  int rank;
  MPI_Datatype datatype;
};

typedef struct Xt_exchanger_ *Xt_exchanger;

Xt_exchanger
xt_exchanger_irecv_isend_new(int nsend, int nrecv,
                             struct Xt_redist_msg * send_msgs,
                             struct Xt_redist_msg * recv_msgs,
                             MPI_Comm comm, int tag_offset);
Xt_exchanger
xt_exchanger_irecv_isend_packed_new(int nsend, int nrecv,
                                    struct Xt_redist_msg * send_msgs,
                                    struct Xt_redist_msg * recv_msgs,
                                    MPI_Comm comm, int tag_offset);
Xt_exchanger
xt_exchanger_irecv_send_new(int nsend, int nrecv,
                            struct Xt_redist_msg * send_msgs,
                            struct Xt_redist_msg * recv_msgs,
                            MPI_Comm comm, int tag_offset);
Xt_exchanger
xt_exchanger_mix_isend_irecv_new(int nsend, int nrecv,
                                 struct Xt_redist_msg * send_msgs,
                                 struct Xt_redist_msg * recv_msgs,
                                 MPI_Comm comm, int tag_offset);

void xt_exchanger_s_exchange(Xt_exchanger exchanger, const void * src_data,
                             void * dst_data);

void xt_exchanger_delete(Xt_exchanger);

typedef Xt_exchanger (*exchanger_new_func) (int nsend, int nrecv,
                                            struct Xt_redist_msg * send_msgs,
                                            struct Xt_redist_msg * recv_msgs,
                                            MPI_Comm comm, int tag_offset);

static exchanger_new_func parse_options(int *argc, char ***argv);

static void
test_bcast(MPI_Comm comm, exchanger_new_func exchanger_new);

int main(int argc, char **argv)
{

  // init mpi
  xt_mpi_call(MPI_Init(&argc, &argv), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  int my_rank, comm_size;

  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &my_rank), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  if (comm_size != 3) {

    xt_finalize();
    MPI_Finalize();
    return 77;
  }

  exchanger_new_func exchanger_new = parse_options(&argc, &argv);

  test_bcast(MPI_COMM_WORLD, exchanger_new);

  // gather pattern
  for (int i = 0; i < 3; ++i) {

    // setup

    int nsend[3] = {0, 1, 1};
    int nrecv[3] = {2, 0, 0};

    MPI_Datatype MPI_INT_CUSTOM;

    MPI_Type_indexed(1, (int[]){1}, (int[]){1}, MPI_INT, &MPI_INT_CUSTOM);
    MPI_Type_commit(&MPI_INT_CUSTOM);

    struct Xt_redist_msg send_msgs[3][1] =
      {{{.rank=-1, .datatype=MPI_DATATYPE_NULL}},
       {{.rank=(0+i)%3, .datatype=MPI_INT}},
       {{.rank=(0+i)%3, .datatype=MPI_INT}}};
    struct Xt_redist_msg recv_msgs[2][3][2] =
      {{{{.rank=(1+i)%3, .datatype=MPI_INT},
         {.rank=(2+i)%3, .datatype=MPI_INT_CUSTOM}},
        {{.rank=-1, .datatype=MPI_DATATYPE_NULL},
         {.rank=-1, .datatype=MPI_DATATYPE_NULL}},
        {{.rank=-1, .datatype=MPI_DATATYPE_NULL},
         {.rank=-1, .datatype=MPI_DATATYPE_NULL}}},
       {{{.rank=(2+i)%3, .datatype=MPI_INT},
         {.rank=(1+i)%3, .datatype=MPI_INT_CUSTOM}},
        {{.rank=-1, .datatype=MPI_DATATYPE_NULL},
         {.rank=-1, .datatype=MPI_DATATYPE_NULL}},
        {{.rank=-1, .datatype=MPI_DATATYPE_NULL},
         {.rank=-1, .datatype=MPI_DATATYPE_NULL}}}};

    Xt_exchanger exchanger[2] = {exchanger_new(nsend[(my_rank+3-i)%3],
                                               nrecv[(my_rank+3-i)%3],
                                               send_msgs[(my_rank+3-i)%3],
                                               recv_msgs[0][(my_rank+3-i)%3],
                                               MPI_COMM_WORLD, 0),
                                 exchanger_new(nsend[(my_rank+3-i)%3],
                                               nrecv[(my_rank+3-i)%3],
                                               send_msgs[(my_rank+3-i)%3],
                                               recv_msgs[1][(my_rank+3-i)%3],
                                               MPI_COMM_WORLD, 0)};

    // test

    int src_data[1] = {(my_rank+3-i)%3};
    int dst_data[2][2] = {{-1, -1}, {-1, -1}};

    xt_exchanger_s_exchange(exchanger[0], (void*)src_data, (void*)dst_data[0]);
    xt_exchanger_s_exchange(exchanger[1], (void*)src_data, (void*)dst_data[1]);

    if ((my_rank+3-i)%3 == 0) {
      if ((dst_data[0][0]!=1) || (dst_data[0][1]!=2)) PUT_ERR("invalid data\n");
      if ((dst_data[1][0]!=2) || (dst_data[1][1]!=1)) PUT_ERR("invalid data\n");
    }

    // cleanup

    MPI_Type_free(&MPI_INT_CUSTOM);
    xt_exchanger_delete(exchanger[0]);
    xt_exchanger_delete(exchanger[1]);
  }

  // all-to-all pattern
  {

    // setup

    int nsend = 2;
    int nrecv = 2;

    MPI_Datatype MPI_INT_CUSTOM_1, MPI_INT_CUSTOM_2;

    MPI_Type_indexed(1, (int[]){1}, (int[]){1}, MPI_INT, &MPI_INT_CUSTOM_1);
    MPI_Type_indexed(1, (int[]){1}, (int[]){2}, MPI_INT, &MPI_INT_CUSTOM_2);
    MPI_Type_commit(&MPI_INT_CUSTOM_1);
    MPI_Type_commit(&MPI_INT_CUSTOM_2);

    struct Xt_redist_msg send_msgs[2][3][2] =
      {{{{.rank=1, .datatype=MPI_INT}, {.rank=2, .datatype=MPI_INT}},
        {{.rank=0, .datatype=MPI_INT}, {.rank=2, .datatype=MPI_INT}},
        {{.rank=0, .datatype=MPI_INT}, {.rank=1, .datatype=MPI_INT}}},
       {{{.rank=2, .datatype=MPI_INT}, {.rank=1, .datatype=MPI_INT}},
        {{.rank=2, .datatype=MPI_INT}, {.rank=0, .datatype=MPI_INT}},
        {{.rank=1, .datatype=MPI_INT}, {.rank=0, .datatype=MPI_INT}}}};
    struct Xt_redist_msg recv_msgs[2][3][2] =
      {{{{.rank=1, .datatype=MPI_INT_CUSTOM_1},
         {.rank=2, .datatype=MPI_INT_CUSTOM_2}},
        {{.rank=0, .datatype=MPI_INT},
         {.rank=2, .datatype=MPI_INT_CUSTOM_2}},
        {{.rank=0, .datatype=MPI_INT},
         {.rank=1, .datatype=MPI_INT_CUSTOM_1}}},
       {{{.rank=2, .datatype=MPI_INT_CUSTOM_2},
         {.rank=1, .datatype=MPI_INT_CUSTOM_1}},
        {{.rank=2, .datatype=MPI_INT_CUSTOM_2},
         {.rank=0, .datatype=MPI_INT}},
        {{.rank=1, .datatype=MPI_INT_CUSTOM_1},
         {.rank=0, .datatype=MPI_INT}}}};

    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        Xt_exchanger exchanger = exchanger_new(nsend, nrecv,
                                               send_msgs[i][my_rank],
                                               recv_msgs[j][my_rank],
                                               MPI_COMM_WORLD, 0);

        // test

        int src_data[1] = {my_rank};
        int dst_data[3] = {my_rank, my_rank, my_rank};

        xt_exchanger_s_exchange(exchanger, (void*)src_data, (void*)dst_data);

        for (int k = 0; k < 3; ++k)
          if (dst_data[k] != k) PUT_ERR("invalid data\n");

        // cleanup

        xt_exchanger_delete(exchanger);
      }
    }
    MPI_Type_free(&MPI_INT_CUSTOM_1);
    MPI_Type_free(&MPI_INT_CUSTOM_2);
  }

  // round robin pattern
  for (int i = 0; i < 2; ++i) {

    // setup

    int nsend = 1;
    int nrecv = 1;

    struct Xt_redist_msg send_msgs[2][1] =
      {{{.rank=(my_rank+2)%3, .datatype=MPI_INT}},
       {{.rank=(my_rank+1)%3, .datatype=MPI_INT}}};
    struct Xt_redist_msg recv_msgs[2][1] =
      {{{.rank=(my_rank+1)%3, .datatype=MPI_INT}},
       {{.rank=(my_rank+2)%3, .datatype=MPI_INT}}};

    Xt_exchanger exchanger = exchanger_new(nsend, nrecv, send_msgs[i],
                                           recv_msgs[i], MPI_COMM_WORLD, 0);

    // test

    int src_data[1] = {my_rank};
    int dst_data[1] = {-1};

    xt_exchanger_s_exchange(exchanger, (void*)(src_data), (void*)(dst_data));

    if (dst_data[0] != (my_rank+1+i)%3) PUT_ERR("invalid data\n");

    // cleanup

    xt_exchanger_delete(exchanger);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static exchanger_new_func parse_options(int *argc, char ***argv)
{
  exchanger_new_func exchanger_new = xt_exchanger_mix_isend_irecv_new;
  int opt;
  while ((opt = getopt(*argc, *argv, "m:")) != -1) {
    switch (opt) {
    case 'm':
      if (!strcmp(optarg, "irecv_isend"))
        exchanger_new = xt_exchanger_irecv_isend_new;
      else if (!strcmp(optarg, "irecv_isend_packed"))
        exchanger_new = xt_exchanger_irecv_isend_packed_new;
      else if (!strcmp(optarg, "irecv_send"))
        exchanger_new = xt_exchanger_irecv_send_new;
      else if (!strcmp(optarg, "mix_irecv_isend"))
        exchanger_new = xt_exchanger_mix_isend_irecv_new;
      else
      {
        fprintf(stderr, "Unknown exchanger constructor requested %s\n",
                optarg);
        exit(EXIT_FAILURE);
      }
    }
  }
  return exchanger_new;
}

static void
test_bcast(MPI_Comm comm, exchanger_new_func exchanger_new)
{
  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // bcast pattern
  for (int i = 0; i < comm_size; ++i) {

    // setup

    struct Xt_redist_msg *send_msgs = NULL;
    int nsend = 0;
    int nrecv = my_rank != i;

    if (my_rank == i) {
      nsend = comm_size - 1;
      send_msgs = malloc((size_t)nsend * sizeof (*send_msgs));
      if (send_msgs == 0) {
        perror("Failed to allocate message meta-data");
        abort();
      }
      for (size_t j = 0; j < (size_t)i; ++j)
        send_msgs[j] = (struct Xt_redist_msg){.rank=(int)j,
                                              .datatype=MPI_INT};
      for (size_t j = (size_t)i; j < (size_t)nsend; ++j)
        send_msgs[j] = (struct Xt_redist_msg){.rank=(1+(int)j)%comm_size,
                                              .datatype=MPI_INT};
    }
    struct Xt_redist_msg recv_msgs[2] =
      {{.rank=-1, .datatype=MPI_DATATYPE_NULL},
       {.rank=i, .datatype=MPI_INT}};

    Xt_exchanger exchanger = exchanger_new(nsend, nrecv,
                                           send_msgs,
                                           recv_msgs + (my_rank != i),
                                           comm, 0);

    // test

    int src_data[1] = { my_rank == i ? 4711 : -1 };
    int dst_data[1] = { my_rank == i ? 4711 : -1 };

    xt_exchanger_s_exchange(exchanger, (void*)(src_data), (void*)(dst_data));

    if (dst_data[0] != 4711) PUT_ERR("invalid data\n");

    // cleanup
    free(send_msgs);
    xt_exchanger_delete(exchanger);
  }
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
