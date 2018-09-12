/**
 * @file test_exchanger_parallel.c
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>
#include <yaxt.h>

#include "../src/xt_exchanger.h"
#include "../src/xt_redist_internal.h"

#include "tests.h"

struct test_message {

  int rank;       // rank of communication partner
  const int *pos;   // positions to be sent/received
  int num_pos; // number of positions
};

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

Xt_exchanger
xt_exchanger_neigh_alltoall_new(int nsend, int nrecv,
                                struct Xt_redist_msg *send_msgs,
                                struct Xt_redist_msg *recv_msgs,
                                MPI_Comm comm, int tag_offset);

typedef Xt_exchanger (*exchanger_new_func) (int nsend, int nrecv,
                                            struct Xt_redist_msg * send_msgs,
                                            struct Xt_redist_msg * recv_msgs,
                                            MPI_Comm comm, int tag_offset);

static exchanger_new_func *parse_options(int *argc, char ***argv);

static void
test_bcast(MPI_Comm comm, exchanger_new_func exchanger_new);
static void
test_gather(MPI_Comm comm, exchanger_new_func exchanger_new);
static void
test_all2all(MPI_Comm comm, exchanger_new_func exchanger_new);
static void
test_rr(MPI_Comm comm, exchanger_new_func exchanger_new);

int main(int argc, char **argv)
{

  // init mpi
  xt_mpi_call(MPI_Init(&argc, &argv), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  int my_rank, comm_size;

  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &my_rank), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  exchanger_new_func *exchangers_new = parse_options(&argc, &argv);

  for (size_t i = 0; exchangers_new[i] != (exchanger_new_func)0; ++i) {
    exchanger_new_func exchanger_new = exchangers_new[i];

    test_bcast(MPI_COMM_WORLD, exchanger_new);

    test_gather(MPI_COMM_WORLD, exchanger_new);

    test_all2all(MPI_COMM_WORLD, exchanger_new);

    test_rr(MPI_COMM_WORLD, exchanger_new);
  }
  free(exchangers_new);
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static exchanger_new_func *parse_options(int *argc, char ***argv)
{
  exchanger_new_func *exchangers_new = malloc(2 * sizeof (*exchangers_new));
  exchangers_new[0] = xt_exchanger_mix_isend_irecv_new;
  exchangers_new[1] = (exchanger_new_func)0;
  size_t cur_ex = 0;
  int opt;
  while ((opt = getopt(*argc, *argv, "m:")) != -1) {
    switch (opt) {
    case 'm':
      if (!strcmp(optarg, "irecv_isend"))
        exchangers_new[cur_ex] = xt_exchanger_irecv_isend_new;
      else if (!strcmp(optarg, "irecv_isend_packed"))
        exchangers_new[cur_ex] = xt_exchanger_irecv_isend_packed_new;
      else if (!strcmp(optarg, "irecv_send"))
        exchangers_new[cur_ex] = xt_exchanger_irecv_send_new;
      else if (!strcmp(optarg, "mix_irecv_isend"))
        exchangers_new[cur_ex] = xt_exchanger_mix_isend_irecv_new;
      else if (!strcmp(optarg, "neigh_alltoall"))
#if MPI_VERSION >= 3
        exchangers_new[cur_ex] = xt_exchanger_neigh_alltoall_new;
#else
      {
        fputs("xt_exchanger_neigh_alltoall_new requires MPI version 3.0 or "
              "higher\n", stderr);
        continue;
      }
#endif
      else {
        fprintf(stderr, "Unknown exchanger constructor requested %s\n",
                optarg);
        exit(EXIT_FAILURE);
      }
      ++cur_ex;
      void *temp = realloc(exchangers_new,
                           sizeof (*exchangers_new) * (cur_ex + 1));
      if (!temp) {
        perror("failed reallocation");
        abort();
      }
      exchangers_new = temp;
      exchangers_new[cur_ex] = (exchanger_new_func)0;
    }
  }
  return exchangers_new;
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
    int test_async = (exchanger_new != xt_exchanger_irecv_send_new);
    for (int async = 0; async < 1 + test_async; ++async) {

      int src_data[1] = { my_rank == i ? 4711 : -1 };
      int dst_data[1] = { my_rank == i ? 4711 : -1 };

      if (async) {
        Xt_request request;
        int flag;
        xt_exchanger_a_exchange(exchanger, (void*)(src_data),
                                (void*)(dst_data), &request);
        xt_request_test(&request, &flag);
        xt_request_wait(&request);
        xt_request_test(&request, &flag);
        if (!flag) PUT_ERR("invalid flag result\n");
      } else {
        xt_exchanger_s_exchange(exchanger, (void*)(src_data),
                                (void*)(dst_data));
      }

      if (dst_data[0] != 4711) PUT_ERR("invalid data\n");
    }

    // cleanup
    free(send_msgs);
    xt_exchanger_delete(exchanger);
  }
}

static void
test_gather(MPI_Comm comm, exchanger_new_func exchanger_new)
{
  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  // gather pattern
  // prepare datatypes outside of loop for load-balance
  MPI_Datatype *dt_by_ofs = malloc((size_t)comm_size * sizeof (*dt_by_ofs));
  if (!dt_by_ofs) {
    perror("Failed to allocate receive datatypes");
    abort();
  }
  dt_by_ofs[0] = MPI_INT;
  for (size_t j = 1; j < (size_t)comm_size; ++j)
  {
    MPI_Type_indexed(1, (int[]){1}, (int[]){(int)j}, MPI_INT, dt_by_ofs + j);
    MPI_Type_commit(dt_by_ofs + j);
  }
  int *dst_data = malloc(((size_t)comm_size - 1) * sizeof (*dst_data) * 2);
  if (comm_size - 1 && !dst_data) {
    perror("Failed to allocate message receive buffer");
    abort();
  }
  for (int i = 0; i < comm_size; ++i) {

    // setup
    int nsend = i != my_rank;
    size_t nrecv = 0;


    struct Xt_redist_msg send_msgs[1] = {{.rank=i, .datatype=MPI_INT}};
    struct Xt_redist_msg *recv_msgs = NULL;
    if (my_rank == i) {
      nrecv = (size_t)comm_size - 1;
      recv_msgs = malloc(nrecv * sizeof (*recv_msgs) * 2);
      if (nrecv && !recv_msgs) {
        perror("Failed to allocate message meta-data");
        abort();
      }
      for (size_t j = 0; j < nrecv; ++j) {
        recv_msgs[j].rank = (i + (int)j + 1)%comm_size;
        recv_msgs[j].datatype = dt_by_ofs[j];

        recv_msgs[j + nrecv].rank
          = (comm_size - (int)j - 1 + i)%comm_size;
        recv_msgs[j + nrecv].datatype = dt_by_ofs[j];
      }
    }

    enum { exchange_fwd, exchange_rev, num_exchanges };
    Xt_exchanger exchanger[num_exchanges];
    for (size_t j = 0; j < num_exchanges; ++j)
      exchanger[j] = exchanger_new(nsend, (int)nrecv,
                                   send_msgs, recv_msgs + j * nrecv, comm, 0);

    // test
    int test_async = (exchanger_new != xt_exchanger_irecv_send_new);
    for (int async = 0; async < 1 + test_async; ++async) {
      int src_data[1] = {(my_rank+comm_size-i)%comm_size};
      for (size_t j = 0; j < 2 * nrecv; ++j)
        dst_data[j] = -1;

      for (size_t j = 0; j < num_exchanges; ++j) {
        if (async) {
          Xt_request request;
          int flag;
          xt_exchanger_a_exchange(exchanger[j], src_data,
                                  dst_data + j * nrecv, &request);
          xt_request_test(&request, &flag);
          xt_request_wait(&request);
          xt_request_test(&request, &flag);
          if (!flag) PUT_ERR("invalid flag result\n");
        } else {
          xt_exchanger_s_exchange(exchanger[j], src_data, dst_data + j * nrecv);
        }
      }

      for (size_t j = 0; j < nrecv; ++j)
        if (((size_t)dst_data[j] != j + 1)
            || (dst_data[nrecv + j] != comm_size - (int)j - 1))
          PUT_ERR("invalid data\n");
    }

    // cleanup
    free(recv_msgs);
    xt_exchanger_delete(exchanger[0]);
    xt_exchanger_delete(exchanger[1]);
  }
  free(dst_data);
  for (size_t j = 1; j < (size_t)comm_size; ++j)
    MPI_Type_free(dt_by_ofs + j);
  free(dt_by_ofs);
}

static void
test_all2all(MPI_Comm comm, exchanger_new_func exchanger_new)
{
  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  // all-to-all pattern
  // setup
  int nsend = comm_size - 1;
  int nrecv = comm_size - 1;

  struct Xt_redist_msg *send_msgs = malloc((size_t)nsend * sizeof (*send_msgs));
  if (nsend && !send_msgs) {
    perror("Failed to allocate tx message meta-data");
    abort();
  }
  for (size_t j = 0; j < (size_t)nsend; ++j)
    send_msgs[j].datatype = MPI_INT;
  struct Xt_redist_msg *recv_msgs = malloc((size_t)nrecv * sizeof (*recv_msgs));
  if (nrecv && !recv_msgs) {
    perror("Failed to allocate rx message meta-data");
    abort();
  }
  int *dst_data = malloc((size_t)comm_size * sizeof (*dst_data));
  if (!dst_data) {
    perror("Failed to allocate recv data");
    abort();
  }
  MPI_Datatype *dt_by_ofs = malloc((size_t)comm_size * sizeof (*dt_by_ofs));
  if (!dt_by_ofs) {
    perror("Failed to allocate datatypes");
    abort();
  }
  dt_by_ofs[0] = my_rank != 0 ? MPI_INT : MPI_DATATYPE_NULL;
  for (size_t i = 1; i < (size_t)comm_size; ++i)
    if (i != (size_t)my_rank) {
      MPI_Type_indexed(1, (int[]){1}, (int[]){(int)i}, MPI_INT, dt_by_ofs + i);
      MPI_Type_commit(dt_by_ofs + i);
    } else
      dt_by_ofs[i] = MPI_DATATYPE_NULL;
  for (size_t i = 0; i < (size_t)comm_size - 1; ++i) {
    for (size_t j = 0; j < (size_t)nsend; ++j) {
      int ofs = my_rank + 1 + (int)i + (int)j;
      send_msgs[j].rank = (ofs + (ofs >= comm_size + my_rank))%comm_size;
    }
    for (size_t j = 0; j < (size_t)comm_size - 1; ++j) {
      for (size_t k = 0; k < (size_t)nrecv; ++k) {
        int ofs = ((int)i + (int)j + (int)k)%(comm_size - 1);
        ofs += ofs >= my_rank;
        recv_msgs[k].rank = ofs;
        recv_msgs[k].datatype = dt_by_ofs[ofs];
      }
      Xt_exchanger exchanger = exchanger_new(nsend, nrecv,
                                             send_msgs,
                                             recv_msgs,
                                             MPI_COMM_WORLD, 0);
      // test
      int test_async = (exchanger_new != xt_exchanger_irecv_send_new);
      for (int async = 0; async < 1 + test_async; ++async) {
        int src_data[1] = {my_rank};
        for (size_t k = 0; k < (size_t)comm_size; ++k)
          dst_data[k] = my_rank;

        if (async) {
          Xt_request request;
          int flag;
          xt_exchanger_a_exchange(exchanger, (void*)src_data, (void*)dst_data,
                                  &request);
          xt_request_test(&request, &flag);
          xt_request_wait(&request);
          xt_request_test(&request, &flag);
          if (!flag) PUT_ERR("invalid flag result\n");
        } else {
          xt_exchanger_s_exchange(exchanger, (void*)src_data, (void*)dst_data);
        }

        for (int k = 0; k < comm_size; ++k)
          if (dst_data[k] != k) PUT_ERR("invalid data\n");
      }

      // cleanup
      xt_exchanger_delete(exchanger);
    }
  }
  for (size_t i = 1; i < (size_t)comm_size; ++i)
    if (i != (size_t)my_rank)
      MPI_Type_free(dt_by_ofs + i);
  free(dt_by_ofs);
  free(dst_data);
  free(recv_msgs);
  free(send_msgs);
}

static void
test_rr(MPI_Comm comm, exchanger_new_func exchanger_new)
{
  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  // round robin pattern
  for (int i = 1; i < comm_size; ++i) {

    // setup
    enum { nsend = 1, nrecv = 1 };
    struct Xt_redist_msg send_msgs[nsend]
      = {{.rank=(my_rank + i)%comm_size, .datatype=MPI_INT}};
    struct Xt_redist_msg recv_msgs[nrecv]
      = {{.rank=(my_rank + comm_size - i)%comm_size, .datatype=MPI_INT}};

    Xt_exchanger exchanger = exchanger_new(nsend, nrecv, send_msgs,
                                           recv_msgs, comm, 0);

    // test
    int test_async = (exchanger_new != xt_exchanger_irecv_send_new);
    for (int async = 0; async < 1 + test_async; ++async) {
      int src_data[1] = {my_rank};
      int dst_data[1] = {-1};

      if (async) {
        Xt_request request;
        int flag;
        xt_exchanger_a_exchange(exchanger, src_data, dst_data, &request);
        xt_request_test(&request, &flag);
        xt_request_wait(&request);
        xt_request_test(&request, &flag);
        if (!flag) PUT_ERR("invalid flag result\n");
      } else {
        xt_exchanger_s_exchange(exchanger, src_data, dst_data);
      }

      if (dst_data[0] != (my_rank + comm_size - i)%comm_size)
        PUT_ERR("invalid data\n");
    }

    // cleanup
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
