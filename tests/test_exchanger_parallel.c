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
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>

#include "yaxt.h"
#include "../src/xt_exchanger.h"

#include "ctest_common.h"
#include "test_redist_common.h"
#include "core/ppm_xfuncs.h"
#include "tests.h"

struct test_message {

  int rank;       // rank of communication partner
  const int *pos;   // positions to be sent/received
  int num_pos; // number of positions
};

static int *parse_options(int *argc, char ***argv);

static void
test_bcast(MPI_Comm comm, Xt_config config);
static void
test_gather(MPI_Comm comm, Xt_config config);
static void
test_all2all(MPI_Comm comm, Xt_config config);
static void
test_rr(MPI_Comm comm, Xt_config config);
static void
test_intercomm_all2all(MPI_Comm comm, Xt_config config);

static int test_freq = 3;

int main(int argc, char **argv)
{

  test_init_mpi(&argc, &argv, MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  // check that non-existent exchanger results in -1
  {
    int no_exch_id = xt_exchanger_id_by_name("non-existant fantasy exchanger");
    if (no_exch_id != -1)
      PUT_ERR("invalid exchanger");
  }
  int *exchanger_code = parse_options(&argc, &argv);
  Xt_config config = xt_config_new();
  static const int mt_configs_2_test[] = {
    XT_MT_NONE,
#ifdef _OPENMP
    XT_MT_OPENMP,
#endif
  };
  static const int send_comm_modes_2_test[] = {
    xt_exch_send_standard,
    xt_exch_send_buffered,
    xt_exch_send_synchronous,
  };
  enum {
    num_mt_modes = sizeof (mt_configs_2_test) / sizeof (mt_configs_2_test[0]),
    num_send_modes = (sizeof (send_comm_modes_2_test)
                      / sizeof (send_comm_modes_2_test[0])),
  };
  int world_size, world_rank;
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &world_size), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &world_rank), MPI_COMM_WORLD);
  size_t buf_size = (sizeof (int) + (size_t)MPI_BSEND_OVERHEAD)
    * (size_t)world_size;
  void *buf = xmalloc(buf_size);
  MPI_Comm *split_comm
    = xmalloc(sizeof(*split_comm) * (size_t)(world_size/2 + 1));
  split_comm[0] = MPI_COMM_WORLD;
  for (int split_rank = 1, world_size_2 = world_size/2;
       split_rank <= world_size_2; ++split_rank)
    xt_mpi_call(MPI_Comm_split(MPI_COMM_WORLD, world_rank < split_rank,
                               world_rank, split_comm+split_rank),
                MPI_COMM_WORLD);
  xt_mpi_call(MPI_Buffer_attach(buf, (int)buf_size), MPI_COMM_WORLD);
  for (int split_rank = 0, world_size_2 = world_size/2;
       split_rank <= world_size_2; ++split_rank) {
    MPI_Comm comm = split_comm[split_rank];
    for (size_t m = 0; m < num_mt_modes; ++m) {
#ifdef _OPENMP
      int thread_support_provided;
      xt_mpi_call(MPI_Query_thread(&thread_support_provided), comm);
      if (mt_configs_2_test[m] == XT_MT_OPENMP
          && thread_support_provided != MPI_THREAD_MULTIPLE)
        continue;
#endif
      xt_config_set_redist_mthread_mode(config, mt_configs_2_test[m]);
#if XT_ENABLE_ALTERNATIVE_EXCHANGER_SEND_MODES
      for (size_t n = 0; n < num_send_modes; ++n) {
        xt_config_set_exchange_send_comm_mode(config, send_comm_modes_2_test[n]);
#endif
        for (size_t i = 0; exchanger_code[i] != -1; ++i) {
          xt_config_set_exchange_method(config, exchanger_code[i]);
          test_bcast(comm, config);

          test_gather(comm, config);

          test_all2all(comm, config);

          test_rr(comm, config);

          test_intercomm_all2all(comm, config);
        }
#if XT_ENABLE_ALTERNATIVE_EXCHANGER_SEND_MODES
      }
#endif
    }
  }
  for (int split_rank = 1, world_size_2 = world_size/2;
       split_rank <= world_size_2; ++split_rank)
    xt_mpi_call(MPI_Comm_free(split_comm+split_rank), MPI_COMM_WORLD);
  free(split_comm);
  {
    int buf_size_;
    void *buf_;
    xt_mpi_call(MPI_Buffer_detach(&buf_, &buf_size_), MPI_COMM_WORLD);
    assert(buf == buf_ && buf_size_ >= 0 && (size_t)buf_size_ == buf_size);
  }
  free(buf);
  xt_config_delete(config);
  free(exchanger_code);

  xt_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static int *parse_options(int *argc, char ***argv)
{
  int *exchangers_new = xmalloc(2 * sizeof (*exchangers_new));
  exchangers_new[0] = xt_exchanger_mix_isend_irecv;
  exchangers_new[1] = -1;
  size_t cur_ex = 0;
  int opt;
  while ((opt = getopt(*argc, *argv, "m:s:")) != -1) {
    switch (opt) {
    case 'm':
      {
        int exchanger_new_id = xt_exchanger_id_by_name(optarg);
        if (exchanger_new_id == -1)
        {
          fprintf(stderr, "Unknown exchanger constructor requested: %s\n",
                  optarg);
          exit(EXIT_FAILURE);
        }
#ifndef XT_CAN_USE_MPI_NEIGHBOR_ALLTOALL
        else if (exchanger_new_id == xt_exchanger_neigh_alltoall)
        {
          fputs("xt_exchanger_neigh_alltoall_new requires MPI version 3.0 or "
                "higher\n", stderr);
          continue;
        }
#endif
#ifndef XT_ENABLE_DDT_EXCHANGER
        else if (exchanger_new_id == xt_exchanger_irecv_isend_ddt_packed)
        {
          fputs("xt_exchanger_irecv_isend_ddt_packed disabled.\n", stderr);
          continue;
        }
#endif
        exchangers_new[cur_ex] = exchanger_new_id;
        ++cur_ex;
        exchangers_new = xrealloc(exchangers_new,
                                  sizeof (*exchangers_new) * (cur_ex + 1));
        exchangers_new[cur_ex] = -1;
      }
      break;
    case 's':
      {
        char *endptr;
        errno = 0;
        long v = strtol(optarg, &endptr, 0);
        if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN))
            || (errno != 0 && v == 0)) {
          perror("failed to parse argument to -s option");
          exit(EXIT_FAILURE);
        }
        if (endptr == optarg) {
          fputs("malformed argument to -s option, no digits were found\n",
                stderr);
          exit(EXIT_FAILURE);
        }
        if (v < 1 || v > INT_MAX) {
          fprintf(stderr, "value of -s option (%ld) out of range [1,%d]\n",
                  v, INT_MAX);
          exit(EXIT_FAILURE);
        }
        test_freq = (int)v;
      }
      break;
    }
  }
  return exchangers_new;
}

static void
data_exchange(Xt_exchanger exchanger,
              const void *src_data, void *dst_data, bool async)
{
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
}


static void
test_bcast(MPI_Comm comm, Xt_config config)
{
  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // bcast pattern
  int incr = comm_size/test_freq + (comm_size < test_freq);
  int test_async = (xt_config_get_exchange_method(config)
                    != xt_exchanger_irecv_send);
  for (int i = 0; i < comm_size; i += incr) {

    // setup

    struct Xt_msg_param *send_msgs = NULL;
    int nsend = 0;
    int nrecv = my_rank != i;

    if (my_rank == i) {
      nsend = comm_size - 1;
      send_msgs = xmalloc((size_t)nsend * sizeof (*send_msgs));
      for (size_t j = 0; j < (size_t)i; ++j)
        send_msgs[j] = (struct Xt_msg_param){.rank=(int)j,
                                             .datatype=MPI_INT,
                                             .displacement=0, .count=1};
      for (size_t j = (size_t)i; j < (size_t)nsend; ++j)
        send_msgs[j] = (struct Xt_msg_param){.rank=(1+(int)j)%comm_size,
                                             .datatype=MPI_INT,
                                             .displacement=0, .count=1};
    }
    struct Xt_msg_param recv_msgs[2] =
      {{.rank=-1, .datatype=MPI_DATATYPE_NULL, .displacement=0, .count=1 },
       {.rank=i, .datatype=MPI_INT, .displacement=0, .count=1}};

    Xt_exchanger exchanger = xt_exchanger_new(nsend, nrecv,
                                              send_msgs,
                                              recv_msgs + (my_rank != i),
                                              comm, 0, config);

    // test
    for (int async = 0; async < 1 + test_async; ++async) {

      int src_data[1] = { my_rank == i ? 4711 : -1 };
      int dst_data[1] = { my_rank == i ? 4711 : -1 };

      data_exchange(exchanger, src_data, dst_data, async);

      if (dst_data[0] != 4711) PUT_ERR("invalid data\n");
    }

    // cleanup
    free(send_msgs);
    xt_exchanger_delete(exchanger);
  }
}

static void
test_gather(MPI_Comm comm, Xt_config config)
{
  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  // gather pattern
  // prepare datatypes outside of loop for load-balance
  MPI_Datatype *dt_by_ofs = xmalloc((size_t)comm_size * sizeof (*dt_by_ofs));
  dt_by_ofs[0] = MPI_INT;
  for (size_t j = 1; j < (size_t)comm_size; ++j)
  {
    MPI_Type_indexed(1, (int[]){1}, (int[]){(int)j}, MPI_INT, dt_by_ofs + j);
    MPI_Type_commit(dt_by_ofs + j);
  }
  int *dst_data = xmalloc(((size_t)comm_size - 1) * sizeof (*dst_data) * 2);
  int incr = comm_size/test_freq + (comm_size < test_freq);
  int test_async = (xt_config_get_exchange_method(config)
                    != xt_exchanger_irecv_send);
  for (int i = 0; i < comm_size; i += incr) {

    // setup
    int nsend = i != my_rank;
    size_t nrecv = 0;


    struct Xt_msg_param send_msgs[1]
      = {{ .rank=i, .datatype=MPI_INT, .displacement=0, .count=1}};
    struct Xt_msg_param *recv_msgs = NULL;
    if (my_rank == i) {
      nrecv = (size_t)comm_size - 1;
      recv_msgs = xmalloc(nrecv * sizeof (*recv_msgs) * 2);
      for (size_t j = 0; j < nrecv; ++j) {
        recv_msgs[j].rank = (i + (int)j + 1)%comm_size;
        recv_msgs[j].datatype = dt_by_ofs[j];
        recv_msgs[j].displacement = 0;
        recv_msgs[j].count=1;
        recv_msgs[j + nrecv].rank
          = (comm_size - (int)j - 1 + i)%comm_size;
        recv_msgs[j + nrecv].datatype = dt_by_ofs[j];
        recv_msgs[j + nrecv].displacement = 0;
        recv_msgs[j + nrecv].count=1;
      }
    }

    enum { exchange_fwd, exchange_rev, num_exchanges };
    Xt_exchanger exchanger[num_exchanges];
    for (size_t j = 0; j < num_exchanges; ++j)
      exchanger[j] = xt_exchanger_new(nsend, (int)nrecv, send_msgs,
                                      recv_msgs + j * nrecv, comm, 0, config);

    // test
    for (int async = 0; async < 1 + test_async; ++async) {
      int src_data[1] = {(my_rank+comm_size-i)%comm_size};
      for (size_t j = 0; j < 2 * nrecv; ++j)
        dst_data[j] = -1;

      for (size_t j = 0; j < num_exchanges; ++j)
        data_exchange(exchanger[j], src_data, dst_data + j * nrecv,
                      async);

      bool mismatch = false;
      for (size_t j = 0; j < nrecv; ++j)
        mismatch |= (((size_t)dst_data[j] != j + 1)
                     || (dst_data[nrecv + j] != comm_size - (int)j - 1));
      if (mismatch)
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
test_all2all(MPI_Comm comm, Xt_config config)
{
  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  // all-to-all pattern
  // setup
  int nsend = comm_size - 1;
  int nrecv = comm_size - 1;

  struct Xt_msg_param *send_msgs = xmalloc((size_t)nsend * sizeof (*send_msgs));
  for (size_t j = 0; j < (size_t)nsend; ++j) {
    send_msgs[j].datatype = MPI_INT;
    send_msgs[j].displacement = 0;
    send_msgs[j].count = 1;
  }
  struct Xt_msg_param *recv_msgs
    = xmalloc((size_t)nrecv * sizeof (*recv_msgs));
  int *dst_data = xmalloc((size_t)comm_size * sizeof (*dst_data));
  MPI_Datatype *dt_by_ofs = xmalloc((size_t)comm_size * sizeof (*dt_by_ofs));
  dt_by_ofs[0] = my_rank != 0 ? MPI_INT : MPI_DATATYPE_NULL;
  for (size_t i = 1; i < (size_t)comm_size; ++i)
    if (i != (size_t)my_rank) {
      MPI_Type_indexed(1, (int[]){1}, (int[]){(int)i}, MPI_INT, dt_by_ofs + i);
      MPI_Type_commit(dt_by_ofs + i);
    } else
      dt_by_ofs[i] = MPI_DATATYPE_NULL;
  size_t comm_size_ = (size_t)comm_size,
    incr = (size_t)(comm_size/test_freq + (comm_size < test_freq));
  int test_async = (xt_config_get_exchange_method(config)
                    != xt_exchanger_irecv_send);
  for (size_t i = 0; i < comm_size_ - 1; i += incr) {
    for (size_t j = 0; j < (size_t)nsend; ++j) {
      int ofs = my_rank + 1 + (int)i + (int)j;
      send_msgs[j].rank = (ofs + (ofs >= comm_size + my_rank))%comm_size;
      send_msgs[j].displacement = 0;
    }
    for (size_t j = 0; j < comm_size_ - 1; j += incr) {
      for (size_t k = 0; k < (size_t)nrecv; ++k) {
        int ofs = ((int)i + (int)j + (int)k)%(comm_size - 1);
        ofs += ofs >= my_rank;
        recv_msgs[k].rank = ofs;
        recv_msgs[k].datatype = dt_by_ofs[ofs]; // replace by MPI_INT
        recv_msgs[k].displacement = 0; // replace by sizeof (int) * ofs
        recv_msgs[k].count = 1;
      }
      Xt_exchanger exchanger
        = xt_exchanger_new(nsend, nrecv, send_msgs, recv_msgs, comm, 0, config);
      // test
      for (int async = 0; async < 1 + test_async; ++async) {
        int src_data[1] = {my_rank};
        for (size_t k = 0; k < comm_size_; ++k)
          dst_data[k] = my_rank;

        data_exchange(exchanger, src_data, dst_data, async);

        bool mismatch = false;
        for (int k = 0; k < comm_size; ++k)
          mismatch |= (dst_data[k] != k);
        if (mismatch)
          PUT_ERR("invalid data\n");
      }

      // cleanup
      xt_exchanger_delete(exchanger);
    }
  }
  for (size_t i = 1; i < comm_size_; ++i)
    if (i != (size_t)my_rank)
      MPI_Type_free(dt_by_ofs + i);
  free(dt_by_ofs);
  free(dst_data);
  free(recv_msgs);
  free(send_msgs);
}

static void
test_rr(MPI_Comm comm, Xt_config config)
{
  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);
  // round robin pattern
  int incr = comm_size/test_freq + (comm_size < test_freq);
  int test_async = (xt_config_get_exchange_method(config)
                    != xt_exchanger_irecv_send);
  for (int i = 1; i < comm_size; i += incr) {

    // setup
    enum { nsend = 1, nrecv = 1 };
    struct Xt_msg_param send_msgs[nsend]
      = {{.rank=(my_rank + i)%comm_size, .datatype=MPI_INT,
            .displacement = 0, .count=1}};
    struct Xt_msg_param recv_msgs[nrecv]
      = {{.rank=(my_rank + comm_size - i)%comm_size, .datatype=MPI_INT,
            .displacement=0, .count=1}};

    Xt_exchanger exchanger
      = xt_exchanger_new(nsend, nrecv, send_msgs, recv_msgs, comm, 0, config);

    // test
    for (int async = 0; async < 1 + test_async; ++async) {
      int src_data[1] = {my_rank};
      int dst_data[1] = {-1};

      data_exchange(exchanger, src_data, dst_data, async);

      if (dst_data[0] != (my_rank + comm_size - i)%comm_size)
        PUT_ERR("invalid data\n");
    }

    // test queries
    if (i == 1) {
      Xt_exchanger_iter msg_iter
        = xt_exchanger_get_iterator(exchanger, config);
      if (msg_iter)
        do {
          MPI_Datatype dt
            = xt_exchanger_iter_get_MPI_Datatype(msg_iter);
          if (dt != MPI_INT)
            PUT_ERR("unexpected datatype\n");
          int msg_rank
            = xt_exchanger_iter_get_rank(msg_iter);
          int msg_count
            = xt_exchanger_iter_get_count(msg_iter);
          if (msg_count != 1)
            PUT_ERR("unexpected count\n");
          int direction
            = xt_exchanger_iter_get_direction(msg_iter);
          if (direction != 0 && direction != 1)
            PUT_ERR("unexpected direction\n");
          if ((direction == 0 && msg_rank != (my_rank + i)%comm_size)
              || (direction == 1
                  && msg_rank != (my_rank + comm_size - i)%comm_size))
            PUT_ERR("unexpected communication\n");
        } while (xt_exchanger_iter_next(msg_iter));
      if (xt_exchanger_iter_get_direction(msg_iter) != -1)
        PUT_ERR("unexpected direction\n");
      xt_exchanger_iter_delete(msg_iter);
    }
    // cleanup
    xt_exchanger_delete(exchanger);
  }
}

static void
test_intercomm_all2all(MPI_Comm comm, Xt_config config)
{
  // inter-communicator's are not defined for virtual topologies, which are
  // used by xt_exchanger_neigh_alltoall_new
  if (xt_config_get_exchange_method(config)
      == xt_exchanger_neigh_alltoall) return;

  int my_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(comm, &my_rank), comm);
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  if (comm_size == 1) return;

  // generate intercomm
  int splitRank = (comm_size * 2) / 3;
  int group = my_rank >= splitRank;
  MPI_Comm intra_group_comm;
  MPI_Comm inter_comm; // split communicator with 2 to 1 ratio
  xt_mpi_call(MPI_Comm_split(comm, group, 0, &intra_group_comm), comm);
  xt_mpi_call(MPI_Intercomm_create(intra_group_comm, 0, comm,
                                   group ? 0 : splitRank,
                                   0, &inter_comm), comm);

  int intra_rank;
  int local_size, remote_size;
  xt_mpi_call(MPI_Comm_rank(inter_comm, &intra_rank), comm);
  xt_mpi_call(MPI_Comm_size(inter_comm, &local_size), comm);
  xt_mpi_call(MPI_Comm_remote_size(inter_comm, &remote_size), comm);

  // all-to-all pattern
  // setup
  int nsend = remote_size;
  int nrecv = remote_size;

  struct Xt_msg_param *send_msgs =
    xmalloc((size_t)nsend * sizeof (*send_msgs));
  for (int i = 0; i < nsend; ++i) {
    send_msgs[i].rank = i;
    send_msgs[i].datatype = MPI_INT;
    send_msgs[i].displacement = 0;
    send_msgs[i].count = 1;
  }

  struct Xt_msg_param *recv_msgs =
    xmalloc((size_t)nrecv * sizeof (*recv_msgs));
  for (int i = 0; i < nrecv; ++i) {
    recv_msgs[i].rank = i;
    recv_msgs[i].datatype = MPI_INT;
    recv_msgs[i].displacement = (MPI_Aint)sizeof (int) * i;
    recv_msgs[i].count = 1;
  }

  int *dst_data = xmalloc((size_t)nrecv * sizeof (*dst_data));

  Xt_exchanger exchanger
    = xt_exchanger_new(nsend, nrecv, send_msgs, recv_msgs,
                       inter_comm, 0, config);
  // test
  int test_async = (xt_config_get_exchange_method(config)
                    != xt_exchanger_irecv_send);
  for (int async = 0; async < 1 + test_async; ++async) {

    int src_data[1] = {my_rank};
    for (int i = 0; i < nrecv; ++i) dst_data[i] = -1;

    data_exchange(exchanger, src_data, dst_data, async);

    int dst_data_offset = (my_rank >= splitRank)?0:splitRank;
    bool mismatch = false;
    for (int i = 0; i < nrecv; ++i)
      mismatch |= (dst_data[i] != i + dst_data_offset);
    if (mismatch)
      PUT_ERR("invalid data\n");
  }

  // cleanup
  xt_exchanger_delete(exchanger);

  free(dst_data);
  free(recv_msgs);
  free(send_msgs);

  xt_mpi_call(MPI_Comm_free(&inter_comm), comm);
  xt_mpi_call(MPI_Comm_free(&intra_group_comm), comm);
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
