/**
 * @file test_id_buckets_parallel.c
 *
 * @copyright Copyright  (C)  2019 Jörg Behrens <behrens@dkrz.de>
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

//#define VERBOSE

#include <mpi.h>
#include <string.h>
#include <yaxt.h>
#include <stdint.h>
#include "core/ppm_xfuncs.h"
#include "../src/xt_id_buckets.h"
#include "../src/xt_memcmp.h"
#include "tests.h"

static void check_buckets(void * ids, size_t id_size, size_t num_ids,
                          Xt_memcmp_func compar_ids,
                          int check_equal_distribution) {

  // sort the ids
  void * sorted_ids = xmalloc(id_size * num_ids);
  memcpy(sorted_ids, ids, id_size * num_ids);
  qsort(sorted_ids, num_ids, id_size, compar_ids);

  int comm_size;
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  // generate buckets
  Xt_id_buckets buckets = xt_id_buckets_new(sorted_ids, id_size, num_ids,
                                            (size_t)comm_size, MPI_COMM_WORLD);

  void * ubounds = xmalloc(id_size * (size_t)comm_size);

  Xt_id_buckets_iter iter = xt_id_buckets_get_iter(buckets);
  void * ubound = xmalloc(id_size);
  memset(ubound, -1, id_size);

  uint64_t * bucket_sizes =
    xmalloc(2 * (size_t)comm_size * sizeof(*bucket_sizes));
  uint64_t * global_bucket_sizes = bucket_sizes + (size_t)comm_size;

  size_t ids_idx = 0, bucket_sum = 0;
  for (int i = 0; i < comm_size; ++i) {

    // get ubound
    xt_id_buckets_iter_get_next_ubound(iter, ubound);

    // get copy of ubound
    memcpy((char*)ubounds + (size_t)i * id_size, ubound, id_size);

    // count number of ids in current bucket
    size_t curr_num_ids = 0;
    while ((ids_idx < num_ids) &&
           (memcmp((char*)sorted_ids + ids_idx * id_size,
                   ubound, id_size) <= 0)) {
      ++curr_num_ids;
      ++ids_idx;
    }
    bucket_sum += curr_num_ids;
    bucket_sizes[i] = (uint64_t)curr_num_ids;
  }
  if (bucket_sum != num_ids) PUT_ERR("unmatched ids");
  xt_mpi_call(
    MPI_Allreduce(bucket_sizes, global_bucket_sizes, comm_size,
                  MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD), MPI_COMM_WORLD);

  // ensure that all processes have the same ubounds
  {
    void * temp_ubounds = xmalloc(id_size * (size_t)comm_size);
    xt_mpi_call(
      MPI_Allreduce(ubounds, temp_ubounds, (int)id_size * comm_size,
                    MPI_CHAR, MPI_MAX, MPI_COMM_WORLD), MPI_COMM_WORLD);
    for (size_t i = 0; i < id_size * (size_t)comm_size; ++i)
      if (((char*)temp_ubounds)[i] != ((char*)ubounds)[i])
        PUT_ERR("ubounds do not match between processes");
    free(temp_ubounds);
  }

  // ensure that no process gets too many points
  if (check_equal_distribution) {

    uint64_t global_num_ids = 0;
    for (int i = 0; i < comm_size; ++i)
      global_num_ids += global_bucket_sizes[i];
    size_t prev_split_point = 0;
    for (int i = 0; i < comm_size; ++i) {

      size_t curr_split_point =
        (size_t)(((long long)global_num_ids * (long long)(i + 1)) /
                 (long long)comm_size);

      size_t max_bucket_size = ((curr_split_point - prev_split_point) * 5) / 4;
      prev_split_point = curr_split_point;

      if ((size_t)(global_bucket_sizes[i]) > max_bucket_size)
        PUT_ERR("bucket is too big");
    }
  }

  free(bucket_sizes);
  free(ubound);
  free(ubounds);
  xt_id_buckets_iter_delete(iter);
  xt_id_buckets_delete(buckets);
  free(sorted_ids);
}

int main(void) {

  // init mpi
  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &comm_size), MPI_COMM_WORLD);

  { // empty buckets on all processes
    char ids[] = {-1};
    size_t id_size = sizeof(ids[0]);
    size_t num_ids = 0;
    check_buckets(ids, id_size, num_ids, xt_memcmp_get_func(id_size), 1);
  }

  { // 10 ids on all processes
    double ids[10];
    size_t id_size = sizeof(ids[0]);
    size_t num_ids = sizeof(ids)/id_size;
    for (int i = 0; i < num_ids; ++i)
      ids[i] = (double)(comm_rank * num_ids + i);
    check_buckets(ids, id_size, num_ids, xt_memcmp_get_func(id_size), 1);
  }

  { // 10 ids on all processes (all ids are identical
    double ids[10];
    size_t id_size = sizeof(ids[0]);
    size_t num_ids = sizeof(ids)/id_size;
    for (int i = 0; i < num_ids; ++i) ids[i] = 0.0d;
    check_buckets(ids, id_size, num_ids, xt_memcmp_get_func(id_size), 0);
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
