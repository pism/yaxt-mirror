/**
 * @file test_redist_collection.c
 *
 * @copyright Copyright  (C)  2012 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "test_redist_common.h"
#include "core/ppm_xfuncs.h"

static void
test_repeated_redist(int cache_size);

static void
test_displacement_variations(void);


int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // general test with one redist
    // set up data
    enum { nvalues = 5, nselect = (nvalues + 1)/2 };
    Xt_xmap xmap = build_odd_selection_xmap(nvalues);

    Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    xt_xmap_delete(xmap);

    // generate redist_collection

    Xt_redist redist_coll
      = xt_redist_collection_new(&redist, 1, -1, MPI_COMM_WORLD);

    // test communicator of redist

    if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_coll),
                                     MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    xt_redist_delete(redist);

    // test exchange

    static const double src_data[nvalues] = {1,2,3,4,5};
    double dst_data[nselect];
    static const double ref_dst_data[nselect] = {1,3,5};

    check_redist(redist_coll, src_data, nselect,
                 dst_data, fill_array_double, NULL,
                 ref_dst_data, MPI_DOUBLE, MPI_DOUBLE);

    Xt_redist redist_coll_copy = xt_redist_copy(redist_coll);

    check_redist(redist_coll_copy, src_data, nselect,
                 dst_data, fill_array_double, NULL,
                 ref_dst_data, MPI_DOUBLE, MPI_DOUBLE);
    // clean up

    xt_redist_delete(redist_coll_copy);
    xt_redist_delete(redist_coll);
  }

  { // test empty redist
    Xt_idxlist src_idxlist = xt_idxempty_new();
    Xt_idxlist dst_idxlist = xt_idxempty_new();
    Xt_xmap xmap =
      xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);

    Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    xt_xmap_delete(xmap);

    // generate redist_collection

    Xt_redist redist_coll
      = xt_redist_collection_new(&redist, 1, -1, MPI_COMM_WORLD);

    // test communicator of redist
    if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_coll),
                                     MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    xt_redist_delete(redist);

    // test exchange

    static const double src_data[1] = {-1};
    double dst_data[1] = {-2};

    xt_redist_s_exchange1(redist_coll, src_data, dst_data);

    static const double ref_dst_data[1] = {-2};

    if (ref_dst_data[0] != dst_data[0])
      PUT_ERR("error in xt_redist_s_exchange\n");

    Xt_redist redist_coll_copy = xt_redist_copy(redist_coll);
    dst_data[0] = -2;

    xt_redist_s_exchange1(redist_coll_copy, src_data, dst_data);

    if (ref_dst_data[0] != dst_data[0])
      PUT_ERR("error in xt_redist_s_exchange\n");

    // clean up
    xt_redist_delete(redist_coll_copy);
    xt_redist_delete(redist_coll);
  }

  // test with one redist used three times (with two different input data
  // displacements -> test of cache) (with default cache size)
  // set up data
  test_repeated_redist(-1);

  // test with one redist used three times (with two different input data
  // displacements -> test of cache) (with cache size == 0)
  // set up data
  test_repeated_redist(0);

  test_displacement_variations();

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
test_repeated_redist(int cache_size)
{
  Xt_xmap xmap = build_odd_selection_xmap(5);

  Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

  xt_xmap_delete(xmap);

  // generate redist_collection
  Xt_redist redists[3] = {redist, redist, redist};
  Xt_redist redist_coll
    = xt_redist_collection_new(redists, 3, cache_size, MPI_COMM_WORLD);

  // test communicator of redist

  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_coll),
                                   MPI_COMM_WORLD))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");

  xt_redist_delete(redist);

  // test exchange
  for (int sync_mode = 0; sync_mode < 2; ++sync_mode)
  {
    static const double src_data[3][5]
      = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

    static const void *const src_data_p[3]
      = {src_data[0],src_data[1],src_data[2]};
    void *dst_data_p[3] = {dst_data[0],dst_data[1],dst_data[2]};

    exchange_func_ptr exchange_func
      = sync_mode == 0 ? xt_redist_s_exchange : wrap_a_exchange;
    exchange_func(redist_coll, 3, (const void **)src_data_p, dst_data_p);

    static const double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        if (ref_dst_data[i][j] != dst_data[i][j])
          PUT_ERR("error in xt_redist_s_exchange\n");
  }

  // test exchange with changed displacements
  for (int sync_mode = 0; sync_mode < 2; ++sync_mode)
  {
    static const double src_data[3][5]
      = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

    static const void *const src_data_p[3]
      = {src_data[1],src_data[0],src_data[2]};
    void *dst_data_p[3] = {dst_data[1],dst_data[0],dst_data[2]};

    exchange_func_ptr exchange_func
      = sync_mode == 0 ? xt_redist_s_exchange : wrap_a_exchange;
    exchange_func(redist_coll, 3, (const void **)src_data_p, dst_data_p);

    static const double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        if (ref_dst_data[i][j] != dst_data[i][j])
          PUT_ERR("error in xt_redist_s_exchange\n");
  }

  // test exchange with original displacements
  for (int sync_mode = 0; sync_mode < 2; ++sync_mode)
  {
    static const double src_data[3][5]
      = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    double dst_data[3][3] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

    static const void *const src_data_p[3]
      = {src_data[0],src_data[1],src_data[2]};
    void *dst_data_p[3] = {dst_data[0],dst_data[1],dst_data[2]};

    exchange_func_ptr exchange_func
      = sync_mode == 0 ? xt_redist_s_exchange : wrap_a_exchange;
    exchange_func(redist_coll, 3, (const void **)src_data_p, dst_data_p);

    static const double ref_dst_data[3][3] = {{1,3,5},{6,8,10},{11,13,15}};

    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        if (ref_dst_data[i][j] != dst_data[i][j])
          PUT_ERR("error in xt_redist_s_exchange\n");
  }

  // clean up

  xt_redist_delete(redist_coll);
}

enum { num_redists = 3 };
enum { nvalues = 5, nselect = nvalues/2+(nvalues&1) };

static void
run_displacement_check(Xt_redist redist_coll, int sync)
{
  static const double src_data[num_redists][nvalues]
    = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
  double dst_data[num_redists][nselect] = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};

  enum { cache_size = 16, cache_overrun = 2 };

  double src_data_[nvalues + cache_size + cache_overrun],
    dst_data_[nselect + cache_size + cache_overrun];

  const void *src_data_p[num_redists] = {src_data[0],src_data[1],NULL};
  void *dst_data_p[num_redists] = {dst_data[0],dst_data[1],NULL};

  for (size_t k = 0; k < cache_size + cache_overrun; ++k) {

    memcpy(src_data_+k, src_data[2], 5 * sizeof(*src_data_));
    for (size_t i = 0; i < num_redists; ++i)
      for (size_t j = 0; j < nselect; ++j)
        dst_data[i][j] = -1;
    memcpy(dst_data_+k, dst_data[2], 3 * sizeof(*dst_data_));

    src_data_p[2] = src_data_+k;
    dst_data_p[2] = dst_data_+k;

    exchange_func_ptr exchange_func
      = sync ? xt_redist_s_exchange : wrap_a_exchange;
    exchange_func(redist_coll, num_redists, src_data_p, dst_data_p);

    static const double ref_dst_data[num_redists][nselect]
      = {{1,3,5},{6,8,10},{11,13,15}};

    for (size_t i = 0; i < num_redists; ++i)
      for (size_t j = 0; j < nselect; ++j)
        if (ref_dst_data[i][j] != ((double *)dst_data_p[i])[j])
          PUT_ERR("error in xt_redist_s_exchange\n");
  }
}


static void
test_displacement_variations(void)
{
  // test with one redist used three times (with different input
  // data displacements until the cache is full)
  // set up data
  Xt_xmap xmap = build_odd_selection_xmap(nvalues);

  Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

  xt_xmap_delete(xmap);

  // generate redist_collection

  Xt_redist redists[num_redists] = {redist, redist, redist};

  Xt_redist redist_coll
    = xt_redist_collection_new(redists, num_redists, -1, MPI_COMM_WORLD);

  // test communicator of redist

  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_coll),
                                   MPI_COMM_WORLD))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");

  xt_redist_delete(redist);


  // test exchange
  run_displacement_check(redist_coll, 0);
  run_displacement_check(redist_coll, 1);
  Xt_redist redist_coll_copy = xt_redist_copy(redist_coll);
  run_displacement_check(redist_coll_copy, 0);
  run_displacement_check(redist_coll_copy, 1);

  // clean up
  xt_redist_delete(redist_coll_copy);
  xt_redist_delete(redist_coll);
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
