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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "ctest_common.h"
#include "test_redist_common.h"
#include "core/ppm_xfuncs.h"

static void
simple_test(MPI_Comm comm, Xt_config config);

static void
test_empty_redist(MPI_Comm comm, Xt_config config);

static void
test_repeated_redist(MPI_Comm comm, Xt_config config, int cache_size);

static void
test_displacement_variations(MPI_Comm comm, Xt_config config,
                             size_t num_redists);

static void
test_collection_of_collection(MPI_Comm comm, Xt_config config);


int main(int argc, char **argv) {

  test_init_mpi(&argc, &argv, MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  Xt_config config = redist_exchanger_option(&argc, &argv);

  simple_test(MPI_COMM_WORLD, config);

  test_empty_redist(MPI_COMM_WORLD, config);

  // test with one redist used three times (with two different input data
  // displacements -> test of cache) (with default cache size)
  // set up data
  test_repeated_redist(MPI_COMM_WORLD, config, -1);

  // test with one redist used three times (with two different input data
  // displacements -> test of cache) (with cache size == 0)
  // set up data
  test_repeated_redist(MPI_COMM_WORLD, config, 0);

  test_displacement_variations(MPI_COMM_WORLD, config, 3);
  test_displacement_variations(MPI_COMM_WORLD, config, 40);

  test_collection_of_collection(MPI_COMM_WORLD, config);

  xt_config_delete(config);
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
simple_test(MPI_Comm comm, Xt_config config)
{ // general test with one redist
  // set up data
  enum { nvalues = 5, nselect = (nvalues + 1)/2 };
  Xt_xmap xmap = build_odd_selection_xmap(nvalues, comm);

  Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

  xt_xmap_delete(xmap);

  // generate redist_collection

  Xt_redist redist_coll
    = xt_redist_collection_custom_new(&redist, 1, -1, comm, config);

  // test communicator of redist
  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_coll), comm))
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

static void
test_empty_redist(MPI_Comm comm, Xt_config config)
{ // test empty redist
  Xt_idxlist src_idxlist = xt_idxempty_new();
  Xt_idxlist dst_idxlist = xt_idxempty_new();
  Xt_xmap xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, comm);

  xt_idxlist_delete(src_idxlist);
  xt_idxlist_delete(dst_idxlist);

  Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

  xt_xmap_delete(xmap);

  // generate redist_collection
  Xt_redist redist_coll
    = xt_redist_collection_custom_new(&redist, 1, -1, comm, config);

  // test communicator of redist
  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_coll), comm))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");

  xt_redist_delete(redist);

  // test exchange
  static const double src_data[1] = {0};
  double dst_data[1];
  static const double ref_dst_data[1] = {-2};
  check_redist(redist_coll, src_data, 1, dst_data, fill_array_double,
               &(double){-2}, ref_dst_data, MPI_DOUBLE, MPI_DOUBLE);


  Xt_redist redist_coll_copy = xt_redist_copy(redist_coll);
  check_redist(redist_coll_copy, src_data, 1, dst_data, fill_array_double,
               &(double){-2}, ref_dst_data, MPI_DOUBLE, MPI_DOUBLE);

  // clean up
  xt_redist_delete(redist_coll_copy);
  xt_redist_delete(redist_coll);
}


static void
test_repeated_redist(MPI_Comm comm, Xt_config config, int cache_size)
{
  enum { num_slice = 3,
         src_slice_len = 5, dst_slice_len = (src_slice_len+1)/2 };
  Xt_xmap xmap = build_odd_selection_xmap(src_slice_len, comm);

  MPI_Datatype dt[3] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
  Xt_redist redist = xt_redist_p2p_new(xmap, dt[0]);

  xt_xmap_delete(xmap);

  // generate redist_collection
  Xt_redist redists[num_slice] = {redist, redist, redist};
  Xt_redist redist_coll
    = xt_redist_collection_custom_new(redists, num_slice, cache_size, comm,
                                      config);

  // test communicator of redist
  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_coll), comm))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");

  xt_redist_delete(redist);

  static const size_t dst_num_elem[num_slice]
    = { dst_slice_len, dst_slice_len, dst_slice_len };
  static const prepare_dst prep[num_slice]
    = { fill_array_double, fill_array_double, fill_array_double };

  // test exchange
  {
    static const double src_data[num_slice][src_slice_len]
      = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    static const double ref_dst_data[num_slice][dst_slice_len]
      = {{1,3,5},{6,8,10},{11,13,15}};
    double dst_data[num_slice][dst_slice_len];

    static const void *const src_data_p[num_slice]
      = { src_data[0], src_data[1], src_data[2]};
    void *const dst_data_p[num_slice]
      = { dst_data[0], dst_data[1], dst_data[2] };
    const void *const ref_dst_data_p[num_slice]
      = { ref_dst_data[0], ref_dst_data[1], ref_dst_data[2] };
    check_redist_coll(redist_coll, sync_mode_test_all,
                      num_slice, src_data_p,
                      dst_num_elem, dst_data_p,
                      prep, NULL, ref_dst_data_p, dt, dt);
  }
  // test exchange with changed displacements
  {
    static const double src_data[num_slice][src_slice_len]
      = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    static const double ref_dst_data[num_slice][dst_slice_len]
      = {{1,3,5},{6,8,10},{11,13,15}};
    double dst_data[num_slice][dst_slice_len];
    static const void *const src_data_p[num_slice]
      = {src_data[1],src_data[0],src_data[2]};
    void *const dst_data_p[num_slice]
      = { dst_data[1], dst_data[0], dst_data[2] };
    const void *const ref_dst_data_p[num_slice]
      = { ref_dst_data[0], ref_dst_data[1], ref_dst_data[2] };

    check_redist_coll(redist_coll, sync_mode_test_all,
                      num_slice, src_data_p,
                      dst_num_elem, dst_data_p,
                      prep, NULL, ref_dst_data_p, dt, dt);
  }

  // test exchange with original displacements
  {
    static const double src_data[num_slice][src_slice_len]
      = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    static const double ref_dst_data[num_slice][dst_slice_len]
      = {{1,3,5},{6,8,10},{11,13,15}};
    double dst_data[num_slice][3];
    static const void *const src_data_p[num_slice]
      = { src_data[0], src_data[1], src_data[2] };
    void *const dst_data_p[num_slice]
      = { dst_data[0], dst_data[1], dst_data[2] };
    const void *const ref_dst_data_p[num_slice]
      = { ref_dst_data[0], ref_dst_data[1], ref_dst_data[2] };

    check_redist_coll(redist_coll, sync_mode_test_all,
                      num_slice, src_data_p,
                      dst_num_elem, dst_data_p,
                      prep, NULL, ref_dst_data_p, dt, dt);
  }

  // clean up

  xt_redist_delete(redist_coll);
}

enum {
  nvalues = 5,
  nselect = nvalues/2+(nvalues&1),
  cache_size = 16,
  cache_overrun = 2,
  num_src_elems = nvalues + cache_size + cache_overrun,
};

static void
run_displacement_check(Xt_redist redist_coll, int sync,
                       size_t num_redists)
{
  size_t num_dst_elems = num_redists * nselect + cache_size + cache_overrun;
  double (*src_data)[nvalues]
    = xmalloc(num_redists * sizeof (*src_data));
  for (size_t j = 0; j < num_redists; ++j)
    for (size_t i = 0; i < nvalues; ++i)
      src_data[j][i] = (double)(1 + j * nvalues + i);

  MPI_Datatype dt[num_redists];
  for (size_t j = 0; j < num_redists; ++j)
    dt[j] = MPI_DOUBLE;
  double src_data_[num_src_elems];

  double *dst_data
    = xmalloc(num_dst_elems * sizeof (*dst_data));

  const void *src_data_p[num_redists];
  void *dst_data_p[num_redists];
  for (size_t j = 0; j < num_redists-1; ++j) {
    src_data_p[j] = src_data[j];
    dst_data_p[j] = dst_data+nselect*j;
  }

  double *ref_dst_data
    = xmalloc(num_dst_elems * sizeof (*ref_dst_data));
  const void *ref_dst_data_p[num_redists];
  for (size_t j = 0; j < num_redists-1; ++j) {
    ref_dst_data_p[j] = ref_dst_data + j*nselect;
    for (size_t i = 0; i < nselect; ++i)
      ref_dst_data[j*nselect + i] = src_data[j][i*2];
  }

  for (size_t k = 0; k < cache_size + cache_overrun; ++k) {

    memcpy(src_data_+k, src_data[num_redists-1],
           nvalues * sizeof(*src_data_));
    for (size_t i = 0; i < num_dst_elems; ++i)
      dst_data[i] = -1;

    src_data_p[num_redists-1] = src_data_ + k;
    dst_data_p[num_redists-1] = dst_data + nselect*(num_redists-1) + k;
    ref_dst_data_p[num_redists-1] = ref_dst_data + nselect*(num_redists-1) + k;

    /* put every second value from src_data[3] into ref_dst_data
     * starting at (num_redists-1)*nselect + k, i.e. with offset k
     * vs. a contiguous transformation */
    for (size_t i = 0; i < k; ++i)
      ref_dst_data[(num_redists-1)*nselect + i] = -1;
    for (size_t i = 0; i < nselect; ++i)
      ref_dst_data[(num_redists-1)*nselect + k + i]
        = src_data[num_redists-1][i*2];
    for (size_t i = num_redists*nselect + k; i < num_dst_elems; ++i)
      ref_dst_data[i] = -1;
    size_t num_dst_elems_by_redist[num_redists];
    prepare_dst prep[num_redists];
    for (size_t i = 0; i < num_redists; ++i) {
      num_dst_elems_by_redist[i] = nselect;
      prep[i] = fill_array_double;
    }
    check_redist_coll(redist_coll,
                      sync ? sync_mode_test_s : sync_mode_test_a,
                      (int)num_redists, src_data_p, num_dst_elems_by_redist,
                      dst_data_p,
                      prep, NULL, ref_dst_data_p, dt, dt);
  }
  free(ref_dst_data);
  free(dst_data);
  free(src_data);
}


static void
test_displacement_variations(MPI_Comm comm, Xt_config config,
                             size_t num_redists)
{
  // test with one redist used three times (with different input
  // data displacements until the cache is full)
  // set up data
  Xt_xmap xmap = build_odd_selection_xmap(nvalues, comm);

  Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

  xt_xmap_delete(xmap);

  // generate redist_collection

  Xt_redist redists[num_redists];
  for (size_t j = 0; j < num_redists; ++j)
    redists[j] = redist;

  Xt_redist redist_coll
    = xt_redist_collection_custom_new(redists, (int)num_redists, -1,
                                      comm, config);

  // test communicator of redist

  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_coll), comm))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");

  xt_redist_delete(redist);


  // test exchange
  run_displacement_check(redist_coll, 0, num_redists);
  run_displacement_check(redist_coll, 1, num_redists);
  Xt_redist redist_coll_copy = xt_redist_copy(redist_coll);
  run_displacement_check(redist_coll_copy, 0, num_redists);
  run_displacement_check(redist_coll_copy, 1, num_redists);

  // clean up
  xt_redist_delete(redist_coll_copy);
  xt_redist_delete(redist_coll);
}


static void
test_collection_of_collection(MPI_Comm comm, Xt_config config)
{
  // general test with one redist
  // set up data
  enum {
    num_redist_coll = 2,
    num_p2p_redist = 2,
  };
  Xt_xmap xmap = build_odd_selection_xmap(nvalues, comm);

  Xt_redist redists[num_p2p_redist];
  redists[0]
    = xt_redist_p2p_new(xmap, MPI_DOUBLE);
  redists[1]
    = xt_redist_p2p_new(xmap, MPI_INT);

  xt_xmap_delete(xmap);

  // generate redist_collection
  Xt_redist redist_coll[num_redist_coll];
  redist_coll[0]
    = xt_redist_collection_custom_new(redists, 2, -1, comm, config);
  for (size_t i = 1; i < num_redist_coll; ++i)
    redist_coll[i] = redist_coll[0];

  xt_redist_delete(redists[0]);
  xt_redist_delete(redists[1]);

  Xt_redist redist_coll_coll
    = xt_redist_collection_custom_new(redist_coll, 2, 2, comm, config);

  // test communicator of redist
  if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_coll_coll), comm))
    PUT_ERR("error in xt_redist_get_MPI_Comm\n");

  xt_redist_delete(redist_coll[0]);

  // test exchange
  static const struct src_databuf {
    double data_dbl[nvalues];
    int data_int[nvalues];
  } src_data = { .data_dbl = {1,2,3,4,5}, .data_int = {1,2,3,4,5} };
  struct dst_databuf {
    double data_dbl[nselect];
    int data_int[nselect];
  } dst_data[num_redist_coll];
  static const struct dst_databuf ref_dst_data
    = { {1,3,5}, {1,3,5} };
  const void *const src_data_p[num_redist_coll*num_p2p_redist]
    = { src_data.data_dbl, src_data.data_int,
        src_data.data_dbl, src_data.data_int };
  static const size_t dst_num_elem[num_redist_coll*num_p2p_redist]
    = { nselect, nselect, nselect, nselect };
  void *const dst_data_p[num_redist_coll*num_p2p_redist]
    = { dst_data[0].data_dbl, dst_data[0].data_int,
        dst_data[1].data_dbl, dst_data[1].data_int };
  static const prepare_dst prep[num_redist_coll*num_p2p_redist]
    = { fill_array_double, fill_array_int,
        fill_array_double, fill_array_int };
  const void *ref_dst_data_p[num_redist_coll*num_p2p_redist]
    = { ref_dst_data.data_dbl, ref_dst_data.data_int,
        ref_dst_data.data_dbl, ref_dst_data.data_int };
  MPI_Datatype dt[num_redist_coll*num_p2p_redist]
    = { MPI_DOUBLE, MPI_INT, MPI_DOUBLE, MPI_INT };
  check_redist_coll(redist_coll_coll, sync_mode_test_all,
                    num_redist_coll * num_p2p_redist,
                    src_data_p, dst_num_elem,
                    dst_data_p, prep, NULL,
                    ref_dst_data_p, dt, dt);

  /* check the same with the same data at different relative offset */
  static const struct src_databuf_displaced {
    double data_dbl[nvalues+1];
    int data_int[nvalues];
  } src_data_displaced = { .data_dbl = {1,2,3,4,5,-2},
                           .data_int = {1,2,3,4,5} };
  const void *const src_data_displaced_p[num_redist_coll*num_p2p_redist]
    = { src_data_displaced.data_dbl, src_data_displaced.data_int,
        src_data_displaced.data_dbl, src_data_displaced.data_int };
  check_redist_coll(redist_coll_coll, sync_mode_test_all,
                    num_redist_coll * num_p2p_redist,
                    src_data_displaced_p, dst_num_elem,
                    dst_data_p, prep, NULL,
                    ref_dst_data_p, dt, dt);

  Xt_redist redist_coll_coll_copy = xt_redist_copy(redist_coll_coll);

  check_redist_coll(redist_coll_coll_copy, sync_mode_test_all,
                    num_redist_coll * num_p2p_redist,
                    src_data_p, dst_num_elem,
                    dst_data_p, prep, NULL,
                    ref_dst_data_p, dt, dt);
  // clean up
  xt_redist_delete(redist_coll_coll_copy);
  xt_redist_delete(redist_coll_coll);
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
