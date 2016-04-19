/**
 * @file test_redist_p2p.c
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

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "test_redist_common.h"

static void
run_redist_test(Xt_redist redist, const void *src_data, size_t elem_size,
                size_t num_values, const void *restrict ref_values);

int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  // offset-free test:
  {
    // source index list
    Xt_int src_index_list[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
    enum { src_num_indices
           = sizeof(src_index_list) / sizeof(src_index_list[0]) };

    Xt_idxlist src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    // destination index list
    static const Xt_int dst_index_list[] = {5,4,3,96,1,5,4,5,4,3,13,2,1};
    enum { dst_num_indices
           = sizeof(dst_index_list) / sizeof(dst_index_list[0]) };

    Xt_idxlist dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

    // xmap
    Xt_xmap xmap
      = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // redist_p2p
    Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    // test communicator of redist

    if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist),
                                     MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    // test exchange
    static const double src_data[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
    static const double ref_dst_data[] = {0,2,13,9,7,0,2,0,2,13,4,6,7};
    enum { num_ref_values = sizeof(ref_dst_data) / sizeof(ref_dst_data[0]) };
    run_redist_test(redist, src_data, sizeof (src_data[0]), num_ref_values,
                    ref_dst_data);
    Xt_redist redist_copy = xt_redist_copy(redist);
    xt_redist_delete(redist);
    run_redist_test(redist_copy, src_data, sizeof (src_data[0]), num_ref_values,
                    ref_dst_data);
    // clean up

    xt_redist_delete(redist_copy);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);
  }

  // check offsets
  {
    // source index list
    static const Xt_int src_index_list[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
    enum { src_num = sizeof(src_index_list) / sizeof(src_index_list[0]) };

    Xt_idxlist src_idxlist = xt_idxvec_new(src_index_list, src_num);

    // destination index list
    static const Xt_int dst_index_list[] = {5,4,3,96,1,5,4,5,4,3,13,2,1};
    enum { dst_num = sizeof(dst_index_list) / sizeof(dst_index_list[0]) };

    Xt_idxlist dst_idxlist = xt_idxvec_new(dst_index_list, dst_num);

    // xmap
    Xt_xmap xmap
      = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // redist_p2p with offsets

    int src_pos[src_num];
    int dst_pos[dst_num];

    for (int i = 0; i < src_num; i++)
      src_pos[i] = i;
    for (int i = 0; i < dst_num; i++)
      dst_pos[i] = dst_num-1-i;

    Xt_redist redist
      = xt_redist_p2p_off_new(xmap, src_pos, dst_pos, MPI_DOUBLE);

    // test communicator of redist

    if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist),
                                     MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    // test exchange

    static const double src_data[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
    assert(sizeof(src_data)/sizeof(src_data[0]) == src_num);
    static const double ref_dst_data[] = {7,6,4,13,2,0,2,0,7,9,13,2,0};
    run_redist_test(redist, src_data, sizeof (src_data[0]), dst_num,
                    ref_dst_data);
    Xt_redist redist_copy = xt_redist_copy(redist);
    xt_redist_delete(redist);
    run_redist_test(redist_copy, src_data, sizeof (src_data[0]), dst_num,
                    ref_dst_data);
    // clean up

    xt_redist_delete(redist_copy);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);
  }

  // check offset extents
  {
    // source index list
    static const Xt_int src_index_list[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
    enum { src_num = sizeof(src_index_list) / sizeof(src_index_list[0]) };

    Xt_idxlist src_idxlist = xt_idxvec_new(src_index_list, src_num);

    // destination index list
    static const Xt_int dst_index_list[] = {5,4,3,96,1,5,4,5,4,3,13,2,1};
    enum { dst_num = sizeof(dst_index_list) / sizeof(dst_index_list[0]) };

    Xt_idxlist dst_idxlist = xt_idxvec_new(dst_index_list, dst_num);

    // xmap
    Xt_xmap xmap
      = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // redist_p2p with extents of offsets
    struct Xt_offset_ext
      src_pos[1] = { { .start = 0, .size = src_num, .stride =  1 } },
      dst_pos[1] = { { .start = dst_num - 1, .size = dst_num, .stride = -1 } };

    Xt_redist redist
      = xt_redist_p2p_ext_new(xmap,
                              sizeof (src_pos)/sizeof (src_pos[0]), src_pos,
                              sizeof (dst_pos)/sizeof (dst_pos[0]), dst_pos,
                              MPI_LONG);

    // test communicator of redist

    if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist),
                                     MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    // test exchange
    static const long src_data[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
    assert(sizeof(src_data)/sizeof(src_data[0]) == src_num);
    static const long ref_dst_data[] = {7,6,4,13,2,0,2,0,7,9,13,2,0};
    run_redist_test(redist, src_data, sizeof (src_data[0]), dst_num,
                    ref_dst_data);
    Xt_redist redist_copy = xt_redist_copy(redist);
    xt_redist_delete(redist);
    run_redist_test(redist_copy, src_data, sizeof (src_data[0]), dst_num,
                    ref_dst_data);

    // clean up
    xt_redist_delete(redist_copy);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
run_redist_test(Xt_redist redist, const void *src_data, size_t elem_size,
                size_t num_values, const void *restrict ref_values)
{
  size_t dst_size = num_values * elem_size;
  void *restrict dst_data = malloc(dst_size);
  if (dst_size && !dst_data) {
    perror("could not allocate destination buffer");
    abort();
  }
  memset(dst_data, 0xff, dst_size);
  const void *src_data_p[1] = { src_data };
  void *dst_data_p[1] = { dst_data };
  xt_redist_s_exchange(redist, 1, src_data_p, dst_data_p);
  if (memcmp(ref_values, dst_data, dst_size))
    PUT_ERR("error in xt_redist_s_exchange\n");

  free(dst_data);
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
