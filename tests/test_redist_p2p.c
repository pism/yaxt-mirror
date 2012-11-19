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
#include <assert.h>
#include <mpi.h>

#include <yaxt.h>

#include "tests.h"

int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  // offset-free test:
  {
    // source index list
    Xt_int src_index_list[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
    Xt_count src_num_indices
      = sizeof(src_index_list) / sizeof(src_index_list[0]);

    Xt_idxlist src_idxlist;

    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    // destination index list
    Xt_int dst_index_list[] = {41,5,4,7,3,96,1,5,48,54,4,5,4,6,3,13,7,2,1};
    Xt_count dst_num_indices
      = sizeof(dst_index_list) / sizeof(dst_index_list[0]);

    Xt_idxlist dst_idxlist;

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

    // xmap
    Xt_xmap xmap;

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // redist_p2p
    Xt_redist redist;

    redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    // test exchange

    double src_data[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
    double dst_data[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                         -1,-1,-1,-1,-1,-1,-1,-1,-1};

    double * src_data_p = &src_data[0];
    double * dst_data_p = &dst_data[0];

    xt_redist_s_exchange(redist, (void**)(&src_data_p), 1,
                         (void**)(&dst_data_p), 1);

    double ref_dst_data[] = {-1,0,2,-1,13,9,7,0,-1,-1,2,0,2,-1,13,4,-1,6,7};

    unsigned i;

    for (i = 0; i < sizeof(dst_data) / sizeof(dst_data[0]); ++i)
      if (ref_dst_data[i] != dst_data[i])
        PUT_ERR("error in xt_redist_s_exchange\n");

    // clean up

    xt_redist_delete(redist);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);
  }

  // check offsets
  {
    // source index list
    Xt_int src_index_list[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
    unsigned src_num = sizeof(src_index_list) / sizeof(src_index_list[0]);

    Xt_idxlist src_idxlist;

    src_idxlist = xt_idxvec_new(src_index_list, src_num);

    // destination index list
    Xt_int dst_index_list[] = {41,5,4,7,3,96,1,5,48,54,4,5,4,6,3,13,7,2,1};
    unsigned dst_num = sizeof(dst_index_list) / sizeof(dst_index_list[0]);

    Xt_idxlist dst_idxlist;

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num);

    // xmap
    Xt_xmap xmap;

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // redist_p2p with offsets

    unsigned i;

    int src_pos[src_num];
    int dst_pos[dst_num];

    for (i = 0; i < src_num; i++)
      src_pos[i] = i;
    for (i = 0; i < dst_num; i++)
      dst_pos[i] = dst_num-1-i;

    Xt_redist redist;
    redist = xt_redist_p2p_off_new(xmap, src_pos, dst_pos, MPI_DOUBLE);

    // test exchange

    double src_data[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
    double dst_data[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                         -1,-1,-1,-1,-1,-1,-1,-1,-1};
    assert(sizeof(src_data)/sizeof(src_data[0]) == src_num);
    assert(sizeof(dst_data)/sizeof(dst_data[0]) == dst_num);

    double * src_data_p = &src_data[0];
    double * dst_data_p = &dst_data[0];

    xt_redist_s_exchange(redist, (void**)(&src_data_p), 1,
                         (void**)(&dst_data_p), 1);

    double ref_dst_data[] = {-1,0,2,-1,13,9,7,0,-1,-1,
                             2,0,2,-1,13,4,-1,6,7};

    for (i = 0; i < dst_num; ++i)
      if (ref_dst_data[dst_pos[i]] != dst_data[i])
        PUT_ERR("error in xt_redist_s_exchange\n");

    // clean up

    xt_redist_delete(redist);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
