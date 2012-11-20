/**
 * @file test_xmap_all2all.c
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

#include <mpi.h>

#include <yaxt.h>

#define VERBOSE
#include "tests.h"

static void
shift_idx(Xt_int idx[], int num, int offset)  {
  int rank, i;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for (i=0; i<num; i++) {
    idx[i] = idx[i] + rank*offset;
  }
}

int main(void) {

   // init mpi
   xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

   xt_initialize(MPI_COMM_WORLD);

   int my_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   {
      // source index list
      Xt_int src_index_list[] = {1,2,3,4,5,6,7};
      Xt_int src_num_indices
        = sizeof(src_index_list) / sizeof(src_index_list[0]);
      shift_idx(src_index_list, src_num_indices, 7);

      Xt_idxlist src_idxlist;

      src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

      // destination index list
      Xt_int dst_index_list[] = {7,6,5,4,3,2,1};
      Xt_int dst_num_indices
        = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
      shift_idx(dst_index_list, dst_num_indices, 7);

      Xt_idxlist dst_idxlist;

      dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

      // test of exchange map
      Xt_xmap xmap;

      xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

      // test results

      if (xt_xmap_get_num_destinations(xmap) != 1)
        PUT_ERR("error in xt_xmap_get_num_destinations\n");

      if (xt_xmap_get_num_sources(xmap) != 1)
        PUT_ERR("error in xt_xmap_get_num_destinations\n");

      int rank;

      xt_xmap_get_destination_ranks(xmap, &rank);
      if (rank != my_rank)
        PUT_ERR("error in xt_xmap_get_destination_ranks\n");

      xt_xmap_get_source_ranks(xmap, &rank);
      if (rank != my_rank)
        PUT_ERR("error in xt_xmap_get_source_ranks\n");

      // clean up
      xt_xmap_delete(xmap);
      xt_idxlist_delete(src_idxlist);
      xt_idxlist_delete(dst_idxlist);
   }

   {
      // source index list
      Xt_int src_index_list[] = {5,67,4,5,13,9,2,1,0,96,13,12,1,3};
      Xt_int src_num_indices
        = sizeof(src_index_list) / sizeof(src_index_list[0]);
      shift_idx(src_index_list, src_num_indices, 100);

      Xt_idxlist src_idxlist;

      src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

      // destination index list
      Xt_int dst_index_list[] = {41,5,4,7,3,96,1,5,48,54,4,5,4,6,3,13,7,2,1};
      Xt_int dst_num_indices
        = sizeof(dst_index_list) / sizeof(dst_index_list[0]);
      shift_idx(dst_index_list, dst_num_indices, 100);

      Xt_idxlist dst_idxlist;

      dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

      Xt_xmap xmap;

      xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

      // test results

      // the global id intersection should look as follows:
      // 1,1,2,3,3,4,4,4,5,5,5,13,96

      if (xt_xmap_get_num_destinations(xmap) != 1)
        PUT_ERR("error in xt_xmap_get_num_destinations\n");

      if (xt_xmap_get_num_sources(xmap) != 1)
        PUT_ERR("error in xt_xmap_get_num_destinations\n");

      int rank;

      xt_xmap_get_destination_ranks(xmap, &rank);
      if (rank != my_rank)
        PUT_ERR("error in xt_xmap_get_destination_ranks\n");

      xt_xmap_get_source_ranks(xmap, &rank);
      if (rank != my_rank)
        PUT_ERR("error in xt_xmap_get_source_ranks\n");

      // clean up
      xt_xmap_delete(xmap);
      xt_idxlist_delete(src_idxlist);
      xt_idxlist_delete(dst_idxlist);
   }

   MPI_Finalize();

   return TEST_EXIT_CODE;
}
