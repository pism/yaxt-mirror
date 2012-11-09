/**
 * @file test_redist_p2p_parallel.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
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

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"

static inline int imin(int a, int b) {
  return  (((a)<(b))?(a):(b));
}

static inline int imax(int a, int b) {
  return  (((a)>(b))?(a):(b));
}

int main(void) {

  // init mpi

  int rank, size;

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);
  xt_mpi_call(MPI_Comm_size(MPI_COMM_WORLD, &size), MPI_COMM_WORLD);

  {
    Xt_count i;
    // source index list
    Xt_int src_index_list[10];
    Xt_count src_num_indices = 10;
    for (i = 0; i < src_num_indices; ++i)
      src_index_list[i] = rank * 10 + i;

    Xt_idxlist src_idxlist;

    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    // destination index list
    Xt_int dst_index_list[10];
    Xt_count dst_num_indices = 10;
    for (i = 0; i < dst_num_indices; ++i)
      dst_index_list[i] = (rank * 10 + i + 2) % (size * 10);

    Xt_idxlist dst_idxlist;

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

    // xmap
    Xt_xmap xmap;

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // redist_p2p
    Xt_redist redist;

    redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    // test exchange

    double src_data[10];
    double dst_data[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

    for (i = 0; i < src_num_indices; ++i)
      src_data[i] = rank * 10 + i;

    double * src_data_p = &src_data[0];
    double * dst_data_p = &dst_data[0];

    xt_redist_s_exchange(redist, (void**)(&src_data_p), 1,
                         (void**)(&dst_data_p), 1);

    //unsigned i;

    for (i = 0; i < sizeof(dst_data) / sizeof(dst_data[0]); ++i)
      if (dst_index_list[i] != dst_data[i])
        PUT_ERR("error in xt_redist_s_exchange\n");

    // clean up

    xt_redist_delete(redist);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);
  }

  // test nonuniform numbers of send and receive partners

  {
    Xt_count i;
    // source index list
    Xt_int src_index_list[size];
    Xt_count src_num_indices = 0;

    if (rank == 0) src_num_indices = size;
    for (i = 0; i < src_num_indices; ++i)
      src_index_list[i] = i;

    Xt_idxlist src_idxlist;

    src_idxlist = xt_idxvec_new(src_index_list, src_num_indices);

    // destination index list
    Xt_int dst_index_list[size];
    Xt_count dst_num_indices = size;
    for (i = 0; i < dst_num_indices; ++i)
      dst_index_list[i] = i;

    Xt_idxlist dst_idxlist;

    dst_idxlist = xt_idxvec_new(dst_index_list, dst_num_indices);

    // xmap
    Xt_xmap xmap;

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, MPI_COMM_WORLD);

    // redist_p2p
    Xt_redist redist;

    redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);

    // test exchange

    double src_data[size];
    double dst_data[size];
    for (i = 0; i < size; ++i) {
      if (rank == 0)
        src_data[i] = i;
      else
        src_data[i] = -2;
      dst_data[i] = -1;
    }

    double * src_data_p = &src_data[0];
    double * dst_data_p = &dst_data[0];

    xt_redist_s_exchange(redist, (void**)(&src_data_p), 1,
                         (void**)(&dst_data_p), 1);

    //unsigned i;

    for (i = 0; i < size; ++i)
      if (dst_data[i] != i)
        PUT_ERR("error in xt_redist_s_exchange\n");

    // clean up

    xt_redist_delete(redist);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(src_idxlist);
    xt_idxlist_delete(dst_idxlist);
  }

  // test redist with blocks:
  {
    // the global index domain (1dim problem):
    int ngdom = 2*size;
    int gdoma[ngdom]; // start state (index distribution) of global domain
    int gdomb[ngdom]; // end state ""
    int gsurfdata[ngdom];
    int gdepth[ngdom]; // think: ocean depth of an one dim. ocean
    int gvol_size; // volume of deep ocean
    int ig2col_off[ngdom]; // offset of surface data within vol

    gvol_size = 0;
    for (int i = 0; i < ngdom; i++) {
      gdoma[i] = i;
      gdomb[i] = ngdom-1-i;
      gsurfdata[i] = 100+i;
      gdepth[i] = i+1;
      ig2col_off[i] = gvol_size;
      gvol_size += gdepth[i];
    }

    int nwin = ngdom/size; // my local window size of the global surface domain
    // start of my window within global index domain (== global offset)
    int ig0 = rank*nwin;
    if (nwin*size != ngdom) PUT_ERR("internal error\n");

    // local index
    Xt_int iveca[nwin], ivecb[nwin];
    for (int i = 0; i < nwin; i++) {
      int ig = ig0+i;
      iveca[i]= gdoma[ig];
      ivecb[i]= gdomb[ig];
    }

    Xt_idxlist idxlist_a = xt_idxvec_new(iveca, nwin);
    Xt_idxlist idxlist_b = xt_idxvec_new(ivecb, nwin);

    Xt_xmap xmap = xt_xmap_all2all_new(idxlist_a, idxlist_b, MPI_COMM_WORLD);

    // simple redist:
    Xt_redist redist = xt_redist_p2p_new(xmap, MPI_INT);
    int a_surfdata[nwin];
    int b_surfdata[nwin];
    int b_surfdata_ref[nwin];
    for (int i = 0; i < nwin; i++) {
      a_surfdata[i] = gsurfdata[iveca[i]];
      b_surfdata[i] = -1;
      b_surfdata_ref[i] = gsurfdata[ivecb[i]];
    }

    xt_redist_s_exchange1(redist, a_surfdata, b_surfdata);
    for (int i = 0; i < nwin; ++i) {
      if (b_surfdata[i] != b_surfdata_ref[i])
        PUT_ERR("error in xt_redist_s_exchange\n");
    }

    xt_redist_delete(redist);

    // generate global volume data
    int gvoldata[gvol_size];
    for (int i = 0; i < ngdom; i++) {
      for (int j = 0; j <  gdepth[i]; j++) {
        int p = ig2col_off[i] + j;
        gvoldata[p] = i*100 + j;
      }
    }

    // generate blocks

    int src_block_offsets[nwin];
    int src_block_sizes[nwin];
    int dst_block_offsets[nwin];
    int dst_block_sizes[nwin];
    int a_vol_size = 0; // state a volume of my proc
    int b_vol_size = 0; // state b volume of my proc
    // we only need local size but simply oversize here
    int a_voldata[gvol_size];
    int b_voldata[gvol_size]; // ..
    int b_voldata_ref[gvol_size]; // ..

    for (int i = 0; i < gvol_size; i++) {
      a_voldata[i] = -1;
      b_voldata[i] = -1;
      b_voldata_ref[i] = -1;
    }

    {
      int qa=0;
      for (int i = 0; i < nwin; i++) {
        int ia = iveca[i];
        if (i==0)
          src_block_offsets[i] = 0;
        else
          src_block_offsets[i] = src_block_offsets[i-1]+src_block_sizes[i-1];
        src_block_sizes[i]   = gdepth[ia];
        for (int j = 0; j < gdepth[ia]; j++) {
          int p = ig2col_off[ia] + j;
          a_voldata[qa] = gvoldata[p];
          qa++;
        }
        a_vol_size += src_block_sizes[i];
      }
    }

    {
      int qb=0;
      for (int i = 0; i < nwin; i++) {
        int ib = ivecb[i];
        if (i==0)
          dst_block_offsets[i] = 0;
        else
          dst_block_offsets[i] = dst_block_offsets[i-1]+dst_block_sizes[i-1];
        dst_block_sizes[i]   = gdepth[ib];
        for (int j = 0; j < gdepth[ib]; j++) {
          int p = ig2col_off[ib] + j;
          b_voldata_ref[qb] = gvoldata[p];
          qb++;
        }
        b_vol_size += dst_block_sizes[i];
      }
    }

    // redist with blocks:
    Xt_redist block_redist;
    block_redist
      = xt_redist_p2p_blocks_off_new(xmap,
                                     src_block_offsets, src_block_sizes, nwin,
                                     dst_block_offsets, dst_block_sizes, nwin,
                                     MPI_INT);

    xt_redist_s_exchange1(block_redist, a_voldata, b_voldata);

    for (int i = 0; i < gvol_size; i++) {
      if (b_voldata[i] != b_voldata_ref[i])
        PUT_ERR("error in xt_redist_s_exchange (1) for volume data\n");

    }

    // redist with blocks but without explicit offsets:
    Xt_redist block_redist2;
    block_redist2 = xt_redist_p2p_blocks_new(xmap,
                                             src_block_sizes, nwin,
                                             dst_block_sizes, nwin,
                                             MPI_INT);

    for (int i = 0; i < gvol_size; i++) {
      b_voldata[i] = -1;
    }

    xt_redist_s_exchange1(block_redist2, a_voldata, b_voldata);

    for (int i = 0; i < gvol_size; i++) {
      if (b_voldata[i] != b_voldata_ref[i])
        PUT_ERR("error in xt_redist_s_exchange (2) for volume data\n");

    }

    xt_redist_delete(block_redist);
    xt_redist_delete(block_redist2);
    xt_xmap_delete(xmap);
    xt_idxlist_delete(idxlist_a);
    xt_idxlist_delete(idxlist_b);

    // end of test
  }

  MPI_Finalize();

 return TEST_EXIT_CODE;
}
