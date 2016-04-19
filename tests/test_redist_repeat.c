/**
 * @file test_redist_repeat.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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

#include <string.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "test_redist_common.h"


int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // general test with one redist
    // set up data
    enum { src_slice_len = 5, dst_slice_len = 3 };
    Xt_xmap xmap = build_odd_selection_xmap(src_slice_len);
    Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);
    xt_xmap_delete(xmap);

    // generate redist_repeatection

    static const double src_data[] = {1,2,3,4,5};
    double dst_data[] = {-1,-1,-1};
    Xt_redist redist_repeat;
    MPI_Aint src_extent = (MPI_Aint)(sizeof(src_data)),
      dst_extent = (MPI_Aint)(sizeof(dst_data));
    int displacements[1] = {0};

    redist_repeat = xt_redist_repeat_new(redist, src_extent, dst_extent, 1,
                                         displacements);

    // test communicator of redist

    if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_repeat),
                                     MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    xt_redist_delete(redist);

    // test exchange
    static const double ref_dst_data[] = {1,3,5};

    check_redist(redist_repeat, src_data,
                 sizeof(dst_data), dst_data, ref_dst_data);

    // clean up
    xt_redist_delete(redist_repeat);
  }

  {
    // test with one redist used three times (two exchanges)
    // set up data
    enum { src_slice_len = 5, dst_slice_len = (src_slice_len+1)/2 };
    Xt_xmap xmap = build_odd_selection_xmap(src_slice_len);
    Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);
    xt_xmap_delete(xmap);

    // generate redist_repeatection

    enum { num_repetitions = 3 };
    Xt_redist redist_repeat;
    static const double src_data[num_repetitions][src_slice_len]
      = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
    double dst_data[num_repetitions][dst_slice_len];
    MPI_Aint src_extent =
      (MPI_Aint)((size_t)(src_data[1] - src_data[0]) * sizeof (double));
    MPI_Aint dst_extent =
      (MPI_Aint)((size_t)(dst_data[1] - dst_data[0]) * sizeof (double));
    int displacements[num_repetitions] = {0, 1, 2};

    redist_repeat = xt_redist_repeat_new(redist, src_extent, dst_extent,
                                         num_repetitions, displacements);

    // test communicator of redist

    if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_repeat),
                                     MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    xt_redist_delete(redist);

    // test exchange
    {
      for (size_t j = 0; j < num_repetitions; ++j)
        for (size_t i = 0; i < dst_slice_len; ++i)
          dst_data[j][i] = -1.0;

      static const double ref_dst_data[num_repetitions][dst_slice_len]
        = {{1,3,5},{6,8,10},{11,13,15}};
      check_redist(redist_repeat, src_data,
                   sizeof (dst_data), dst_data, ref_dst_data);
    }

    {
      for (size_t j = 0; j < num_repetitions; ++j)
        for (size_t i = 0; i < dst_slice_len; ++i)
          dst_data[j][i] = -1.0;
      static const double ref_dst_data[num_repetitions][dst_slice_len]
        = {{1,3,5},{6,8,10},{11,13,15}};
      check_redist(redist_repeat, src_data,
                   sizeof (dst_data), dst_data, ref_dst_data);
    }

    // clean up

    xt_redist_delete(redist_repeat);
  }

  { // test with one redist used three times with gaps between the levels
    enum { src_slice_len = 5, dst_slice_len = (src_slice_len+1)/2 };
    Xt_xmap xmap = build_odd_selection_xmap(src_slice_len);
    Xt_redist redist = xt_redist_p2p_new(xmap, MPI_DOUBLE);
    xt_xmap_delete(xmap);

    // generate redist_repeatection
    MPI_Aint src_extent = (MPI_Aint)(5 * sizeof (double));
    MPI_Aint dst_extent = (MPI_Aint)(3 * sizeof (double));
    int displacements[3] = {0, 2, 4};
    Xt_redist redist_repeat
      = xt_redist_repeat_new(redist, src_extent, dst_extent, 3, displacements);

    // test communicator of redist

    if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_repeat),
                                     MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    xt_redist_delete(redist);

    // test exchange
    {
      static const double src_data[5][5]
        = {{1,2,3,4,5},{-1},{6,7,8,9,10},{-1},{11,12,13,14,15}};
      double dst_data[5][3]
        = {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};
      static const double ref_dst_data[5][3]
        = {{1,3,5},{-1,-1,-1},{6,8,10},{-1,-1,-1},{11,13,15}};
      check_redist(redist_repeat, src_data,
                   sizeof (dst_data), dst_data, ref_dst_data);
    }

    // clean up

    xt_redist_delete(redist_repeat);
  }


  { // test of redist_repeat generated with overlapping repetition of input redist
    enum { npt = 9 };
    Xt_xmap xmap = build_odd_selection_xmap(6);
    static const int src_pos[npt] = {1,2,3,4,5,6,7,8,9};
    static const int dst_pos[npt] = {0,2,4,6,8,10,12,14,16};
    Xt_redist redist = xt_redist_p2p_off_new(xmap, src_pos, dst_pos, MPI_DOUBLE);

    xt_xmap_delete(xmap);

    // init data
    double src_data[npt], dst_data[npt];
    static const double ref_dst_data[npt]
        = {102.0, 103.0, 104.0, 105.0, 106.0, 107.0, -1.0, -1.0, -1.0};
    for(int i = 0; i<npt; i++) {
      src_data[i] = 100+i+1;
      dst_data[i] = -1;
    }

    // generate reference destination data: {102, 103, ..., 107, -1, -1, ... }
    xt_redist_s_exchange1(redist, src_data, dst_data);
    xt_redist_s_exchange1(redist, src_data+1, dst_data+1);

    // explicit check of ref data:
    for(int i=0; i<npt; i++)
      if (dst_data[i] != ref_dst_data[i]) {
        PUT_ERR("error in xt_redist_repeat\n");
        break;
      }
    for(int i = 0; i<npt; i++)
      dst_data[i] = -1;

    // generate redist_repeat
    MPI_Aint src_extent = (MPI_Aint)(1 * sizeof (double));
    MPI_Aint dst_extent = (MPI_Aint)(1 * sizeof (double));
    int displacements[2] = {0,1};
    Xt_redist redist_repeat
      = xt_redist_repeat_new(redist, src_extent, dst_extent, 2, displacements);

    // test communicator of redist
    if (!communicators_are_congruent(xt_redist_get_MPI_Comm(redist_repeat),
                                     MPI_COMM_WORLD))
      PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    // test exchange
    check_redist(redist_repeat, src_data,
                 sizeof (dst_data), dst_data, ref_dst_data);

    // clean up
    xt_redist_delete(redist);
    xt_redist_delete(redist_repeat);

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
