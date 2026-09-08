/**
 * @file test_idxempty.c
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

#include <stdlib.h>
#include <limits.h>
#include <assert.h>

#include <mpi.h>

#include <yaxt.h>

#include "tests.h"
#include "ctest_common.h"
#include "test_idxlist_utils.h"
#include "core/ppm_xfuncs.h"

int main(int argc, char **argv) {

  test_init_mpi(&argc, &argv, MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // test packing and unpacking

    Xt_idxlist idxempty = xt_idxempty_new();

    check_idxlist(idxempty, NULL, 0);

    Xt_idxlist idxempty_copy
      = idxlist_pack_unpack_copy(idxempty);

    // check the received index vector

    check_idxlist(idxempty_copy, NULL, 0);

    // compute intersection between the two index vectors

    Xt_idxlist intersection
      = xt_idxlist_get_intersection(idxempty, idxempty_copy);

    // check the computed intersection (should be identically to the original
    // index vector)

    check_idxlist(intersection, NULL, 0);

    // check the conversion to stripes
    {
      struct Xt_stripe *stripes;
      int num_stripes;

      xt_idxlist_get_index_stripes(idxempty, &stripes, &num_stripes);

      if (stripes != NULL || num_stripes != 0)
        PUT_ERR("ERROR: xt_idxlist_get_index_stripes\n");
      free(stripes);
    }

    // check the computation of a bounding box

    enum { ndim = 3 };
    Xt_int global_size[ndim];
    Xt_int global_start_index = 0;
    struct Xt_bounds bounds[ndim];

    for (size_t i = 0; i < ndim; ++i)
      global_size[i] = 10;

    xt_idxlist_get_bounding_box(idxempty, ndim, global_size,
                                global_start_index, bounds);

    for (size_t i = 0; i < ndim; ++i)
      if (bounds[i].size != 0)
        PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

    // check that position and index queries return the expected
    // (missing) values
    {
      Xt_int idx = (Xt_int)0xdead;
      int pos, rc;
      if ((rc = xt_idxlist_get_index_at_position(idxempty, 0, &idx)) != 1
          || idx != (Xt_int)0xdead)
        PUT_ERR("retrieved index position from empty list, rc=%d, idx=%lld!\n",
                rc, (long long)idx);
      int pos_vec[] = { -1, 0, 2 };
      enum {
        num_pos = sizeof (pos_vec) / sizeof (pos_vec[0]),
        undef_val = -5,
      };
      Xt_int idx_vec[num_pos];
      for (size_t i = 0; i < num_pos; ++i)
        idx_vec[i] = 7;
      if ((rc = xt_idxlist_get_indices_at_positions(
             idxempty, pos_vec, num_pos, idx_vec, undef_val)) != num_pos)
        PUT_ERR("retrieved index position from empty list, rc=%d!\n", rc);
      for (size_t i = 0; i < num_pos; ++i)
        if (idx_vec[i] != undef_val)
          PUT_ERR("incorrectly initialized index at"
                  " position %zu from empty list to %lld!\n",
                  i, (long long)idx_vec[i]);
      pos = 0xbeef;
      if ((rc = xt_idxlist_get_position_of_index(idxempty, 0, &pos)) != 1
          || pos != -1)
        PUT_ERR("incorrect position query result: rc=%d, pos=%d\n", rc, pos);
      for (size_t i = 0; i < num_pos; ++i) {
        idx_vec[i] = (Xt_int)(i^(size_t)0xbad);
        pos_vec[i] = -((int)(i+1)^0xdead);
      }
      if ((rc = xt_idxlist_get_positions_of_indices(
             idxempty, idx_vec, num_pos, pos_vec, 0)) != num_pos)
        PUT_ERR("empty index list does not reject all indices, rc=%d\n",
                rc);
      for (size_t i = 0; i < num_pos; ++i)
        if (idx_vec[i] != (Xt_int)(i^(size_t)0xbad)
            || pos_vec[i] != -((int)(i+1)^0xdead))
          PUT_ERR("encountered unexpected result at idx_vec[%zu]=%lld,"
                  " pos_vec[%zu]=%d\n", i, (long long)idx_vec[i],
                  i, pos_vec[i]);
    }
    // clean up
    xt_idxlist_delete(idxempty);
    xt_idxlist_delete(idxempty_copy);
    xt_idxlist_delete(intersection);
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
