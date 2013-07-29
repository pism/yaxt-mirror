/**
 * @file test_idxstripes.c
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

#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#include <mpi.h>

#include <yaxt.h>

#include "core/ppm_xfuncs.h"
#include "tests.h"
#include "test_idxlist_utils.h"

static void
do_tests(Xt_idxlist idxlist, Xt_int * ref_indices, int num_indices);

int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // general tests
    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[3] = {{0,5,1},{10,5,1},{20,5,1}};
    int num_stripes = 3;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    Xt_int ref_indices[15] = {0,1,2,3,4, 10,11,12,13,14, 20,21,22,23,24};
    int ref_num_indices = 15;

    do_tests(idxstripes, ref_indices, ref_num_indices);

    xt_idxlist_delete(idxstripes);
  }

  { // general tests
    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[3] = {{0,5,1},{10,5,2},{20,5,3}};
    int num_stripes = 3;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    Xt_int ref_indices[15] = {0,1,2,3,4, 10,12,14,16,18, 20,23,26,29,32};
    int ref_num_indices = 15;

    do_tests(idxstripes, ref_indices, ref_num_indices);

    xt_idxlist_delete(idxstripes);
  }

  { // general tests
    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[2] = {{0,5,6},{1,5,3}};
    int num_stripes = 2;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    Xt_int ref_indices[10] = {0,6,12,18,24, 1,4,7,10,13};
    int ref_num_indices = 10;

    do_tests(idxstripes, ref_indices, ref_num_indices);

    xt_idxlist_delete(idxstripes);
  }

  { // general tests
    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[2] = {{0,5,-1},{1,5,1}};
    int num_stripes = 2;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    Xt_int ref_indices[10] = {0,-1,-2,-3,-4, 1,2,3,4,5};
    int ref_num_indices = 10;

    do_tests(idxstripes, ref_indices, ref_num_indices);

    xt_idxlist_delete(idxstripes);
  }

  { // general tests
    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[2] = {{9,5,-2},{0,5,2}};
    int num_stripes = 2;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    Xt_int ref_indices[10] = {9,7,5,3,1, 0,2,4,6,8};
    int ref_num_indices = 10;

    do_tests(idxstripes, ref_indices, ref_num_indices);

    xt_idxlist_delete(idxstripes);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[2] = {{0,4,1},{6,4,1}};
    struct Xt_stripe stripes_b[1] = {{1,8,1}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 2);
    idxstripes_b = xt_idxstripes_new(stripes_b, 1);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_int ref_indices[6] = {1,2,3, 6,7,8};
    int ref_num_indices = 6;

    do_tests(intersection[0], ref_indices, ref_num_indices);
    do_tests(intersection[1], ref_indices, ref_num_indices);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[3] = {{0,4,1},{6,4,1}, {11,4,1}};
    struct Xt_stripe stripes_b[2] = {{1,7,1}, {9,5,1}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 3);
    idxstripes_b = xt_idxstripes_new(stripes_b, 2);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_int ref_indices[9] = {1,2,3, 6,7, 9, 11,12,13};
    int ref_num_indices = 9;

    do_tests(intersection[0], ref_indices, ref_num_indices);
    do_tests(intersection[1], ref_indices, ref_num_indices);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[2] = {{0,3,1},{8,3,1}};
    struct Xt_stripe stripes_b[2] = {{3,5,1}, {11,3,1}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 2);
    idxstripes_b = xt_idxstripes_new(stripes_b, 2);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_int *ref_indices = NULL;
    int ref_num_indices = 0;

    do_tests(intersection[0], ref_indices, ref_num_indices);
    do_tests(intersection[1], ref_indices, ref_num_indices);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[1] = {{0,10,1}};
    struct Xt_stripe stripes_b[2] = {{0,5,2}, {9,5,-2}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 1);
    idxstripes_b = xt_idxstripes_new(stripes_b, 2);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_int ref_indices[10] = {0,1,2,3,4,5,6,7,8,9};
    int ref_num_indices = 10;

    do_tests(intersection[0], ref_indices, ref_num_indices);
    do_tests(intersection[1], ref_indices, ref_num_indices);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[2] = {{0,5,3},{1,5,7}};
    struct Xt_stripe stripes_b[2] = {{0,7,2}, {24,10,-1}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 2);
    idxstripes_b = xt_idxstripes_new(stripes_b, 2);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_int ref_indices[6] = {0,6,8,12,15,22};
    int ref_num_indices = 6;

    do_tests(intersection[0], ref_indices, ref_num_indices);
    do_tests(intersection[1], ref_indices, ref_num_indices);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[1] = {{0,10,1}};
    struct Xt_stripe stripes_b[2] = {{5,5,1}, {4,5,-1}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 1);
    idxstripes_b = xt_idxstripes_new(stripes_b, 2);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_int ref_indices[10] = {0,1,2,3,4,5,6,7,8,9};
    int ref_num_indices = 10;

    do_tests(intersection[0], ref_indices, ref_num_indices);
    do_tests(intersection[1], ref_indices, ref_num_indices);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[2] = {{0,10,1}, {20,5,1}};
    struct Xt_stripe stripes_b[2] = {{3,5,1},  {17,5,1}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 2);
    idxstripes_b = xt_idxstripes_new(stripes_b, 2);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_int ref_indices[7] = {3,4,5,6,7,20,21};
    int ref_num_indices = 7;

    do_tests(intersection[0], ref_indices, ref_num_indices);
    do_tests(intersection[1], ref_indices, ref_num_indices);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[10] = {{0,2,1},{3,2,1},{5,2,1},{8,2,1},{10,2,1},
                                      {14,2,1},{17,2,1},{20,2,1},{23,2,1},
                                      {25,2,1}};
    struct Xt_stripe stripes_b[5] = {{5,3,1},{8,2,1},{19,1,1},{20,2,1},
                                     {30,2,1}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 10);
    idxstripes_b = xt_idxstripes_new(stripes_b, 5);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_int ref_indices[6] = {5,6,8,9,20,21};
    int ref_num_indices = 6;

    do_tests(intersection[0], ref_indices, ref_num_indices);
    do_tests(intersection[1], ref_indices, ref_num_indices);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[3] = {{0,5,1},{1,5,1},{2,5,1}};
    struct Xt_stripe stripes_b[1] = {{-2,10,1}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 3);
    idxstripes_b = xt_idxstripes_new(stripes_b, 1);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_int ref_indices[2][15] = {{0,1,2,3,4,5,6},
                                 {0,1,1,2,2,2,3,3,3,4,4,4,5,5,6}};
    int ref_num_indices[2] = {7,15};

    do_tests(intersection[0], ref_indices[0], ref_num_indices[0]);
    do_tests(intersection[1], ref_indices[1], ref_num_indices[1]);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[1] = {{0,5,2}};
    struct Xt_stripe stripes_b[1] = {{1,5,2}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 1);
    idxstripes_b = xt_idxstripes_new(stripes_b, 1);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_int ref_indices[1] = {0};
    int ref_num_indices = 0;

    do_tests(intersection[0], ref_indices, ref_num_indices);
    do_tests(intersection[1], ref_indices, ref_num_indices);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  { // intersection test

    Xt_idxlist idxstripes_a, idxstripes_b;

    struct Xt_stripe stripes_a[1] = {{0,20,5}};
    struct Xt_stripe stripes_b[1] = {{1,15,7}};

    idxstripes_a = xt_idxstripes_new(stripes_a, 1);
    idxstripes_b = xt_idxstripes_new(stripes_b, 1);

    // compute intersection

    Xt_idxlist intersection[2];

    intersection[0] = xt_idxlist_get_intersection(idxstripes_a, idxstripes_b);
    intersection[1] = xt_idxlist_get_intersection(idxstripes_b, idxstripes_a);

    // check intersection

    Xt_int ref_indices[3] = {15,50,85};
    int ref_num_indices = 3;

    do_tests(intersection[0], ref_indices, ref_num_indices);
    do_tests(intersection[1], ref_indices, ref_num_indices);

    // clean up

    xt_idxlist_delete(idxstripes_a);
    xt_idxlist_delete(idxstripes_b);
    xt_idxlist_delete(intersection[0]);
    xt_idxlist_delete(intersection[1]);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: mixed valid and invalid positions
    const Xt_int undef_idx = XT_INT_MIN;
    struct Xt_stripe stripes[] = {{.start=0, .nstrides=5, .stride=1},
                                  {.start=10, .nstrides=5, .stride=1},
                                  {.start=20, .nstrides=5, .stride=-1}};
    int num_stripes = sizeof(stripes) / sizeof(stripes[0]);
    Xt_idxlist idxlist = xt_idxstripes_new(stripes, num_stripes);

    int pos[] = {0,2,7,9,11,100,11,200,9,300,18,400,5};
    int num_pos = sizeof(pos) / sizeof(pos[0]);
    Xt_int ref_sel_idx[num_pos];
    int ref_undef_count = 0;
    for (int ip=0; ip<num_pos; ip++) {
      int p = pos[ip];
      if (xt_idxlist_get_index_at_position(idxlist, p, &ref_sel_idx[ip]) != 0) {
        ref_sel_idx[ip] = undef_idx;
        ref_undef_count++;
      }
    }
    Xt_int sel_idx[num_pos];
    int undef_count = xt_idxlist_get_indices_at_positions(idxlist, pos, num_pos,
                                                          sel_idx, undef_idx);
    if (undef_count != ref_undef_count)
      PUT_ERR("test_idxvec.c: (undef_count != ref_undef_count)\n");
    for (int ip=0; ip<num_pos; ip++) {
      if (sel_idx[ip] != ref_sel_idx[ip])
        PUT_ERR("test_idxvec.c: (sel_idx[ip] != ref_sel_idx[ip])\n");
    }
    xt_idxlist_delete(idxlist);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: only valid positions
    const Xt_int undef_idx = XT_INT_MIN;
    struct Xt_stripe stripes[] = {{.start= 0, .nstrides=3, .stride= 1},
                                  {.start=10, .nstrides=2, .stride= 1},
                                  {.start=20, .nstrides=6, .stride=-1},
                                  {.start=30, .nstrides=7, .stride=-1}};
    int num_stripes = sizeof(stripes) / sizeof(stripes[0]);
    Xt_idxlist idxlist = xt_idxstripes_new(stripes, num_stripes);

    int pos[] = {-1,0,1,2,3,4,23,5,6,7,8,9,10,11,12, 0,2,100,2000};
    int num_pos = sizeof(pos) / sizeof(pos[0]);
    Xt_int ref_sel_idx[num_pos];
    int ref_undef_count = 0;
    for (int ip=0; ip<num_pos; ip++) {
      int p = pos[ip];
      if (xt_idxlist_get_index_at_position(idxlist, p, &ref_sel_idx[ip]) != 0) {
        ref_sel_idx[ip] = undef_idx;
        ref_undef_count++;
      }
    }
    Xt_int sel_idx[num_pos];
    int undef_count = xt_idxlist_get_indices_at_positions(idxlist, pos, num_pos,
                                                          sel_idx, undef_idx);
    if (undef_count != ref_undef_count)
      PUT_ERR("test_idxvec.c: (undef_count != ref_undef_count)\n");
    for (int ip=0; ip<num_pos; ip++) {
      if (sel_idx[ip] != ref_sel_idx[ip])
        PUT_ERR("test_idxvec.c: (sel_idx[ip] != ref_sel_idx[ip])\n");
    }
    xt_idxlist_delete(idxlist);

  }

  {
    // check idxvec_get_indices_at_positions
    // case: complete permutation
    const Xt_int undef_idx = XT_INT_MIN;
    struct Xt_stripe stripes[] = {{.start= 0, .nstrides=3, .stride= 1},
                                  {.start=10, .nstrides=2, .stride= 1},
                                  {.start=20, .nstrides=6, .stride=-1},
                                  {.start=30, .nstrides=7, .stride=-1}};
    int num_stripes = sizeof(stripes) / sizeof(stripes[0]);
    Xt_idxlist idxlist = xt_idxstripes_new(stripes, num_stripes);

    int pos[] = {4,7,2,5,9,0,10,6,11,8,12,1,3};
    int num_pos = sizeof(pos) / sizeof(pos[0]);
    Xt_int ref_sel_idx[num_pos];
    int ref_undef_count = 0;
    for (int ip=0; ip<num_pos; ip++) {
      int p = pos[ip];
      if (xt_idxlist_get_index_at_position(idxlist, p, &ref_sel_idx[ip]) != 0) {
        ref_sel_idx[ip] = undef_idx;
        ref_undef_count++;
      }
    }
    Xt_int sel_idx[num_pos];
    int undef_count = xt_idxlist_get_indices_at_positions(idxlist, pos,
                                                          num_pos, sel_idx,
                                                          undef_idx);
    if (undef_count != ref_undef_count)
      PUT_ERR("test_idxvec.c: (undef_count != ref_undef_count)\n");
    for (int ip=0; ip<num_pos; ip++) {
      if (sel_idx[ip] != ref_sel_idx[ip])
        PUT_ERR("test_idxvec.c: (sel_idx[ip] != ref_sel_idx[ip])\n");
    }
    xt_idxlist_delete(idxlist);

  }

  {
    // check idxvec_get_indices_at_positions
    // case: only invalid positions
    const Xt_int undef_idx = XT_INT_MIN;
    struct Xt_stripe stripes[] = {{.start=0, .nstrides=5, .stride=1},
                                  {.start=10, .nstrides=5, .stride=1},
                                  {.start=20, .nstrides=5, .stride=-1}};
    int num_stripes = sizeof(stripes) / sizeof(stripes[0]);
    Xt_idxlist idxlist = xt_idxstripes_new(stripes, num_stripes);

    int pos[] = {-10,200,700,90,90,18,141};
    int num_pos = sizeof(pos) / sizeof(pos[0]);
    Xt_int ref_sel_idx[num_pos];
    int ref_undef_count = 0;
    for (int ip=0; ip<num_pos; ip++) {
      int p = pos[ip];
      if (xt_idxlist_get_index_at_position(idxlist, p, &ref_sel_idx[ip]) != 0) {
        ref_sel_idx[ip] = undef_idx;
        ref_undef_count++;
      }
    }

    Xt_int sel_idx[num_pos];
    int undef_count
      = xt_idxlist_get_indices_at_positions(idxlist, pos, num_pos,
                                            sel_idx, undef_idx);
    if (undef_count != ref_undef_count)
      PUT_ERR("test_idxvec.c: (undef_count != ref_undef_count)\n");
    for (int ip=0; ip<num_pos; ip++) {
      if (sel_idx[ip] != ref_sel_idx[ip])
        PUT_ERR("test_idxvec.c: (sel_idx[ip] != ref_sel_idx[ip])\n");
    }
    xt_idxlist_delete(idxlist);
  }

  { // test with overlapping stripes
    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[2] = {{0,5,1},{1,5,1}};
    int num_stripes = 2;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    Xt_int ref_indices[10] = {0,1,2,3,4, 1,2,3,4,5};
    int ref_num_indices = 10;

    do_tests(idxstripes, ref_indices, ref_num_indices);

    xt_idxlist_delete(idxstripes);
  }

  { // check get_bounding_box

    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[1] = {{-1,-1,-1}};
    int num_stripes = 0;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    unsigned ndim = 3;
    Xt_int global_size[ndim];
    Xt_int global_start_index = 0;
    struct Xt_bounds bounds[ndim];

    for (unsigned i = 0; i < ndim; ++i)
      global_size[i] = 4;

    xt_idxlist_get_bounding_box(idxstripes, ndim, global_size,
                                global_start_index, bounds);

    for (unsigned i = 0; i < ndim; ++i)
      if (bounds[i].size != 0)
        PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

    xt_idxlist_delete(idxstripes);
  }

  { // check get_bounding_box

    Xt_idxlist idxstripes;

    struct Xt_stripe stripes[3] = {{47,2,-12}, {32,2,12}, {36,2,12}};
    int num_stripes = 3;

    idxstripes = xt_idxstripes_new(stripes, num_stripes);

    unsigned ndim = 3;
    Xt_int global_size[ndim];
    Xt_int global_start_index = 1;
    struct Xt_bounds bounds[ndim];

    global_size[0] = 5;
    global_size[1] = 4;
    global_size[2] = 3;

    xt_idxlist_get_bounding_box(idxstripes, ndim, global_size,
                                global_start_index, bounds);

    Xt_int ref_start[3] = {2,2,1};

    for (unsigned i = 0; i < ndim; ++i)
      if (bounds[i].size != 2 || bounds[i].start != ref_start[i])
        PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

    xt_idxlist_delete(idxstripes);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
do_tests(Xt_idxlist idxlist, Xt_int * ref_indices, int num_indices) {

  check_idxlist(idxlist, ref_indices, num_indices);

  struct Xt_stripe * stripes;
  int num_stripes;
  Xt_idxlist temp_idxlist;

  xt_idxlist_get_index_stripes(idxlist, &stripes, &num_stripes);
  temp_idxlist = xt_idxvec_from_stripes_new(stripes, num_stripes);

  check_idxlist(temp_idxlist, ref_indices, num_indices);

  xt_idxlist_delete(temp_idxlist);

  free(stripes);

  { // test packing and unpacking

    size_t buffer_size;

    buffer_size = xt_idxlist_get_pack_size(idxlist, MPI_COMM_WORLD);

    // allocate send buffers

    void* send_buffer, * recv_buffer;

    send_buffer = xmalloc(buffer_size);
    recv_buffer = xmalloc(buffer_size);

    // pack the index vector

    int position = 0;

    assert(buffer_size <= INT_MAX);
    xt_idxlist_pack(idxlist, send_buffer, (int)buffer_size, &position,
                    MPI_COMM_WORLD);

    // send the buffer

    int rank;

    xt_mpi_call(MPI_Comm_rank(MPI_COMM_WORLD, &rank), MPI_COMM_WORLD);

    MPI_Request request;

    xt_mpi_call(MPI_Isend(send_buffer, (int)buffer_size, MPI_PACKED, rank, 0,
                              MPI_COMM_WORLD, &request), MPI_COMM_WORLD);
    xt_mpi_call(MPI_Request_free(&request), MPI_COMM_WORLD);

    // receive the buffer

    xt_mpi_call(MPI_Recv(recv_buffer, (int)buffer_size, MPI_PACKED, rank, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE),
                    MPI_COMM_WORLD);

    // unpack the buffer

    Xt_idxlist idxlist_copy;

    position = 0;
    idxlist_copy = xt_idxlist_unpack(recv_buffer, (int)buffer_size, &position,
                                     MPI_COMM_WORLD);

    // check copy

    check_idxlist(idxlist_copy, ref_indices, num_indices);

    // clean up

    free(send_buffer);
    free(recv_buffer);
    xt_idxlist_delete(idxlist_copy);
  }

  { // test copying

    Xt_idxlist idxlist_copy;

    idxlist_copy = xt_idxlist_copy(idxlist);

    // check copy

    check_idxlist(idxlist_copy, ref_indices, num_indices);

    // clean up

    xt_idxlist_delete(idxlist_copy);
  }
}
