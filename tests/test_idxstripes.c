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
#include <limits.h>
#include <stdlib.h>

#include <mpi.h>

#include <yaxt.h>

#include "core/ppm_xfuncs.h"
#include "tests.h"
#include "test_idxlist_utils.h"

static void
do_tests(Xt_idxlist idxlist, const Xt_int *ref_indices, int num_indices);

static void
stripe_test_general(int num_stripes, const struct Xt_stripe *stripes,
                    int num_ref_indices, const Xt_int *ref_indices);

static void
stripe_test_intersection(int num_stripes_a, const struct Xt_stripe *stripes_a,
                         int num_stripes_b, const struct Xt_stripe *stripes_b,
                         int num_ref_indices, const Xt_int *ref_indices);

static void
stripe_test_asymmetric_intersection(
  int num_stripes_a, const struct Xt_stripe *stripes_a,
  int num_stripes_b, const struct Xt_stripe *stripes_b,
  int num_ref_indices_a, const Xt_int *ref_indices_a,
  int num_ref_indices_b, const Xt_int *ref_indices_b);

static void
check_idxvec_get_indices_at_positions(int num_stripes,
                                      const struct Xt_stripe *stripes,
                                      int num_pos, const int *pos);

static void
check_bb(int num_stripes, const struct Xt_stripe *stripes,
         unsigned ndim, const Xt_int global_size[ndim],
         const struct Xt_bounds ref_bounds[ndim],
         Xt_int global_start_index);

static void
check_pos_ext(size_t num_stripes, const struct Xt_stripe stripes[num_stripes],
              size_t num_search_stripes,
              const struct Xt_stripe search_stripes[num_search_stripes],
              size_t ref_num_ext,
              const struct Xt_pos_ext ref_pos_ext[ref_num_ext],
              int single_match_only, int ref_unmatched,
              const char *test_desc);

int main(void) {

  // init mpi

  xt_mpi_call(MPI_Init(NULL, NULL), MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // general tests
    enum {
      num_stripes = 3,
      num_ref_indices = 15,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start =  0, .stride = 1, .nstrides = 5},
         {.start = 10, .stride = 1, .nstrides = 5},
         {.start = 20, .stride = 1, .nstrides = 5}};

    static const Xt_int ref_indices[num_ref_indices]
      = {0,1,2,3,4, 10,11,12,13,14, 20,21,22,23,24};

    stripe_test_general(num_stripes, stripes, num_ref_indices, ref_indices);
  }

  { // general tests
    enum {
      num_stripes = 3,
      num_ref_indices = 15,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start =  0, .stride = 1, .nstrides = 5},
         {.start = 10, .stride = 2, .nstrides = 5},
         {.start = 20, .stride = 3, .nstrides = 5}};

    static const Xt_int ref_indices[num_ref_indices]
      = {0,1,2,3,4, 10,12,14,16,18, 20,23,26,29,32};

    stripe_test_general(num_stripes, stripes, num_ref_indices, ref_indices);
  }

  { // general tests
    enum {
      num_stripes = 2,
      num_ref_indices = 10,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 0, .stride = 6, .nstrides = 5},
         {.start = 1, .stride = 3, .nstrides = 5}};

    static const Xt_int ref_indices[num_ref_indices]
      = {0,6,12,18,24, 1,4,7,10,13};

    stripe_test_general(num_stripes, stripes, num_ref_indices, ref_indices);
  }

  { // general tests
    enum {
      num_stripes = 2,
      num_ref_indices = 10,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 0, .stride = -1, .nstrides = 5},
         {.start = 1, .stride =  1, .nstrides = 5}};

    static const Xt_int ref_indices[num_ref_indices]
      = {0,-1,-2,-3,-4, 1,2,3,4,5};

    stripe_test_general(num_stripes, stripes, num_ref_indices, ref_indices);
  }

  { // general tests
    enum {
      num_stripes = 2,
      num_ref_indices = 10,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 9, .stride = -2, .nstrides = 5},
         {.start = 0, .stride =  2, .nstrides = 5}};

    static const Xt_int ref_indices[num_ref_indices]
      = {9,7,5,3,1, 0,2,4,6,8};

    stripe_test_general(num_stripes, stripes, num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 2,
      num_stripes_b = 1,
      num_ref_indices = 6,
    };
    static const struct Xt_stripe
      stripes_a[num_stripes_a] = {{.start = 0, .stride = 1, .nstrides = 4},
                                  {.start = 6, .stride = 1, .nstrides = 4}},
      stripes_b[num_stripes_b] = {{.start = 1, .stride = 1, .nstrides = 8}};
    static const Xt_int ref_indices[num_ref_indices] = {1,2,3, 6,7,8};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 3,
      num_stripes_b = 2,
      num_ref_indices = 9,
    };
    static const struct Xt_stripe
      stripes_a[num_stripes_a] = {{.start =  0, .stride = 1, .nstrides = 4},
                                  {.start =  6, .stride = 1, .nstrides = 4},
                                  {.start = 11, .stride = 1, .nstrides = 4}},
      stripes_b[num_stripes_b] = {{.start =  1, .stride = 1, .nstrides = 7},
                                  {.start =  9, .stride = 1, .nstrides = 5}};
    static const Xt_int ref_indices[num_ref_indices]
      = {1,2,3, 6,7, 9, 11,12,13};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 2,
      num_stripes_b = 2,
      num_ref_indices = 0,
    };
    static const struct Xt_stripe
      stripes_a[num_stripes_a] = {{.start =  0, .stride = 1, .nstrides = 3},
                                  {.start =  8, .stride = 1, .nstrides = 3}},
      stripes_b[num_stripes_b] = {{.start =  3, .stride = 1, .nstrides = 5},
                                  {.start = 11, .stride = 1, .nstrides = 3}};
    static const Xt_int *ref_indices = NULL;

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 1,
      num_stripes_b = 2,
      num_ref_indices = 10,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 0, .nstrides = 10, .stride =  1}},
      stripes_b[num_stripes_b] = {{.start = 0, .nstrides =  5, .stride =  2},
                                  {.start = 9, .nstrides =  5, .stride = -2}};
    static const Xt_int ref_indices[num_ref_indices] = {0,1,2,3,4,5,6,7,8,9};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 2,
      num_stripes_b = 2,
      num_ref_indices = 6,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start =  0, .stride =  3, .nstrides =  5},
         {.start =  1, .stride =  7, .nstrides =  5}},
      stripes_b[num_stripes_b] = {{.start =  0, .stride =  2, .nstrides =  7},
                                  {.start = 24, .stride = -1, .nstrides = 10}};
    static const Xt_int ref_indices[6] = {0,6,8,12,15,22};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 1,
      num_stripes_b = 2,
      num_ref_indices = 10,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 0, .stride =  1, .nstrides = 10}},
      stripes_b[num_stripes_b] = {{.start = 5, .stride =  1, .nstrides =  5},
                                  {.start = 4, .stride = -1, .nstrides =  5}};
    static const Xt_int ref_indices[10] = {0,1,2,3,4,5,6,7,8,9};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 2,
      num_stripes_b = 2,
      num_ref_indices = 7,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start =  0, .stride = 1, .nstrides = 10},
         {.start = 20, .stride = 1, .nstrides =  5}},
      stripes_b[num_stripes_b] = {{.start =  3, .stride = 1, .nstrides =  5},
                                  {.start = 17, .stride = 1, .nstrides =  5}};
    static const Xt_int ref_indices[7] = {3,4,5,6,7,20,21};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 10,
      num_stripes_b = 5,
      num_ref_indices = 6,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start =  0, .stride = 1, .nstrides = 2},
         {.start =  3, .stride = 1, .nstrides = 2},
         {.start =  5, .stride = 1, .nstrides = 2},
         {.start =  8, .stride = 1, .nstrides = 2},
         {.start = 10, .stride = 1, .nstrides = 2},
         {.start = 14, .stride = 1, .nstrides = 2},
         {.start = 17, .stride = 1, .nstrides = 2},
         {.start = 20, .stride = 1, .nstrides = 2},
         {.start = 23, .stride = 1, .nstrides = 2},
         {.start = 25, .stride = 1, .nstrides = 2}},
      stripes_b[num_stripes_b] =  {{.start =  5, .stride = 1, .nstrides = 3},
                                   {.start =  8, .stride = 1, .nstrides = 2},
                                   {.start = 19, .stride = 1, .nstrides = 1},
                                   {.start = 20, .stride = 1, .nstrides = 2},
                                   {.start = 30, .stride = 1, .nstrides = 2}};
    static const Xt_int ref_indices[6] = {5,6,8,9,20,21};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 3,
      num_stripes_b = 1,
      num_ref_indices_a = 7,
      num_ref_indices_b = 15,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start =  0, .stride = 1, .nstrides =  5},
         {.start =  1, .stride = 1, .nstrides =  5},
         {.start =  2, .stride = 1, .nstrides =  5}},
      stripes_b[num_stripes_b] = {{.start = -2, .stride = 1, .nstrides = 10}};
    static const Xt_int ref_indices_a[num_ref_indices_a] = {0,1,2,3,4,5,6},
      ref_indices_b[num_ref_indices_b] = {0,1,1,2,2,2,3,3,3,4,4,4,5,5,6};

    stripe_test_asymmetric_intersection(num_stripes_a, stripes_a,
                                        num_stripes_b, stripes_b,
                                        num_ref_indices_a, ref_indices_a,
                                        num_ref_indices_b, ref_indices_b);
  }

  { // intersection test
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 0,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 0, .stride = 2, .nstrides = 5}},
      stripes_b[num_stripes_b] = {{.start = 1, .stride = 2, .nstrides = 5}};

    static const Xt_int ref_indices[1] = {0};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 3,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 0, .stride = 5, .nstrides = 20}},
      stripes_b[num_stripes_b] = {{.start = 1, .stride = 7, .nstrides = 15}};
    static const Xt_int ref_indices[3] = {15,50,85};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  { // intersection test, both ranges overlap in range but have no
    // indices in common because of stride
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 0,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 34, .stride = 29, .nstrides = 12}},
      stripes_b[num_stripes_b] = {{.start = 36, .stride = 7, .nstrides = 2}};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, NULL);
  }

  { // intersection test, same as before but with negative stride
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 0,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 353, .stride = -29, .nstrides = 12}},
      stripes_b[num_stripes_b] = {{.start = 36, .stride = 7, .nstrides = 2}};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, NULL);
  }

  {
    /* test case where overlap equals start of one and end of other stripe */
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 1,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 95, .stride = -29, .nstrides = 2}},
      stripes_b[num_stripes_b] = {{.start = 81, .stride = 14, .nstrides = 2}};
    static const Xt_int ref_indices[num_ref_indices] = {95};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  {
    /* test case where overlap equals end of both stripes */
    enum {
      num_stripes_a = 1,
      num_stripes_b = 1,
      num_ref_indices = 1,
    };
    static const struct Xt_stripe stripes_a[num_stripes_a]
      = {{.start = 546, .stride = 14, .nstrides = 2}},
      stripes_b[num_stripes_b] = {{.start = 354, .stride = 206, .nstrides = 2}};
    static const Xt_int ref_indices[num_ref_indices] = {560};

    stripe_test_intersection(num_stripes_a, stripes_a,
                             num_stripes_b, stripes_b,
                             num_ref_indices, ref_indices);
  }

  {
    // generate idxvec from stripes
    static const struct Xt_stripe stripes[]
      = {{.start = 4, .stride = 1, .nstrides = 1},
         {.start = 5, .stride = 1, .nstrides = 1},
         {.start = 10, .stride = -10, .nstrides = 2}};
    int num_stripes = sizeof(stripes) / sizeof(stripes[0]);
    Xt_idxlist idxvec_a = xt_idxvec_from_stripes_new(stripes, num_stripes);

    // generate idxvec containing single elemente that is also with idxvec_a
    Xt_int index_vector[] = {5};
    int num_indices = sizeof(index_vector) / sizeof(index_vector[0]);
    Xt_idxlist idxvec_b = xt_idxvec_new(index_vector, num_indices);

    // intersect both index vectors
    Xt_idxlist intersection =
      xt_idxlist_get_intersection(idxvec_a, idxvec_b);

    // test result
    if (xt_idxlist_get_num_indices(intersection) != 1) {
      PUT_ERR("wrong number of indices in intersection\n");
    } else {
      Xt_int idx;
      xt_idxlist_get_index_at_position(intersection, 0, &idx);
      if (idx != index_vector[0])
        PUT_ERR("wrong index in intersection\n");
    }
    xt_idxlist_delete(intersection);
    xt_idxlist_delete(idxvec_b);
    xt_idxlist_delete(idxvec_a);
  }

  {
    Xt_int index_vector[] = {
      3375, 3376,
      3379, 3380, 3381,
      3387, 3388, 3389, 3390, 3391, 3392, 3393,
      3421, 3422, 3423, 3424, 3425, 3426, 3427,
      3444,
      3458, 3459,
      3461, 3462, 3463, 3464, 3465, 3466, 3467, 3468, 3469, 3470,
        3471, 3472, 3473, 3474, 3475, 3476, 3477, 3478, 3479, 3480,
      3529,
      3606, 3607, 3608,
      3611, 3612, 3613, 3614,
      3617,
      3620, 3621, 3622, 3623, 3624, 3625, 3626, 3627, 3628, 3629,
        3630, 3631,
      3684, 3685, 3686, 3687, 3688, 3689, 3690, 3691, 3692, 3693,
        3694, 3695, 3696, 3697, 3698, 3699, 3700, 3701, 3702, 3703,
        3704, 3705, 3706, 3707, 3708, 3709,
      3713, 3714, 3715, 3716, 3717, 3718, 3719, 3720, 3721, 3722,
        3723, 3724, 3725, 3726, 3727, 3728, 3729, 3730, 3731,
      3741, 3742,
      3931, 3932,
      3374,
      3382,
      3385,
      3394,
      3404,
      3408,
      3412,
      3440,
      3443,
      3457,
      3481,
      3483,
      3527,
      3619,
      3735,
      3743,
      3925,
      3930,
      3377, 3378,
      3383, 3384,
      3386,
      3395,
      3397, 3398,
      3400,
      3402, 3403,
      3407,
      3409, 3410,
      3413,
      3420,
      3441, 3442,
      3445,
      3448, 3449,
      3451,
      3460,
      3482,
      3519, 3520,
      3526,
      3528,
      3530,
      3592, 3593,
      3595, 3596, 3597,
      3609, 3610,
      3615, 3616,
      3618,
      3644,
      3710, 3711, 3712,
      3732, 3733,
      3736, 3737,
      3748, 3749,
      3753, 3754,
      3759, 3760,
      3766, 3767,
      3919, 3920,
      3924,
      3926,
      3933, 3934,
      2589,
      2602,
      2680,
      3326,
      3340, 3341,
      3396,
      3401,
      3411,
      3414,
      3418,
      3446, 3447,
      3450,
      3515,
      3521,
      3525,
      3582,
      3590, 3591,
      3594,
      3642,
      3734,
      3738,
      3747,
      3750,
      3761,
      3765,
      3865,
      3918,
      3923,
      3935
    };
    int num_indices = sizeof(index_vector) / sizeof(index_vector[0]);
    Xt_idxlist idxlist = xt_idxvec_new(index_vector, num_indices);

    struct  Xt_stripe stripes[] = {
      {.start = 3326, .stride = 14, .nstrides = 2 },
      {.start = 3341, .stride = 33, .nstrides = 1 },
      {.start = 3374, .stride = 1 , .nstrides = 25},
      {.start = 3400, .stride = 1 , .nstrides = 5 },
      {.start = 3407, .stride = 1 , .nstrides = 8 },
      {.start = 3418, .stride = 2 , .nstrides = 1 },
      {.start = 3420, .stride = 1 , .nstrides = 8 },
      {.start = 3440, .stride = 1 , .nstrides = 12},
      {.start = 3457, .stride = 1 , .nstrides = 27},
      {.start = 3515, .stride = 4 , .nstrides = 1 },
      {.start = 3519, .stride = 1 , .nstrides = 3 },
      {.start = 3525, .stride = 1 , .nstrides = 6 },
      {.start = 3582, .stride = 8 , .nstrides = 1 },
      {.start = 3590, .stride = 1 , .nstrides = 8 },
      {.start = 3606, .stride = 1 , .nstrides = 26},
      {.start = 3642, .stride = 2 , .nstrides = 2 },
      {.start = 3684, .stride = 1 , .nstrides = 55},
      {.start = 3741, .stride = 1 , .nstrides = 3 },
      {.start = 3747, .stride = 1 , .nstrides = 4 },
      {.start = 3753, .stride = 1 , .nstrides = 2 },
      {.start = 3759, .stride = 1 , .nstrides = 3 },
      {.start = 3765, .stride = 1 , .nstrides = 3 },
      {.start = 3865, .stride = 53, .nstrides = 1 },
      {.start = 3918, .stride = 1 , .nstrides = 3 },
      {.start = 3923, .stride = 1 , .nstrides = 4 },
      {.start = 3930, .stride = 1 , .nstrides = 6 }
    };
    int num_stripes = sizeof(stripes) / sizeof(stripes[0]);

    int num_ext;
    struct Xt_pos_ext * pos_ext;
    int single_match_only = 1;

    xt_idxlist_get_pos_exts_of_index_stripes(idxlist, num_stripes, stripes,
                                             &num_ext, &pos_ext,
                                             single_match_only);
    { // testing of results
      Xt_idxlist intersection = xt_idxvec_from_stripes_new(stripes, num_stripes);
      for (int i = 0, k = 0; i < num_ext; ++i) {
        int abs_pos_ext_size = abs(pos_ext[i].size),
          jsign = pos_ext[i].size >= 0 ? 1 : -1;
        for (int j = 0; j < abs_pos_ext_size; ++j, ++k) {
          Xt_int intersection_index;
          xt_idxlist_get_index_at_position(intersection, k, &intersection_index);
          int send_pos = pos_ext[i].start + jsign * j;
          Xt_int orig_index;
          xt_idxlist_get_index_at_position(idxlist, send_pos, &orig_index);
          if (intersection_index != orig_index) {
            PUT_ERR("error in xt_idxlist_get_pos_exts_of_index_stripes\n");
            fprintf(stderr, "intersection pos %d index %"XT_INT_FMT" orig pos %d index %"XT_INT_FMT"\n",
                    k, intersection_index, send_pos, orig_index);
          }
        }
      }
      xt_idxlist_delete(intersection);
    }

    free(pos_ext);
    xt_idxlist_delete(idxlist);
  }

  {
    enum {
      NUM_ITERATIONS=128,
      MAX_NUM_INDICES=1024,
      MAX_INDEX=1024,
    };

    /* demonstrates problem of invalid intersection
     * get_stripe_intersection({916,-287,2}, {603,300,2})
     * gives an intersection with nstrides == 1 even though we'd need
     * nstrides == 0
     */

    srand(6544668);

    for (int iteration = 0; iteration < NUM_ITERATIONS; ++iteration) {

      int num_indices = rand()%MAX_NUM_INDICES;

      Xt_int * indices = malloc((size_t)num_indices * sizeof(*indices));

      for (int i = 0; i < num_indices; ++i)
        indices[i] = (Xt_int)(rand()%(2*MAX_INDEX)-MAX_INDEX);

      Xt_idxlist idxlist = xt_idxvec_new(indices, num_indices);

      struct Xt_stripe * stripes;
      int num_stripes;

      xt_idxlist_get_index_stripes(idxlist, &stripes, &num_stripes);

      int num_ext;
      struct Xt_pos_ext * pos_ext;
      int single_match_only = 1;

      xt_idxlist_get_pos_exts_of_index_stripes(idxlist, num_stripes, stripes,
                                               &num_ext, &pos_ext,
                                               single_match_only);
      { // testing of results
        Xt_idxlist intersection =
          xt_idxvec_from_stripes_new(stripes, num_stripes);
        for (int i = 0, k = 0; i < num_ext; ++i) {
          int abs_pos_ext_size = abs(pos_ext[i].size),
            jsign = pos_ext[i].size >= 0 ? 1 : -1;
          for (int j = 0; j < abs_pos_ext_size; ++j, ++k) {
            Xt_int intersection_index;
            xt_idxlist_get_index_at_position(intersection, k, &intersection_index);
            int send_pos = pos_ext[i].start + jsign * j;
            Xt_int orig_index;
            xt_idxlist_get_index_at_position(idxlist, send_pos, &orig_index);
            if (intersection_index != orig_index) {
              PUT_ERR("error in xt_idxlist_get_pos_exts_of_index_stripes\n");
              fprintf(stderr, "intersection pos %d index %"XT_INT_FMT
                      " orig pos %d index %"XT_INT_FMT"\n", k,
                      intersection_index, send_pos, orig_index);
            }
          }
        }
        xt_idxlist_delete(intersection);
      }

      xt_idxlist_delete(idxlist);
      free(indices);
      free(stripes);
      free(pos_ext);
    }
  }

  {
    static const Xt_int index_vector[] = {
      178, 179, 180, 181, 182, 183, 184,
      186, 187, 188, 189, 190,
      194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
        208, 209, 210, 211, 212,
      217,
      223,
      426,
      428, 429, 430,
      434, 435, 436, 437, 438, 439, 440,
      442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455,
        456, 457, 458,
      670, 671, 672, 673, 674, 675, 676, 677,
      682,
      684, 685, 686, 687, 688, 689, 690,
      692,
      695,
      703, 704, 705, 706, 707,
      894, 895, 896, 897, 898, 899, 900, 901,
      906, 907, 908,
      913,
      915,
      921, 922, 923, 924, 925, 926, 927,
      1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103,
      1107, 1108, 1109, 1110, 1111,
      1113, 1114,
      1119, 1120, 1121,
      2095, 2096, 2097, 2098,
      2100,
      2102, 2103, 2104, 2105,
      2107, 2108, 2109, 2110,
      2112,
      2118,
      2120, 2121, 2122, 2123, 2124, 2125,
      2127, 2128, 2129, 2130,
      2134,
      2140, 2141, 2142, 2143,
      2145,
      2148, 2149,
      2151, 2152, 2153, 2154, 2155, 2156,
      683,
      691,
      903,
      914,
      1105,
      1115,
      2099,
      2106,
      2111,
      2115,
      2126,
      2132,
      2139,
      2144,
      2147,
      2150,
      2305,
      427,
      465, 466,
      678,
      693,
      902,
      909,
      1104,
      1112,
      2101,
      2113, 2114,
      2116, 2117,
      2119,
      2131,
      2136,
      2138,
      2146,
      2297,
      2302,
      2304,
      2307
    };
    enum { num_indices = sizeof(index_vector) / sizeof(index_vector[0]) };
    Xt_idxlist idxlist = xt_idxvec_new(index_vector, num_indices);

    static const struct  Xt_stripe stripes[] = {
      {.start = 670, .stride = 1, .nstrides = 9},
      {.start = 682, .stride = 1, .nstrides = 12},
      {.start = 695, .stride = 8, .nstrides = 1},
      {.start = 703, .stride = 1, .nstrides = 5},
      {.start = 894, .stride = 1, .nstrides = 10},
      {.start = 906, .stride = 1, .nstrides = 4},
      {.start = 913, .stride = 1, .nstrides = 3},
      {.start = 921, .stride = 1, .nstrides = 7}
    };
    enum { num_stripes = sizeof(stripes) / sizeof(stripes[0]) };

    int num_ext;
    struct Xt_pos_ext *restrict pos_ext;
    int single_match_only = 1;

    int retval = xt_idxlist_get_pos_exts_of_index_stripes(
      idxlist, num_stripes, stripes, &num_ext, (struct Xt_pos_ext **)&pos_ext,
      single_match_only);

    { // testing of results
      if (retval != 0)
          PUT_ERR("error in xt_idxlist_get_pos_exts_of_index_stripes\n");
      Xt_idxlist intersection
        = xt_idxvec_from_stripes_new(stripes, num_stripes);
      for (int i = 0, k = 0; i < num_ext; ++i) {
        int abs_pos_ext_size = abs(pos_ext[i].size),
          jsign = pos_ext[i].size >= 0 ? 1 : -1;
        for (int j = 0; j < abs_pos_ext_size; ++j, ++k) {
          Xt_int intersection_index;
          xt_idxlist_get_index_at_position(intersection, k,
                                           &intersection_index);
          int send_pos = pos_ext[i].start + jsign * j;
          Xt_int orig_index;
          xt_idxlist_get_index_at_position(idxlist, send_pos, &orig_index);
          if (intersection_index != orig_index) {
            PUT_ERR("error in xt_idxlist_get_pos_exts_of_index_stripes\n");
            fprintf(stderr, "intersection pos %d index %"XT_INT_FMT
                    " orig pos %d index %"XT_INT_FMT"\n",
                    k, intersection_index, send_pos, orig_index);
          }
        }
      }
      xt_idxlist_delete(intersection);
    }

    free(pos_ext);
    xt_idxlist_delete(idxlist);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: mixed valid and invalid positions
    static const struct Xt_stripe stripes[]
      = {{.start=0, .nstrides=5, .stride=1},
         {.start=10, .nstrides=5, .stride=1},
         {.start=20, .nstrides=5, .stride=-1}};
    int pos[] = {0,2,7,9,11,100,11,200,9,300,18,400,5};
    enum {
      num_pos = sizeof(pos) / sizeof(pos[0]),
      num_stripes = sizeof(stripes) / sizeof(stripes[0]),
    };
    check_idxvec_get_indices_at_positions(num_stripes, stripes, num_pos, pos);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: only valid positions
    static const struct Xt_stripe stripes[]
      = {{.start= 0, .nstrides=3, .stride= 1},
         {.start=10, .nstrides=2, .stride= 1},
         {.start=20, .nstrides=6, .stride=-1},
         {.start=30, .nstrides=7, .stride=-1}};
    int pos[] = {-1,0,1,2,3,4,23,5,6,7,8,9,10,11,12, 0,2,100,2000};
    enum {
      num_stripes = sizeof(stripes) / sizeof(stripes[0]),
      num_pos = sizeof(pos) / sizeof(pos[0]),
    };
    check_idxvec_get_indices_at_positions(num_stripes, stripes, num_pos, pos);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: complete permutation
    static const struct Xt_stripe stripes[]
      = {{.start= 0, .nstrides=3, .stride= 1},
         {.start=10, .nstrides=2, .stride= 1},
         {.start=20, .nstrides=6, .stride=-1},
         {.start=30, .nstrides=7, .stride=-1}};
    int pos[] = {4,7,2,5,9,0,10,6,11,8,12,1,3};
    enum {
      num_stripes = sizeof(stripes) / sizeof(stripes[0]),
      num_pos = sizeof(pos) / sizeof(pos[0]),
    };
    check_idxvec_get_indices_at_positions(num_stripes, stripes, num_pos, pos);
  }

  {
    // check idxvec_get_indices_at_positions
    // case: only invalid positions
    static const struct Xt_stripe stripes[]
      = {{.start=0, .nstrides=5, .stride=1},
         {.start=10, .nstrides=5, .stride=1},
         {.start=20, .nstrides=5, .stride=-1}};
    int pos[] = {-10,200,700,90,90,18,141};
    enum {
      num_stripes = sizeof(stripes) / sizeof(stripes[0]),
      num_pos = sizeof(pos) / sizeof(pos[0]),
    };
    check_idxvec_get_indices_at_positions(num_stripes, stripes, num_pos, pos);
  }

  { // test with overlapping stripes
    enum {
      num_stripes = 2,
      num_ref_indices = 10,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 0, .stride = 1, .nstrides = 5},
         {.start = 1, .stride = 1, .nstrides = 5}};
    static const Xt_int ref_indices[num_ref_indices] = {0,1,2,3,4, 1,2,3,4,5};

    stripe_test_general(num_stripes, stripes, num_ref_indices, ref_indices);
  }

  { // check get_bounding_box
    enum {
      num_stripes = 0,
      num_ref_indices = 10,
      ndim = 3,
    };
    static const struct Xt_stripe stripes[1]
      = {{.start = -1, .stride = -1, .nstrides = -1}};
    static const Xt_int global_size[ndim] = { 4, 4, 4 };
    static const struct Xt_bounds ref_bounds[ndim]
      = { { .start = -1, .size = 0 },
          { .start = -1, .size = 0 },
          { .start = -1, .size = 0 } };
    Xt_int global_start_index = 0;

    check_bb(num_stripes, stripes,
             ndim, global_size, ref_bounds, global_start_index);
  }

  { // check get_bounding_box
    enum {
      num_stripes = 3,
      ndim = 3,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 47, .stride = -12, .nstrides = 2},
         {.start = 32, .stride =  12, .nstrides = 2},
         {.start = 36, .stride =  12, .nstrides = 2}};
    static const Xt_int global_size[ndim] = { 5, 4, 3 };
    static const struct Xt_bounds ref_bounds[ndim]
      = { { .start = 2, .size = 2 },
          { .start = 2, .size = 2 },
          { .start = 1, .size = 2 } };
    Xt_int global_start_index = 1;

    check_bb(num_stripes, stripes,
             ndim, global_size, ref_bounds, global_start_index);
  }

  {
    enum {
      num_stripes = 1,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 0,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 10}},
      search_stripe = {.start = 10, .stride = -1, .nstrides = 5 };

    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
      = {{.start = 9, .size = -5}};

    check_pos_ext(num_stripes, stripes,
                  1, &search_stripe, num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "simple inverted stripe");
  }

  {
    enum {
      num_stripes = 1,
      num_search_stripes = 2,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 5,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 10}},
      search_stripes[num_search_stripes]
        = {{.start = 10, .stride = -1, .nstrides = 5 },
           {.start = 10, .stride = -1, .nstrides = 5 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
                     = {{.start = 9, .size = -5}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "simple inverted stripe");
  }

  {
    enum {
      num_stripes = 2,
      num_search_stripes = 1,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 4,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 10},
         {.start = 15, .stride = 1, .nstrides = 10}},
      search_stripes[num_search_stripes]
        = {{.start = 10, .stride = 1, .nstrides = 6 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
                     = {{.start = 9, .size = 2}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "search inc stripe over inc gap");
  }

  {
    enum {
      num_stripes = 2,
      num_search_stripes = 1,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 4,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 25, .stride = -1, .nstrides = 11},
         {.start = 10, .stride = -1, .nstrides = 10}},
      search_stripes[num_search_stripes]
        = {{.start = 10, .stride = 1, .nstrides = 6 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
                     = {{.start = 11, .size = -2}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "search inc stripe over dec gap");
  }

  {
    enum {
      num_stripes = 2,
      num_search_stripes = 1,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 4,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 25, .stride = -1, .nstrides = 11},
         {.start = 10, .stride = -1, .nstrides = 10}},
      search_stripes[num_search_stripes]
        = {{.start = 15, .stride = -1, .nstrides = 6 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
                     = {{.start = 10, .size = 2}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "search dec stripe over dec gap");
  }

  {
    enum {
      num_stripes = 2,
      num_search_stripes = 1,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 4,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 10},
         {.start = 15, .stride = 1, .nstrides = 10}},
      search_stripes[num_search_stripes]
        = {{.start = 15, .stride = -1, .nstrides = 6 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
                     = {{.start = 10, .size = -2}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "search dec stripe over inc gap");
  }

  {
    enum {
      num_stripes = 3,
      num_search_stripes = 1,
      num_ref_pos_ext = 1,
      num_ref_unmatched = 8,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 10},
         {.start = 15, .stride = 1, .nstrides = 10},
         {.start = 29, .stride = 1, .nstrides = 10}},
      search_stripes[num_search_stripes]
        = {{.start = 32, .stride = -1, .nstrides = 30 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
                     = {{.start = 23, .size = -22}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "search dec stripe over 2 inc gap");
  }

  {
    enum {
      num_stripes = 5,
      num_search_stripes = 1,
      num_ref_pos_ext = 5,
      num_ref_unmatched = 0,
    };
    static const struct Xt_stripe stripes[num_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 10},
         {.start = 15, .stride = 1, .nstrides = 10},
         {.start = 29, .stride = 1, .nstrides = 10},
         {.start = 14, .stride = -1, .nstrides = 4},
         {.start = 28, .stride = -1, .nstrides = 4}},
      search_stripes[num_search_stripes]
        = {{.start = 32, .stride = -1, .nstrides = 30 }};
    static const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext]
      = {{.start = 23, .size = -4},{.start = 34, .size = 4},
         {.start = 19, .size = -10},{.start = 30, .size = 4},
         {.start = 9, .size = -8}};

    check_pos_ext(num_stripes, stripes,
                  num_search_stripes, search_stripes,
                  num_ref_pos_ext, ref_pos_ext, 1,
                  num_ref_unmatched, "search dec stripe over jumbled stripes");
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
do_tests(Xt_idxlist idxlist, const Xt_int *ref_indices, int num_indices) {

  check_idxlist(idxlist, ref_indices, num_indices);

  struct Xt_stripe * stripes;
  int num_stripes;
  Xt_idxlist temp_idxlist;

  xt_idxlist_get_index_stripes(idxlist, &stripes, &num_stripes);
  temp_idxlist = xt_idxvec_from_stripes_new(stripes, num_stripes);

  check_idxlist(temp_idxlist, ref_indices, num_indices);

  xt_idxlist_delete(temp_idxlist);

  free(stripes);

  {
    // test packing and unpacking
    Xt_idxlist idxlist_copy
      = idxlist_pack_unpack_copy(idxlist);

    // check copy
    check_idxlist(idxlist_copy, ref_indices, num_indices);

    // clean up
    xt_idxlist_delete(idxlist_copy);
  }

  {
    // test copying
    Xt_idxlist idxlist_copy = xt_idxlist_copy(idxlist);

    // check copy
    check_idxlist(idxlist_copy, ref_indices, num_indices);

    // clean up
    xt_idxlist_delete(idxlist_copy);
  }
}

static void
stripe_test_general(int num_stripes, const struct Xt_stripe *stripes,
                    int num_ref_indices, const Xt_int *ref_indices)
{
  Xt_idxlist idxstripes = xt_idxstripes_new(stripes, num_stripes);
  do_tests(idxstripes, ref_indices, num_ref_indices);
  int num_ext;
  struct Xt_pos_ext *pos_ext;
  xt_idxlist_get_pos_exts_of_index_stripes(idxstripes, num_stripes, stripes,
                                           &num_ext, &pos_ext, 1);
  size_t num_pos = 0;
  for (int i = 0; i < num_ext; ++i)
  {
    size_t ext_size = (size_t)pos_ext[i].size;
    if (num_pos != (size_t)pos_ext[i].start)
      PUT_ERR("position/start mismatch at extent %d\n", i);
    num_pos += ext_size;
  }
  if (num_pos != (size_t)xt_idxlist_get_num_indices(idxstripes))
    PUT_ERR("index list length/positions overlap mismatch\n");

  free(pos_ext);
  xt_idxlist_delete(idxstripes);
  Xt_idxlist idxvec = xt_idxvec_new(ref_indices, num_ref_indices);
  idxstripes = xt_idxstripes_from_idxlist_new(idxvec);
  check_idxlist(idxstripes, ref_indices, num_ref_indices);
  xt_idxlist_delete(idxvec);
  xt_idxlist_delete(idxstripes);
}

static void
stripe_test_intersection(int num_stripes_a, const struct Xt_stripe *stripes_a,
                         int num_stripes_b, const struct Xt_stripe *stripes_b,
                         int num_ref_indices, const Xt_int *ref_indices)
{
  Xt_idxlist idxstripes_a = xt_idxstripes_new(stripes_a, num_stripes_a),
    idxstripes_b = xt_idxstripes_new(stripes_b, num_stripes_b);

  // compute intersections
  Xt_idxlist intersection[2]
    = { xt_idxlist_get_intersection(idxstripes_a, idxstripes_b),
        xt_idxlist_get_intersection(idxstripes_b, idxstripes_a) };

  // check intersection
  do_tests(intersection[0], ref_indices, num_ref_indices);
  do_tests(intersection[1], ref_indices, num_ref_indices);

  // clean up
  xt_idxlist_delete(idxstripes_a);
  xt_idxlist_delete(idxstripes_b);
  xt_idxlist_delete(intersection[0]);
  xt_idxlist_delete(intersection[1]);
}

static void
stripe_test_asymmetric_intersection(
  int num_stripes_a, const struct Xt_stripe *stripes_a,
  int num_stripes_b, const struct Xt_stripe *stripes_b,
  int num_ref_indices_a, const Xt_int *ref_indices_a,
  int num_ref_indices_b, const Xt_int *ref_indices_b)
{
  Xt_idxlist idxstripes_a = xt_idxstripes_new(stripes_a, num_stripes_a),
    idxstripes_b = xt_idxstripes_new(stripes_b, num_stripes_b);

  // compute intersection
  Xt_idxlist intersection[2]
    = { xt_idxlist_get_intersection(idxstripes_a, idxstripes_b),
        xt_idxlist_get_intersection(idxstripes_b, idxstripes_a) };

  // check intersection
  do_tests(intersection[0], ref_indices_a, num_ref_indices_a);
  do_tests(intersection[1], ref_indices_b, num_ref_indices_b);

  // clean up
  xt_idxlist_delete(idxstripes_a);
  xt_idxlist_delete(idxstripes_b);
  xt_idxlist_delete(intersection[0]);
  xt_idxlist_delete(intersection[1]);
}

/* test whether
 *  xt_idxlist_get_index_at_position and
 *  xt_idxlist_get_indices_at_positions
 * give the same results
 */
static void
check_idxvec_get_indices_at_positions(int num_stripes,
                                      const struct Xt_stripe *stripes,
                                      int num_pos, const int *pos)
{
  static const Xt_int undef_idx = XT_INT_MIN;
  Xt_idxlist idxlist = xt_idxstripes_new(stripes, num_stripes);
  Xt_int ref_sel_idx[num_pos];
  int ref_undef_count = 0;
  for (int i=0; i<num_pos; i++) {
    int p = pos[i];
    if (xt_idxlist_get_index_at_position(idxlist, p, &ref_sel_idx[i]) != 0) {
      ref_sel_idx[i] = undef_idx;
      ref_undef_count++;
    }
  }
  Xt_int sel_idx[num_pos];
  int undef_count = xt_idxlist_get_indices_at_positions(idxlist, pos, num_pos,
                                                        sel_idx, undef_idx);
  if (undef_count != ref_undef_count)
    PUT_ERR("test_idxvec.c: (undef_count != ref_undef_count)\n");
  for (int i=0; i<num_pos; i++) {
    if (sel_idx[i] != ref_sel_idx[i])
      PUT_ERR("test_idxvec.c: (sel_idx[%d] != ref_sel_idx[%d])\n", i, i);
  }
  xt_idxlist_delete(idxlist);
}

static void
check_bb(int num_stripes, const struct Xt_stripe *stripes,
         unsigned ndim, const Xt_int global_size[ndim],
         const struct Xt_bounds ref_bounds[ndim],
         Xt_int global_start_index)
{
  Xt_idxlist idxstripes = xt_idxstripes_new(stripes, num_stripes);
  struct Xt_bounds bounds[ndim];

  xt_idxlist_get_bounding_box(idxstripes, ndim, global_size,
                              global_start_index, bounds);

  for (unsigned i = 0; i < ndim; ++i)
    if ((bounds[i].size != ref_bounds[i].size) |
        ((ref_bounds[i].size != 0) & (bounds[i].start != ref_bounds[i].start)))
      PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

  xt_idxlist_delete(idxstripes);
}

static void
check_pos_ext(size_t num_stripes, const struct Xt_stripe stripes[num_stripes],
              size_t num_search_stripes,
              const struct Xt_stripe search_stripes[num_search_stripes],
              size_t num_ref_pos_ext,
              const struct Xt_pos_ext ref_pos_ext[num_ref_pos_ext],
              int single_match_only, int ref_unmatched,
              const char *test_desc)
{
  Xt_idxlist idxstripes = xt_idxstripes_new(stripes, (int)num_stripes);
  int num_ext;
  struct Xt_pos_ext *pos_ext;
  int unmatched
    = xt_idxlist_get_pos_exts_of_index_stripes(idxstripes,
                                               (int)num_search_stripes,
                                               search_stripes,
                                               &num_ext, &pos_ext,
                                               single_match_only);
  (void)test_desc;
  if (unmatched != ref_unmatched)
    PUT_ERR("error in number of unmatched indices for %s", test_desc);
  if (num_ext < 0 || (size_t)num_ext != num_ref_pos_ext)
    PUT_ERR("error finding %s\n", test_desc);
  for (size_t i = 0; i < num_ref_pos_ext; ++i) {
    if (pos_ext[i].start != ref_pos_ext[i].start)
      PUT_ERR("incorrect starting position found in %s\n", test_desc);
    if (pos_ext[i].size != ref_pos_ext[i].size)
      PUT_ERR("incorrect position extent length found in %s\n", test_desc);
  }
  free(pos_ext);
  xt_idxlist_delete(idxstripes);
}
