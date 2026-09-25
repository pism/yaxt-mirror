/**
 * @file test_idxsection.c
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

#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#include <mpi.h>

#include <yaxt.h>

#include "core/ppm_xfuncs.h"

#define VERBOSE
#include "tests.h"
#include "ctest_common.h"
#include "test_idxlist_utils.h"


static void
do_tests(Xt_idxlist idxlist, const Xt_int *ref_indices, int num_indices,
         const struct Xt_stripe *ref_stripes, int ref_num_stripes);

static void
test_1d_section(void);
static void
test_2d_section(void);
static void
test_3d_section1(void);
static void
test_3d_section2(void);
static void
test_4d_section(void);
static void
test_2d_simple(void);
static void
test_1d_intersection1(void);
static void
test_1d_intersection2(void);
static void
test_2d_intersection1(void);
static void
test_2d_1(void);
static void
test_2d_2(void);

static void
test_get_positions1(void);
static void
test_get_positions2(void);
static void
test_get_positions3(void);
static void
test_get_positions4(void);
static void
test_get_positions5(void);
static void
test_get_positions6(void);
static void
test_other_intersection(void);
static void
test_stripe_intersection(void);
static void
test_sec4d_stripe_intersection(void);
static void
test_signed_sizes1(void);
static void
test_signed_sizes2(void);
static void
test_signed_size_intersections(void);
static void
test_signed_size_positions(void);
static void
test_section_with_stride1(void);
static void
test_section_with_stride2(void);
static void
test_bounding_box1(void);
static void
test_bounding_box2(void);
static void
test_bounding_box3(void);

int main(int argc, char **argv) {

  test_init_mpi(&argc, &argv, MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  test_1d_section();
  test_2d_section();
  test_3d_section1();
  test_3d_section2();
  test_4d_section();
  test_2d_simple();
  test_1d_intersection1();
  test_1d_intersection2();
  test_2d_intersection1();
  test_2d_1();
  test_2d_2();

  test_get_positions1();
  test_get_positions2();
  test_get_positions3();
  test_get_positions4();
  test_get_positions5();
  test_get_positions6();

  test_other_intersection();
  test_stripe_intersection();
  test_sec4d_stripe_intersection();

  test_signed_sizes1();
  test_signed_sizes2();

  test_signed_size_intersections();
  test_signed_size_positions();
  test_section_with_stride1();
  test_section_with_stride2();
  test_bounding_box1();
  test_bounding_box2();
  test_bounding_box3();

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void
do_tests(Xt_idxlist idxlist, const Xt_int *ref_indices, int num_indices,
         const struct Xt_stripe *ref_stripes, int ref_num_stripes) {

  check_idxlist(idxlist, ref_indices, num_indices);

  struct Xt_stripe * stripes;
  int num_stripes;

  xt_idxlist_get_index_stripes(idxlist, &stripes, &num_stripes);

  check_stripes(stripes, num_stripes, ref_stripes, ref_num_stripes);

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

  { // test copying

    Xt_idxlist idxlist_copy;

    idxlist_copy = xt_idxlist_copy(idxlist);

    // check copy

    check_idxlist(idxlist_copy, ref_indices, num_indices);

    // clean up

    xt_idxlist_delete(idxlist_copy);
  }
}

static void
test_1d_section(void)
{
  Xt_int start = 0;
  enum { num_dimensions = 1 };
  static const Xt_int global_size[num_dimensions] = {10};
  static const int local_size [num_dimensions] = {5};
  static const Xt_int local_start[num_dimensions] = {3};

  // create index section
  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  // testing
  static const Xt_int ref_indices[5] = {3,4,5,6,7};
  static const struct Xt_stripe ref_stripes[1]
    = {{.start = 3, .stride = 1, .nstrides = 5}};
  do_tests(idxsection, ref_indices, 5, ref_stripes, 1);

  // clean up
  xt_idxlist_delete(idxsection);
}

static void
test_2d_section(void)
{
  Xt_int start = 0;
  enum { num_dimensions = 2 };
  static const Xt_int global_size[num_dimensions] = {5,6};
  static const int local_size [num_dimensions] = {3,2};
  static const Xt_int local_start[num_dimensions] = {1,2};

  // create index section
  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  // testing
  static const Xt_int ref_indices[6] = {8,9,14,15,20,21};
  static const struct Xt_stripe ref_stripes[3]
    = {{.start = 8,  .stride = 1, .nstrides = 2},
       {.start = 14, .stride = 1, .nstrides = 2},
       {.start = 20, .stride = 1, .nstrides = 2}};
  do_tests(idxsection, ref_indices, 6, ref_stripes, 3);

  // clean up
  xt_idxlist_delete(idxsection);
}

static void
test_3d_section1(void)
{
  Xt_int start = 0;
  enum { num_dimensions = 3 };
  static const Xt_int global_size[num_dimensions] = {4,4,4};
  static const int local_size[num_dimensions] = {4,2,2};
  static const Xt_int local_start[num_dimensions] = {0,1,1};

  // create index section
  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  // testing
  static const Xt_int ref_indices[16]
    = {5,6,9,10, 21,22,25,26, 37,38,41,42, 53,54,57,58};
  static const struct Xt_stripe ref_stripes[8]
    = {{.start =  5, .stride = 1, .nstrides = 2},
       {.start =  9, .stride = 1, .nstrides = 2},
       {.start = 21, .stride = 1, .nstrides = 2},
       {.start = 25, .stride = 1, .nstrides = 2},
       {.start = 37, .stride = 1, .nstrides = 2},
       {.start = 41, .stride = 1, .nstrides = 2},
       {.start = 53, .stride = 1, .nstrides = 2},
       {.start = 57, .stride = 1, .nstrides = 2}};
  do_tests(idxsection, ref_indices, 16, ref_stripes, 8);

  // clean up
  xt_idxlist_delete(idxsection);
}

static void
test_3d_section2(void)
{
  Xt_int start = 0;
  enum { num_dimensions = 3 };
  static const Xt_int global_size[num_dimensions] = {3,4,5};
  static const int local_size[num_dimensions] = {2,2,3};
  static const Xt_int local_start[num_dimensions] = {1,1,1};

  // create index section
  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  // testing
  static const Xt_int ref_indices[12] = {26,27,28,31,32,33,46,47,48,51,52,53};
  static const struct Xt_stripe ref_stripes[4]
    = {{.start = 26, .stride = 1, .nstrides = 3},
       {.start = 31, .stride = 1, .nstrides = 3},
       {.start = 46, .stride = 1, .nstrides = 3},
       {.start = 51, .stride = 1, .nstrides = 3}};
  do_tests(idxsection, ref_indices, 12, ref_stripes, 4);

  // clean up
  xt_idxlist_delete(idxsection);
}

static void
test_4d_section(void)
{
  Xt_int start = 0;
  enum { num_dimensions = 4 };
  static const Xt_int global_size[4] = {3,4,4,3};
  static const int local_size[4] = {2,3,3,2};
  static const Xt_int local_start[4] = {0,1,1,1};

  // create index section

  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  // testing

  static const Xt_int ref_indices[36]
    = {16,17,19,20,22,23, 28,29,31,32,34,35, 40,41,43,44,46,47,
       64,65,67,68,70,71, 76,77,79,80,82,83, 88,89,91,92,94,95};
  static const struct Xt_stripe ref_stripes[18]
    = {{.start = 16, .stride = 1, .nstrides = 2},
       {.start = 19, .stride = 1, .nstrides = 2},
       {.start = 22, .stride = 1, .nstrides = 2},
       {.start = 28, .stride = 1, .nstrides = 2},
       {.start = 31, .stride = 1, .nstrides = 2},
       {.start = 34, .stride = 1, .nstrides = 2},
       {.start = 40, .stride = 1, .nstrides = 2},
       {.start = 43, .stride = 1, .nstrides = 2},
       {.start = 46, .stride = 1, .nstrides = 2},
       {.start = 64, .stride = 1, .nstrides = 2},
       {.start = 67, .stride = 1, .nstrides = 2},
       {.start = 70, .stride = 1, .nstrides = 2},
       {.start = 76, .stride = 1, .nstrides = 2},
       {.start = 79, .stride = 1, .nstrides = 2},
       {.start = 82, .stride = 1, .nstrides = 2},
       {.start = 88, .stride = 1, .nstrides = 2},
       {.start = 91, .stride = 1, .nstrides = 2},
       {.start = 94, .stride = 1, .nstrides = 2}};
  do_tests(idxsection, ref_indices, 36, ref_stripes, 18);

  // clean up

  xt_idxlist_delete(idxsection);
}

static void
test_2d_simple(void)
{
  Xt_int start = 0;
  enum { num_dimensions = 2 };
  static const Xt_int global_size[num_dimensions] = {5,10};
  static const int local_size [num_dimensions] = {3,4};
  static const Xt_int local_start[num_dimensions] = {1,2};

  // create index section

  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  // testing

  static const Xt_int ref_indices[12] = {12,13,14,15,22,23,24,25,32,33,34,35};

  check_idxlist(idxsection, ref_indices, 12);

  // clean up

  xt_idxlist_delete(idxsection);
}

static void
test_intersection(
  Xt_int start_a, int num_dimensions_a,
  const Xt_int *global_size_a, const int *local_size_a,
  const Xt_int *local_start_a,
  Xt_int start_b, int num_dimensions_b,
  const Xt_int *global_size_b, const int *local_size_b,
  const Xt_int *local_start_b,
  size_t num_ref_indices, const Xt_int ref_indices[num_ref_indices],
  size_t num_ref_stripes, const struct Xt_stripe ref_stripes[num_ref_stripes])
{
  Xt_idxlist idxsection_a
    = xt_idxsection_new(start_a, num_dimensions_a, global_size_a,
                        local_size_a, local_start_a),
    idxsection_b
    = xt_idxsection_new(start_b, num_dimensions_b, global_size_b,
                        local_size_b, local_start_b);
  // compute intersection
  Xt_idxlist intersection
    = xt_idxlist_get_intersection(idxsection_a, idxsection_b);

  // check intersection
  do_tests(intersection, ref_indices, (int)num_ref_indices,
           ref_stripes, (int)num_ref_stripes);

  xt_idxlist_delete(intersection);
  xt_idxlist_delete(idxsection_a);
  xt_idxlist_delete(idxsection_b);
}

static void
test_1d_intersection1(void)
{
  Xt_int start = 0;
  enum { num_dimensions = 1 };
  static const Xt_int global_size_a[num_dimensions] = {10};
  static const int local_size_a [num_dimensions] = {5};
  static const Xt_int local_start_a[num_dimensions] = {4};

  static const Xt_int global_size_b[num_dimensions] = {15};
  static const int local_size_b [num_dimensions] = {6};
  static const Xt_int local_start_b[num_dimensions] = {7};

  // create index sections
  enum {
    num_ref_indices = 2,
    num_ref_stripes = 1,
  };
  static const Xt_int ref_indices[num_ref_indices] = {7,8};
  static const struct Xt_stripe ref_stripes[num_ref_stripes]
    = {{.start = 7, .stride = 1, .nstrides = 2}};
  test_intersection(
    start, num_dimensions, global_size_a, local_size_a, local_start_a,
    start, num_dimensions, global_size_b, local_size_b, local_start_b,
    num_ref_indices, ref_indices,
    num_ref_stripes, ref_stripes);
}

static void
test_1d_intersection2(void)
{
  // 1D intersection test (with empty intersection)
  Xt_int start = 0;
  enum { num_dimensions = 1 };
  static const Xt_int global_size_a[num_dimensions] = {10};
  static const int local_size_a [num_dimensions] = {1};
  static const Xt_int local_start_a[num_dimensions] = {3};

  static const Xt_int global_size_b[num_dimensions] = {10};
  static const int local_size_b [num_dimensions] = {5};
  static const Xt_int local_start_b[num_dimensions] = {4};

  enum {
    num_ref_indices = 0,
    num_ref_stripes = 0,
  };
  static const Xt_int ref_indices[1] = {0};
  static const struct Xt_stripe ref_stripes[1]
      = {{.start = 0, .stride = 1, .nstrides = 1}};

  test_intersection(
    start, num_dimensions, global_size_a, local_size_a, local_start_a,
    start, num_dimensions, global_size_b, local_size_b, local_start_b,
    num_ref_indices, ref_indices,
    num_ref_stripes, ref_stripes);
}

static void
test_2d_intersection1(void)
{
  Xt_int start = 0;
  enum { num_dimensions = 2 };
  static const Xt_int global_size_a[num_dimensions] = {6,6};
  static const int local_size_a [num_dimensions] = {4,2};
  static const Xt_int local_start_a[num_dimensions] = {1,1};

  static const Xt_int global_size_b[num_dimensions] = {6,6};
  static const int local_size_b [num_dimensions] = {3,3};
  static const Xt_int local_start_b[num_dimensions] = {3,2};

  // check intersection
  enum {
    num_ref_indices = 2,
    num_ref_stripes = 1,
  };
  static const Xt_int ref_indices[2] = {20,26};
  static const struct Xt_stripe ref_stripes[1]
    = {{.start = 20, .stride = 6, .nstrides = 2}};
  test_intersection(
    start, num_dimensions, global_size_a, local_size_a, local_start_a,
    start, num_dimensions, global_size_b, local_size_b, local_start_b,
    num_ref_indices, ref_indices,
    num_ref_stripes, ref_stripes);
}

static void
test_2d_1(void)
{
  Xt_int start = 0;
  enum { num_dimensions = 2 };
  static const Xt_int global_size[num_dimensions] = {4,4};
  static const int local_size [num_dimensions] = {2,2};
  static const Xt_int local_start[num_dimensions] = {0,2};

  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);
  enum { num_ref_indices = 4 };
  static const Xt_int ref_indices[num_ref_indices] = {2,3,6,7};

  check_idxlist(idxsection, ref_indices, num_ref_indices);

  xt_idxlist_delete(idxsection);
}

static void
test_2d_2(void)
{
  Xt_int start = 1;
  enum { num_dimensions = 2 };
  static const Xt_int global_size[num_dimensions] = {4,4};
  static const int local_size [num_dimensions] = {2,2};
  static const Xt_int local_start[num_dimensions] = {0,2};

  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  enum { num_ref_indices = 4 };
  static const Xt_int ref_indices[num_ref_indices] = {3,4,7,8};

  check_idxlist(idxsection, ref_indices, num_ref_indices);

  xt_idxlist_delete(idxsection);
}

static void
test_get_positions_common(
  int num_dimensions,
  Xt_int start,
  const Xt_int global_size[num_dimensions],
  const int local_size[num_dimensions],
  const Xt_int local_start[num_dimensions],
  bool single_match_only,
  int num_selection,
  const Xt_int selection[num_selection],
  const int ref_positions[num_selection])
{
  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  int positions[num_selection];
  int count_missing_pos_ref = 0;
  for (int i = 0; i < num_selection; ++i)
    count_missing_pos_ref += ref_positions[i] == -1;
  if (xt_idxlist_get_positions_of_indices(
        idxsection, selection,
        num_selection, positions, single_match_only) != count_missing_pos_ref)
    PUT_ERR("xt_idxlist_get_positions_of_indices returned incorrect"
            " num_unmatched\n");

  bool mismatch = false;
  for (int i=0; i<num_selection; i++)
    mismatch |= (positions[i] != ref_positions[i]);
  if (mismatch)
    PUT_ERR("xt_idxlist_get_positions_of_indices returned incorrect"
            " position\n");

  if (!single_match_only)
    for (int i=0; i<num_selection; i++) {
      int p;
      xt_idxlist_get_position_of_index(idxsection, selection[i], &p);
      if (p != ref_positions[i])
        PUT_ERR("xt_idxlist_get_position_of_index returned incorrect"
                " position\n");
    }

  xt_idxlist_delete(idxsection);
}


static void
test_get_positions1(void)
{
  // check get_positions_of_indices
  Xt_int start = 0;
  enum { num_dimensions = 2 };
  static const Xt_int global_size[num_dimensions] = {4,4};
  static const int local_size [num_dimensions] = {2,2};
  static const Xt_int local_start[num_dimensions] = {0,2};
  // we have indices = {2,3,6,7}
  static const Xt_int selection[] = {1,2,5,6,7,8};
  enum { num_selection = sizeof(selection) / sizeof(*selection) };

  static const int ref_positions[]
    = {1*0-1, 2*0+0, 5*0-1, 6*0+2, 7*0+3, 8*0-1};

  test_get_positions_common(num_dimensions, start, global_size,
                            local_size, local_start, 0, num_selection,
                            selection, ref_positions);
}

static void
test_get_positions2(void)
{
  // check get_positions_of_indices with single_match_only = 0
  static const Xt_int start = 0;
  enum { num_dimensions = 2 };
  static const Xt_int global_size[num_dimensions] = {4,4};
  static const int local_size[num_dimensions] = {2,2};
  static const Xt_int local_start[num_dimensions] = {0,2};
  // we have indices = {2,3,6,7}
  static const Xt_int selection[] = {2,1,5,7,6,7,7,6,8};
  enum { num_selection = sizeof(selection) / sizeof(*selection) };

  static const int ref_positions[]
    = {2*0+0, 1*0-1, 5*0-1, 7*0+3, 6*0+2, 7*0+3, 7*0+3, 6*0+2, 8*0-1};

  test_get_positions_common(num_dimensions, start, global_size,
                            local_size, local_start, 0, num_selection,
                            selection, ref_positions);
}

static void
test_get_positions3(void)
{
  // check get_positions_of_indices with single_match_only = 1
  static const Xt_int start = 0;
  enum { num_dimensions = 2 };
  static const Xt_int global_size[num_dimensions] = {4,4};
  static const int local_size [num_dimensions] = {2,2};
  static const Xt_int local_start[num_dimensions] = {0,2};

  // we have indices = {2,3,6,7}
  static const Xt_int selection[] = {2,1,5,7,6,7,7,6,8};
  enum { num_selection = sizeof(selection) / sizeof(*selection) };
  static const int ref_positions[]
    = {2*0+0, 1*0-1, 5*0-1, 7*0+3, 6*0+2, 7*0-1, 7*0-1, 6*0-1, 8*0-1};

  test_get_positions_common(num_dimensions, start, global_size,
                            local_size, local_start, 1, num_selection,
                            selection, ref_positions);
}

static void
test_get_positions4(void)
{
  // check get_positions_of_indices with single_match_only = 1
  // for mesh with one dimension enumerated in descending fashion
  Xt_int start = 0;
  enum { num_dimensions = 2 };
  static const Xt_int global_size[num_dimensions] = {4,4};
  static const int local_size [num_dimensions] = {2,-2};
  static const Xt_int local_start[num_dimensions] = {0,2};

  // we have indices = {3,2,7,6}
  static const Xt_int selection[] = {2,1,5,7,6,7,7,6,8};
  enum { num_selection = sizeof(selection) / sizeof(*selection) };

  static const int ref_positions[]
    = {2*0+1, 1*0-1, 5*0-1, 7*0+2, 6*0+3, 7*0-1, 7*0-1, 6*0-1, 8*0-1};

  test_get_positions_common(num_dimensions, start, global_size,
                            local_size, local_start, 1, num_selection,
                            selection, ref_positions);
}

static void
test_get_positions5(void)
{
  // check get_positions_of_indices with single_match_only = 1
  // for mesh with dimension enumerated in alternating fashion
  // to create an unsorted case
  // then lookup positions outside the valid ranges
  Xt_int start = 3;
  enum { num_dimensions = 3 };
  static const Xt_int global_size[num_dimensions] = {-3,5,-1};
  static const int local_size [num_dimensions] = {-3,-3,1};
  static const Xt_int local_start[num_dimensions] = {0,1,0};

  // we have indices = {...}
  static const Xt_int selection[] = { 12, 264 };
  enum { num_selection = sizeof(selection) / sizeof(*selection) };

  static const int ref_positions[]
    = {12*0-1, 264*0-1};

  test_get_positions_common(num_dimensions, start, global_size,
                            local_size, local_start, 1, num_selection,
                            selection, ref_positions);
}

static void
test_get_positions6(void)
{
  Xt_int start = 3;
  enum { num_dimensions = 2 };
  static const Xt_int global_size[2] = {3,-3};
  static const Xt_int local_start[2] = {0,0};
  static const int local_size[2] = {-3,-2};
  /* the above mesh represents indices: 10 11 7 8 4 5 */
  static const Xt_int selection[] = { 6 };
  static const int ref_positions[] = { -1 };
  test_get_positions_common(num_dimensions, start, global_size,
                            local_size, local_start, 1, 1,
                            selection, ref_positions);
}


static void
test_other_intersection(void)
{
  // check idxsection_get_intersection_with_other_idxlist
  Xt_int start = 0;
  enum { num_dimensions = 2 };
  static const Xt_int global_size[num_dimensions] = {4,4};
  static const int local_size [num_dimensions] = {2,2};
  static const Xt_int local_start[num_dimensions] = {0,2};

  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  // we have indices = {2,3,6,7}
  static const Xt_int sel_idx[] = {2,1,5,7,6,7,7,6,8};
  enum { num_sel_idx = sizeof(sel_idx) / sizeof(*sel_idx) };

  Xt_idxlist sel_idxlist = xt_idxvec_new(sel_idx, num_sel_idx);

  Xt_idxlist inter_idxlist
    = xt_idxlist_get_intersection(idxsection, sel_idxlist);

  static const Xt_int ref_inter_idx[] = {2,6,6,7,7,7};
  enum { num_ref_inter_idx = sizeof(ref_inter_idx) / sizeof(*ref_inter_idx) };

  xt_idxlist_delete(sel_idxlist);
  xt_idxlist_delete(idxsection);
  check_idxlist(inter_idxlist, ref_inter_idx, num_ref_inter_idx);
  xt_idxlist_delete(inter_idxlist);
}

static void
test_stripe_intersection(void)
{
  // check idxsection_get_idxstripes_intersection
  Xt_int start = 0;
  enum { num_dimensions = 2 };
  static const Xt_int global_size[num_dimensions] = {4,4};
  static const int local_size [num_dimensions] = {2,2};
  static const Xt_int local_start[num_dimensions] = {0,2};

  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  // we have indices = {2,1,5,7,6,7,7,6,8}
  static const struct Xt_stripe stripes[]
    = { { .start = 2, .stride = -1, .nstrides = 2 },
        { .start = 5, .stride = 2, .nstrides = 2 },
        { .start = 6, .stride = 1, .nstrides = 2 },
        { .start = 7, .stride = -1, .nstrides = 2 },
        { .start = 8, .stride = 1, .nstrides = 1 },
  };
  enum { num_stripes = sizeof(stripes) / sizeof(stripes[0]) };

  Xt_idxlist sel_idxlist = xt_idxstripes_new(stripes, num_stripes);

  Xt_idxlist inter_idxlist
    = xt_idxlist_get_intersection(idxsection, sel_idxlist);

  static const Xt_int ref_inter_idx[] = {2,6,6,7,7,7};
  enum { num_ref_inter_idx = sizeof(ref_inter_idx) / sizeof(*ref_inter_idx) };

  check_idxlist(inter_idxlist, ref_inter_idx, num_ref_inter_idx);
  xt_idxlist_delete(inter_idxlist);

  Xt_idxlist inter_idxlist_r
    = xt_idxlist_get_intersection(sel_idxlist, idxsection);
  xt_idxlist_delete(sel_idxlist);
  xt_idxlist_delete(idxsection);

  static const Xt_int ref_inter_idx_r[] = {2,6,7};
  enum { num_ref_inter_idx_r = sizeof(ref_inter_idx_r) / sizeof(*ref_inter_idx_r) };

  check_idxlist(inter_idxlist_r, ref_inter_idx_r, num_ref_inter_idx_r);

  xt_idxlist_delete(inter_idxlist_r);
}

static void
test_sec4d_stripe_intersection(void)
{
  // test intersection of idxsection with non-zero offset and
  // stripes
  enum { irank = 4, num_ref_indices=160, bias=1};
  static const Xt_int gsz[irank] = { 1, 1, 16, 19 };
  static const int lsz[irank] = { 1, 1, 16, 10 };
  static const Xt_int lst[irank] = { 0, 0, 0, 0 };
  Xt_int ref_inter_idx[num_ref_indices];

  for (Xt_int j = 0; j < lsz[2]; ++j)
    for (Xt_int i = 0; i < lsz[3]; ++i)
      ref_inter_idx[i+j*lsz[3]] = (Xt_int)(i + bias + j*19);

  Xt_idxlist section = xt_idxsection_new(bias, irank, gsz, lsz, lst),
    stripes = xt_idxstripes_new((struct Xt_stripe[1]){ { 0, 1, 304 } }, 1);

  Xt_idxlist intersect
    = xt_idxlist_get_intersection(section, stripes);
  xt_idxlist_delete(stripes);
  xt_idxlist_delete(section);
  check_idxlist(intersect, ref_inter_idx, num_ref_indices);
  xt_idxlist_delete(intersect);
}


static void
test_signed_sizes1(void)
{
  // test 2D section with negative global size
  // iterate through all sign combinations of -/+ for local and global
  // and for x and y, giving 2^2^2 combinations
  for (int i = 0; i < 16; ++i) {

    Xt_int start = 0;
    enum { num_dimensions = 2 };
    static const Xt_int global_size[4][num_dimensions]
      = {{7,10},{7,-10},{-7,10},{-7,-10}};
    static const int local_size[4][num_dimensions]
      = {{3,4},{3,-4},{-3,4},{-3,-4}};
    static const Xt_int local_start[num_dimensions] = {1,2};

    // create index section

    Xt_idxlist idxsection
      = xt_idxsection_new(start, num_dimensions, global_size[i >> 2],
                          local_size[i & 3], local_start);

    // testing

    static const Xt_int ref_indices[16][12] =
      {{12, 13, 14, 15, 22, 23, 24, 25, 32, 33, 34, 35},
       {15, 14, 13, 12, 25, 24, 23, 22, 35, 34, 33, 32},
       {32, 33, 34, 35, 22, 23, 24, 25, 12, 13, 14, 15},
       {35, 34, 33, 32, 25, 24, 23, 22, 15, 14, 13, 12},
       {17, 16, 15, 14, 27, 26, 25, 24, 37, 36, 35, 34},
       {14, 15, 16, 17, 24, 25, 26, 27, 34, 35, 36, 37},
       {37, 36, 35, 34, 27, 26, 25, 24, 17, 16, 15, 14},
       {34, 35, 36, 37, 24, 25, 26, 27, 14, 15, 16, 17},
       {52, 53, 54, 55, 42, 43, 44, 45, 32, 33, 34, 35},
       {55, 54, 53, 52, 45, 44, 43, 42, 35, 34, 33, 32},
       {32, 33, 34, 35, 42, 43, 44, 45, 52, 53, 54, 55},
       {35, 34, 33, 32, 45, 44, 43, 42, 55, 54, 53, 52},
       {57, 56, 55, 54, 47, 46, 45, 44, 37, 36, 35, 34},
       {54, 55, 56, 57, 44, 45, 46, 47, 34, 35, 36, 37},
       {37, 36, 35, 34, 47, 46, 45, 44, 57, 56, 55, 54},
       {34, 35, 36, 37, 44, 45, 46, 47, 54, 55, 56, 57}};

    check_idxlist(idxsection, ref_indices[i], 12);

    // clean up

    xt_idxlist_delete(idxsection);
  }
}

static void
test_signed_sizes2(void)
{
  // test 2D section with negative global size
  for (int i = 0; i < 16; ++i) {

    Xt_int start = 0;
    enum { num_dimensions = 2 };
    static const Xt_int global_size[4][num_dimensions]
      = {{5,6},{5,-6},{-5,6},{-5,-6}};
    static const int local_size [4][num_dimensions]
      = {{2,3},{2,-3},{-2,3},{-2,-3}};
    static const Xt_int local_start[num_dimensions] = {1,2};

    // create index section

    Xt_idxlist idxsection
      = xt_idxsection_new(start, num_dimensions, global_size[i >> 2],
                          local_size[i & 3], local_start);

    // testing

    static const Xt_int ref_indices[16][6] =
      {{8,9,10,14,15,16},
       {10,9,8,16,15,14},
       {14,15,16,8,9,10},
       {16,15,14,10,9,8},
       {9,8,7,15,14,13},
       {7,8,9,13,14,15},
       {15,14,13,9,8,7},
       {13,14,15,7,8,9},
       {20,21,22,14,15,16},
       {22,21,20,16,15,14},
       {14,15,16,20,21,22},
       {16,15,14,22,21,20},
       {21,20,19,15,14,13},
       {19,20,21,13,14,15},
       {15,14,13,21,20,19},
       {13,14,15,19,20,21}};

    check_idxlist(idxsection, ref_indices[i], 6);

    // clean up
    xt_idxlist_delete(idxsection);
  }
}

static void
test_signed_size_intersections(void)
{
  // test intersection of 2D section with negative global size
  for (int i = 0; i < 16; ++i) {
    for (int j = 0; j < 16; ++j) {

      Xt_int start = 0;
      enum { num_dimensions = 2 };
      static const Xt_int global_size[4][num_dimensions]
        = {{5,10},{5,-10},{-5,10},{-5,-10}};
      static const int local_size [4][num_dimensions]
        = {{3,4},{3,-4},{-3,4},{-3,-4}};
      static const Xt_int local_start[num_dimensions] = {1,2};

      static const Xt_int indices[16][12] =
        {{12, 13, 14, 15, 22, 23, 24, 25, 32, 33, 34, 35},
         {15, 14, 13, 12, 25, 24, 23, 22, 35, 34, 33, 32},
         {32, 33, 34, 35, 22, 23, 24, 25, 12, 13, 14, 15},
         {35, 34, 33, 32, 25, 24, 23, 22, 15, 14, 13, 12},
         {17, 16, 15, 14, 27, 26, 25, 24, 37, 36, 35, 34},
         {14, 15, 16, 17, 24, 25, 26, 27, 34, 35, 36, 37},
         {37, 36, 35, 34, 27, 26, 25, 24, 17, 16, 15, 14},
         {34, 35, 36, 37, 24, 25, 26, 27, 14, 15, 16, 17},
         {32, 33, 34, 35, 22, 23, 24, 25, 12, 13, 14, 15},
         {35, 34, 33, 32, 25, 24, 23, 22, 15, 14, 13, 12},
         {12, 13, 14, 15, 22, 23, 24, 25, 32, 33, 34, 35},
         {15, 14, 13, 12, 25, 24, 23, 22, 35, 34, 33, 32},
         {37, 36, 35, 34, 27, 26, 25, 24, 17, 16, 15, 14},
         {34, 35, 36, 37, 24, 25, 26, 27, 14, 15, 16, 17},
         {17, 16, 15, 14, 27, 26, 25, 24, 37, 36, 35, 34},
         {14, 15, 16, 17, 24, 25, 26, 27, 34, 35, 36, 37}};

      // create index section
      Xt_idxlist idxsection_a
        = xt_idxsection_new(start, num_dimensions, global_size[i >> 2],
                            local_size[i & 3], local_start),
        idxsection_b
        = xt_idxsection_new(start, num_dimensions, global_size[j >> 2],
                            local_size[j & 3], local_start);

      // create reference index vectors

      Xt_idxlist idxvec_a = xt_idxvec_new(indices[i], 12),
        idxvec_b = xt_idxvec_new(indices[j], 12);

      // also test with stripes now that those use specific
      // intersection routine
      Xt_idxlist idxstripes_b = xt_idxstripes_from_idxlist_new(idxvec_b);

      // testing

      Xt_idxlist idxsection_intersection
        = xt_idxlist_get_intersection(idxsection_a, idxsection_b),
        idxsection_intersection_stripes
        = xt_idxlist_get_intersection(idxsection_a, idxstripes_b),
        idxsection_intersection_other
        = xt_idxlist_get_intersection(idxsection_a, idxvec_b),
        idxvec_intersection = xt_idxlist_get_intersection(idxvec_a, idxvec_b);

      check_idxlist(idxsection_intersection,
                    xt_idxlist_get_indices_const(idxvec_intersection),
                    xt_idxlist_get_num_indices(idxvec_intersection));
      check_idxlist(idxsection_intersection_stripes,
                    xt_idxlist_get_indices_const(idxvec_intersection),
                    xt_idxlist_get_num_indices(idxvec_intersection));
      check_idxlist(idxsection_intersection_other,
                    xt_idxlist_get_indices_const(idxvec_intersection),
                    xt_idxlist_get_num_indices(idxvec_intersection));

      // clean up
      xt_idxlist_delete(idxvec_intersection);
      xt_idxlist_delete(idxsection_intersection_other);
      xt_idxlist_delete(idxsection_intersection_stripes);
      xt_idxlist_delete(idxsection_intersection);
      xt_idxlist_delete(idxstripes_b);
      xt_idxlist_delete(idxvec_b);
      xt_idxlist_delete(idxvec_a);
      xt_idxlist_delete(idxsection_b);
      xt_idxlist_delete(idxsection_a);
    }
  }
}

static void
test_signed_size_positions(void)
{
  // test 2D section with negative global size
  Xt_int start = 0;
  enum { num_dimensions = 2 };
  static const Xt_int global_size[num_dimensions] = {-5,6};
  static const int local_size [num_dimensions] = {-2,-3};
  static const Xt_int local_start[num_dimensions] = {1,2};

  // create index section
  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  // testing

  static const Xt_int ref_indices[6] = {16,15,14,22,21,20};

  check_idxlist(idxsection, ref_indices, 6);

  // check get_positions_of_indices
  enum {
    num_idxpos_q = 34,
  };
  static const Xt_int indices[num_idxpos_q]
    = {-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,14,16,17,18,19,
       20,20,21,22,23,24,25,26,27,28,29,30};
  static const int ref_positions[num_idxpos_q]
    = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,1,-1,0,-1,-1,-1,
       5,-1,4,3,-1,-1,-1,-1,-1,-1,-1,-1};
  int positions[num_idxpos_q];

  if (xt_idxlist_get_positions_of_indices(idxsection, indices, num_idxpos_q,
                                          positions, 1) != 28)
    PUT_ERR("error in xt_idxlist_get_positions_of_indices"
            " (wrong number of unmatched indices)\n");

  {
    bool mismatch = false;
    for (int i = 0; i < num_idxpos_q; ++i)
      mismatch |= (ref_positions[i] != positions[i]);
    if (mismatch)
      PUT_ERR("error in xt_idxlist_get_positions_of_indices"
              " (wrong position)\n");
  }
  {
    // replace not-found positions with other out-of-range values
    Xt_int answer_ref[num_idxpos_q];
    for (size_t i = 0; i < num_idxpos_q; ++i)
      if (positions[i] == -1) {
        positions[i] = 6+(int)i;
        answer_ref[i] = -2;
      } else
        answer_ref[i] = indices[i];

    Xt_int answer_idx[num_idxpos_q];
    xt_idxlist_get_indices_at_positions(idxsection, positions,
                                        num_idxpos_q, answer_idx,
                                        -2);
    bool mismatch = false;
    for (int i = 0; i < num_idxpos_q; ++i)
      mismatch |= (answer_idx[i] != answer_ref[i]);
    if (mismatch)
      PUT_ERR("error in xt_idxlist_get_indices_at_positions"
              " (wrong index)\n");
  }
  // clean up

  xt_idxlist_delete(idxsection);
}



static void
test_section_with_stride1(void)
{
  Xt_idxlist idxsection;
  Xt_int start = 0;
  enum { num_dimensions = 3 };
  static const Xt_int global_size[num_dimensions] = {5,5,2};
  static const int local_size [num_dimensions] = {3,4,1};
  static const Xt_int local_start[num_dimensions] = {2,0,1};

  // create index section
  idxsection = xt_idxsection_new(start, num_dimensions, global_size,
                                 local_size, local_start);

  // testing
  enum {
    num_ref_indices = 12,
    num_ref_stripes = 3,
  };
  static const Xt_int ref_indices[num_ref_indices]
    = {21, 23, 25, 27, 31, 33, 35, 37, 41, 43, 45, 47};
  static const struct Xt_stripe ref_stripes[num_ref_stripes]
    = { { 21, 2, 4 },
        { 31, 2, 4 },
        { 41, 2, 4 } };
  do_tests(idxsection, ref_indices, num_ref_indices,
           ref_stripes, num_ref_stripes);

  // clean up
  xt_idxlist_delete(idxsection);
}

static void
test_section_with_stride2(void)
{
  // test 2D section with stride in x and y
  Xt_int start = 0;
  enum { num_dimensions = 4 };
  static const Xt_int global_size[4] = {3,2,5,2};
  static const int local_size [4] = {3,1,4,1};
  static const Xt_int local_start[4] = {0,1,1,0};

  // create index section
  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  // testing

  static const Xt_int ref_indices[12]
    = {12, 14, 16, 18, 32, 34, 36, 38, 52, 54, 56, 58};

  check_idxlist(idxsection, ref_indices, 12);

  // clean up
  xt_idxlist_delete(idxsection);
}

static void
test_bounding_box1(void)
{
  Xt_int start = 0;
  enum { num_dimensions = 3 };
  static const Xt_int global_size[num_dimensions] = {4,4,4};
  static const int local_size [num_dimensions] = {0,0,0};
  static const Xt_int local_start[num_dimensions] = {2,0,1};

  // create index section
  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  enum { ndim = 3 };
  static const Xt_int global_size_bb[ndim] = { 4, 4, 4 };
  Xt_int global_start_index = 0;
  struct Xt_bounds bounds[ndim];

  xt_idxlist_get_bounding_box(idxsection, ndim, global_size_bb,
                              global_start_index, bounds);

  for (unsigned i = 0; i < ndim; ++i)
    if (bounds[i].size != 0)
      PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

  xt_idxlist_delete(idxsection);
}

static void
test_bounding_box2(void)
{
  Xt_int start = 1;
  enum { num_dimensions = 3 };
  static const Xt_int global_size[num_dimensions] = {5,4,3};
  static const int local_size [num_dimensions] = {2,2,2};
  static const Xt_int local_start[num_dimensions] = {2,2,1};

  // create index section
  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  enum { ndim = 3 };
  static const Xt_int global_size_bb[ndim] = { 5, 4, 3 };
  Xt_int global_start_index = 1;
  struct Xt_bounds bounds[ndim];

  xt_idxlist_get_bounding_box(idxsection, ndim, global_size_bb,
                              global_start_index, bounds);

  static const Xt_int ref_start[3] = {2,2,1};

  bool mismatch = false;
  for (int i = 0; i < ndim; ++i)
    mismatch |= (bounds[i].size != 2 || bounds[i].start != ref_start[i]);
  if (mismatch)
    PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

  xt_idxlist_delete(idxsection);
}

static void
test_bounding_box3(void)
{
  Xt_int start = 1;
  enum { num_dimensions = 4 };
  static const Xt_int global_size[num_dimensions] = {5,2,2,3};
  static const int local_size [num_dimensions] = {2,2,1,2};
  static const Xt_int local_start[num_dimensions] = {2,0,1,1};

  // create index section
  Xt_idxlist idxsection
    = xt_idxsection_new(start, num_dimensions, global_size,
                        local_size, local_start);

  enum { ndim = 3 };
  static const Xt_int global_size_bb[ndim] = { 5, 4, 3 };
  Xt_int global_start_index = 1;
  struct Xt_bounds bounds[ndim];

  xt_idxlist_get_bounding_box(idxsection, ndim, global_size_bb,
                              global_start_index, bounds);

  static const Xt_int ref_start[3] = {2,1,1};
  static const Xt_int ref_size[3] = {2,3,2};

  bool mismatch = false;
  for (int i = 0; i < ndim; ++i)
    mismatch |= (bounds[i].size != ref_size[i]
                 || bounds[i].start != ref_start[i]);
  if (mismatch)
    PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

  xt_idxlist_delete(idxsection);
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
