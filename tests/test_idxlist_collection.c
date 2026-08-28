/**
 * @file test_idxlist_collection.c
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
#include <string.h>

#include <mpi.h>
#include <yaxt.h>

#define VERBOSE

#include "tests.h"
#include "ctest_common.h"
#include "test_idxlist_utils.h"
#include "core/ppm_xfuncs.h"

static void
test_coll_of_coll(Xt_idxlist collectionlist);

static void
test_idxsections(void);

int main(int argc, char **argv) {

  test_init_mpi(&argc, &argv, MPI_COMM_WORLD);

  xt_initialize(MPI_COMM_WORLD);

  { // test packing and unpacking of collection index lists

    enum {
      num_indices = 7,
      num_ref_stripes = 2,
      num_lists = 2,
    };
    static const Xt_int index_list[num_lists][num_indices]
      = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1}};

    // create two index lists
    Xt_idxlist idxlists[num_lists];
    for (size_t i = 0; i < num_lists; ++i)
      idxlists[i] = xt_idxvec_new(index_list[i], num_indices);

    // generate a collection index list
    Xt_idxlist collectionlist = xt_idxlist_collection_new(idxlists, num_lists);
    for (size_t i = 0; i < num_lists; ++i)
      xt_idxlist_delete(idxlists[i]);

    // test generated collection list
    check_idxlist(collectionlist, index_list[0], (int)(num_lists*num_indices));
    int pos = -2, rc
      = xt_idxlist_get_position_of_index_off(collectionlist, 8, &pos, 5);
    if (rc != 1 || pos != -1)
      PUT_ERR("invalid return for index not to be found");

    rc = xt_idxlist_get_position_of_index_off(collectionlist, 1, &pos, num_indices*2);
    if (rc != 1 || pos != -1)
      PUT_ERR("invalid return for index not to be found");

    static const struct Xt_stripe ref_stripes[num_ref_stripes]
      = {{.start = 1, .stride = 1, .nstrides = 7},
         {.start = 7, .stride = -1, .nstrides = 7}};
    Xt_idxlist collectionlist_copy
      = idxlist_pack_unpack_copy(collectionlist);
    check_idxlist_copy(collectionlist, collectionlist_copy,
                       2 * (size_t)num_indices, index_list[0],
                       (size_t)num_ref_stripes, ref_stripes);
    xt_idxlist_delete(collectionlist_copy);
    xt_idxlist_delete(collectionlist);
  }

  { // test copying of collection index lists

    enum {
      num_lists = 2,
      num_indices = 7,
      num_ref_stripes = 2,
    };

    static const Xt_int index_list[num_lists][num_indices]
      = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1}};

    Xt_idxlist idxlists[num_lists];

    // create two index lists
    for (size_t i = 0; i < num_lists; ++i)
      idxlists[i] = xt_idxvec_new(index_list[i], num_indices);

    // generate a collection index list
    Xt_idxlist collectionlist = xt_idxlist_collection_new(idxlists, 2);

    for (size_t i = 0; i < num_lists; ++i)
      xt_idxlist_delete(idxlists[i]);

    // test generated collection list
    check_idxlist(collectionlist, &index_list[0][0], 2 * num_indices);

    // copy the index list
    Xt_idxlist collectionlist_copy = xt_idxlist_copy(collectionlist);

    // check copy of collection list
    check_idxlist(collectionlist_copy, &(index_list[0][0]),
                  2 * num_indices);

    // check the conversion to stripes

    struct Xt_stripe *stripes;
    int num_stripes;

    xt_idxlist_get_index_stripes(collectionlist, &stripes, &num_stripes);

    static const struct Xt_stripe ref_stripes[2]
      = {{.start = 1, .stride = 1, .nstrides = 7},
         {.start = 7, .stride = -1, .nstrides = 7}};

    check_stripes(stripes, num_stripes, ref_stripes, num_ref_stripes);

    // clean up
    free(stripes);

    xt_idxlist_delete(collectionlist);
    xt_idxlist_delete(collectionlist_copy);
  }

  { // test computing intersection of index list collections

    enum {
      num_indices = 7,
      num_vecs = 3
    };
    static const Xt_int index_list[num_vecs][num_indices]
      = {{1,2,3,4,5,6,7},{7,6,5,4,3,2,1},{2,6,1,4,7,3,0}},
      sorted_index_list[num_vecs * num_indices]
        = {0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,6,6,6,7,7,7};

    Xt_idxlist idxlists[num_vecs];

    // create two index lists
    for (size_t i = 0; i < num_vecs; ++i)
      idxlists[i] = xt_idxvec_new(index_list[i], num_indices);

    // generate a collection index list
    Xt_idxlist collectionlist = xt_idxlist_collection_new(idxlists, 3);

    for (size_t i = 0; i < num_vecs; ++i)
      xt_idxlist_delete(idxlists[i]);

    // test generated collection list
    check_idxlist(collectionlist, &(index_list[0][0]), (int)(3*num_indices));

    // test computing of intersection
    Xt_idxlist ref_idxvec = xt_idxvec_new(&(index_list[0][0]), 3*num_indices),
      intersection = xt_idxlist_get_intersection(ref_idxvec, collectionlist);

    check_idxlist(intersection, sorted_index_list, (int)(num_vecs*num_indices));

    xt_idxlist_delete(intersection);

    intersection = xt_idxlist_get_intersection(collectionlist, ref_idxvec);

    check_idxlist(intersection, sorted_index_list, (int)(num_vecs*num_indices));

    xt_idxlist_delete(intersection);
    xt_idxlist_delete(ref_idxvec);
    xt_idxlist_delete(collectionlist);
  }

  { // test using various tpyes of index lists

    enum {
      num_lists = 3,
      section_ndims = 2,
    };
    static const Xt_int index_list[] = {1,3,5,7,9,11};
    static const struct Xt_stripe stripes[]
      = {{.start = 0, .stride = 2, .nstrides = 5},
         {.start = 1, .stride = 2, .nstrides = 5}};
    static const Xt_int global_size[section_ndims] = {10,10};
    static const int local_size[section_ndims] = {5,5};
    static const Xt_int local_start[section_ndims] = {2,2};

    static const Xt_int ref_index_list[] = {1,3,5,7,9,11,
                                            0,2,4,6,8,1,3,5,7,9,
                                            22,23,24,25,26,
                                            32,33,34,35,36,
                                            42,43,44,45,46,
                                            52,53,54,55,56,
                                            62,63,64,65,66};

    Xt_idxlist idxlists[num_lists]
      = {xt_idxvec_new(index_list, sizeof(index_list)/sizeof(index_list[0])),
         xt_idxstripes_new(stripes, sizeof(stripes)/sizeof(stripes[0])),
         xt_idxsection_new(0, section_ndims, global_size, local_size, local_start)};

    // generate a collection index list
    Xt_idxlist collectionlist = xt_idxlist_collection_new(idxlists, num_lists);

    for (size_t i = 0; i < num_lists; ++i)
      xt_idxlist_delete(idxlists[i]);

    // test generated collection list
    check_idxlist(collectionlist, ref_index_list,
                  (int)(sizeof(ref_index_list)/sizeof(ref_index_list[0])));

    xt_idxlist_delete(collectionlist);
  }

  { // check get_bounding_box
    enum { num_lists = 2 };
    Xt_idxlist idxlists[num_lists] = {xt_idxempty_new(),
                                      xt_idxempty_new()};

    // generate a collection index list
    Xt_idxlist collectionlist = xt_idxlist_collection_new(idxlists, num_lists);

    for (size_t i = 0; i < num_lists; ++i)
      xt_idxlist_delete(idxlists[i]);

    enum { ndim = 3 };
    Xt_int global_size_bb[ndim];
    static const Xt_int global_start_index = 0;
    struct Xt_bounds bounds[ndim];

    for (size_t i = 0; i < ndim; ++i)
      global_size_bb[i] = 4;

    xt_idxlist_get_bounding_box(collectionlist, ndim, global_size_bb,
                                global_start_index, bounds);

    bool mismatch = false;
    for (size_t i = 0; i < ndim; ++i)
      mismatch |= (bounds[i].size != 0);

    if (mismatch)
      PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

    xt_idxlist_delete(collectionlist);
  }

  {
    enum {
      num_lists = 2,
      num_indices = 3,
    };
    // check get_bounding_box
    static const Xt_int indices[num_lists][num_indices]
      = {{45,35,32},{32,48,33}};

    Xt_idxlist idxlists[num_lists];

    // create two index lists
    for (size_t i = 0; i < num_lists; ++i)
      idxlists[i] = xt_idxvec_new(indices[i], num_indices);

    // generate a collection index list
    Xt_idxlist collectionlist = xt_idxlist_collection_new(idxlists, 2);

    for (size_t i = 0; i < num_lists; ++i)
      xt_idxlist_delete(idxlists[i]);

    enum { ndim = 3 };
    static const Xt_int global_size[ndim] = { 5, 4, 3 };
    static const Xt_int global_start_index = 1;
    struct Xt_bounds bounds[ndim];

    xt_idxlist_get_bounding_box(collectionlist, ndim, global_size,
                                global_start_index, bounds);

    static const Xt_int ref_start[3] = { 2, 2, 1 };

    bool mismatch = false;
    for (size_t i = 0; i < ndim; ++i)
      mismatch |= (bounds[i].size != 2 || bounds[i].start != ref_start[i]);

    if (mismatch)
      PUT_ERR("ERROR: xt_idxlist_get_bounding_box\n");

    test_coll_of_coll(collectionlist);

    xt_idxlist_delete(collectionlist);
  }

  test_idxsections();

  xt_finalize();
  xt_mpi_call(MPI_Finalize(), MPI_COMM_WORLD);

  return TEST_EXIT_CODE;
}

static void
test_idxsections(void)
{
  enum { irank = 4, num_regions=5, global_index_start=-1 };

  /* global size */
  static const Xt_int gsz[irank]
    = { 1, 1, 122, 132 };
  static const int gsz_i[irank]
    = { 1, 1, 122, 132 };
  /* slice size */
  static const int lsz[num_regions][irank]
    = { { 1, 1, 15,  8 }, // core region
        { 1, 1,  3, 14 }, // 1st halo
        { 1, 1,  3, 14 }, // 2nd halo
        { 1, 1, 15,  3 }, // 3rd halo
        { 1, 1, 15,  3 } }; // 4th halo

  /* slice start coordinate */
  static const Xt_int lst[num_regions][irank]
    = { { 0, 0, 46, 42 },
        { 0, 0, 43, 39 },
        { 0, 0, 61, 39 },
        { 0, 0, 46, 50 },
        { 0, 0, 46, 39 } };

  Xt_idxlist component_idxlists[num_regions];

  size_t num_indices = 0;
  for (size_t i = 0; i < num_regions; ++i) {
    size_t product = 1;
    for (size_t j = 0; j < irank; ++j)
      product *= (size_t)lsz[i][j];
    num_indices += product;
  }
  assert(num_indices < INT_MAX);
  Xt_int *restrict indices = xmalloc(sizeof (*indices) * num_indices);

  size_t idx_ofs = 0;

  /* sub-slices or regions */
  for (size_t i = 0; i < num_regions; ++i) {
    component_idxlists[i] = xt_idxsection_new(
      global_index_start, irank, gsz, lsz[i], lst[i]);

    Xt_int idx_stride4 = (Xt_int)(gsz[3] * gsz[2] * gsz[0]),
      idx_stride3 = (Xt_int)(gsz[3] * gsz[2]),
      idx_stride2 = gsz[3];
    for (size_t i4 = 0; i4 < (size_t)lsz[i][0]; ++i4) {
      Xt_int dim_ofs4 = (Xt_int)((lst[i][0]+(Xt_int)i4) * idx_stride4);
      for (size_t i3 = 0; i3 < (size_t)lsz[i][1]; ++i3) {
        Xt_int dim_ofs3
          = (Xt_int)((lst[i][1]+(Xt_int)i3) * idx_stride3 + dim_ofs4);
        for (size_t i2 = 0; i2 < (size_t)lsz[i][2]; ++i2) {
          Xt_int dim_ofs2
            = (Xt_int)((lst[i][2]+(Xt_int)i2) * idx_stride2 + dim_ofs3);
          for (size_t i1 = 0; i1 < (size_t)lsz[i][3]; ++i1)
            indices[idx_ofs+i1] = (Xt_int)
              (lst[i][3] + (Xt_int)i1 + global_index_start + dim_ofs2);
          idx_ofs += (size_t)lsz[i][3];
        }
      }
    }
  }

  Xt_idxlist collection_idxlist
    = xt_idxlist_collection_new(component_idxlists, num_regions);
  check_idxlist(collection_idxlist, indices, (int)num_indices);

  // global view
  Xt_idxlist global_section = xt_idxsection_new(
    global_index_start, irank, gsz, gsz_i, (Xt_int [irank]){ 0, 0, 0, 0 });

  Xt_idxlist intersection
    = xt_idxlist_get_intersection(global_section, collection_idxlist);
  xt_sort_xt_int(indices, num_indices);
  check_idxlist(intersection, indices, (int)num_indices);
  xt_idxlist_delete(intersection);
  xt_idxlist_delete(global_section);
  xt_idxlist_delete(collection_idxlist);
  for (size_t i = 0; i < num_regions; ++i)
    xt_idxlist_delete(component_idxlists[i]);
  free(indices);
}

static void
test_coll_of_coll(Xt_idxlist collectionlist)
{
  enum {
    num_vec = 2,
    num_lists = 3,
    num_indices = 3,
  };
  static const Xt_int indices[num_vec][num_indices]
    = { { 0, 2, 1 }, { 9, 8, 7 } };

  Xt_idxlist idxvecs[num_vec], idxlists[num_lists];

  // create another pair of index vectors
  for (size_t i = 0; i < num_vec; ++i)
    idxvecs[i] = xt_idxvec_new(indices[i], num_indices);

  // generate a collection index list
  idxlists[0] = xt_idxlist_collection_new(idxvecs, num_vec);
  idxlists[1] = xt_idxempty_new();
  idxlists[2] = collectionlist;

  // generate a collection index list containing another collection
  Xt_idxlist coll_of_coll = xt_idxlist_collection_new(idxlists, num_lists);

  for (size_t i = 0; i < num_vec; ++i) {
    Xt_idxlist intersection
      = xt_idxlist_get_intersection(coll_of_coll, idxvecs[i]);
    Xt_int ref_idx[num_indices];
    memcpy(ref_idx, indices[i], sizeof (indices[i]));
    xt_sort_xt_int(ref_idx, num_indices);
    check_idxlist(intersection, ref_idx, num_indices);
    xt_idxlist_delete(intersection);
  }

  for (size_t i = 0; i < num_vec; ++i)
    xt_idxlist_delete(idxvecs[i]);

  xt_idxlist_delete(idxlists[0]);
  xt_idxlist_delete(idxlists[1]);
  xt_idxlist_delete(coll_of_coll);
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
