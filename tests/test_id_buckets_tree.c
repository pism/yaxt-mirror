/**
 * @file test_id_buckets_tree.c
 *
 * @copyright Copyright  (C)  2019 Jörg Behrens <behrens@dkrz.de>
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

//#define VERBOSE

#include <mpi.h>
#include <string.h>
#include <stdint.h>
#include <yaxt.h>
#include "core/ppm_xfuncs.h"
#include "../src/xt_id_buckets.h"
#include "tests.h"

static void check_leaves_char(
  char * tree, size_t num_leaves, char * ref_ubounds) {

  // test leaves of tree
  Xt_id_buckets_iter iter = xt_id_buckets_iter_new(tree);
  char ubound = -1;
  for (size_t i = 0; i < num_leaves; ++i) {
    xt_id_buckets_iter_get_next_ubound(iter, &ubound);
    if (ubound != ref_ubounds[i]) PUT_ERR("wrong ubound");
  }
  xt_id_buckets_iter_delete(iter);
}

static void check_leaf_removing_char(
  char * tree, size_t tree_size, size_t num_leaves, char * ref_ubounds) {

  check_leaves_char(tree, num_leaves, ref_ubounds);

  if (num_leaves == 1) return;

  size_t * remove_idx = xmalloc((num_leaves - 1) * sizeof(remove_idx));
  uint64_t * keep_flags = xmalloc(num_leaves * sizeof(*keep_flags));
  size_t * dummy = xmalloc(num_leaves * sizeof(*dummy));
  char * temp_ref_ubounds = xmalloc(num_leaves);
  size_t temp_num_leaves;

  for (size_t num_leaves_to_remove = 1; num_leaves_to_remove < num_leaves;
       ++num_leaves_to_remove) {

    for (size_t i = 0; i < num_leaves_to_remove; ++i) remove_idx[i] = i;

    while(1) {

      for (uint64_t i = 0; i < num_leaves; ++i) keep_flags[i] = 1;
      for (size_t i = 0; i < num_leaves_to_remove; ++i)
        keep_flags[remove_idx[i]] = 0;
      for (size_t i = 0, j = 0; i < num_leaves; ++i)
        if (keep_flags[i]) temp_ref_ubounds[j++] = ref_ubounds[i];
      char * temp_tree = xmalloc(tree_size);
      memcpy(temp_tree, tree, tree_size);
      size_t temp_tree_size = tree_size;
      temp_num_leaves = num_leaves;

      xt_id_buckets_tree_remove_leaves(
        &temp_tree, &temp_tree_size, &temp_num_leaves, keep_flags, dummy);
      if (temp_num_leaves != num_leaves - num_leaves_to_remove)
        PUT_ERR("wrong number of leaves after xt_id_buckets_tree_remove_leaves");

      check_leaf_removing_char(
        temp_tree, temp_tree_size, temp_num_leaves, temp_ref_ubounds);

      free(temp_tree);

      {
        size_t i = num_leaves_to_remove;
        for (i = num_leaves_to_remove; i > 0; --i)
          if (remove_idx[i-1] < num_leaves + i-1 - num_leaves_to_remove) break;
        if (i == 0) break;
        remove_idx[i-1]++;
        for (; i < num_leaves_to_remove; ++i)
          remove_idx[i] = remove_idx[i-1] + 1;
      }
    }
  }
  free(temp_ref_ubounds);
  free(dummy);
  free(keep_flags);
  free(remove_idx);
}

int main(void) {

  { // check removing of tree leaves
    // tree:
    //         / \
    //        /   \
    //       /     \
    //      /       \
    //     /         \
    //    / \       / \
    //   /   \     /   \
    //  /\   /\   /\   /\
    // X  X X  X X  X X  X
    char tree[] =        {0B01010100,
                          0B01101011,
                          0B10001101,
                          0B10111010,
                          0B00110101,
                          0B11000110};
    char ref_ubounds[] = {0B00011111,
                          0B00111111,
                          0B01011111,
                          0B01111111,
                          0B10011111,
                          0B10111111,
                          0B11011111,
                          0B11111111};
    check_leaf_removing_char(tree, sizeof(tree), 8, ref_ubounds);
  }

  { // check splitting of tree leaves
    // tree:
    //  /\
    // X  X
    char tree[] = {0B01000110};

    check_leaves_char(tree, 2, (char[]){0B01111111, 0B11111111});

    size_t tree_size = sizeof(tree);
    size_t id_size = sizeof(char);

    char * tree_copy = xmalloc(tree_size);
    memcpy(tree_copy, tree, tree_size);
    size_t num_leaves = 2;
    size_t split_flag[] = {1,1};

    xt_id_buckets_tree_split_leaves(
      &tree_copy, &tree_size, &num_leaves, id_size, split_flag);

    if (num_leaves != 4) PUT_ERR("wrong number of leaves");

    char ref_ubounds[] = {0B00111111,
                          0B01111111,
                          0B10111111,
                          0B11111111};
    check_leaves_char(tree_copy, num_leaves, ref_ubounds);
  }

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
