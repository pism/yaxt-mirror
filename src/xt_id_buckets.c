/**
 * @file xt_id_buckets.c
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

#include <string.h>
#include <mpi.h>
#include <stdint.h>

#include "core/ppm_xfuncs.h"
#include "core/core.h"
#include "xt/xt_core.h"
#include "xt/xt_mpi.h"
#include "ensure_array_size.h"
#include "xt_id_buckets.h"
#include "xt_memcmp.h"

struct Xt_id_buckets_ {
  size_t tree_size;
  size_t id_size;
  size_t num_buckets;
  char tree[];
};

// Each leaf in the tree is described by a number of bits. These bits define
// how the current leaf is to be adjusted to become the next leaf. The bits for
// the leaves are stored continuously in the tree. Leaves have id_size bytes. In
// their initial state all bytes are set to one.
// The tree level is a pointer to the i'the bit in a leaf. Its initial state is
// 0. The tree bit index is a pointer to the j'th bit in the tree. To
// turn the current leaf into the next one, you have to perform the following
// actions.
//
// do {leaf_bits[--i] = 1;} while(tree_bits[j++]);
// while(tree_bits[j++]) leaf_bits[i++] = tree_bits[j++];

struct Xt_id_buckets_iter_ {

  size_t curr_tree_level,
         curr_tree_bit_idx;
  char const * tree;
};

static void Xt_id_buckets_iter_init(
  Xt_id_buckets_iter iter, char const * tree) {

  iter->curr_tree_level = 1;
  iter->curr_tree_bit_idx = 0;
  iter->tree = tree;
}

Xt_id_buckets_iter xt_id_buckets_iter_new(char const * tree) {
  Xt_id_buckets_iter iter = xmalloc(1 * sizeof(*iter));

  Xt_id_buckets_iter_init(iter, tree);

  return iter;
}

void xt_id_buckets_iter_delete(Xt_id_buckets_iter iter) {

  free(iter);
}

static char get_bit(char const * data, size_t bit_idx) {

  size_t byte_idx = bit_idx / 8;
  char mask = (char)(((unsigned char)128) >> (bit_idx % 8));
  return data[byte_idx] & mask;
}

static void unset_bit(char * data, size_t bit_idx) {

  size_t byte_idx = bit_idx / 8;
  char mask = (char)(~(((unsigned char)128) >> (bit_idx % 8)));
  data[byte_idx] &= mask;
}

static void set_bit(char * data, size_t bit_idx) {

  size_t byte_idx = bit_idx / 8;
  char mask = (char)(((unsigned char)128) >> (bit_idx % 8));
  data[byte_idx] |= mask;
}

void xt_id_buckets_iter_get_next_ubound(
  Xt_id_buckets_iter iter, char * ubound) {

  size_t curr_tree_level = iter->curr_tree_level;
  size_t curr_tree_bit_idx = iter->curr_tree_bit_idx;
  char const * tree = iter->tree;

  // go up the tree
  do {
    set_bit(ubound, --curr_tree_level);
  } while (get_bit(tree, curr_tree_bit_idx++));

  // go down until we reach a next leaf
  while(get_bit(tree, curr_tree_bit_idx++)) {

    if (get_bit(tree, curr_tree_bit_idx++)) set_bit(ubound, curr_tree_level++);
    else                                  unset_bit(ubound, curr_tree_level++);
  }

  iter->curr_tree_level = curr_tree_level;
  iter->curr_tree_bit_idx = curr_tree_bit_idx;
}

static void get_basic_tree(
  char ** tree, size_t * tree_size, size_t * num_leaves) {

  // 1st bit == 0 => go up once
  // 2nd bit == 1 => go down once
  // 3rd bit == 0 => set lowest bit to zero
  // 4th bit == 0 => mark end of first leaf
  // 5th bit == 0 => go up once
  // 6th bit == 1 => go down once
  // 7th bit == 1 => set lowest bit to one
  // 8th bit == 0 => mark end of second leaf

  ENSURE_ARRAY_SIZE(*tree, *tree_size, 1);
  (*tree)[0] = 0B01000110;
  *num_leaves = 2;
}

void xt_id_buckets_tree_remove_leaves(
  char ** tree_, size_t * tree_size_, size_t * num_leaves_,
  uint64_t * keep_flags, size_t * bucket_ids) {

  char * tree = *tree_;
  size_t tree_size = *tree_size_;
  size_t num_leaves = *num_leaves_;

  int nothing_to_do = 1;

  for (size_t i = 0; i < num_leaves; ++i) {
    if (keep_flags[i] == 0) {
      nothing_to_do = 0;
      break;
    }
  }

  if (nothing_to_do) return;

  char * new_tree = xmalloc(tree_size);

  size_t prev_tree_bit_idx = 0;
  size_t curr_tree_bit_idx = 0;
  size_t new_tree_bit_idx = 0;
  size_t new_num_leaves = 0;

  for (size_t i = 0; i < num_leaves; ++i) {

    prev_tree_bit_idx = curr_tree_bit_idx;

    if (keep_flags[i]) {

      // go up until we reach a branch
      while (get_bit(tree, curr_tree_bit_idx++)) {
        set_bit(new_tree, new_tree_bit_idx++);
      }
      unset_bit(new_tree, new_tree_bit_idx++);

      // go down until we reach a next leaf
      while (get_bit(tree, curr_tree_bit_idx++)) {
        set_bit(new_tree, new_tree_bit_idx++);
        if (get_bit(tree, curr_tree_bit_idx++))
          set_bit(new_tree, new_tree_bit_idx++);
        else
          unset_bit(new_tree, new_tree_bit_idx++);
      };
      unset_bit(new_tree, new_tree_bit_idx++);

      keep_flags[new_num_leaves] = keep_flags[i];
      bucket_ids[new_num_leaves] = bucket_ids[i];
      new_num_leaves++;

    } else {

      // if this is the last leaf, just skip it
      if (i == num_leaves - 1) continue;

      // go up until we reach a branch
      while (get_bit(tree, curr_tree_bit_idx++));

      size_t temp_tree_bit_idx = curr_tree_bit_idx;

      // go down until we reach a next leaf
      while (get_bit(tree, curr_tree_bit_idx++)) ++curr_tree_bit_idx;

      size_t skip_start_bit_idx = curr_tree_bit_idx - 3;

      // down moves of current leaf and up moves of next leaf cancel each
      // other out, determine how many bits we need to skip
      while ((skip_start_bit_idx > temp_tree_bit_idx) &&
             get_bit(tree, curr_tree_bit_idx++)) {
        skip_start_bit_idx -= 2;
      };

      // if all down moves have been canceled out
      if (skip_start_bit_idx == temp_tree_bit_idx)
        skip_start_bit_idx--; // skip the end mark for the up moves
                              // of the current leaf

      while (skip_start_bit_idx > prev_tree_bit_idx) {
        if (get_bit(tree, --skip_start_bit_idx))
          set_bit(tree, --curr_tree_bit_idx);
        else
          unset_bit(tree, --curr_tree_bit_idx);
      }
    }
  }

  free(tree);
  size_t new_tree_size = (new_tree_bit_idx + 7) / 8;
  *tree_ = xrealloc(new_tree, new_tree_size);
  *tree_size_ = new_tree_size;
  *num_leaves_ = new_num_leaves;
}

int xt_id_buckets_tree_split_leaves(
  char ** tree_, size_t * tree_size_, size_t * num_leaves_, size_t id_size,
  size_t * split_flags) {

  char * tree = *tree_;
  size_t tree_size = *tree_size_;
  size_t num_leaves = *num_leaves_;

  // count number of potential splits
  size_t num_splits = 0;
  for (size_t i = 0; i < num_leaves; ++i) if (split_flags[i]) ++num_splits;

  if (num_splits == 0) return 0;

  size_t max_tree_level = 8 * id_size;

  // now we generate a new tree from the old one and try to incorporate
  // appropriate bucket splits

  // 7 additional bits per potential split
  char * new_tree = xmalloc(tree_size + (num_splits * 7 + 7) / 8);

  size_t curr_new_tree_bit_idx = 0;
  size_t new_num_leaves = num_leaves;

  size_t curr_tree_level = 1;
  size_t prev_tree_bit_idx = 0;
  size_t curr_tree_bit_idx = 0;

  int split_flag = 0;

  for (size_t i = 0; i < num_leaves; ++i) {

    // parse tree until current leaf and copy respective bits
    {
      // go up until we reach a branch
      --curr_tree_level;
      while (get_bit(tree, curr_tree_bit_idx++)) {
        --curr_tree_level;
        set_bit(new_tree, curr_new_tree_bit_idx++);
      }
      unset_bit(new_tree, curr_new_tree_bit_idx++);

      // go down until we reach a next leaf
      while (get_bit(tree, curr_tree_bit_idx++)) {
        set_bit(new_tree, curr_new_tree_bit_idx++);
        if (get_bit(tree, curr_tree_bit_idx++))
          set_bit(new_tree, curr_new_tree_bit_idx++);
        else
          unset_bit(new_tree, curr_new_tree_bit_idx++);
        ++curr_tree_level;
      };
    }

    prev_tree_bit_idx = curr_tree_bit_idx;

    // if we need to split the current bucket and it is possible
    if (split_flags[i] && (curr_tree_level < max_tree_level)) {

      // generate split

      // new left leaf
      set_bit(new_tree, curr_new_tree_bit_idx++);   // mark change to leaf
      unset_bit(new_tree, curr_new_tree_bit_idx++); // left
      unset_bit(new_tree, curr_new_tree_bit_idx++); // mark end of leaf

      // go up only once and reset original bit
      unset_bit(new_tree, curr_new_tree_bit_idx++);

      // go down once and generate new right leaf
      set_bit(new_tree, curr_new_tree_bit_idx++); // mark change to leaf
      set_bit(new_tree, curr_new_tree_bit_idx++); // right
      unset_bit(new_tree, curr_new_tree_bit_idx++); // make end of leaf

      // go up once to reach original level
      set_bit(new_tree, curr_new_tree_bit_idx++);

      ++new_num_leaves;
    } else {
      // no split -> mark end of original leaf
      unset_bit(new_tree, curr_new_tree_bit_idx++);
    }
  }

  if (num_leaves != new_num_leaves) {
    free(tree);
    *num_leaves_ = new_num_leaves;
    size_t new_tree_size = (curr_new_tree_bit_idx + 7) / 8;
    *tree_size_ = new_tree_size ;
    *tree_ = xrealloc(new_tree, new_tree_size );
  } else {
    free(new_tree);
  }

  return new_num_leaves != num_leaves;
}

static void distribute_buckets_recursive(
  uint64_t * bucket_sizes, size_t num_buckets, size_t bucket_sum_offset,
  size_t num_ids, size_t start_idx, size_t end_idx, size_t num_req_buckets,
  size_t * bucket_ids) {

  if ((start_idx >= end_idx) || (num_buckets <= 1)) {
    for (size_t i = 0; i < num_buckets; ++i) bucket_ids[i] = start_idx;
    return;
  }

  size_t middle_idx = (start_idx + end_idx) / 2;

  size_t bucket_split_point =
    (size_t)(((long long)num_ids * (long long)(middle_idx + 1)) /
             (long long)num_req_buckets);

  size_t bucket_split_idx = 0;
  size_t bucket_sum = bucket_sum_offset;
  size_t curr_bucket_size = (size_t)*bucket_sizes;

  while ((bucket_split_idx < num_buckets) &&
         (bucket_sum + curr_bucket_size < bucket_split_point)) {

    bucket_sum += curr_bucket_size;
    curr_bucket_size = (size_t)bucket_sizes[++bucket_split_idx];
  }

  // determine which rank group gets the current bucket
  if ((bucket_split_point - bucket_sum) > (curr_bucket_size / 2)) {
    bucket_sum += curr_bucket_size;
    ++bucket_split_idx;
  }

  distribute_buckets_recursive(
    bucket_sizes, bucket_split_idx, bucket_sum_offset, num_ids,
    start_idx, middle_idx, num_req_buckets, bucket_ids);
  distribute_buckets_recursive(
    bucket_sizes + bucket_split_idx, num_buckets - bucket_split_idx,
    bucket_sum, num_ids, middle_idx + 1, end_idx, num_req_buckets,
    bucket_ids + bucket_split_idx);
}

static void distribute_buckets(
  uint64_t * bucket_sizes, size_t num_buckets, size_t num_ids,
  size_t num_req_buckets, size_t * bucket_ids) {

  distribute_buckets_recursive(
    bucket_sizes, num_buckets, 0, num_ids, 0, num_req_buckets-1,
    num_req_buckets, bucket_ids);
}

static int split_buckets(uint64_t * bucket_sizes, size_t * bucket_ids,
                         size_t num_ids, size_t num_req_buckets,
                         size_t id_size, char ** tree, size_t * tree_size,
                         size_t * num_leaves) {

  size_t bucket_sum = 0;
  size_t prev_bucket_split_point = 0;
  size_t bucket_idx = 0;
  size_t prev_bucket_idx = 0;

  // merge zero size buckets
  xt_id_buckets_tree_remove_leaves(
    tree, tree_size, num_leaves, bucket_sizes, bucket_ids);

  // determine which buckets need to be split
  for (size_t idx = 0; idx < num_req_buckets; ++idx) {

    size_t bucket_split_point =
      (size_t)(((long long)num_ids * (long long)(idx + 1)) /
               (long long)num_req_buckets);

    size_t curr_ids_per_rank = 0;

    while ((bucket_idx < *num_leaves) && (bucket_ids[bucket_idx] == idx)) {
      curr_ids_per_rank += (size_t)bucket_sizes[bucket_idx];
      bucket_ids[bucket_idx] = 0;
      bucket_idx++;
    }

    // if there are currently too many ids
    if ((((bucket_split_point - prev_bucket_split_point) * 5) / 4) <
        curr_ids_per_rank) {

      // if we should split the first bucket
      if (bucket_sum < prev_bucket_split_point)
        bucket_ids[prev_bucket_idx] = 1;

      // if we should split the last bucket
      if (bucket_sum + curr_ids_per_rank > bucket_split_point)
        bucket_ids[bucket_idx-1] = 1;
    }

    bucket_sum += curr_ids_per_rank;
    prev_bucket_split_point = bucket_split_point;
    prev_bucket_idx = bucket_idx;
  }

  // at this point the bucket_ids have been overwritten by 1 if the respective
  // bucket is to be split or 0 otherwise
  return
    xt_id_buckets_tree_split_leaves(
      tree, tree_size, num_leaves, id_size, bucket_ids);
}

static int check_tree(uint64_t * bucket_sizes, size_t id_size,
                      size_t num_req_buckets, char ** tree, size_t * tree_size,
                      size_t * num_leaves_) {

  size_t num_leaves = *num_leaves_;
  uint64_t num_ids = 0;
  for (size_t i = 0; i < num_leaves; ++i) num_ids += bucket_sizes[i];

  size_t * bucket_ids = xmalloc(num_leaves * sizeof(*bucket_ids));

  // distribute buckets into num_req_buckets equal parts
  distribute_buckets(
    bucket_sizes, num_leaves, (size_t)num_ids, num_req_buckets, bucket_ids);

  // check which final bucket get to much data and split current buckets
  // accordingly
  int split_flag =
    split_buckets(bucket_sizes, bucket_ids, (size_t)num_ids, num_req_buckets,
                  id_size, tree, tree_size, num_leaves_);

  free(bucket_ids);

  return split_flag;
}

static void finalise_tree(
  char ** tree, size_t * tree_size, size_t num_leaves,
  uint64_t * bucket_sizes, size_t num_req_buckets, MPI_Comm comm) {

  size_t num_ids = 0;
  for (size_t i = 0; i < num_leaves; ++i) num_ids += bucket_sizes[i];

  size_t * bucket_ids = xmalloc(num_leaves * sizeof(*bucket_ids));

  // generate finale bucket indices for all current buckets
  distribute_buckets(
    bucket_sizes, num_leaves, num_ids, num_req_buckets, bucket_ids);

  if (num_leaves > 0) {
    // merge buckets that have the same final bucket id
    uint64_t * keep_flags = bucket_sizes;
    keep_flags[num_leaves - 1] = 1;
    for (size_t i = 0; i < num_leaves - 1; ++i)
      // zero if identical
      keep_flags[i] = (uint64_t)(bucket_ids[i+1] - bucket_ids[i]);
    xt_id_buckets_tree_remove_leaves(
      tree, tree_size, &num_leaves, keep_flags, bucket_ids);
  }

  // if we have too many leave
  if (num_leaves > num_req_buckets)
    Xt_abort(comm, "INTERNAL_ERROR(finalise_tree) more buckets than ranks",
             __FILE__, __LINE__);
  // if we have too less leaves
  else if (num_leaves < num_req_buckets) {

    // 4 additional bits per empty leaf
    char * new_tree =
      xmalloc(*tree_size + (num_req_buckets - num_leaves) * 4);
    char * old_tree = *tree;
    size_t new_tree_bit_idx = 0;
    size_t curr_tree_bit_idx = 0;

    // ensure that the first rank gets the first bucket
    if (num_leaves > 0) bucket_ids[0] = 0;

    // add some empty leaves
    for (size_t i = 0, j = 0, last_down_bit = 1; i < num_req_buckets; ++i) {

      // if there is a bucket for the current rank
      if ((j < num_leaves) && (bucket_ids[j] == i)) {

        // copy leaf

        // go up until we reach a branch
        while (get_bit(old_tree, curr_tree_bit_idx++)) {
          set_bit(new_tree, new_tree_bit_idx++);
        }
        unset_bit(new_tree, new_tree_bit_idx++);

        // go down until we reach a next leaf
        while (get_bit(old_tree, curr_tree_bit_idx++)) {
          set_bit(new_tree, new_tree_bit_idx++);
          if (get_bit(old_tree, curr_tree_bit_idx++)) {
            set_bit(new_tree, new_tree_bit_idx++);
            last_down_bit = 1;
          } else {
            unset_bit(new_tree, new_tree_bit_idx++);
            last_down_bit = 0;
          }
        };
        unset_bit(new_tree, new_tree_bit_idx++);

        ++j;

      } else {

        // add empty leaf

        // go up only once
        unset_bit(new_tree, new_tree_bit_idx++);

        // go down once
        set_bit(new_tree, new_tree_bit_idx++);
        if (last_down_bit) set_bit(new_tree, new_tree_bit_idx++);
        else unset_bit(new_tree, new_tree_bit_idx++);
        unset_bit(new_tree, new_tree_bit_idx++);
      }
    }

    free(old_tree);
    size_t new_tree_size = (new_tree_bit_idx + 7) / 8;
    *tree = xrealloc(new_tree, new_tree_size);
    *tree_size = new_tree_size;
  }

  free(bucket_ids);
}

void xt_generate_id_buckets(
  void const * ids, size_t id_size, size_t num_ids, size_t num_req_buckets,
  MPI_Comm comm, void * ubounds) {

  if (id_size < 1)
    Xt_abort(comm, "ERROR(xt_generate_id_buckets): invalid id_size "
             "(id_size < 1)", __FILE__, __LINE__);

  int is_inter;
  xt_mpi_call(MPI_Comm_test_inter(comm, &is_inter), comm);
  if (id_size < 1)
    Xt_abort(comm, "ERROR(xt_generate_id_buckets): intercomm is not supported",
             __FILE__, __LINE__);

  unsigned char * ubound = xmalloc(id_size);

  uint64_t * bucket_counts = NULL;
  size_t bucket_counts_size = 0;

  char * tree = NULL;
  size_t tree_size = 0;
  size_t num_leaves;
  get_basic_tree(&tree, &tree_size, &num_leaves);

  uint64_t * global_bucket_counts;

  Xt_memcmp_func compare_ids = xt_memcmp_get_func(id_size);

  // sort the ids
  void const * sorted_ids = ids;
  for (size_t i = 1; i < num_ids; ++i) {

    if (compare_ids((char*)ids + (i - 1) * id_size,
                    (char*)ids + i * id_size) > 0) {

      void * temp_ids = xmalloc(id_size * num_ids);
      memcpy(temp_ids, ids, id_size * num_ids);
      qsort(temp_ids, num_ids, id_size, compare_ids);
      sorted_ids = temp_ids;
      break;
    }
  }

  do {

    ENSURE_ARRAY_SIZE(bucket_counts, bucket_counts_size, 2 * num_leaves);
    uint64_t * local_bucket_counts = bucket_counts;
    global_bucket_counts = bucket_counts + num_leaves;
    uint64_t remain_count = (uint64_t)num_ids;

    struct Xt_id_buckets_iter_ iter;
    Xt_id_buckets_iter_init(&iter, tree);
    memset(ubound, 255, id_size);

    char const * curr_ids = sorted_ids;

    for (size_t i = 0; i < num_leaves; ++i) {

      xt_id_buckets_iter_get_next_ubound(&iter, ubound);

      uint64_t curr_bucket_count = 0;

      // count current bucket
      while ((curr_bucket_count < remain_count) &&
             (compare_ids(curr_ids + curr_bucket_count * id_size, ubound) <= 0))
        ++curr_bucket_count;

      local_bucket_counts[i] = curr_bucket_count;
      remain_count -= curr_bucket_count;
      curr_ids += curr_bucket_count * id_size;
    }

    // allgather bucket counts
    xt_mpi_call(
      MPI_Allreduce(
        local_bucket_counts, global_bucket_counts, (int)num_leaves,
        MPI_UINT64_T, MPI_SUM, comm), comm);

    // determine splits and adjust tree if necessary
  } while (check_tree(global_bucket_counts, id_size, num_req_buckets,
                      &tree, &tree_size, &num_leaves));

  if (sorted_ids != ids) free((void*)sorted_ids);

  // ensure that the number of leaves in the tree matches the number of
  // processes
  finalise_tree(&tree, &tree_size, num_leaves, global_bucket_counts,
                num_req_buckets, comm);

  free(bucket_counts);

  // generate upper bounds
  Xt_id_buckets_iter iter = xt_id_buckets_iter_new(tree);
  memset(ubound, -1, id_size);
  for (size_t i = 0; i < num_req_buckets; ++i) {
    xt_id_buckets_iter_get_next_ubound(iter, ubound);
    memcpy(((char*)ubounds) + i * id_size, ubound, id_size);
  }
  free(ubound);
  xt_id_buckets_iter_delete(iter);
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
