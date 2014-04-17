/**
 * @file xt_xmap_intersection.c
 *
 * @copyright Copyright  (C)  2012 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <mpi.h>

#include "xt/xt_idxlist.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_xmap.h"
#include "xt_xmap_internal.h"
#include "xt/xt_mpi.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_xmap_intersection.h"
#include "ensure_array_size.h"

static MPI_Comm     xmap_intersection_get_communicator(Xt_xmap xmap);
static int          xmap_intersection_get_num_destinations(Xt_xmap xmap);
static int          xmap_intersection_get_num_sources(Xt_xmap xmap);
static void
xmap_intersection_get_destination_ranks(Xt_xmap xmap, int * ranks);
static void
xmap_intersection_get_source_ranks(Xt_xmap xmap, int * ranks);
static Xt_xmap_iter xmap_intersection_get_in_iterator(Xt_xmap xmap);
static Xt_xmap_iter xmap_intersection_get_out_iterator(Xt_xmap xmap);
static void         xmap_intersection_delete(Xt_xmap xmap);
static int          xmap_iterator_intersection_next(Xt_xmap_iter iter);
static int          xmap_intersection_iterator_get_rank(Xt_xmap_iter iter);
static int const *
xmap_intersection_iterator_get_transfer_pos(Xt_xmap_iter iter);
static int
xmap_intersection_iterator_get_num_transfer_pos(Xt_xmap_iter iter);
static void         xmap_intersection_iterator_delete(Xt_xmap_iter iter);
static int          xmap_intersection_get_max_src_pos(Xt_xmap xmap);
static int          xmap_intersection_get_max_dst_pos(Xt_xmap xmap);


static const struct Xt_xmap_iter_vtable
xmap_iterator_intersection_vtable = {
  .next                 = xmap_iterator_intersection_next,
  .get_rank             = xmap_intersection_iterator_get_rank,
  .get_transfer_pos     = xmap_intersection_iterator_get_transfer_pos,
  .get_num_transfer_pos = xmap_intersection_iterator_get_num_transfer_pos,
  .delete               = xmap_intersection_iterator_delete};

typedef struct Xt_xmap_iter_intersection_ *Xt_xmap_iter_intersection;

struct Xt_xmap_iter_intersection_ {

  const struct Xt_xmap_iter_vtable * vtable;

  struct exchange_data * msg;
  int msgs_left;
};

static inline Xt_xmap_iter_intersection
xmii(void *iter)
{
  return (Xt_xmap_iter_intersection)iter;
}


static const struct Xt_xmap_vtable xmap_intersection_vtable = {
        .get_communicator      = xmap_intersection_get_communicator,
        .get_num_destinations  = xmap_intersection_get_num_destinations,
        .get_num_sources       = xmap_intersection_get_num_sources,
        .get_destination_ranks = xmap_intersection_get_destination_ranks,
        .get_source_ranks      = xmap_intersection_get_source_ranks,
        .get_out_iterator      = xmap_intersection_get_out_iterator,
        .get_in_iterator       = xmap_intersection_get_in_iterator,
        .delete                = xmap_intersection_delete,
        .get_max_src_pos       = xmap_intersection_get_max_src_pos,
        .get_max_dst_pos       = xmap_intersection_get_max_dst_pos};

struct exchange_data {
  // list of relative positions in memory to send or receive
  int * transfer_pos;
  int num_transfer_pos;
  int rank;
};

struct Xt_xmap_intersection_ {

  const struct Xt_xmap_vtable * vtable;

  struct exchange_data *in_msg, *out_msg;
  int n_in, n_out;

  // we need the max position in order to enable quick range-checks
  // for xmap-users like redist
  int max_src_pos; // max possible pos over all src transfer_pos (always >= 0)
  int max_dst_pos; // same for dst

  MPI_Comm comm;
};

typedef struct Xt_xmap_intersection_ *Xt_xmap_intersection;

static inline Xt_xmap_intersection
xmi(void *xmap)
{
  return (Xt_xmap_intersection)xmap;
}

static MPI_Comm xmap_intersection_get_communicator(Xt_xmap xmap) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  return xmap_intersection->comm;
}

static int xmap_intersection_get_num_destinations(Xt_xmap xmap) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  // the number of destination equals the number of source messages
  return xmap_intersection->n_out;
}

static int xmap_intersection_get_num_sources(Xt_xmap xmap) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  // the number of source equals the number of destination messages
  return xmap_intersection->n_in;
}

static void xmap_intersection_get_destination_ranks(Xt_xmap xmap, int * ranks) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  for (int i = 0; i < xmap_intersection->n_out; ++i)
    ranks[i] = xmap_intersection->out_msg[i].rank;
}

static void xmap_intersection_get_source_ranks(Xt_xmap xmap, int * ranks) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  for (int i = 0; i < xmap_intersection->n_in; ++i)
    ranks[i] = xmap_intersection->in_msg[i].rank;
}

/* compute list positions for recv direction */
static void
generate_dir_transfer_pos_dst(int num_intersections,
                              const struct Xt_com_list
                              intersections[num_intersections],
                              Xt_idxlist mypart_idxlist,
                              int *resCount,
                              struct exchange_data **resSets,
                              Xt_int ** indices_to_remove,
                              int *num_indices_to_remove_per_intersection)
{
  struct exchange_data *restrict resSets_
    = *resSets = xmalloc((size_t)num_intersections * sizeof(**resSets));

  int *restrict position_mask
    = xcalloc((size_t)xt_idxlist_get_num_indices(mypart_idxlist),
              sizeof(*position_mask));

  int new_num_intersections = 0;
  size_t total_num_indices_to_remove = 0;
  size_t curr_indices_to_remove_size = 0;
  Xt_int *restrict indices_to_remove_ = *indices_to_remove;
  int *restrict intersection_pos = NULL;

  for (int i = 0; i < num_intersections; ++i) {

    const Xt_int *restrict intersection_idxvec
      = xt_idxlist_get_indices_const(intersections[i].list);
    int max_intersection_size
      = xt_idxlist_get_num_indices(intersections[i].list);
    intersection_pos
      = xrealloc(intersection_pos,
                 (size_t)max_intersection_size * sizeof(*intersection_pos));

    int retval = xt_idxlist_get_positions_of_indices(
      mypart_idxlist, intersection_idxvec, max_intersection_size,
      intersection_pos, 1);
    assert(retval != 1);

    // we have to enforce single_match_only not only within a single
    // intersection, but also between all intersections

    int intersection_size = 0;
    int num_indices_to_remove_isect = 0;

    /* at most intersection_size many indices need to be removed */
    ENSURE_ARRAY_SIZE(indices_to_remove_, curr_indices_to_remove_size,
                      total_num_indices_to_remove
                      + (size_t)max_intersection_size);

    for (int j = 0; j < max_intersection_size; ++j) {

      int pos = intersection_pos[j];
      /* the predicate effectively conditionalizes adding of either
       * the position to intersection_pos
       *   if the current value was NOT already in another intersection
       * or
       * the index to indices_to_remove_
       *   if the current value was already in another intersection
       */
      int predicate = position_mask[pos];
      intersection_pos[intersection_size] = pos;
      indices_to_remove_[total_num_indices_to_remove
                         + (size_t)num_indices_to_remove_isect]
        = intersection_idxvec[j];
      intersection_size += predicate ^ 1;
      num_indices_to_remove_isect += predicate;
      position_mask[pos] = 1;
    }

    total_num_indices_to_remove += (size_t)num_indices_to_remove_isect;
    num_indices_to_remove_per_intersection[i] = num_indices_to_remove_isect;

    if (intersection_size > 0) {
      resSets_[new_num_intersections].transfer_pos = intersection_pos;
      resSets_[new_num_intersections].num_transfer_pos = intersection_size;
      resSets_[new_num_intersections].rank = intersections[i].rank;
      new_num_intersections++;
      intersection_pos = NULL;
    }
  }
  free(intersection_pos);

  *indices_to_remove = xrealloc(indices_to_remove_,
                                (size_t)total_num_indices_to_remove
                                * sizeof (**indices_to_remove));
  *resCount = new_num_intersections;
  if (num_intersections != new_num_intersections)
    *resSets = xrealloc(resSets_,
                        (size_t)new_num_intersections * sizeof(**resSets));

  free(position_mask);
}

/* compute list positions for send direction */
static void
generate_dir_transfer_pos_src(int num_intersections,
                              const struct Xt_com_list
                              intersections[num_intersections],
                              Xt_idxlist mypart_idxlist,
                              int *resCount,
                              struct exchange_data **resSets,
                              Xt_int * indices_to_remove,
                              int * num_indices_to_remove_per_intersection)
{

  struct exchange_data *restrict resSets_ = *resSets
    = xmalloc((size_t)num_intersections * sizeof(**resSets));

  int new_num_intersections = 0;
  int offset = 0;

  Xt_int * new_intersection_idxvec = NULL;
  size_t curr_new_intersection_idxvec_size = 0;
  int *restrict intersection_pos = NULL;

  for (int i = 0; i < num_intersections; ++i) {

    const Xt_int *restrict intersection_idxvec
      = xt_idxlist_get_indices_const(intersections[i].list);
    int intersection_size
      = xt_idxlist_get_num_indices(intersections[i].list);
    intersection_pos = xrealloc(intersection_pos,
                                (size_t)intersection_size
                                * sizeof(*intersection_pos));

    int num_indices_to_remove = num_indices_to_remove_per_intersection[i];
    if (num_indices_to_remove > 0) {

      ENSURE_ARRAY_SIZE(
        new_intersection_idxvec, curr_new_intersection_idxvec_size,
        intersection_size - num_indices_to_remove + 1);
      int new_intersection_size = 0;

      for (int j = 0; j < intersection_size; ++j) {

        int discard = 0;

        Xt_int idx = intersection_idxvec[j];
        /* could be improved with BLOOM-filter if
         * num_indices_to_remove was sufficiently large */
        for (int k = 0; k < num_indices_to_remove; ++k)
          discard |= (idx == indices_to_remove[offset + k]);

        new_intersection_idxvec[new_intersection_size] = idx;
        new_intersection_size += !discard;
      }

      intersection_idxvec = new_intersection_idxvec;
      intersection_size = new_intersection_size;
      offset = offset + num_indices_to_remove;
    }

    int retval;
    retval = xt_idxlist_get_positions_of_indices(
      mypart_idxlist, intersection_idxvec, intersection_size,
      intersection_pos, 0);
    assert(retval != 1);

    if (intersection_size > 0) {
      resSets_[new_num_intersections].transfer_pos = intersection_pos;
      resSets_[new_num_intersections].num_transfer_pos = intersection_size;
      resSets_[new_num_intersections].rank = intersections[i].rank;
      new_num_intersections++;
      intersection_pos = NULL;
    }
  }

  free(new_intersection_idxvec);

  *resCount = new_num_intersections;
  if (num_intersections != new_num_intersections)
    *resSets = xrealloc(resSets_,
                        (size_t)new_num_intersections * sizeof(**resSets));
}

static void
exchange_points_to_remove(int num_src_intersections,
                          const struct Xt_com_list
                          src_com[num_src_intersections],
                          int num_dst_intersections,
                          const struct Xt_com_list
                          dst_com[num_dst_intersections],
                          Xt_int ** src_indices_to_remove,
                          int *restrict num_src_indices_to_remove_per_intersection,
                          Xt_int * dst_indices_to_remove,
                          int *restrict num_dst_indices_to_remove_per_intersection,
                          MPI_Comm comm) {

  enum {
    HEADER_EXCHANGE_TAG = 0,
    DATA_EXCHANGE_TAG = 1,
  };

  MPI_Request * requests
    = xmalloc((size_t)(num_src_intersections + 2 * num_dst_intersections) *
              sizeof(*requests));
  MPI_Request *restrict recv_requests = requests,
    *restrict send_header_requests = requests + num_src_intersections,
    *restrict send_data_requests = send_header_requests + num_dst_intersections;

  // set up receives for indices that need to be removed from the send messages
  for (int i = 0; i < num_src_intersections; ++i)
    xt_mpi_call(MPI_Irecv(num_src_indices_to_remove_per_intersection + i,
                          1, MPI_INT, src_com[i].rank, HEADER_EXCHANGE_TAG,
                          comm, recv_requests+i), comm);

  // send indices that need to be removed on the target side due to duplicated
  // receives
  int offset = 0;
  unsigned num_nonempty_dst_intersections = 0;
  for (int i = 0; i < num_dst_intersections; ++i) {
    xt_mpi_call(MPI_Isend(num_dst_indices_to_remove_per_intersection + i,
                          1, MPI_INT, dst_com[i].rank, HEADER_EXCHANGE_TAG,
                          comm, send_header_requests + i), comm);

    if (num_dst_indices_to_remove_per_intersection[i] > 0) {

      xt_mpi_call(MPI_Isend(dst_indices_to_remove + offset,
                            num_dst_indices_to_remove_per_intersection[i],
                            Xt_int_dt, dst_com[i].rank, DATA_EXCHANGE_TAG,
                            comm,
                            send_data_requests + num_nonempty_dst_intersections),
                  comm);
      offset = offset
        + num_dst_indices_to_remove_per_intersection[i];
      ++num_nonempty_dst_intersections;
    }
  }

  // wait for the receiving of headers to complete
  xt_mpi_call(MPI_Waitall(num_src_intersections + num_dst_intersections,
                          recv_requests, MPI_STATUSES_IGNORE), comm);

  size_t total_num_src_indices_to_recv = 0;

  for (int i = 0; i < num_src_intersections; ++i)
    total_num_src_indices_to_recv
      += (size_t)num_src_indices_to_remove_per_intersection[i];

  unsigned num_nonempty_src_intersections = 0;
  if (total_num_src_indices_to_recv > 0) {

    *src_indices_to_remove = xmalloc(total_num_src_indices_to_recv
                                     * sizeof(**src_indices_to_remove));

    // set up receive for indices that need to be removed
    offset = 0;
    for (int i = 0; i < num_src_intersections; ++i)
      if (num_src_indices_to_remove_per_intersection[i] > 0) {
        xt_mpi_call(MPI_Irecv((*src_indices_to_remove) + offset,
                              num_src_indices_to_remove_per_intersection[i],
                              Xt_int_dt, src_com[i].rank, DATA_EXCHANGE_TAG,
                              comm,
                              recv_requests + num_nonempty_src_intersections),
                    comm);

        offset = offset
          + num_src_indices_to_remove_per_intersection[i];
        ++num_nonempty_src_intersections;
      }

  } else {
    *src_indices_to_remove = NULL;
  }

  /* move data request handles for compact wait */
  memcpy(recv_requests + num_nonempty_src_intersections, send_data_requests,
         num_nonempty_dst_intersections * sizeof (recv_requests[0]));

  // wait until all communication is completed
  xt_mpi_call(MPI_Waitall((int)num_nonempty_src_intersections
                          + (int)num_nonempty_dst_intersections,
                          requests, MPI_STATUSES_IGNORE), comm);

  free(requests);
}

static void
generate_transfer_pos(struct Xt_xmap_intersection_ *xmap,
                      int num_src_intersections,
                      const struct Xt_com_list src_com[num_src_intersections],
                      int num_dst_intersections,
                      const struct Xt_com_list dst_com[num_dst_intersections],
                      Xt_idxlist src_idxlist_local,
                      Xt_idxlist dst_idxlist_local,
                      MPI_Comm comm) {

  int * num_src_indices_to_remove_per_intersection =
    xmalloc((size_t)num_src_intersections
            * sizeof(*num_src_indices_to_remove_per_intersection));
  int * num_dst_indices_to_remove_per_intersection =
    xmalloc((size_t)num_dst_intersections
            * sizeof(*num_dst_indices_to_remove_per_intersection));
  Xt_int * src_indices_to_remove = NULL, * dst_indices_to_remove = NULL;

  generate_dir_transfer_pos_dst(
    num_dst_intersections, dst_com, dst_idxlist_local,
    &(xmap->n_in), &(xmap->in_msg), &dst_indices_to_remove,
    num_dst_indices_to_remove_per_intersection);

  // exchange the points that need to be removed
  exchange_points_to_remove(
    num_src_intersections, src_com, num_dst_intersections, dst_com,
    &src_indices_to_remove, num_src_indices_to_remove_per_intersection,
    dst_indices_to_remove, num_dst_indices_to_remove_per_intersection, comm);

  free(dst_indices_to_remove);
  free(num_dst_indices_to_remove_per_intersection);

  generate_dir_transfer_pos_src(
    num_src_intersections, src_com, src_idxlist_local,
    &(xmap->n_out), &(xmap->out_msg),
    src_indices_to_remove, num_src_indices_to_remove_per_intersection);

  free(src_indices_to_remove);
  free(num_src_indices_to_remove_per_intersection);
}

static int check_destination_coverage(struct exchange_data * in_msgs,
                                      int num_msgs, int num_indices,
                                      MPI_Comm comm) {

  if (num_indices == 0)
    return 0;

  size_t *bit_map;
  int num_bits = (int)(sizeof (*bit_map) * CHAR_BIT);
  size_t num_ints = (size_t)((num_indices + num_bits - 1) / num_bits);

  // create bitmap for destination indices
  bit_map = xcalloc(num_ints, sizeof (*bit_map));
  bit_map[num_ints-1] = ~(((size_t)1 << (num_indices%num_bits)) - 1);

  // for all destination messages
  for (int i = 0; i < num_msgs; ++i) {

    // check all transfer positions
    for (int j = 0; j < in_msgs[i].num_transfer_pos; ++j) {

      int curr_pos = (int)(in_msgs[i].transfer_pos[j]);

      if (curr_pos >= num_indices)
        Xt_abort(comm, "ERROR: invalid transer position (pos >= num_indices)",
                 __FILE__, __LINE__);

      bit_map[curr_pos/num_bits] |= (size_t)1 << (curr_pos%num_bits);
    }
  }

  // check resulting bit map
  size_t all_bits_set = 1;
  for (size_t i = 0; i < num_ints; ++i)
    all_bits_set &= bit_map[i] == SIZE_MAX;

  free(bit_map);

  return !all_bits_set;
}

Xt_xmap
xt_xmap_intersection_new(int num_src_intersections,
                         const struct Xt_com_list
                         src_com[num_src_intersections],
                         int num_dst_intersections,
                         const struct Xt_com_list
                         dst_com[num_dst_intersections],
                         Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                         MPI_Comm comm) {

  Xt_xmap_intersection xmap = xmalloc(sizeof (*xmap));

  xmap->vtable = &xmap_intersection_vtable;

  xt_mpi_call(MPI_Comm_dup(comm, &(xmap->comm)), comm);

  // generate exchange lists
  generate_transfer_pos(xmap,
                        num_src_intersections, src_com,
                        num_dst_intersections, dst_com,
                        src_idxlist, dst_idxlist, comm);

  // we could also calculate the (more precise) max pos using only xmap data
  // but using this simple estimate we are still okay for usage checks
  xmap->max_src_pos = xt_idxlist_get_num_indices(src_idxlist);
  xmap->max_dst_pos = xt_idxlist_get_num_indices(dst_idxlist);

  // check if all indices in the destination index list have been found on
  // the source processes
  if (check_destination_coverage(xmap->in_msg, xmap->n_in,
                                 xt_idxlist_get_num_indices(dst_idxlist), comm))
    Xt_abort(comm, "ERROR: destination intersections do not match with "
             "destination index list", __FILE__, __LINE__);

  return (Xt_xmap)xmap;
}

static int xmap_intersection_get_max_src_pos(Xt_xmap xmap) {
  return xmi(xmap)->max_src_pos;
}

static int xmap_intersection_get_max_dst_pos(Xt_xmap xmap) {
  return xmi(xmap)->max_dst_pos;
}


static void xmap_intersection_delete(Xt_xmap xmap) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  for (int i = 0; i < xmap_intersection->n_in; ++i)
    free(xmap_intersection->in_msg[i].transfer_pos);

  for (int i = 0; i < xmap_intersection->n_out; ++i)
    free(xmap_intersection->out_msg[i].transfer_pos);

  free(xmap_intersection->in_msg);
  free(xmap_intersection->out_msg);

  xt_mpi_call(MPI_Comm_free(&(xmap_intersection->comm)), MPI_COMM_WORLD);
  free(xmap_intersection);
}

static Xt_xmap_iter xmap_intersection_get_in_iterator(Xt_xmap xmap) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  if (xmap_intersection->n_in == 0)
    return NULL;

  Xt_xmap_iter_intersection iter = xmalloc(sizeof (*iter));

  iter->vtable = &xmap_iterator_intersection_vtable;
  iter->msg = xmap_intersection->in_msg;
  iter->msgs_left = xmap_intersection->n_in - 1;

  return (Xt_xmap_iter)iter;
}

static Xt_xmap_iter xmap_intersection_get_out_iterator(Xt_xmap xmap) {

  Xt_xmap_intersection xmap_intersection = xmi(xmap);

  if (xmap_intersection->n_out == 0)
    return NULL;

  Xt_xmap_iter_intersection iter = xmalloc(sizeof (*iter));

  iter->vtable = &xmap_iterator_intersection_vtable;
  iter->msg = xmap_intersection->out_msg;
  iter->msgs_left = xmap_intersection->n_out - 1;

  return (Xt_xmap_iter)iter;
}

static int xmap_iterator_intersection_next(Xt_xmap_iter iter) {

  Xt_xmap_iter_intersection iter_intersection = xmii(iter);

  if (iter_intersection == NULL || iter_intersection->msgs_left == 0)
    return 0;

  iter_intersection->msg++;
  iter_intersection->msgs_left--;

  return 1;
}

static int xmap_intersection_iterator_get_rank(Xt_xmap_iter iter) {

  assert(iter != NULL);
  return xmii(iter)->msg->rank;
}

static int const *
xmap_intersection_iterator_get_transfer_pos(Xt_xmap_iter iter) {

  assert(iter != NULL);
  return xmii(iter)->msg->transfer_pos;
}

static int
xmap_intersection_iterator_get_num_transfer_pos(Xt_xmap_iter iter) {

  assert(iter != NULL);
  return xmii(iter)->msg->num_transfer_pos;
}

static void xmap_intersection_iterator_delete(Xt_xmap_iter iter) {

  free(iter);
}
