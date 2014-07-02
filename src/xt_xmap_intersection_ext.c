/**
 * @file xt_xmap_intersection_ext.c
 *
 * @copyright Copyright  (C)  2012 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "xt/xt_idxlist.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_xmap.h"
#include "xt_xmap_internal.h"
#include "xt/xt_mpi.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_xmap_intersection.h"
#include "xt_arithmetic_util.h"
#include "ensure_array_size.h"

static MPI_Comm     xmap_intersection_ext_get_communicator(Xt_xmap xmap);
static int          xmap_intersection_ext_get_num_destinations(Xt_xmap xmap);
static int          xmap_intersection_ext_get_num_sources(Xt_xmap xmap);
static void
xmap_intersection_ext_get_destination_ranks(Xt_xmap xmap, int * ranks);
static void
xmap_intersection_ext_get_source_ranks(Xt_xmap xmap, int * ranks);
static Xt_xmap_iter xmap_intersection_ext_get_in_iterator(Xt_xmap xmap);
static Xt_xmap_iter xmap_intersection_ext_get_out_iterator(Xt_xmap xmap);
static void xmap_intersection_ext_delete(Xt_xmap xmap);
static int xmap_intersection_ext_get_max_src_pos(Xt_xmap xmap);
static int xmap_intersection_ext_get_max_dst_pos(Xt_xmap xmap);


static const struct Xt_xmap_vtable xmap_intersection_vtable = {
        .get_communicator      = xmap_intersection_ext_get_communicator,
        .get_num_destinations  = xmap_intersection_ext_get_num_destinations,
        .get_num_sources       = xmap_intersection_ext_get_num_sources,
        .get_destination_ranks = xmap_intersection_ext_get_destination_ranks,
        .get_source_ranks      = xmap_intersection_ext_get_source_ranks,
        .get_out_iterator      = xmap_intersection_ext_get_out_iterator,
        .get_in_iterator       = xmap_intersection_ext_get_in_iterator,
        .delete                = xmap_intersection_ext_delete,
        .get_max_src_pos       = xmap_intersection_ext_get_max_src_pos,
        .get_max_dst_pos       = xmap_intersection_ext_get_max_dst_pos};

struct exchange_ext {
  // list of relative position extents in index list to send or receive
  struct Xt_pos_ext *transfer_pos_ext;
  /* generated on-demand */
  int *transfer_pos;
  int num_transfer_pos, num_transfer_pos_ext;
  int rank;
};

struct Xt_xmap_intersection_ext_ {

  const struct Xt_xmap_vtable * vtable;

  struct exchange_ext *in_msg, *out_msg;
  int n_in, n_out;

  // we need the max position in order to enable quick range-checks
  // for xmap-users like redist
  int max_src_pos; // max possible pos over all src transfer_pos (always >= 0)
  int max_dst_pos; // same for dst

  MPI_Comm comm;
};

typedef struct Xt_xmap_intersection_ext_ *Xt_xmap_intersection_ext;

static inline Xt_xmap_intersection_ext xmie(void *xmap)
{
  return (Xt_xmap_intersection_ext)xmap;
}

static MPI_Comm xmap_intersection_ext_get_communicator(Xt_xmap xmap)
{
  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);
  return xmap_intersection_ext->comm;
}

static int xmap_intersection_ext_get_num_destinations(Xt_xmap xmap)
{
  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);
  // the number of destination equals the number of source messages
  return xmap_intersection_ext->n_out;
}

static int xmap_intersection_ext_get_num_sources(Xt_xmap xmap)
{
  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);
  // the number of sources equals the number of destination messages
  return xmap_intersection_ext->n_in;
}

static void
xmap_intersection_ext_get_destination_ranks(Xt_xmap xmap, int *restrict ranks)
{
  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);
  size_t n_out = (size_t)xmap_intersection_ext->n_out;
  struct exchange_ext *restrict out_msg = xmap_intersection_ext->out_msg;
  for (size_t i = 0; i < n_out; ++i)
    ranks[i] = out_msg[i].rank;
}

static void
xmap_intersection_ext_get_source_ranks(Xt_xmap xmap, int *restrict ranks)
{
  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);
  size_t n_in = (size_t)xmap_intersection_ext->n_in;
  struct exchange_ext *restrict in_msg = xmap_intersection_ext->in_msg;
  for (size_t i = 0; i < n_in; ++i)
    ranks[i] = in_msg[i].rank;
}

static int xmap_intersection_ext_get_max_src_pos(Xt_xmap xmap) {
  return xmie(xmap)->max_src_pos;
}

static int xmap_intersection_ext_get_max_dst_pos(Xt_xmap xmap) {
  return xmie(xmap)->max_dst_pos;
}

static void
xt_free_exchange_ext(size_t num_msg, struct exchange_ext *restrict msg)
{
  for (size_t i = 0; i < num_msg; ++i) {
    free(msg[i].transfer_pos);
    free(msg[i].transfer_pos_ext);
  }
  free(msg);
}


static void xmap_intersection_ext_delete(Xt_xmap xmap) {

  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);

  xt_free_exchange_ext((size_t)xmap_intersection_ext->n_in,
                       xmap_intersection_ext->in_msg);

  xt_free_exchange_ext((size_t)xmap_intersection_ext->n_out,
                       xmap_intersection_ext->out_msg);

  xt_mpi_call(MPI_Comm_free(&(xmap_intersection_ext->comm)), MPI_COMM_WORLD);
  free(xmap_intersection_ext);
}

static void
generate_transfer_ext(struct Xt_xmap_intersection_ext_ *xmap,
                      int num_src_intersections,
                      const struct Xt_com_list src_com[num_src_intersections],
                      int num_dst_intersections,
                      const struct Xt_com_list dst_com[num_dst_intersections],
                      Xt_idxlist src_idxlist_local,
                      Xt_idxlist dst_idxlist_local,
                      MPI_Comm comm);

Xt_xmap
xt_xmap_intersection_ext_new(int num_src_intersections,
                             const struct Xt_com_list
                             src_com[num_src_intersections],
                             int num_dst_intersections,
                             const struct Xt_com_list
                             dst_com[num_dst_intersections],
                             Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                             MPI_Comm comm) {

  Xt_xmap_intersection_ext xmap = xmalloc(sizeof (*xmap));

  xmap->vtable = &xmap_intersection_vtable;

  /* this dup could be saved if we knew comm was already internal to
   * the library somehow */
  xt_mpi_call(MPI_Comm_dup(comm, &(xmap->comm)), comm);
  comm = xmap->comm;

  // generate exchange lists
  generate_transfer_ext(xmap,
                        num_src_intersections, src_com,
                        num_dst_intersections, dst_com,
                        src_idxlist, dst_idxlist, comm);

  // we could also calculate the (more precise) max pos using only xmap data
  // but using this simple estimate we are still okay for usage checks
  xmap->max_src_pos = xt_idxlist_get_num_indices(src_idxlist);
  xmap->max_dst_pos = xt_idxlist_get_num_indices(dst_idxlist);

  return (Xt_xmap)xmap;
}

struct Xt_pos_ext_vec {
  size_t num_pos_ext;
  struct Xt_pos_ext *pos_ext;
};

static struct Xt_pos_ext_vec
generate_dir_transfer_ext_dst(
  int num_intersections,
  const struct Xt_com_list intersections[num_intersections],
  Xt_idxlist mypart_idxlist,
  int *resCount,
  struct exchange_ext **resSets,
  int (*restrict dst_removals_per_intersection)[2]);

static bool
xt_idxlist_pos_ext_cover_full(Xt_idxlist idxlist,
                              struct Xt_pos_ext_vec coverage);

static struct Xt_pos_ext *
exchange_pos_ext_modifications(
  int num_src_intersections,
  const struct Xt_com_list src_com[num_src_intersections],
  int num_dst_intersections,
  const struct Xt_com_list dst_com[num_dst_intersections],
  struct exchange_ext dst_ext[num_dst_intersections],
  int (*restrict src_removals_per_intersection)[2],
  int (*restrict dst_removals_per_intersection)[2],
  MPI_Comm comm);

static void
remap_dst_intersections(int num_dst_intersections,
                        const struct Xt_com_list dst_com[num_dst_intersections],
                        Xt_idxlist mypart_idxlist,
                        int resCount,
                        struct exchange_ext resSets[resCount],
                        int (*removals_per_intersection)[2]);

static void
generate_dir_transfer_pos_ext_src(
  int num_intersections,
  const struct Xt_com_list intersections[num_intersections],
  Xt_idxlist mypart_idxlist,
  int *resCount,
  struct exchange_ext **resSets,
  int (*restrict removals_per_intersection)[2],
  struct Xt_pos_ext *pos_updates);

static void
generate_transfer_ext(struct Xt_xmap_intersection_ext_ *xmap,
                      int num_src_intersections,
                      const struct Xt_com_list src_com[num_src_intersections],
                      int num_dst_intersections,
                      const struct Xt_com_list dst_com[num_dst_intersections],
                      Xt_idxlist src_idxlist_local,
                      Xt_idxlist dst_idxlist_local,
                      MPI_Comm comm) {

  /* {dst|src}_removals_per_intersection[i][0] denotes the number of
   * indices to be removed from the intersection with {src|dst}_com[i].rank.
   * {dst|src}_removals_per_intersection[rank][1] denotes the number of
   * pos_ext needed to represent this change (0 if either none or all
   * indices got removed).
   */
  int (*src_removals_per_intersection)[2] =
    xmalloc(((size_t)num_dst_intersections + (size_t)num_src_intersections)
            * sizeof(*src_removals_per_intersection)),
    (*dst_removals_per_intersection)[2]
    = src_removals_per_intersection + num_src_intersections;

  {
    struct Xt_pos_ext_vec coverage
      = generate_dir_transfer_ext_dst(
        num_dst_intersections, dst_com, dst_idxlist_local,
        &(xmap->n_in), &(xmap->in_msg), dst_removals_per_intersection);

    if (!xt_idxlist_pos_ext_cover_full(dst_idxlist_local, coverage))
      Xt_abort(comm, "ERROR: destination intersections do not match with "
               "destination index list", __FILE__, __LINE__);
    free(coverage.pos_ext);
  }

  // exchange pos_ext of lists where additional indices need to be removed
  struct Xt_pos_ext *pos_updates
    = exchange_pos_ext_modifications(
      num_src_intersections, src_com, num_dst_intersections, dst_com,
      xmap->in_msg,
      src_removals_per_intersection,
      dst_removals_per_intersection, comm);

  remap_dst_intersections(num_dst_intersections, dst_com, dst_idxlist_local,
                          xmap->n_in, xmap->in_msg,
                          dst_removals_per_intersection);

  src_removals_per_intersection =
    xrealloc(src_removals_per_intersection, (size_t)num_src_intersections
             * sizeof(*src_removals_per_intersection));

  generate_dir_transfer_pos_ext_src(
    num_src_intersections, src_com, src_idxlist_local,
    &(xmap->n_out), &(xmap->out_msg),
    src_removals_per_intersection, pos_updates);

  free(src_removals_per_intersection);
  free(pos_updates);
}

struct Xt_pos_ext_overlap {
  int skip, overlap, tail;
};

struct Xt_pos_ext_overlap
Xt_get_pos_ext_overlap(struct Xt_pos_ext a, struct Xt_pos_ext b)
{
  /* == 0 if a.size >= 0 ; == ~0 if a.size < 0 */
  int aSizeMaskNeg = isign_mask(a.size),
    /* compute start and end indices of ranges */
    a_s = a.start +  (aSizeMaskNeg & (a.size + 1)),
    a_e   = a.start + (~aSizeMaskNeg & (a.size - 1)),
    bSizeMaskNeg = isign_mask(b.size),
    b_s = b.start +  (bSizeMaskNeg & (b.size + 1)),
    b_e = b.start + (~bSizeMaskNeg & (b.size - 1));
  /* does overlap exist? */
  if ((b_s > a_e) | (a_s > b_e))
    return (struct Xt_pos_ext_overlap){ a.size, 0, 0};
  else {
    /* determine length of overlap parts */
    int lowSkipA = b_s - a_s;
    int lowSkipB = -lowSkipA;
    lowSkipA = (int)((unsigned)(lowSkipA + abs(lowSkipA))/2U);
    lowSkipB = (int)((unsigned)(lowSkipB + abs(lowSkipB))/2U);
    int overlapLen = imin(b_e - b_s - lowSkipB + 1,
                          abs(a.size) - lowSkipA);
    int highSkipA = abs(a.size) - lowSkipA - overlapLen;
    /* then adjust lengths to direction of overlap (from
     * perspective of a */
    int aSkipLen = (~aSizeMaskNeg & lowSkipA)
      | (aSizeMaskNeg & -highSkipA),
      aTailLen = (aSizeMaskNeg & -lowSkipA)
      | (~aSizeMaskNeg & highSkipA);
    return (struct Xt_pos_ext_overlap){ aSkipLen, overlapLen, aTailLen };
  }
}



static void
cut_pos_ext_from_pos_exts(struct Xt_pos_ext pos_ext,
                          size_t *num_pos_exts,
                          size_t *size_pos_exts,
                          struct Xt_pos_ext **pos_exts);

static struct Xt_pos_ext_vec
generate_dir_transfer_ext_dst(
  int num_intersections,
  const struct Xt_com_list intersections[num_intersections],
  Xt_idxlist mypart_idxlist,
  int *resCount, struct exchange_ext **resSets,
  int (*restrict dst_removals_per_intersection)[2])
{
  struct exchange_ext *restrict resSets_
    = *resSets = xmalloc((size_t)num_intersections * sizeof(**resSets));

  int new_num_intersections = 0;

  // we have to enforce single_match_only not only within a single
  // intersection, but also between all intersections
  /* ranges already covered from previous intersections, i.e. which
   * must not be transmitted twice */
  size_t num_mypart_covered_pos_ext = 0, size_mypart_covered_pos_ext = 8;
  struct Xt_pos_ext *restrict mypart_covered_pos_ext
    = xmalloc(sizeof (*mypart_covered_pos_ext) * size_mypart_covered_pos_ext);


  struct Xt_pos_ext *restrict isect_transfer_pos_ext = NULL;

  for (int i = 0; i < num_intersections; ++i) {

    int num_stripes, num_indices_to_remove = 0;
    struct Xt_stripe *restrict intersection_idxstripes;
    xt_idxlist_get_index_stripes(intersections[i].list,
                                 (struct Xt_stripe **)&intersection_idxstripes,
                                 &num_stripes);
    struct Xt_pos_ext *restrict isect_pos_exts = NULL;
    int num_isect_pos_exts;
    int retval = xt_idxlist_get_pos_exts_of_index_stripes(
      mypart_idxlist, num_stripes, intersection_idxstripes,
      &num_isect_pos_exts, (struct Xt_pos_ext **)&isect_pos_exts, 1);
    assert(retval == 0);
    int isect_pos_exts_size_psum = 0;
    int intersection_size = xt_idxlist_get_num_indices(intersections[i].list);
    /* start with all indices from intersection as used,
       later split ranges, if overlaps are found */
    size_t num_isect_transfer_pos_ext = 1, size_isect_transfer_pos_ext = 8;
    isect_transfer_pos_ext
      = xrealloc(intersection_idxstripes, sizeof (*isect_transfer_pos_ext)
                 * size_isect_transfer_pos_ext);
    intersection_idxstripes = NULL;
    isect_transfer_pos_ext[0]
      = (struct Xt_pos_ext){ .start = 0, .size = intersection_size };
    /* find overlap(s) with previously found ranges for all
     * stripes mapped to position extents */
    for (size_t j = 0; j < (size_t)num_isect_pos_exts; ++j) {
      struct Xt_pos_ext isect_pos_ext = isect_pos_exts[j];
      /* ensure isect_pos_ext is oriented with ascending positions */
      int isign_mask_isect_pos_ext_size = isign_mask(isect_pos_ext.size);
      isect_pos_ext.start
        += isign_mask_isect_pos_ext_size & (isect_pos_ext.size + 1);
      int isect_pos_ext_orig_size = isect_pos_ext.size;
      isect_pos_ext.size = abs(isect_pos_ext.size);
      /* if isect_pos_ext is of negative orientation, add nothing at
       * the end but subtract from already incremented psum */
      isect_pos_exts_size_psum += isect_pos_ext.size;
      int progress = -isect_pos_ext.size;
      for (size_t k = 0; k < num_mypart_covered_pos_ext; ++k) {
        struct Xt_pos_ext_overlap overlap_desc
          = Xt_get_pos_ext_overlap(isect_pos_ext, mypart_covered_pos_ext[k]);
        if (overlap_desc.overlap) {
          /* insert overlap into updates
           * by ...*/
          /* ...first inserting the skipped part into
           * mypart_covered_pos_ext, since that is sorted
           * and obviously precedes mypart_covered_pos_ext[k],
           * mypart_covered_pos_ext[k] can be seemlessly extended...
           */
          mypart_covered_pos_ext[k].start -= overlap_desc.skip;
          mypart_covered_pos_ext[k].size += overlap_desc.skip;
          /* ...and optionally merged with its predecessor, if the
           * intervening range becomes zero, ... */
          if (k > 0
              && (mypart_covered_pos_ext[k].start
                  == (mypart_covered_pos_ext[k - 1].start
                      + mypart_covered_pos_ext[k - 1].size)))
          {
            mypart_covered_pos_ext[k - 1].size
              += mypart_covered_pos_ext[k].size;
            memmove(mypart_covered_pos_ext + k, mypart_covered_pos_ext + k + 1,
                    (--num_mypart_covered_pos_ext - k)
                    * sizeof (*mypart_covered_pos_ext));
            --k;
          }
          progress = (~isign_mask_isect_pos_ext_size
                      & (progress + overlap_desc.skip))
            | (isign_mask_isect_pos_ext_size
               & (isect_pos_ext_orig_size + overlap_desc.tail));
          /* ... then splitting isect_transfer_pos accordingly, ... */
          num_indices_to_remove += overlap_desc.overlap;
          cut_pos_ext_from_pos_exts((struct Xt_pos_ext){
              .start = isect_pos_exts_size_psum + progress,
                .size = overlap_desc.overlap },
            &num_isect_transfer_pos_ext, &size_isect_transfer_pos_ext,
            (struct Xt_pos_ext **)&isect_transfer_pos_ext);
          progress += overlap_desc.overlap;
          /* ... lastly the search can continue with the tail ... */
          isect_pos_ext.start += overlap_desc.skip + overlap_desc.overlap;
          /* ... if there is any */
          if ((isect_pos_ext.size = overlap_desc.tail))
            continue;
          /* ... or start the next intersection range */
          else
            goto next_isect_pos_ext;
        }
        /* isect_pos_ext does not overlap current part,
         * can there be future overlaps, i.e. is already
         * isect_pos_ext.start + isect_pos_ext.size
         *        < mypart_covered_pos_ext[k].start ?
         */
        if (isect_pos_ext.start + isect_pos_ext.size
            < mypart_covered_pos_ext[k].start)
        {
          /* isect precedes mypart_covered_pos_ext[k] with a hole
           * but might be a seemless extension of
           * mypart_covered_pos_ext[k-1] */
          if (k > 0
              && (mypart_covered_pos_ext[k - 1].start
                  + mypart_covered_pos_ext[k - 1].size == isect_pos_ext.start))
            mypart_covered_pos_ext[k - 1].size += isect_pos_ext.size;
          else
          {
            ENSURE_ARRAY_SIZE(mypart_covered_pos_ext,
                              size_mypart_covered_pos_ext,
                              num_mypart_covered_pos_ext + 1);
            memmove(mypart_covered_pos_ext + k + 1, mypart_covered_pos_ext + k,
                    sizeof (*mypart_covered_pos_ext)
                    * (num_mypart_covered_pos_ext - k));
            mypart_covered_pos_ext[k] = isect_pos_ext;
            ++num_mypart_covered_pos_ext;
          }
          goto next_isect_pos_ext;
        }
        else if (isect_pos_ext.start + isect_pos_ext.size
                 == mypart_covered_pos_ext[k].start)
        {
          mypart_covered_pos_ext[k].start = isect_pos_ext.start;
          mypart_covered_pos_ext[k].size += isect_pos_ext.size;
          if (k > 0
              && (mypart_covered_pos_ext[k].start
                  == (mypart_covered_pos_ext[k - 1].start
                      + mypart_covered_pos_ext[k - 1].size)))
          {
            mypart_covered_pos_ext[k - 1].size
              += mypart_covered_pos_ext[k].size;
            memmove(mypart_covered_pos_ext + k, mypart_covered_pos_ext + k + 1,
                    (--num_mypart_covered_pos_ext - k)
                    * sizeof (*mypart_covered_pos_ext));
          }
          goto next_isect_pos_ext;
        }
      }
      /* isect_pos_ext was not found -> append */
      if (num_mypart_covered_pos_ext > 0
          && (mypart_covered_pos_ext[num_mypart_covered_pos_ext - 1].start
              + mypart_covered_pos_ext[num_mypart_covered_pos_ext - 1].size
              == isect_pos_ext.start))
      {
        mypart_covered_pos_ext[num_mypart_covered_pos_ext - 1].size
          += isect_pos_ext.size;
      }
      else
      {
        ENSURE_ARRAY_SIZE(mypart_covered_pos_ext, size_mypart_covered_pos_ext,
                          num_mypart_covered_pos_ext + 1);
        mypart_covered_pos_ext[num_mypart_covered_pos_ext] = isect_pos_ext;
        ++num_mypart_covered_pos_ext;
      }
    next_isect_pos_ext:
      ;
    }

    if (intersection_size > num_indices_to_remove) {
      resSets_[new_num_intersections].transfer_pos_ext
        = xrealloc(isect_transfer_pos_ext, sizeof (*isect_transfer_pos_ext)
                   * num_isect_transfer_pos_ext);
      /* start with empty cache of positions to transfer */
      resSets_[new_num_intersections].transfer_pos = NULL;
      resSets_[new_num_intersections].num_transfer_pos
        = intersection_size - num_indices_to_remove;
      resSets_[new_num_intersections].num_transfer_pos_ext
        = (int)num_isect_transfer_pos_ext;
      resSets_[new_num_intersections].rank = intersections[i].rank;
      ++new_num_intersections;
      isect_transfer_pos_ext = NULL;
    }
    dst_removals_per_intersection[i][0] = num_indices_to_remove;
    dst_removals_per_intersection[i][1]
      = ((num_indices_to_remove == intersection_size)
         | (num_indices_to_remove == 0))?0:(int)num_isect_transfer_pos_ext;
    free(isect_transfer_pos_ext);
    free(isect_pos_exts);
  }
  mypart_covered_pos_ext
    = xrealloc(mypart_covered_pos_ext,
               sizeof (*mypart_covered_pos_ext) * num_mypart_covered_pos_ext);

  *resCount = new_num_intersections;
  if (num_intersections != new_num_intersections)
    *resSets = xrealloc(resSets_,
                        (size_t)new_num_intersections * sizeof(**resSets));
  return (struct Xt_pos_ext_vec){
    .num_pos_ext = num_mypart_covered_pos_ext,
      .pos_ext = mypart_covered_pos_ext};
}

static void
cut_pos_ext_from_pos_exts(struct Xt_pos_ext pos_ext,
                          size_t *num_pos_exts,
                          size_t *size_pos_exts,
                          struct Xt_pos_ext **pos_exts)
{
  struct Xt_pos_ext *restrict pos_exts_ = *pos_exts;
  size_t num_pos_exts_ = *num_pos_exts;
  size_t i = num_pos_exts_;
  while (pos_exts_[--i].start > pos_ext.start)
    ;
  int db_skip = pos_ext.start - pos_exts_[i].start;
  if ((!db_skip) & (pos_ext.size == pos_exts_[i].size))
  {
    /* delete fully overlapped transfer part */
    memmove(pos_exts_ + i, pos_exts_ + i + 1,
            sizeof (*pos_exts_) * (num_pos_exts_ - i - 1));
    *num_pos_exts = --num_pos_exts_;
  }
  else if (db_skip + pos_ext.size == pos_exts_[i].size)
  {
    /* pos_ext overlaps end of pos_exts_[i] */
    pos_exts_[i].size -= pos_ext.size;
  }
  else if (db_skip == 0)
  {
    /* pos_ext overlaps start of pos_exts_[i] */
    pos_exts_[i].start = pos_ext.start + pos_ext.size;
    pos_exts_[i].size -= pos_ext.size;
  }
  else
  {
    struct Xt_pos_ext orig = pos_exts_[i];
    ENSURE_ARRAY_SIZE(*pos_exts, *size_pos_exts, num_pos_exts_ + 1);
    pos_exts_ = *pos_exts;
    memmove(pos_exts_ + i + 1, pos_exts_ + i,
            (num_pos_exts_ - i) * sizeof (*pos_exts_));
    pos_exts_[i] = (struct Xt_pos_ext){.start = orig.start,
                                       .size = db_skip };
    pos_exts_[i + 1] = (struct Xt_pos_ext){
      .start = pos_ext.start + pos_ext.size,
      .size = orig.size - db_skip - pos_ext.size };
    *num_pos_exts = ++num_pos_exts_;
  }
}

static struct Xt_pos_ext *
exchange_pos_ext_modifications(
  int num_src_intersections,
  const struct Xt_com_list src_com[num_src_intersections],
  int num_dst_intersections,
  const struct Xt_com_list dst_com[num_dst_intersections],
  struct exchange_ext dst_ext[num_dst_intersections],
  int (*restrict src_removals_per_intersection)[2],
  int (*restrict dst_removals_per_intersection)[2],
  MPI_Comm comm)
{
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
    xt_mpi_call(MPI_Irecv(src_removals_per_intersection[i],
                          2, MPI_INT, src_com[i].rank, HEADER_EXCHANGE_TAG,
                          comm, recv_requests + i), comm);

  /* send rebuilt pos_ext vectors that needed to be modified on the
   * target side due to duplicated receives
   */
  unsigned num_active_dst = 0, num_dst_changes = 0;
  for (int i = 0; i < num_dst_intersections; ++i) {
    xt_mpi_call(MPI_Isend(dst_removals_per_intersection[i],
                          2, MPI_INT, dst_com[i].rank, HEADER_EXCHANGE_TAG,
                          comm, send_header_requests + i), comm);

    if (dst_removals_per_intersection[i][1] > 0) {

      assert(dst_removals_per_intersection[i][1]
             == dst_ext[num_active_dst].num_transfer_pos_ext);
      assert(dst_com[i].rank
             == dst_ext[num_active_dst].rank);
      xt_mpi_call(MPI_Isend(dst_ext[num_active_dst].transfer_pos_ext,
                            dst_removals_per_intersection[i][1],
                            MPI_2INT, dst_com[i].rank, DATA_EXCHANGE_TAG,
                            comm, send_data_requests + num_dst_changes),
                  comm);
      ++num_dst_changes;
    }
    num_active_dst += (unsigned)((dst_removals_per_intersection[i][0] == 0)
                                 | (dst_removals_per_intersection[i][1] != 0));
  }

  // wait for the receiving of headers to complete
  xt_mpi_call(MPI_Waitall(num_src_intersections + num_dst_intersections,
                          recv_requests, MPI_STATUSES_IGNORE), comm);

  size_t total_num_pos_ext_to_recv = 0;

  for (size_t i = 0; i < (size_t)num_src_intersections; ++i)
    total_num_pos_ext_to_recv += (size_t)src_removals_per_intersection[i][1];

  struct Xt_pos_ext *src_updated_pos_ext;
  unsigned num_src_changes = 0;
  if (total_num_pos_ext_to_recv > 0) {

    src_updated_pos_ext
      = xmalloc(total_num_pos_ext_to_recv * sizeof(*src_updated_pos_ext));

    /* set up receive for pos_ext that need to be modified because
     * indices needed to be removed from the intersection */
    size_t offset = 0;
    for (int i = 0; i < num_src_intersections; ++i)
      if (src_removals_per_intersection[i][1] > 0) {
        xt_mpi_call(MPI_Irecv(src_updated_pos_ext + offset,
                              src_removals_per_intersection[i][1],
                              MPI_2INT, src_com[i].rank, DATA_EXCHANGE_TAG,
                              comm, recv_requests + num_src_changes), comm);

        offset += (size_t)src_removals_per_intersection[i][1];
        ++num_src_changes;
      }
  } else
    src_updated_pos_ext = NULL;

  /* move data request handles for compact wait */
  memcpy(recv_requests + num_src_changes, send_data_requests,
         num_dst_changes * sizeof (recv_requests[0]));

  // wait until all communication is completed
  xt_mpi_call(MPI_Waitall((int)num_src_changes + (int)num_dst_changes,
                          recv_requests, MPI_STATUSES_IGNORE), comm);

  free(requests);
  return src_updated_pos_ext;
}

/* tests if sorted pos_ext in \a coverage do indeed fully cover \a idxlist
 */
static bool
xt_idxlist_pos_ext_cover_full(Xt_idxlist idxlist,
                              struct Xt_pos_ext_vec coverage)
{
  int idxlist_size = xt_idxlist_get_num_indices(idxlist);
  if ((idxlist_size == 0) & (coverage.num_pos_ext == 0))
    return true;
  if ((idxlist_size == 0) ^ (coverage.num_pos_ext == 0))
    return false;
  /* at this point both idxlist and coverage are non-empty */
  int after_last_end_pos = 0;
  int continuations_hold = 1;
  for (size_t i = 0; i < coverage.num_pos_ext; ++i)
  {
    continuations_hold &= (coverage.pos_ext[i].start == after_last_end_pos);
    after_last_end_pos += coverage.pos_ext[i].size;
  }
  continuations_hold &= (after_last_end_pos == idxlist_size);
  return continuations_hold;
}

static void
remap_intersection(Xt_idxlist mypart_idxlist,
                   Xt_idxlist intersection,
                   size_t num_pos_updates,
                   struct Xt_pos_ext pos_updates[num_pos_updates],
                   struct exchange_ext *resSet,
                   int single_match_only);

static void
remap_dst_intersections(
  int num_dst_intersections,
  const struct Xt_com_list intersections[num_dst_intersections],
  Xt_idxlist mypart_idxlist,
  int resCount,
  struct exchange_ext resSets[resCount],
  int (*removals_per_intersection)[2])
{
  size_t resIdx = 0;
  for (size_t i = 0; i < (size_t)num_dst_intersections; ++i)
  {
    int intersection_size
      = xt_idxlist_get_num_indices(intersections[i].list);

    int num_indices_to_remove = removals_per_intersection[i][0];

    if (num_indices_to_remove != intersection_size) {} else
      /* intersection is made redundant */
      continue;

    struct Xt_pos_ext *pos_updates = resSets[resIdx].transfer_pos_ext;
    remap_intersection(mypart_idxlist, intersections[i].list,
                       (size_t)removals_per_intersection[i][1],
                       pos_updates, resSets + resIdx, 1);
    free(pos_updates);
    ++resIdx;
  }
  assert(resIdx == (size_t)resCount);
}


/* compute updated positions for send direction */
static void
generate_dir_transfer_pos_ext_src(
  int num_intersections,
  const struct Xt_com_list intersections[num_intersections],
  Xt_idxlist mypart_idxlist,
  int *resCount,
  struct exchange_ext **resSets,
  int (*restrict removals_per_intersection)[2],
  struct Xt_pos_ext *pos_updates)
{
  struct exchange_ext *restrict resSets_ = *resSets
    = xmalloc((size_t)num_intersections * sizeof(**resSets));

  /* count total number of intersections that results */
  int new_num_intersections = 0;
  /* indexes into pos_updates */
  size_t intersection_pos_ext = 0;

  for (int i = 0; i < num_intersections; ++i) {

    int intersection_size
      = xt_idxlist_get_num_indices(intersections[i].list);

    int num_indices_to_remove = removals_per_intersection[i][0];

    if (num_indices_to_remove != intersection_size) {} else
      /* intersection is made redundant */
      continue;

    remap_intersection(mypart_idxlist, intersections[i].list,
                       (size_t)removals_per_intersection[i][1],
                       pos_updates + intersection_pos_ext,
                       resSets_ + new_num_intersections, 0);

    /* evaluate cache lazily */
    resSets_[new_num_intersections].transfer_pos = NULL;
    resSets_[new_num_intersections].num_transfer_pos
      = intersection_size - num_indices_to_remove;
    resSets_[new_num_intersections].rank = intersections[i].rank;
    new_num_intersections++;
    intersection_pos_ext += (size_t)removals_per_intersection[i][1];
  }

  *resCount = new_num_intersections;
  if (num_intersections != new_num_intersections)
    *resSets = xrealloc(resSets_,
                        (size_t)new_num_intersections * sizeof(**resSets));
}

/* match index stripes of intersection to corresponding positions in
 * partition list, optionally updating the stripes
 * @param num_pos_updates number of position extents describing the
 *                        subset of positions from intersections to use
 * @param pos_updates list of position extents to use from \a intersection
 */
static void
remap_intersection(Xt_idxlist mypart_idxlist,
                   Xt_idxlist intersection,
                   size_t num_pos_updates,
                   struct Xt_pos_ext pos_updates[num_pos_updates],
                   struct exchange_ext *resSet,
                   int single_match_only)
{
  struct Xt_stripe *restrict intersection_idxstripes;
  int num_stripes;
  xt_idxlist_get_index_stripes(intersection,
                               (struct Xt_stripe **)&intersection_idxstripes,
                               &num_stripes);
  if (num_pos_updates) {
    /* trim intersection_idxstripes to those actually used */
    size_t num_refined_intersection_idxstripes = 0,
      size_refined_intersection_idxstripes = num_pos_updates;
    struct Xt_stripe *restrict refined_intersection_idxstripes
      = xmalloc(size_refined_intersection_idxstripes
                * sizeof (*refined_intersection_idxstripes));
    size_t i_stripe = 0;
    int nstrides_psum = 0;
    for (size_t i_pos_ext = 0; i_pos_ext < num_pos_updates; ++i_pos_ext)
    {
      int pos = pos_updates[i_pos_ext].start;
      int size = pos_updates[i_pos_ext].size;
      while (nstrides_psum + intersection_idxstripes[i_stripe].nstrides <= pos)
      {
        nstrides_psum += intersection_idxstripes[i_stripe].nstrides;
        ++i_stripe;
      }
      do {
        int instripe_pos = pos - nstrides_psum;
        ENSURE_ARRAY_SIZE(refined_intersection_idxstripes,
                          size_refined_intersection_idxstripes,
                          num_refined_intersection_idxstripes + 1);
        struct Xt_stripe cur_stripe = intersection_idxstripes[i_stripe];
        int cur_stripe_nstrides = cur_stripe.nstrides;
        int overlap = imin(cur_stripe_nstrides - instripe_pos, size);
        cur_stripe.start
          = (Xt_int)(cur_stripe.start
                     + (Xt_int)instripe_pos * cur_stripe.stride);
        cur_stripe.nstrides = overlap;
        refined_intersection_idxstripes[num_refined_intersection_idxstripes]
          = cur_stripe;
        ++num_refined_intersection_idxstripes;
        i_stripe += (instripe_pos + overlap == cur_stripe_nstrides);
        nstrides_psum += (instripe_pos + overlap == cur_stripe_nstrides)
          ? cur_stripe_nstrides : 0;
        pos += overlap;
        size -= overlap;
      } while (size);
    }
    free(intersection_idxstripes);
    intersection_idxstripes = refined_intersection_idxstripes;
    num_stripes = (int)num_refined_intersection_idxstripes;
  }

  /* match back intersection_idxstripes to positions in mypart */
  resSet->transfer_pos_ext = NULL;
  int retval = xt_idxlist_get_pos_exts_of_index_stripes(
    mypart_idxlist, num_stripes, intersection_idxstripes,
    &resSet->num_transfer_pos_ext,
    &resSet->transfer_pos_ext, single_match_only);
  assert(retval == 0);
  free(intersection_idxstripes);
}


/* iterator operations */

static int xmap_intersection_ext_iterator_next(Xt_xmap_iter iter);
static int xmap_intersection_ext_iterator_get_rank(Xt_xmap_iter iter);
static int const *
xmap_intersection_ext_iterator_get_transfer_pos(Xt_xmap_iter iter);
static int
xmap_intersection_ext_iterator_get_num_transfer_pos(Xt_xmap_iter iter);
static const struct Xt_pos_ext *
xmap_intersection_ext_iterator_get_transfer_pos_ext(Xt_xmap_iter iter);
static int
xmap_intersection_ext_iterator_get_num_transfer_pos_ext(Xt_xmap_iter iter);
static void xmap_intersection_ext_iterator_delete(Xt_xmap_iter iter);

static const struct Xt_xmap_iter_vtable
xmap_iterator_intersection_ext_vtable = {
  .next                 = xmap_intersection_ext_iterator_next,
  .get_rank             = xmap_intersection_ext_iterator_get_rank,
  .get_transfer_pos     = xmap_intersection_ext_iterator_get_transfer_pos,
  .get_num_transfer_pos = xmap_intersection_ext_iterator_get_num_transfer_pos,
  .get_transfer_pos_ext = xmap_intersection_ext_iterator_get_transfer_pos_ext,
  .get_num_transfer_pos_ext
  = xmap_intersection_ext_iterator_get_num_transfer_pos_ext,
  .delete               = xmap_intersection_ext_iterator_delete};

typedef struct Xt_xmap_iter_intersection_ext_ *Xt_xmap_iter_intersection_ext;

struct Xt_xmap_iter_intersection_ext_ {

  const struct Xt_xmap_iter_vtable *vtable;

  struct exchange_ext *msg;
  int msgs_left;
};

static Xt_xmap_iter xmap_intersection_ext_get_in_iterator(Xt_xmap xmap) {

  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);

  if (xmap_intersection_ext->n_in == 0)
    return NULL;

  Xt_xmap_iter_intersection_ext iter = xmalloc(sizeof (*iter));

  iter->vtable = &xmap_iterator_intersection_ext_vtable;
  iter->msg = xmap_intersection_ext->in_msg;
  iter->msgs_left = xmap_intersection_ext->n_in - 1;

  return (Xt_xmap_iter)iter;
}

static Xt_xmap_iter xmap_intersection_ext_get_out_iterator(Xt_xmap xmap) {

  Xt_xmap_intersection_ext xmap_intersection_ext = xmie(xmap);

  if (xmap_intersection_ext->n_out == 0)
    return NULL;

  Xt_xmap_iter_intersection_ext iter = xmalloc(sizeof (*iter));

  iter->vtable = &xmap_iterator_intersection_ext_vtable;
  iter->msg = xmap_intersection_ext->out_msg;
  iter->msgs_left = xmap_intersection_ext->n_out - 1;

  return (Xt_xmap_iter)iter;
}

static inline Xt_xmap_iter_intersection_ext
xmiei(void *iter)
{
  return (Xt_xmap_iter_intersection_ext)iter;
}

static int xmap_intersection_ext_iterator_next(Xt_xmap_iter iter) {

  Xt_xmap_iter_intersection_ext iter_intersection = xmiei(iter);

  if (iter_intersection == NULL || iter_intersection->msgs_left == 0)
    return 0;

  iter_intersection->msg++;
  iter_intersection->msgs_left--;

  return 1;
}

static int xmap_intersection_ext_iterator_get_rank(Xt_xmap_iter iter) {

  assert(iter != NULL);
  return xmiei(iter)->msg->rank;
}

static int const *
xmap_intersection_ext_iterator_get_transfer_pos(Xt_xmap_iter iter) {

  assert(iter != NULL);
  struct exchange_ext *restrict msg = xmiei(iter)->msg;
  if ((!msg->num_transfer_pos) | (msg->transfer_pos != NULL)) { } else {
    size_t num_transfer_pos = (size_t)msg->num_transfer_pos;
    int *restrict transfer_pos = msg->transfer_pos
      = xmalloc(num_transfer_pos * sizeof (*msg->transfer_pos));
    size_t num_transfer_pos_ext = (size_t)msg->num_transfer_pos_ext;
    struct Xt_pos_ext *restrict transfer_pos_ext = msg->transfer_pos_ext;
    size_t ofs = 0;
    for (size_t i = 0; i < num_transfer_pos_ext; ++i) {
      int abssize = abs(transfer_pos_ext[i].size);
      int step = isign(transfer_pos_ext[i].size);
      for (int j = 0; j < abssize; ++j)
        transfer_pos[ofs + (size_t)j] = transfer_pos_ext[i].start + j * step;
      ofs += (size_t)abssize;
    }
    assert(ofs == num_transfer_pos);
  }
  return msg->transfer_pos;
}

static int
xmap_intersection_ext_iterator_get_num_transfer_pos(Xt_xmap_iter iter) {
  assert(iter != NULL);
  return xmiei(iter)->msg->num_transfer_pos;
}

static const struct Xt_pos_ext *
xmap_intersection_ext_iterator_get_transfer_pos_ext(Xt_xmap_iter iter) {
  assert(iter != NULL);
  return xmiei(iter)->msg->transfer_pos_ext;
}

static int
xmap_intersection_ext_iterator_get_num_transfer_pos_ext(Xt_xmap_iter iter) {
  assert(iter != NULL);
  return xmiei(iter)->msg->num_transfer_pos_ext;
}

static void xmap_intersection_ext_iterator_delete(Xt_xmap_iter iter) {

  free(iter);
}
