/**
 * @file xt_xmap_dist_dir_bucket_gen.h
 *
 * @copyright Copyright  (C)  2024 Jörg Behrens <behrens@dkrz.de>
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
 *
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

#ifndef XT_XMAP_DIST_DIR_BUCKET_GEN_H
#define XT_XMAP_DIST_DIR_BUCKET_GEN_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "xt/xt_core.h"
#include "xt/xt_config.h"
#include "xt/xt_idxlist.h"

enum Xt_dist_dir_bucket_gen_types {
  Xt_dist_dir_bucket_gen_type_send = 1,
  Xt_dist_dir_bucket_gen_type_recv = 2,
  Xt_dist_dir_bucket_gen_type_sendrecv = 3,
};


typedef int (*Xt_xmdd_bucket_gen_init_state)(
  void *gen_state,
  Xt_idxlist src_idxlist,
  Xt_idxlist dst_idxlist,
  MPI_Comm intra_comm, int tag_offset_intra,
  MPI_Comm inter_comm, int tag_offset_inter,
  Xt_config config);

typedef void (*Xt_xmdd_bucket_gen_destroy_state)(void *gen_state);

typedef int (*Xt_xmdd_bucket_gen_get_intersect_max_num)(
  void *gen_state, int type);

typedef struct Xt_com_list (*Xt_xmdd_bucket_gen_next)(
  void *gen_state, int type);

/**
 * Interface of bucket generator
 * Essentially, the generator needs to be able to enumerate all
 * buckets used to form intersections. This means the generator only
 * needs to produce buckets that actually can intersect and can skip
 * buckets that won't intersect the requested type of list.
 */
struct Xt_xmdd_bucket_gen_desc {
  /** The init function sets up the generator state. */
  Xt_xmdd_bucket_gen_init_state init;
  /** The destroy function clean up the generator state. Can be zero
   * if no cleaning is needed. */
  Xt_xmdd_bucket_gen_destroy_state destroy;
  Xt_xmdd_bucket_gen_get_intersect_max_num get_intersect_max_num;
  /** The next function returns the next bucket and corresponding rank
   * (ranks can be skipped when the intersection will be empty
   * anyway).
   * Any previously returned buckets become invalid. */
  Xt_xmdd_bucket_gen_next next;
  /** gen_state_size is the size of the generator state
   *
   * The distributed directory will provide memory aligned to
   * pointer variables of this size. */
  size_t gen_state_size;
};

#endif // XT_XMAP_DIST_DIR_BUCKET_GEN_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
