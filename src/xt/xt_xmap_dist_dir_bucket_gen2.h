/**
 * @file xt_xmap_dist_dir_bucket_gen2.h
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

#ifndef XT_XMAP_DIST_DIR_BUCKET_GEN2_H
#define XT_XMAP_DIST_DIR_BUCKET_GEN2_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stddef.h>

#include <mpi.h>

#include <xt/xt_xmap_dist_dir_bucket_gen.h>
#include <xt/xt_idxlist.h>
#include <xt/xt_config.h>

enum Xt_dist_dir_bucket_gen_types {
  Xt_dist_dir_bucket_gen_type_send = 1,
  Xt_dist_dir_bucket_gen_type_recv = 2,
  Xt_dist_dir_bucket_gen_type_sendrecv = 3,
};

struct Xt_xmdd_bucket_gen_comms {
  MPI_Comm intra_comm, inter_comm;
  int tag_offset_intra, tag_offset_inter;
};

/**
 * Initialize bucket generator state.
 * @return 0 if index lists should be used as-is, non-zero if
 * conversion to stripes appears profitable.
 */
typedef int (*Xt_xmdd_bucket_gen_init_state)(
  void *gen_state,
  Xt_idxlist src_idxlist,
  Xt_idxlist dst_idxlist,
  Xt_config config,
  const struct Xt_xmdd_bucket_gen_comms *comms,
  void *init_params);

typedef void (*Xt_xmdd_bucket_gen_destroy_state)(void *gen_state);

/**
 * @param type one of the values from enum
 * Xt_dist_dir_bucket_gen_types for a distributed directory
 * constructed on an intra-communicator, this will be always
 * Xt_dist_dir_bucket_gen_type_sendrecv, on an inter-communicator,
 * Xt_dist_dir_bucket_gen_type_send and
 * Xt_dist_dir_bucket_gen_type_recv will be used in succession because
 * the buckets differ for the local and remote group
 */
typedef int (*Xt_xmdd_bucket_gen_get_intersect_max_num)(
  void *gen_state, int type);

/**
 * @param type one of the values from enum
 * Xt_dist_dir_bucket_gen_types for a distributed directory
 * constructed on an intra-communicator, this will be always
 * Xt_dist_dir_bucket_gen_type_sendrecv, on an inter-communicator,
 * Xt_dist_dir_bucket_gen_type_send and
 * Xt_dist_dir_bucket_gen_type_recv will be used in succession because
 * the buckets differ for the local and remote group
 */
typedef struct Xt_com_list (*Xt_xmdd_bucket_gen_next)(
  void *gen_state, int type);


/**
 * Define interface of bucket generator
 *
 * Essentially, the generator needs to be able to enumerate all
 * buckets used to form intersections. Conversely, this also means the
 * generator only needs to produce buckets that actually can intersect
 * and is permitted to skip buckets that won't intersect the requested
 * type of list.
 *
 * @param[in,out] gen generator interface object
 * @param[in] init This function is called to set up the generator
 *                 state.
 * @param[in] destroy The destroy function cleans up the generator
 *                 state. Can be zero, if no cleaning is needed.
 *
 * @param[in] next The next function returns the next bucket and
 * corresponding rank (ranks can be skipped when the intersection will
 * be empty anyway).  Any previously returned buckets become invalid.
 * @param[in] get_intersect_max_num This function returns, for a given
 * state the maximal number of buckets that will be generated
 *
 * @param gen_state_size number of bytes to allocate for each generator state
 * @param init_params global parameters passed to each invocation of
 * the \a init function
 * @note it remains the responsibility of the caller to dispose of
 * resources held by init_params
 * @see xt_xmdd_bucket_gen_set_extended_init_params
 */
void
xt_xmdd_bucket_gen_define_interface(
  Xt_xmdd_bucket_gen gen,
  Xt_xmdd_bucket_gen_init_state init,
  Xt_xmdd_bucket_gen_destroy_state destroy,
  Xt_xmdd_bucket_gen_get_intersect_max_num get_intersect_max_num,
  Xt_xmdd_bucket_gen_next next,
  size_t gen_state_size,
  void *init_params);

/**
 * Inquire state-independent initialization parameters of bucket
 * generator.
 * @param gen bucket generator to inquire for
 * @return pointer to initialization parameters previously defined via
 * \a xt_xmddd_bucket_gen_define_interface
 */
void *
xt_xmdd_bucket_gen_get_init_params(Xt_xmdd_bucket_gen gen);

/**
 * Change state-independent initialization parameters of bucket
 * generator.
 * @param gen bucket generator to update state parameters for
 * @param init_params global parameters passed to each invocation of
 * the \a init function
 * @note previously set initial parameters will not be
 * destroyed automatically
 */
void
xt_xmdd_bucket_gen_set_init_params(Xt_xmdd_bucket_gen gen,
                                   void *init_params);

/**
 * Inquire initialization function of bucket generator.
 *
 * @param gen bucket generator to inquire for
 * @return pointer to initialization function previously defined via
 * \a xt_xmddd_bucket_gen_define_interface
 */
Xt_xmdd_bucket_gen_init_state
xt_xmdd_bucket_gen_get_init(Xt_xmdd_bucket_gen gen);

/**
 * Change state initialization function of bucket generator.
 *
 * @param gen bucket generator to update state initializer for
 * @param init pointer of initialization function
 * @note previously set initial parameters will not be
 * destroyed automatically
 */
void
xt_xmdd_bucket_gen_set_init(Xt_xmdd_bucket_gen gen,
                            Xt_xmdd_bucket_gen_init_state init);



typedef void (*Xt_xmdd_bucket_gen_iparm_destroy)(
  void *init_params, size_t init_params_size);

typedef void (*Xt_xmdd_bucket_gen_iparm_copy)(
  void *restrict dest_init_params,
  const void *restrict src_init_params,
  size_t init_params_size);

/**
 * Add handling of init_params resources to bucket generator.
 *
 * By supplying the following additional information, responsibility
 * for the object backing init_params can be delegated to the
 * bucket generator.
 *
 * @param gen bucket generator to enable resource handling for
 * @param init_params_size If non-zero, ownership of the object pointed to
 * by init_params passes to the bucket generator object it's assigned
 * to. The init_params must be dynamically allocated with malloc,
 * because when xt_xmdd_bucket_gen_delete is called, free will be
 * called for init_params.
 * @param destroy When destroy is non-zero, it will be called
 * during xt_xmdd_bucket_gen_delete (before free).
 * @param cpy Override for memcpy of init_params.
 */
void
xt_xmdd_bucket_gen_set_extended_init_params(
  Xt_xmdd_bucket_gen gen,
  size_t init_params_size,
  Xt_xmdd_bucket_gen_iparm_destroy destroy,
  Xt_xmdd_bucket_gen_iparm_copy cpy);


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
