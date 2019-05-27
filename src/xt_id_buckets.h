/**
 * @file xt_id_buckets.h
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
#ifndef XT_ID_BUCKETS_H
#define XT_ID_BUCKETS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

typedef struct Xt_id_buckets_ *Xt_id_buckets;
typedef struct Xt_id_buckets_iter_ *Xt_id_buckets_iter;

void xt_id_buckets_iter_get_next_ubound(Xt_id_buckets_iter iter, char * ubound);

void xt_id_buckets_tree_remove_leaves(
  char ** tree, size_t * tree_size, size_t * num_leaves, uint64_t * keep_flags,
  size_t * bucket_ids);
int xt_id_buckets_tree_split_leaves(
  char ** tree, size_t * tree_size, size_t * num_leaves, size_t id_size,
  size_t * split_flags);

Xt_id_buckets xt_id_buckets_new(
  void const * ids, size_t id_size, size_t num_ids, size_t num_reg_buckets,
  MPI_Comm comm);
void xt_id_buckets_delete(Xt_id_buckets id_buckets);

void xt_id_buckets_get_ubounts(Xt_id_buckets id_buckets, void * ubounds);

Xt_id_buckets_iter xt_id_buckets_get_iter(Xt_id_buckets id_buckets);
Xt_id_buckets_iter xt_id_buckets_iter_new(char const * tree);
int xt_id_buckets_iter_next(Xt_id_buckets_iter iter, void * ubound);
void xt_id_buckets_iter_delete(Xt_id_buckets_iter iter);

#endif

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
