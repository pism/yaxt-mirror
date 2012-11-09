/**
 * @file yaxt_f2c.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <mpi.h>

#include "core/core.h"

#include "xt/xt_idxlist.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_xmap.h"
#include "xt/xt_xmap_all2all.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_redist_p2p.h"

struct xt_idxlist_f {
   Xt_idxlist cptr;
};

struct xt_xmap_f {
   struct Xt_xmap *cptr;
};

struct xt_redist_f {
  struct xt_redist *cptr;
};

void xt_initialize_f(MPI_Fint *comm_f) {
  MPI_Comm comm_c;
  comm_c = MPI_Comm_f2c(*comm_f);
  xt_initialize(comm_c);
}

void xt_abort_cmsl_f(MPI_Fint *comm_f, const char *msg, const char *source, int line) {
  MPI_Comm comm_c = MPI_Comm_f2c(*comm_f);
  Xt_abort(comm_c, msg, source, line);
}

void xt_abort_msl_f(const char *msg, const char *source, int line) {
  MPI_Comm comm_c = SymPrefix(default_comm);
  Xt_abort(comm_c, msg, source, line);
}

struct xt_idxlist_f xt_idxvec_new_f(Xt_idx *idxvec, Xt_count num_indices) {
  struct xt_idxlist_f idxlist_f;
  idxlist_f.cptr = xt_idxvec_new(idxvec, num_indices);
  return idxlist_f;
}

void xt_idxlist_delete_f(struct xt_idxlist_f *idxlist_f) {
  xt_idxlist_delete(idxlist_f->cptr);
}

struct xt_xmap_f xt_xmap_all2all_new_f(struct xt_idxlist_f *src_idxlist_f, struct xt_idxlist_f *dst_idxlist_f, MPI_Fint *comm_f) {
  MPI_Comm comm_c = MPI_Comm_f2c(*comm_f);
  struct xt_xmap_f xmap_f;
  xmap_f.cptr = xt_xmap_all2all_new(src_idxlist_f->cptr, dst_idxlist_f->cptr, comm_c);
  return xmap_f;
}

void xt_xmap_delete_f(struct xt_xmap_f *xmap_f) {
  xt_xmap_delete(xmap_f->cptr);
}

struct xt_redist_f xt_redist_p2p_blocks_off_new_f(struct xt_xmap_f *xmap_f,
                                              int *src_block_offsets, int *src_block_sizes, int src_block_num,
                                              int *dst_block_offsets, int *dst_block_sizes, int dst_block_num,
                                              MPI_Fint *datatype_f) {
  MPI_Datatype datatype_c = MPI_Type_f2c(*datatype_f);
  struct xt_redist_f redist_f;

  redist_f.cptr =  xt_redist_p2p_blocks_off_new(xmap_f->cptr,
                                                src_block_offsets, src_block_sizes, src_block_num,
                                                dst_block_offsets, dst_block_sizes, dst_block_num,
                                                datatype_c);

  return redist_f;
}

struct xt_redist_f xt_redist_p2p_blocks_new_f(struct xt_xmap_f *xmap_f,
                                              int *src_block_sizes, int src_block_num,
                                              int *dst_block_sizes, int dst_block_num,
                                              MPI_Fint *datatype_f) {
  MPI_Datatype datatype_c = MPI_Type_f2c(*datatype_f);
  struct xt_redist_f redist_f;

  redist_f.cptr =  xt_redist_p2p_blocks_new(xmap_f->cptr,
                                            src_block_sizes, src_block_num,
                                            dst_block_sizes, dst_block_num,
                                            datatype_c);

  return redist_f;
}

struct xt_redist_f
xt_redist_p2p_off_new_f(struct xt_xmap_f *xmap_f,
                        MPI_Fint *src_offsets, MPI_Fint *dst_offsets,
                        MPI_Fint *datatype_f) {
  MPI_Datatype datatype_c = MPI_Type_f2c(*datatype_f);
  struct xt_redist_f redist_f;
  assert(sizeof (MPI_Fint) == sizeof (int));

  redist_f.cptr =  xt_redist_p2p_off_new(xmap_f->cptr, src_offsets, dst_offsets, datatype_c);

  return redist_f;
}

struct xt_redist_f xt_redist_p2p_new_f(struct xt_xmap_f *xmap_f, MPI_Fint *datatype_f) {
  MPI_Datatype datatype_c = MPI_Type_f2c(*datatype_f);
  struct xt_redist_f redist_f;

  redist_f.cptr =  xt_redist_p2p_new(xmap_f->cptr, datatype_c);

  return redist_f;
}

void xt_redist_delete_f(struct xt_redist_f *redist_f) {
  xt_redist_delete(redist_f->cptr);
}

void xt_redist_s_exchange1_f(struct xt_redist_f *redist_f, void **src_data, void **dst_data) {
  xt_redist_s_exchange(redist_f->cptr, src_data, 1, dst_data, 1);
}

struct xt_idxlist_f xt_idxstripes_new_f(struct Xt_stripe const *stripes,  Xt_count num_stripes) {
  struct xt_idxlist_f idxlist_f;
  idxlist_f.cptr = xt_idxstripes_new(stripes, num_stripes);
  return idxlist_f;
}
