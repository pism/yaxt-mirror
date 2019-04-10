/**
 * @file yaxt_f2c.c
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
#include "config.h"
#endif

#include <assert.h>
#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#if defined __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreserved-id-macro"
#pragma GCC diagnostic ignored "-Wstrict-prototypes"
#endif
#include "cfortran.h"
#if defined __clang__
#pragma GCC diagnostic pop
#endif

#include "xt/xt_mpi.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_idxmod.h"
#include "xt/xt_xmap.h"
#include "xt/xt_xmap_intersection.h"
#include "xt/xt_xmap_all2all.h"
#include "xt/xt_xmap_dist_dir.h"
#include "xt/xt_xmap_dist_dir_intercomm.h"
#include "xt/xt_xmap_dist_radix.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_dmap.h"
#include "xt/xt_dmap_offset.h"
#include "xt/xt_redist.h"
#include "xt/xt_redist_p2p.h"
#include "xt/xt_redist_single_array_base.h"
#include "xt/xt_redist_collection_static.h"
#include "xt/xt_redist_collection.h"
#include "xt/xt_sort.h"

struct xt_idxlist_f {
   Xt_idxlist cptr;
};

struct xt_xmap_f {
   Xt_xmap cptr;
};

struct xt_dmap_f {
  Xt_dmap cptr;
};

struct xt_redist_f {
  Xt_redist cptr;
};

struct xt_redist_msg_f {
  MPI_Fint rank, datatype;
};

struct xt_request_f {
  Xt_request cptr;
};

/* These functions are meant to be called from Fortran and don't need
 * an interface declaration in a C header  */
#if (defined __GNUC__ && __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5))\
  || (defined __clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-prototypes"
#endif

void xt_initialize_f(MPI_Fint *comm_f) {
  MPI_Comm comm_c;
  comm_c = MPI_Comm_f2c(*comm_f);
  xt_initialize(comm_c);
}

MPI_Fint xt_get_default_comm_f(void)
{
  return MPI_Comm_c2f(Xt_default_comm);
}

void xt_abort_f(MPI_Fint comm_f, const char *msg, const char *source,
                int line) __attribute__((noreturn));
void xt_abort_f(MPI_Fint comm_f, const char *msg, const char *source,
                int line) {
  MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
  Xt_abort(comm_c, msg, source, line);
}

void xt_sort_int_f2c(int *a, size_t n)
{
  xt_sort_int(a, n);
}

void xt_sort_index_f2c(Xt_int *a, int n, int *idx, int reset_index)
{
  xt_sort_index(a, n, idx, reset_index);
}

void xt_sort_idxpos_f2c(idxpos_type *a, size_t n)
{
  xt_sort_idxpos(a, n);
}

Xt_idxlist xt_idxlist_f2c(struct xt_idxlist_f *p)
{
  return p->cptr;
}

Xt_dmap xt_dmap_f2c(struct xt_dmap_f *p)
{
  return p->cptr;
}

Xt_redist xt_redist_f2c(struct xt_redist_f *p)
{
  return p->cptr;
}

Xt_request xt_request_f2c(struct xt_request_f *p)
{
  return p->cptr;
}

Xt_xmap xt_xmap_f2c(struct xt_xmap_f *p)
{
  return p->cptr;
}


MPI_Fint xt_idxlist_get_pack_size_f2c(struct xt_idxlist_f *idxlist,
                                      MPI_Fint comm_f)
{
  MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
  size_t size = xt_idxlist_get_pack_size(idxlist->cptr, comm_c);
  if (size > XT_MPI_FINT_MAX)
    Xt_abort(comm_c, "pack size too large", __FILE__, __LINE__);
  return (MPI_Fint)size;
}

static void
xt_idxlist_pack_f2c(struct xt_idxlist_f *idxlist, void *buffer,
                    MPI_Fint buffer_size, MPI_Fint *position, MPI_Fint comm_f)
{
  MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
  xt_idxlist_pack(idxlist->cptr, buffer, (int)buffer_size, position,
                  comm_c);
}

FCALLSCSUB5(xt_idxlist_pack_f2c, XT_IDXLIST_PACK_F, xt_idxlist_pack_f,
            PVOID, PVOID, INT, PINT, INT)

static void
xt_idxlist_unpack_f2c(struct xt_idxlist_f *idxlist, void *buffer,
                      MPI_Fint buffer_size, MPI_Fint *position,
                      MPI_Fint comm_f)
{
  MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
  idxlist->cptr = xt_idxlist_unpack(buffer, (int)buffer_size, position, comm_c);
}

FCALLSCSUB5(xt_idxlist_unpack_f2c, XT_IDXLIST_UNPACK_F, xt_idxlist_unpack_f,
            PVOID, PVOID, INT, PINT, INT)

Xt_xmap xt_xmap_all2all_new_f(struct xt_idxlist_f *src_idxlist_f,
                              struct xt_idxlist_f *dst_idxlist_f,
                              MPI_Fint comm_f) {
  MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
  return xt_xmap_all2all_new(src_idxlist_f->cptr, dst_idxlist_f->cptr, comm_c);
}

Xt_xmap xt_xmap_dist_dir_new_f(struct xt_idxlist_f *src_idxlist_f,
                               struct xt_idxlist_f *dst_idxlist_f,
                               MPI_Fint comm_f) {
  MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
  return xt_xmap_dist_dir_new(src_idxlist_f->cptr, dst_idxlist_f->cptr, comm_c);
}

Xt_xmap
xt_xmap_dist_dir_intercomm_new_f(struct xt_idxlist_f *src_idxlist_f,
                                 struct xt_idxlist_f *dst_idxlist_f,
                                 MPI_Fint inter_comm_f, MPI_Fint intra_comm_f)
{
  MPI_Comm inter_comm_c = MPI_Comm_f2c(inter_comm_f),
    intra_comm_c = MPI_Comm_f2c(intra_comm_f);
  return xt_xmap_dist_dir_intercomm_new(src_idxlist_f->cptr,
                                        dst_idxlist_f->cptr,
                                        inter_comm_c, intra_comm_c);
}

Xt_xmap
xt_xmap_dist_radix_new_f(void * src_ids, int num_src_ids,
                         void * dst_ids, int num_dst_ids,
                         int id_size, MPI_Fint comm_f)
{
  MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
  return
    xt_xmap_dist_radix_new(
      src_ids, num_src_ids, dst_ids, num_dst_ids, (size_t)id_size, comm_c);
}

void xt_dmap_reorder_f(struct xt_dmap_f *dmap_f, int reorder_type) {

  xt_dmap_reorder(dmap_f->cptr, (enum xt_reorder_type)reorder_type);
}

void xt_dmap_repeat_f(struct xt_dmap_f *dmap_f, int num_repetitions,
                      MPI_Fint *src_displacements, MPI_Fint *dst_displacements)
{
  assert(sizeof (MPI_Fint) == sizeof (int));

  xt_dmap_repeat(
    dmap_f->cptr, num_repetitions, src_displacements, dst_displacements);
}

void * xt_dmap_generate_redist_f(struct xt_dmap_f *dmap_f, MPI_Fint datatype_f)
{
  return xt_dmap_generate_redist(dmap_f->cptr, MPI_Type_f2c(datatype_f));
}

Xt_dmap xt_dmap_offset_new_f(struct xt_xmap_f *xmap_f, MPI_Fint * src_offsets,
                             MPI_Fint * dst_offsets)
{
  assert(sizeof (MPI_Fint) == sizeof (int));

  return xt_dmap_offset_new(xmap_f->cptr, src_offsets, dst_offsets);
}

Xt_redist xt_redist_p2p_blocks_off_new_f(struct xt_xmap_f *xmap_f,
                                         int *src_block_offsets,
                                         int *src_block_sizes,
                                         int src_block_num,
                                         int *dst_block_offsets,
                                         int *dst_block_sizes,
                                         int dst_block_num,
                                         MPI_Fint datatype_f) {
  MPI_Datatype datatype_c = MPI_Type_f2c(datatype_f);

  return xt_redist_p2p_blocks_off_new(xmap_f->cptr,
                                      src_block_offsets, src_block_sizes,
                                      src_block_num,
                                      dst_block_offsets, dst_block_sizes,
                                      dst_block_num,
                                      datatype_c);
}

Xt_redist xt_redist_p2p_blocks_new_f(struct xt_xmap_f *xmap_f,
                                     int *src_block_sizes, int src_block_num,
                                     int *dst_block_sizes, int dst_block_num,
                                     MPI_Fint datatype_f) {
  MPI_Datatype datatype_c = MPI_Type_f2c(datatype_f);

  return xt_redist_p2p_blocks_new(xmap_f->cptr,
                                  src_block_sizes, src_block_num,
                                  dst_block_sizes, dst_block_num,
                                  datatype_c);
}

Xt_redist
xt_redist_p2p_off_new_f(struct xt_xmap_f *xmap_f,
                        MPI_Fint *src_offsets, MPI_Fint *dst_offsets,
                        MPI_Fint datatype_f) {
  MPI_Datatype datatype_c = MPI_Type_f2c(datatype_f);
  assert(sizeof (MPI_Fint) == sizeof (int));

  return xt_redist_p2p_off_new(xmap_f->cptr, src_offsets, dst_offsets,
                               datatype_c);
}

Xt_redist
xt_redist_p2p_new_f(struct xt_xmap_f *xmap_f, MPI_Fint datatype_f) {
  MPI_Datatype datatype_c = MPI_Type_f2c(datatype_f);

  return xt_redist_p2p_new(xmap_f->cptr, datatype_c);
}

Xt_redist
xt_redist_collection_static_new_f(Xt_redist *redists, MPI_Fint num_redists,
                                  MPI_Aint *src_displacements,
                                  MPI_Aint *dst_displacements,
                                  MPI_Fint comm_f)
{
#if XT_MPI_FINT_MAX != INT_MAX
  assert((long long)num_redists <= (long long)INT_MAX);
#endif
  MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
  if (num_redists < 1)
    Xt_abort(comm_c, "bad case: (num_redists < 1)", __FILE__, __LINE__);

  return xt_redist_collection_static_new(redists, (int)num_redists,
                                         src_displacements,
                                         dst_displacements,
                                         comm_c);
}

Xt_redist
xt_redist_collection_new_f(Xt_redist *redists, MPI_Fint num_redists,
                           MPI_Fint cache_size, MPI_Fint comm_f)
{
#if XT_MPI_FINT_MAX != INT_MAX
  assert((long long)num_redists <= (long long)INT_MAX
         && (long long)cache_size <= (long long)INT_MAX);
#endif
  MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
  if (num_redists < 1)
    Xt_abort(comm_c, "bad case: (num_redists < 1)", __FILE__, __LINE__);
  return xt_redist_collection_new(redists, (int)num_redists,
                                  (int)cache_size, comm_c);
}

static void
xt_slice_c_loc_f2c(void *a, void *p)
{
  *(void **)p = a;
}

FCALLSCSUB2(xt_slice_c_loc_f2c,XT_SLICE_C_LOC,xt_slice_c_loc,PVOID,PVOID)

MPI_Fint
xt_redist_get_mpi_comm_c2f(Xt_redist *redist)
{
  MPI_Comm comm = xt_redist_get_MPI_Comm(*redist);
  return MPI_Comm_c2f(comm);
}

void *
xt_redist_p2p_ext_new_c2f(Xt_xmap *xmap,
                          int num_src_ext, struct Xt_offset_ext src_extents[],
                          int num_dst_ext, struct Xt_offset_ext dst_extents[],
                          MPI_Fint datatype_f)
{
  return xt_redist_p2p_ext_new(*xmap, num_src_ext, src_extents,
                               num_dst_ext, dst_extents,
                               MPI_Type_f2c(datatype_f));
}

void *
xt_redist_single_array_base_new_c2f(int nsend, int nrecv,
                                    const struct xt_redist_msg_f *send_msgs_f,
                                    const struct xt_redist_msg_f *recv_msgs_f,
                                    MPI_Fint comm_f)
{
  assert(nsend >= 0 && nrecv >= 0);
  MPI_Comm comm_c = MPI_Comm_f2c(comm_f);
  struct Xt_redist_msg *send_msgs_c
    = xmalloc((size_t)nsend * sizeof(*send_msgs_c));
  struct Xt_redist_msg *recv_msgs_c
    = xmalloc((size_t)nrecv * sizeof(*recv_msgs_c));
  for (int i = 0; i < nsend; ++i) {
    send_msgs_c[i].rank = send_msgs_f[i].rank;
    send_msgs_c[i].datatype = MPI_Type_f2c(send_msgs_f[i].datatype);
  }
  for (int i = 0; i < nrecv; ++i) {
    recv_msgs_c[i].rank = recv_msgs_f[i].rank;
    recv_msgs_c[i].datatype = MPI_Type_f2c(recv_msgs_f[i].datatype);
  }
  void *redist = xt_redist_single_array_base_new(nsend, nrecv, send_msgs_c,
                                                 recv_msgs_c, comm_c);
  free(recv_msgs_c);
  free(send_msgs_c);
  return redist;
}

void *
xt_xmap_intersection_new_f2c(
  int num_src_intersections,
  const struct Xt_com_list src_com[num_src_intersections],
  int num_dst_intersections,
  const struct Xt_com_list dst_com[num_dst_intersections],
  void *src_idxlist, void *dst_idxlist, MPI_Fint comm)
{
  return xt_xmap_intersection_new(num_src_intersections, src_com,
                                  num_dst_intersections, dst_com,
                                  src_idxlist, dst_idxlist, MPI_Comm_f2c(comm));
}

void *
xt_xmap_intersection_ext_new_f2c(
  int num_src_intersections, const void *src_com,
  int num_dst_intersections, const void *dst_com,
  void *src_idxlist, void *dst_idxlist, MPI_Fint comm)
{
  return xt_xmap_intersection_ext_new(
    num_src_intersections, src_com, num_dst_intersections, dst_com,
    src_idxlist, dst_idxlist, MPI_Comm_f2c(comm));
}

#ifndef HAVE_FC_IS_CONTIGUOUS
int
xt_com_list_contiguous(const struct Xt_com_list *p_com_a,
                       const struct Xt_com_list *p_com_b)
{
  return (p_com_a + 1 == p_com_b);
}

int
xt_redist_msg_contiguous(const struct xt_redist_msg_f *p_msgs_a,
                         const struct xt_redist_msg_f *p_msgs_b)
{
  return (p_msgs_a + 1 == p_msgs_b);
}
#endif

void
xt_mpi_comm_mark_exclusive_f2c(MPI_Fint *comm) {
  xt_mpi_comm_mark_exclusive(MPI_Comm_f2c(*comm));
}
#if (defined __GNUC__ && __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5))\
  || (defined __clang__)
#pragma GCC diagnostic pop
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
