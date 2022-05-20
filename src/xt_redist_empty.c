/**
 * @file xt_redist_empty.c
 *
 * @copyright Copyright  (C)  2022 Jörg Behrens <behrens@dkrz.de>
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
#include <mpi.h>

#include "xt_mpi_internal.h"
#include "xt/xt_redist_empty.h"
#include "xt_redist_internal.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"
#include "xt_config_internal.h"

static Xt_redist
redist_empty_copy(Xt_redist redist);

static void
redist_empty_delete(Xt_redist redist);

static void
redist_empty_s_exchange(Xt_redist redist, int num_arrays,
                      const void **src_data, void **dst_data);

static void
redist_empty_a_exchange(Xt_redist redist, int num_arrays,
                        const void **src_data, void **dst_data,
                        Xt_request *request);

static void
redist_empty_s_exchange1(Xt_redist redist, const void *src_data, void *dst_data);

static void
redist_empty_a_exchange1(Xt_redist redist, const void *src_data, void *dst_data,
                       Xt_request *request);

static int redist_empty_get_num_msg(Xt_redist redist,
                                  enum xt_msg_direction direction);

static MPI_Datatype
redist_empty_get_MPI_Datatype(Xt_redist redist, int rank,
                            enum xt_msg_direction direction,
                            bool do_dup);

static int
redist_empty_get_msg_ranks(Xt_redist redist,
                         enum xt_msg_direction direction,
                         int *restrict *ranks);

static MPI_Comm
redist_empty_get_MPI_Comm(Xt_redist redist);

static const struct xt_redist_vtable redist_empty_vtable = {
  .copy                  = redist_empty_copy,
  .delete                = redist_empty_delete,
  .s_exchange            = redist_empty_s_exchange,
  .a_exchange            = redist_empty_a_exchange,
  .s_exchange1           = redist_empty_s_exchange1,
  .a_exchange1           = redist_empty_a_exchange1,
  .get_num_msg           = redist_empty_get_num_msg,
  .get_msg_MPI_Datatype  = redist_empty_get_MPI_Datatype,
  .get_msg_ranks         = redist_empty_get_msg_ranks,
  .get_MPI_Comm          = redist_empty_get_MPI_Comm
};

typedef struct Xt_redist_empty_ *Xt_redist_empty;

struct Xt_redist_empty_ {

  const struct xt_redist_vtable *vtable;

  MPI_Comm comm;
  int tag_offset;
};

Xt_redist xt_redist_empty_new(MPI_Comm comm) {

  // ensure that yaxt is initialized
  assert(xt_initialized());

  Xt_redist_empty redist = xmalloc(sizeof (*redist));

  redist->vtable = &redist_empty_vtable;
  redist->comm = xt_mpi_comm_smart_dup(comm, &redist->tag_offset);

  return (Xt_redist)redist;
}

static inline Xt_redist_empty xrempty(void *redist) {
  return (Xt_redist_empty)redist;
}

static Xt_redist redist_empty_copy(Xt_redist redist) {

  Xt_redist_empty redist_empty = xrempty(redist);
  return xt_redist_empty_new(redist_empty->comm);
}

static void redist_empty_delete(Xt_redist redist) {

  Xt_redist_empty redist_empty = xrempty(redist);

  xt_mpi_comm_smart_dedup(&redist_empty->comm, redist_empty->tag_offset);

  free(redist_empty);
}

static void redist_empty_s_exchange(Xt_redist redist, int num_arrays,
                                    const void **XT_UNUSED(src_data),
                                    void **XT_UNUSED(dst_data)) {

  if (num_arrays != 0)
    Xt_abort(xrempty(redist)->comm, "ERROR: num_arrays != 0 is not"
             " implemented for this xt_redist type (Xt_redist_empty)",
             __FILE__, __LINE__);
}

static void redist_empty_a_exchange(Xt_redist redist, int num_arrays,
                                    const void **XT_UNUSED(src_data),
                                    void **XT_UNUSED(dst_data),
                                    Xt_request *XT_UNUSED(request)) {

  if (num_arrays != 0)
    Xt_abort(xrempty(redist)->comm, "ERROR: ERROR: num_arrays != 0 is not"
             " implemented for this xt_redist type "
             "(Xt_redist_empty)", __FILE__, __LINE__);
}

static void redist_empty_s_exchange1(Xt_redist XT_UNUSED(redist),
                                     const void *XT_UNUSED(src_data),
                                     void *XT_UNUSED(dst_data)) {

  return;
}

static void redist_empty_a_exchange1(Xt_redist XT_UNUSED(redist),
                                     const void *XT_UNUSED(src_data),
                                     void *XT_UNUSED(dst_data),
                                     Xt_request *XT_UNUSED(request)) {

  return;
}

static int redist_empty_get_num_msg(
  Xt_redist XT_UNUSED(redist), enum xt_msg_direction XT_UNUSED(direction)) {

  return 0;
}

static MPI_Datatype redist_empty_get_MPI_Datatype(
  Xt_redist XT_UNUSED(redist), int XT_UNUSED(rank),
  enum xt_msg_direction XT_UNUSED(direction), bool XT_UNUSED(do_dup)) {

  return MPI_DATATYPE_NULL;
}

static int redist_empty_get_msg_ranks(
  Xt_redist XT_UNUSED(redist), enum xt_msg_direction XT_UNUSED(direction),
  int *restrict *XT_UNUSED(ranks)) {

  return 0;
}

static MPI_Comm redist_empty_get_MPI_Comm(Xt_redist redist) {

  return xrempty(redist)->comm;
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
