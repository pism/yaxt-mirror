/**
 * @file xt_exchanger.c
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

#include <string.h>
#include <stdlib.h>

#include "xt/xt_mpi.h"
#include "xt_exchanger.h"
#include "xt_exchanger_irecv_isend.h"
#include "xt_exchanger_mix_isend_irecv.h"
#include "xt_exchanger_irecv_isend_packed.h"
#include "xt_exchanger_irecv_isend_ddt_packed.h"
#include "xt_exchanger_neigh_alltoall.h"
#include "xt_exchanger_internal.h"

Xt_exchanger
xt_exchanger_new(int nsend, int nrecv,
                 const struct Xt_msg_param *send_msgs,
                 const struct Xt_msg_param *recv_msgs,
                 MPI_Comm comm, int tag_offset,
                 Xt_config config)
{
  Xt_exchanger_new exchanger_new
    = xt_config_get_exchange_new_by_comm(config, comm);
  Xt_exchanger exchanger = exchanger_new(nsend, nrecv, send_msgs, recv_msgs,
                                         comm, tag_offset, config);
  return exchanger;
}



Xt_exchanger
xt_exchanger_copy(Xt_exchanger exchanger, MPI_Comm new_comm, int new_tag_offset)
{
  return exchanger->vtable->copy(exchanger, new_comm, new_tag_offset);
}

void xt_exchanger_delete(Xt_exchanger exchanger) {

  exchanger->vtable->delete(exchanger);
}

void
xt_exchanger_s_exchange(Xt_exchanger exchanger,
                        const void *src_data, void *dst_data)
{
  exchanger->vtable->s_exchange(exchanger, src_data, dst_data);
}

void xt_exchanger_a_exchange(Xt_exchanger exchanger,
                             const void *src_data, void *dst_data,
                             Xt_request *request) {
  if (exchanger->vtable->a_exchange)
    exchanger->vtable->a_exchange(exchanger, src_data, dst_data, request);
  else {
    exchanger->vtable->s_exchange(exchanger, src_data, dst_data);
    *request = XT_REQUEST_NULL;
  }
}

int
xt_exchanger_get_num_msg(Xt_exchanger exchanger,
                         int direction)
{
  return exchanger->vtable->get_num_msg(exchanger, direction);
}

int
xt_exchanger_get_msg_ranks(Xt_exchanger exchanger,
                           enum xt_msg_direction direction,
                           int *restrict *ranks)
{
  return exchanger->vtable->get_msg_ranks(exchanger, direction, ranks);
}

MPI_Datatype
xt_exchanger_get_MPI_Datatype(Xt_exchanger exchanger, int rank,
                              enum xt_msg_direction direction, bool do_dup) {
  return exchanger->vtable->get_MPI_Datatype(exchanger, rank, direction, do_dup);
}

Xt_exchanger_thread_share
xt_exchanger_create_thread_share(Xt_exchanger exchanger, Xt_config config)
{
  if (exchanger->vtable->create_thread_share)
    return exchanger->vtable->create_thread_share(exchanger, config);
  else
    return NULL;
}

void
xt_exchanger_destroy_thread_share(Xt_exchanger exchanger, Xt_exchanger_thread_share thread_share)
{
  if (exchanger->vtable->destroy_thread_share)
    exchanger->vtable->destroy_thread_share(exchanger, thread_share);
}

void
xt_msg_params_copy(size_t n,
                   const struct Xt_msg_param *restrict src,
                   struct Xt_msg_param *restrict dst,
                   MPI_Comm comm, bool dt_dup) {
  for (size_t i = 0; i < n; ++i) {
    dst[i].rank = src[i].rank;
    dst[i].displacement = src[i].displacement;
    dst[i].count = src[i].count;
    if (dt_dup)
      xt_mpi_call(MPI_Type_dup(src[i].datatype, &dst[i].datatype), comm);
    else
      dst[i].datatype = src[i].datatype;
  }
}

void
xt_msg_params_destruct(size_t n, struct Xt_msg_param *msgs,
                       MPI_Comm comm) {

  for (size_t i = 0; i < n; ++i) {
    MPI_Datatype dt = msgs[i].datatype;
    if (dt != MPI_DATATYPE_NULL) {
      xt_mpi_call(MPI_Type_free(&msgs[i].datatype), comm);
      for (size_t j = i + 1; j < n; ++j)
        if (msgs[j].datatype == dt)
          msgs[j].datatype = MPI_DATATYPE_NULL;
    }
  }
}

Xt_exchanger_iter
xt_exchanger_get_iterator(Xt_exchanger exchanger, Xt_config config)
{
  return exchanger->vtable->get_iterator(exchanger, config);
}

void
xt_exchanger_iter_delete(Xt_exchanger_iter iter)
{
  iter->vtable->delete(iter);
}

int
xt_exchanger_iter_next(Xt_exchanger_iter iter)
{
  return iter->vtable->next(iter);
}

MPI_Datatype
xt_exchanger_iter_get_MPI_Datatype(Xt_exchanger_iter iter)
{
  return iter->vtable->get_MPI_Datatype(iter);
}

int
xt_exchanger_iter_get_rank(Xt_exchanger_iter iter)
{
  return iter->vtable->get_rank(iter);
}

int
xt_exchanger_iter_get_count(Xt_exchanger_iter iter)
{
  return iter->vtable->get_count(iter);
}

int
xt_exchanger_iter_get_direction(Xt_exchanger_iter iter)
{
  return iter->vtable->get_direction(iter);
}

MPI_Aint
xt_exchanger_iter_get_displacement(Xt_exchanger_iter iter)
{
  return iter->vtable->get_displacement(iter);
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
