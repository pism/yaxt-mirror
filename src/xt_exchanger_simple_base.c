/**
 * @file xt_exchanger_simple_base.c
 *
 * @copyright Copyright  (C)  2014 Jörg Behrens <behrens@dkrz.de>
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

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_redist_internal.h"
#include "xt_exchanger.h"
#include "xt_exchanger_simple_base.h"

static void xt_exchanger_simple_base_delete(Xt_exchanger exchanger);
static void xt_exchanger_simple_base_s_exchange(Xt_exchanger exchanger,
                                                const void * src_data,
                                                void * dst_data);

static const struct xt_exchanger_vtable exchanger_simple_base_vtable = {
  .delete = xt_exchanger_simple_base_delete,
  .s_exchange = xt_exchanger_simple_base_s_exchange,
};

typedef struct Xt_exchanger_simple_base_ * Xt_exchanger_simple_base;

struct Xt_exchanger_simple_base_ {

  const struct xt_exchanger_vtable * vtable;

  int nsend, nrecv;
  struct Xt_redist_msg * send_msgs, * recv_msgs;
  MPI_Comm comm;

  xt_simple_exchange_func func;
};

Xt_exchanger
xt_exchanger_simple_base_new(int nsend, int nrecv,
                             struct Xt_redist_msg * send_msgs,
                             struct Xt_redist_msg * recv_msgs,
                             MPI_Comm comm,
                             xt_simple_exchange_func func) {

  if (func == NULL)
    Xt_abort(comm, "ERROR(xt_exchanger_simple_base_new): invalid exchange "
             "function pointer", __FILE__, __LINE__);

  Xt_exchanger_simple_base exchanger = xmalloc(1 * sizeof(*exchanger));

  exchanger->vtable = &exchanger_simple_base_vtable;

  xt_mpi_call(MPI_Comm_dup(comm, &(exchanger->comm)), comm);
  exchanger->nsend = nsend;
  exchanger->send_msgs = xt_redist_msgs_copy(nsend, send_msgs, exchanger->comm);
  exchanger->nrecv = nrecv;
  exchanger->recv_msgs = xt_redist_msgs_copy(nrecv, recv_msgs, exchanger->comm);
  exchanger->func = func;

  xt_exchanger_internal_optimize(nsend, exchanger->send_msgs,
                                 sizeof(*(exchanger->send_msgs)),
                                 exchanger->comm);

  xt_exchanger_internal_optimize(nrecv, exchanger->recv_msgs,
                                 sizeof(*(exchanger->recv_msgs)),
                                 exchanger->comm);

  return (Xt_exchanger)exchanger;
}

static void xt_exchanger_simple_base_delete(Xt_exchanger exchanger) {

  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;

  xt_redist_msgs_free(exchanger_sb->nsend, exchanger_sb->send_msgs,
                      exchanger_sb->comm);
  xt_redist_msgs_free(exchanger_sb->nrecv, exchanger_sb->recv_msgs,
                      exchanger_sb->comm);

  xt_mpi_call(MPI_Comm_free(&(exchanger_sb->comm)), MPI_COMM_WORLD);
  free(exchanger_sb);
}

static void xt_exchanger_simple_base_s_exchange(Xt_exchanger exchanger,
                                                const void * src_data,
                                                void * dst_data) {

  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;

  exchanger_sb->func(src_data, dst_data, exchanger_sb->nsend,
                     exchanger_sb->nrecv, exchanger_sb->send_msgs,
                     exchanger_sb->recv_msgs, exchanger_sb->comm);
}
