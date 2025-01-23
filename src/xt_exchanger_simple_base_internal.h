/**
 * @file xt_exchanger_simple_base_internal.h
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
#ifndef XT_EXCHANGER_SIMPLE_BASE_INTERNAL_H
#define XT_EXCHANGER_SIMPLE_BASE_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "core/ppm_visibility.h"
#include "xt/xt_core.h"
#include "xt/xt_config.h"
#include "xt_redist_internal.h"
#include "xt_exchanger.h"

#include "xt_exchanger_simple_base.h"

typedef struct Xt_exchanger_simple_base_ *Xt_exchanger_simple_base;

PPM_DSO_INTERNAL Xt_exchanger_simple_base
xt_exchanger_simple_base_alloc(size_t nmsg, size_t extra);

PPM_DSO_INTERNAL void
xt_exchanger_simple_base_init(
  Xt_exchanger_simple_base exchanger,
  int nsend, int nrecv,
  const struct Xt_msg_param send_msgs[nsend],
  const struct Xt_msg_param recv_msgs[nrecv],
  MPI_Comm comm, int tag_offset,
  xt_simple_s_exchange_func s_func,
  xt_simple_a_exchange_func a_func,
  xt_simple_create_thread_share_func create_thread_share_func,
  xt_simple_destroy_thread_share_func destroy_thread_share_func,
  Xt_config config);

PPM_DSO_INTERNAL extern const struct xt_exchanger_vtable
xt_exchanger_simple_base_vtable;


struct Xt_exchanger_simple_base_ {

  const struct xt_exchanger_vtable * vtable;

  int nmsg[2];
  int config_flags;
  int tag_offset;
  MPI_Comm comm;
  xt_simple_s_exchange_func s_func;
  xt_simple_a_exchange_func a_func;
  xt_simple_create_thread_share_func create_thread_share_func;
  xt_simple_destroy_thread_share_func destroy_thread_share_func;
  struct Xt_msg_param msgs[];
};

#endif // XT_EXCHANGER_SIMPLE_BASE_INTERNAL_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
