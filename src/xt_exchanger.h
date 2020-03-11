/**
 * @file xt_exchanger.h
 * @brief exchanging of data based on information provided by redist's
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
#ifndef XT_EXCHANGER_H
#define XT_EXCHANGER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "xt/xt_core.h"
#include "xt/xt_request.h"
#include "xt_redist_internal.h"

/** \example test_exchanger_parallel.c
 */

typedef struct Xt_exchanger_ *Xt_exchanger;

struct xt_exchanger_vtable {
  Xt_exchanger (*copy)(Xt_exchanger, MPI_Comm, int);
  void (*delete)(Xt_exchanger);
  void (*s_exchange)(Xt_exchanger, const void *, void *);
  void (*a_exchange)(Xt_exchanger, const void *, void *, Xt_request *request);
  int (*get_msg_ranks)(Xt_exchanger, enum xt_msg_direction, int *restrict *);
  MPI_Datatype (*get_MPI_Datatype)(Xt_exchanger, int, enum xt_msg_direction);
};

struct Xt_exchanger_ {
  struct xt_exchanger_vtable * vtable;
};

void xt_exchanger_s_exchange(Xt_exchanger exchanger, const void * src_data,
                             void * dst_data);

void xt_exchanger_a_exchange(Xt_exchanger exchanger,
                             const void * src_data, void * dst_data,
                             Xt_request *request);

Xt_exchanger
xt_exchanger_copy(Xt_exchanger orig, MPI_Comm new_comm, int new_tag_offset);

void xt_exchanger_delete(Xt_exchanger);

void xt_exchanger_internal_optimize(size_t n, void * msgs, size_t msg_type_size,
                                    MPI_Comm comm);

MPI_Datatype
xt_exchanger_get_MPI_Datatype(Xt_exchanger exchanger,
                              int rank, enum xt_msg_direction direction);

int
xt_exchanger_get_msg_ranks(Xt_exchanger exchanger,
                           enum xt_msg_direction direction,
                           int *restrict *ranks);

typedef Xt_exchanger
(*Xt_exchanger_new)(int nsend, int nrecv,
                    const struct Xt_redist_msg *send_msgs,
                    const struct Xt_redist_msg *recv_msgs,
                    MPI_Comm comm, int tag_offset);

#endif // XT_EXCHANGER_H
/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
