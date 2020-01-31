/**
 * @file test_redist_common.h
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
#ifndef TEST_REDIST_COMMON_H
#define TEST_REDIST_COMMON_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>

#include <mpi.h>

#include "yaxt.h"

Xt_xmap
build_odd_selection_xmap(int src_num_indices);

int communicators_are_congruent(MPI_Comm comm1, MPI_Comm comm2);

typedef void (*prepare_dst)(void *dst, const void *dst_prep_info,
                            size_t dst_num_elems);

void
check_redist_(Xt_redist redist, const void *src,
              size_t dst_num_elems,
              void *dst,
              prepare_dst dst_prep,
              const void *dst_prep_info,
              const void *ref_dst_data,
              MPI_Datatype dst_data_dt, MPI_Datatype ref_dst_data_dt,
              const char *file, int line);

#define check_redist(redist, src, dst_size, dst, dst_prep, dst_prep_info, \
                     ref_dst_data, dst_data_dt, ref_dst_data_dt)        \
  check_redist_(redist, src, dst_size, dst, dst_prep, dst_prep_info,    \
                ref_dst_data, dst_data_dt, ref_dst_data_dt, __FILE__, __LINE__)

void
fill_array_double(void *dst, const void *dst_prep_info, size_t dst_num_elems);
void
fill_array_long(void *dst, const void *dst_prep_info, size_t dst_num_elems);
void
fill_array_int(void *dst, const void *dst_prep_info, size_t dst_num_elems);



void
wrap_a_exchange(Xt_redist redist, int num_data_p, const void *src_data_p[],
                void *dst_data_p[]);

typedef void
(*exchange_func_ptr)(Xt_redist redist, int num_data_p, const void *src_data_p[],
                     void *dst_data_p[]);

void
wrap_a_exchange1(Xt_redist redist, const void *src_data_p, void *dst_data_p);

typedef void
(*exchange1_func_ptr)(Xt_redist redist, const void *src_data_p,
                     void *dst_data_p);

void
check_wait_request_(Xt_request *request, const char *file, int line);

#define check_wait_request(request) \
  check_wait_request_(request, __FILE__, __LINE__)
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
