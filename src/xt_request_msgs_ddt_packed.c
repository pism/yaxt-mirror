/**
 * @file xt_request_msgs_ddt_packed.c
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
#include <stdlib.h>
#include <string.h>

#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"

#include "xt/xt_request_msgs_ddt_packed.h"
#include "xt_mpi_internal.h"
#include "xt_request_internal.h"
#include "xt_ddt_internal.h"

static void xt_request_msgs_wait_ddt_packed(Xt_request request);
static int xt_request_msgs_test_ddt_packed(Xt_request request);

static const struct Xt_request_vtable request_msgs_ddt_packed_vtable = {

  .wait = xt_request_msgs_wait_ddt_packed,
  .test = xt_request_msgs_test_ddt_packed,
};

typedef struct Xt_request_msgs_ddt_packed_ * Xt_request_msgs_ddt_packed;

struct Xt_request_msgs_ddt_packed_ {

  const struct Xt_request_vtable * vtable;
  MPI_Request * requests;
  int * ops_completed_buffer;
  MPI_Comm comm;
  int n_requests;
  int n_packed, n_tmp_buffers;
  void * unpacked_data;
  Xt_ddt * ddts;
  enum xt_memtype packed_memtype;
  enum xt_memtype tmp_memtype;
  void *buffers[];
};

Xt_request xt_request_msgs_ddt_packed_new(int n_requests,
                                          const MPI_Request requests[n_requests],
                                          MPI_Comm comm, int n_packed,
                                          int n_tmp_buffers,
                                          Xt_ddt packed_ddts[n_packed],
                                          void *packed_data[n_packed],
                                          void *tmp_buffers[n_tmp_buffers],
                                          void * unpacked_data,
                                          enum xt_memtype packed_memtype,
                                          enum xt_memtype tmp_memtype) {

  assert(n_requests >= 0 && n_packed >= 0 && n_tmp_buffers >= 0);
  size_t hdr_size = sizeof(struct Xt_request_msgs_ddt_packed_),
    bufp_size
    = ((size_t)n_packed + (size_t)n_tmp_buffers) * sizeof(void *);
  Xt_request_msgs_ddt_packed request = xmalloc(hdr_size + bufp_size);

  request->vtable = &request_msgs_ddt_packed_vtable;
  request->n_requests = n_requests;
  request->requests = xmalloc((size_t)n_requests * sizeof(*request->requests) +
                              (size_t)n_requests *
                              sizeof(*request->ops_completed_buffer));
  request->ops_completed_buffer =
    (int *)(request->requests + (size_t)n_requests);
  memcpy(request->requests, requests,
         (size_t)n_requests * sizeof(*request->requests));
  request->comm = comm;
  request->n_packed = n_packed;
  request->n_tmp_buffers = n_tmp_buffers;
  request->ddts = xmalloc((size_t)n_packed * sizeof(*request->ddts));
  for (int i = 0; i < n_packed; ++i) {
    request->ddts[i] = packed_ddts[i];
    xt_ddt_inc_ref_count(request->ddts[i]);
  }
  request->packed_memtype = packed_memtype;
  request->tmp_memtype = tmp_memtype;
  memcpy(request->buffers, packed_data,
         (size_t)n_packed * sizeof(*request->buffers));
  memcpy(request->buffers + n_packed, tmp_buffers,
         (size_t)n_tmp_buffers * sizeof(*request->buffers));
  request->unpacked_data = unpacked_data;

  return (Xt_request)request;
}

static void unpack_data(Xt_request_msgs_ddt_packed request_msgs_ddt_packed) {

//! \TODO merge all unpacking kernels into single kernel call -> less overhead
  int n_packed = request_msgs_ddt_packed->n_packed;
  for (int i = 0; i < n_packed; ++i)
    xt_ddt_unpack_internal(
      request_msgs_ddt_packed->ddts[i],
      request_msgs_ddt_packed->buffers[i],
      request_msgs_ddt_packed->unpacked_data,
      request_msgs_ddt_packed->packed_memtype);
}

static void free_request(Xt_request_msgs_ddt_packed request_msgs_ddt_packed) {

  for (int i = 0; i < request_msgs_ddt_packed->n_packed; ++i)
    xt_gpu_free(
      request_msgs_ddt_packed->buffers[i],
      request_msgs_ddt_packed->packed_memtype);
  for (int i = request_msgs_ddt_packed->n_packed;
       i < request_msgs_ddt_packed->n_packed +
       request_msgs_ddt_packed->n_tmp_buffers; ++i)
    xt_gpu_free(
      request_msgs_ddt_packed->buffers[i],
      request_msgs_ddt_packed->tmp_memtype);
  for (int i = 0; i < request_msgs_ddt_packed->n_packed; ++i)
    xt_ddt_delete(request_msgs_ddt_packed->ddts[i]);
  free(request_msgs_ddt_packed->ddts);
  free(request_msgs_ddt_packed->requests);
  free(request_msgs_ddt_packed);
}

static void xt_request_msgs_wait_ddt_packed(Xt_request request) {

  Xt_request_msgs_ddt_packed request_msgs_ddt_packed =
    (Xt_request_msgs_ddt_packed)request;

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ == 11 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
  xt_mpi_call(MPI_Waitall(request_msgs_ddt_packed->n_requests,
                          request_msgs_ddt_packed->requests,
                          MPI_STATUSES_IGNORE),
              request_msgs_ddt_packed->comm);
#if __GNUC__ == 11 && defined MPICH
#pragma GCC diagnostic pop
#endif
  unpack_data(request_msgs_ddt_packed);
  free_request(request_msgs_ddt_packed);
}

static int xt_request_msgs_test_ddt_packed(Xt_request request) {

  Xt_request_msgs_ddt_packed request_msgs_ddt_packed =
    (Xt_request_msgs_ddt_packed)request;

  int flag =
    xt_mpi_test_some(&(request_msgs_ddt_packed->n_requests),
                     request_msgs_ddt_packed->requests,
                     request_msgs_ddt_packed->ops_completed_buffer,
                     request_msgs_ddt_packed->comm);

  if (flag) {
    unpack_data(request_msgs_ddt_packed);
    free_request(request_msgs_ddt_packed);
  }

  return flag;
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
