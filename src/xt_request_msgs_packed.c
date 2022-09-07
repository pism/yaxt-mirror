/**
 * @file xt_request_msgs_packed.c
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

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"

#include "xt/xt_request_msgs_packed.h"
#include "xt_request_msgs_packed_internal.h"
#include "xt_mpi_internal.h"
#include "xt_request_internal.h"

static void xt_request_msgs_wait_packed(Xt_request request);
static int xt_request_msgs_test_packed(Xt_request request);

static const struct Xt_request_vtable request_msgs_packed_vtable = {

  .wait = xt_request_msgs_wait_packed,
  .test = xt_request_msgs_test_packed,
};

typedef struct Xt_request_msgs_packed_ * Xt_request_msgs_packed;

struct Xt_request_msgs_packed_ {

  const struct Xt_request_vtable * vtable;
  MPI_Request * requests;
  int * ops_completed_buffer;
  MPI_Comm comm;
  int n_requests;
  int n_packed, n_tmp_buffers;
  void * unpacked_data;
  void (*finalizer)(Xt_request, void *);
  size_t extra_buf_ofs;
  void *buffers[];
};

#define no_finalizer ((void (*)(Xt_request, void *))0)

Xt_request
xt_request_msgs_packed_alloc(int n_requests,
                             int n_packed,
                             int n_tmp_buffers,
                             MPI_Comm comm,
                             size_t extra_buf_size)
{
  size_t hdr_size = sizeof(struct Xt_request_msgs_packed_),
    bufp_size
    = ((size_t)n_packed + (size_t)n_tmp_buffers) * sizeof(void *),
    bufr_size = (size_t)n_requests * sizeof(MPI_Request)
    + (size_t)n_requests * sizeof(int);
  /* round allocation up to next multiple of sizeof (void *) to place
   * beginning of user requested extra buffer */
  size_t ofs = ((hdr_size + bufp_size + bufr_size + sizeof (void *) - 1)
                   & -sizeof(void *));
  Xt_request_msgs_packed request = xmalloc(ofs + extra_buf_size);
  request->vtable = &request_msgs_packed_vtable;
  request->comm = comm;
  request->n_requests = n_requests;
  request->n_packed = n_packed;
  request->n_tmp_buffers = n_tmp_buffers;
  request->requests = (MPI_Request *)(request->buffers+n_requests);
  request->finalizer = no_finalizer;
  request->ops_completed_buffer
    = (int *)(request->requests + (size_t)n_requests);
  request->extra_buf_ofs = ofs;
  return (Xt_request)request;
}

static inline void *
extra_buf(Xt_request_msgs_packed request_msgs_packed)
{
  size_t ofs = request_msgs_packed->extra_buf_ofs;
  return (unsigned char *)request_msgs_packed + ofs;
}

Xt_request xt_request_msgs_packed_new(int n_requests,
                                      const MPI_Request requests[n_requests],
                                      MPI_Comm comm, int n_packed,
                                      int n_tmp_buffers,
                                      void *packed_data[n_packed],
                                      void *tmp_buffers[n_tmp_buffers],
                                      void * unpacked_data) {

  assert(n_requests >= 0 && n_packed >= 0 && n_tmp_buffers >= 0);
  Xt_request_msgs_packed request
    = (Xt_request_msgs_packed)xt_request_msgs_packed_alloc(
      n_requests, n_packed, n_tmp_buffers, comm, 0);
  memcpy(request->requests, requests,
         (size_t)n_requests * sizeof(*request->requests));
  memcpy(request->buffers, packed_data,
         (size_t)n_packed * sizeof(*request->buffers));
  memcpy(request->buffers + n_packed, tmp_buffers,
         (size_t)n_tmp_buffers * sizeof(*request->buffers));
  request->unpacked_data = unpacked_data;

  return (Xt_request)request;
}

static void
finish_requests(Xt_request_msgs_packed request_msgs_packed)
{
  if (request_msgs_packed->finalizer != no_finalizer)
    request_msgs_packed->finalizer((Xt_request)request_msgs_packed,
                                   extra_buf(request_msgs_packed));
  free(request_msgs_packed);
}

static void xt_request_msgs_wait_packed(Xt_request request) {

  Xt_request_msgs_packed request_msgs_packed = (Xt_request_msgs_packed)request;

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ == 11 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
  xt_mpi_call(MPI_Waitall(request_msgs_packed->n_requests,
                          request_msgs_packed->requests, MPI_STATUSES_IGNORE),
              request_msgs_packed->comm);
#if __GNUC__ == 11 && defined MPICH
#pragma GCC diagnostic pop
#endif
  finish_requests(request_msgs_packed);
}

static int xt_request_msgs_test_packed(Xt_request request) {

  Xt_request_msgs_packed request_msgs_packed = (Xt_request_msgs_packed)request;

  int flag =
    xt_mpi_test_some(&(request_msgs_packed->n_requests),
                     request_msgs_packed->requests,
                     request_msgs_packed->ops_completed_buffer,
                     request_msgs_packed->comm);

  if (flag)
    finish_requests(request_msgs_packed);

  return flag;
}

MPI_Request *
xt_request_msgs_packed_get_req_ptr(Xt_request request)
{
  Xt_request_msgs_packed request_msgs_packed = (Xt_request_msgs_packed)request;
  assert(request->vtable == &request_msgs_packed_vtable);
  return request_msgs_packed->requests;
}

void **
xt_request_msgs_packed_get_packed_ptr(Xt_request request)
{
  Xt_request_msgs_packed request_msgs_packed = (Xt_request_msgs_packed)request;
  assert(request->vtable == &request_msgs_packed_vtable);
  return request_msgs_packed->buffers;
}

void **
xt_request_msgs_packed_get_unpacked_ptr(Xt_request request)
{
  Xt_request_msgs_packed request_msgs_packed = (Xt_request_msgs_packed)request;
  assert(request->vtable == &request_msgs_packed_vtable);
  return &request_msgs_packed->unpacked_data;
}

void **
xt_request_msgs_packed_get_tmp_buf_ptr(Xt_request request)
{
  Xt_request_msgs_packed request_msgs_packed = (Xt_request_msgs_packed)request;
  assert(request->vtable == &request_msgs_packed_vtable);
  return request_msgs_packed->buffers+request_msgs_packed->n_packed;
}

MPI_Comm
xt_request_msgs_packed_get_comm(Xt_request request)
{
  Xt_request_msgs_packed request_msgs_packed = (Xt_request_msgs_packed)request;
  assert(request->vtable == &request_msgs_packed_vtable);
  return request_msgs_packed->comm;
}

void
xt_request_msgs_packed_set_finalizer(Xt_request request,
                                     void (*finalizer)(Xt_request, void *))
{
  Xt_request_msgs_packed request_msgs_packed = (Xt_request_msgs_packed)request;
  assert(request->vtable == &request_msgs_packed_vtable);
  request_msgs_packed->finalizer = finalizer;
}

void *
xt_request_msgs_packed_get_extra_buf(Xt_request request)
{
  Xt_request_msgs_packed request_msgs_packed = (Xt_request_msgs_packed)request;
  assert(request->vtable == &request_msgs_packed_vtable);
  return extra_buf(request_msgs_packed);
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
