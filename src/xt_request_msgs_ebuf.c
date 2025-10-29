/**
 * @file xt_request_msgs_ebuf.c
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

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"

#include "xt/xt_request_msgs_ebuf.h"
#include "xt_request_msgs_ebuf_internal.h"
#include "xt_config_internal.h"
#include "xt_mpi_internal.h"
#include "xt_request_internal.h"

static void xt_request_msgs_ebuf_wait(Xt_request request);
static int xt_request_msgs_ebuf_test(Xt_request request);

#ifdef _OPENMP
static void
xt_request_msgs_ebuf_wait_omp(Xt_request request);
static int
xt_request_msgs_ebuf_test_omp(Xt_request request);
#endif

static const struct Xt_request_vtable request_msgs_ebuf_vtable = {

  .wait = xt_request_msgs_ebuf_wait,
  .test = xt_request_msgs_ebuf_test,
}
#ifdef _OPENMP
  , request_msgs_ebuf_auto_omp_vtable = {
  .wait = xt_request_msgs_ebuf_wait_omp,
  .test = xt_request_msgs_ebuf_test_omp,
};
#endif
;

typedef struct Xt_request_msgs_ebuf_ * Xt_request_msgs_ebuf;

struct Xt_request_msgs_ebuf_ {
  const struct Xt_request_vtable *vtable;
  MPI_Comm comm;
  int n_requests, auto_free;
  void (*finalizer)(Xt_request, void *);
  size_t extra_buf_ofs;
  MPI_Request requests[];
};

#define no_finalizer ((void (*)(Xt_request, void *))0)

Xt_request
xt_request_msgs_ebuf_alloc(int n_requests,
                           MPI_Comm comm,
                           size_t extra_buf_size,
                           Xt_config config)
{
  size_t hdr_size = sizeof(struct Xt_request_msgs_ebuf_),
    bufr_size = (size_t)n_requests * sizeof(MPI_Request)
    + (size_t)n_requests * sizeof(int);
  /* round allocation up to next multiple of sizeof (void *) to place
   * beginning of user requested extra buffer */
  size_t ofs = ((hdr_size + bufr_size + sizeof (void *) - 1)
                   & -sizeof(void *));
  Xt_request_msgs_ebuf request = xmalloc(ofs + extra_buf_size);
#ifndef _OPENMP
  (void)config;
#else
  int mthread_mode = xt_config_get_redist_mthread_mode(config);
  if (mthread_mode == XT_MT_OPENMP)
    request->vtable = &request_msgs_ebuf_auto_omp_vtable;
  else
#endif
  request->vtable = &request_msgs_ebuf_vtable;
  request->comm = comm;
  request->n_requests = n_requests;
  request->auto_free = 1;
  request->finalizer = no_finalizer;
  request->extra_buf_ofs = ofs;
  return (Xt_request)request;
}

static inline void *
extra_buf(Xt_request_msgs_ebuf request_msgs_ebuf)
{
  size_t ofs = request_msgs_ebuf->extra_buf_ofs;
  return (unsigned char *)request_msgs_ebuf + ofs;
}

Xt_request
xt_request_msgs_ebuf_new(int n_requests,
                         size_t extra_buf_size,
                         Xt_fill_ebuf_requests init,
                         void *data,
                         Xt_request_msgs_ebuf_finalizer finalize,
                         MPI_Comm comm) {
  return xt_request_msgs_ebuf_custom_new(n_requests, extra_buf_size, init, data,
                                         finalize, comm, &xt_default_config);
}

Xt_request
xt_request_msgs_ebuf_custom_new(int n_requests,
                                size_t extra_buf_size,
                                Xt_fill_ebuf_requests init,
                                void *data,
                                Xt_request_msgs_ebuf_finalizer finalize,
                                MPI_Comm comm,
                                Xt_config config)
{
  assert(n_requests >= 0);
  Xt_request request;
  if (n_requests > 0) {
    Xt_request_msgs_ebuf request_me
      = (Xt_request_msgs_ebuf)xt_request_msgs_ebuf_alloc(
        n_requests, comm, extra_buf_size, config);
    request = (Xt_request)request_me;
    xt_request_msgs_ebuf_set_finalizer(request, finalize);
    init(request, request_me->requests, extra_buf(request_me),
         data, comm);
  } else
    request = XT_REQUEST_NULL;
  return (Xt_request)request;
}

static void
finish_requests(Xt_request_msgs_ebuf request_msgs_ebuf)
{
  if (request_msgs_ebuf->finalizer != no_finalizer)
    request_msgs_ebuf->finalizer((Xt_request)request_msgs_ebuf,
                                   extra_buf(request_msgs_ebuf));
  if (request_msgs_ebuf->auto_free)
    free(request_msgs_ebuf);
}

static void xt_request_msgs_ebuf_wait(Xt_request request) {

  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
  xt_mpi_call(MPI_Waitall(request_msgs_ebuf->n_requests,
                          request_msgs_ebuf->requests, MPI_STATUSES_IGNORE),
              request_msgs_ebuf->comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic pop
#endif
  finish_requests(request_msgs_ebuf);
}

static int xt_request_msgs_ebuf_test(Xt_request request) {

  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;

  int flag =
    xt_mpi_test_some(&(request_msgs_ebuf->n_requests),
                     request_msgs_ebuf->requests,
                     (int *)(request_msgs_ebuf->requests
                             + request_msgs_ebuf->n_requests),
                     request_msgs_ebuf->comm);

  if (flag)
    finish_requests(request_msgs_ebuf);

  return flag;
}

MPI_Request *
xt_request_msgs_ebuf_get_req_ptr(Xt_request request)
{
  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;
#ifndef _OPENMP
  assert(request->vtable == &request_msgs_ebuf_vtable);
#else
  assert(request->vtable == &request_msgs_ebuf_vtable
         || request->vtable == &request_msgs_ebuf_auto_omp_vtable);
#endif
  return request_msgs_ebuf->requests;
}

void
xt_request_msgs_ebuf_set_auto_free(Xt_request request, int auto_free)
{
  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;
#ifndef _OPENMP
  assert(request->vtable == &request_msgs_ebuf_vtable);
#else
  assert(request->vtable == &request_msgs_ebuf_vtable
         || request->vtable == &request_msgs_ebuf_auto_omp_vtable);
#endif
  request_msgs_ebuf->auto_free = auto_free;
}

MPI_Comm
xt_request_msgs_ebuf_get_comm(Xt_request request)
{
  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;
#ifndef _OPENMP
  assert(request->vtable == &request_msgs_ebuf_vtable);
#else
  assert(request->vtable == &request_msgs_ebuf_vtable
         || request->vtable == &request_msgs_ebuf_auto_omp_vtable);
#endif
  return request_msgs_ebuf->comm;
}

void
xt_request_msgs_ebuf_set_finalizer(Xt_request request,
                                   Xt_request_msgs_ebuf_finalizer finalizer)
{
  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;
#ifndef _OPENMP
  assert(request->vtable == &request_msgs_ebuf_vtable);
#else
  assert(request->vtable == &request_msgs_ebuf_vtable
         || request->vtable == &request_msgs_ebuf_auto_omp_vtable);
#endif
  request_msgs_ebuf->finalizer = finalizer;
}

void *
xt_request_msgs_ebuf_get_extra_buf(Xt_request request)
{
  Xt_request_msgs_ebuf request_msgs_ebuf = (Xt_request_msgs_ebuf)request;
#ifndef _OPENMP
  assert(request->vtable == &request_msgs_ebuf_vtable);
#else
  assert(request->vtable == &request_msgs_ebuf_vtable
         || request->vtable == &request_msgs_ebuf_auto_omp_vtable);
#endif
  return extra_buf(request_msgs_ebuf);
}


#ifdef _OPENMP
static void
finish_requests_mt(Xt_request_msgs_ebuf request_msgs_ebuf)
{
  if (request_msgs_ebuf->finalizer != no_finalizer)
    request_msgs_ebuf->finalizer((Xt_request)request_msgs_ebuf,
                                   extra_buf(request_msgs_ebuf));
#pragma omp barrier
#pragma omp master
  if (request_msgs_ebuf->auto_free)
    free(request_msgs_ebuf);
}

static size_t
xt_request_msgs_ebuf_get_num_omp_threads(Xt_request request)
{
  size_t num_threads= (size_t)omp_get_max_threads();
  size_t num_req = (size_t)(((Xt_request_msgs_ebuf)request)->n_requests);
  return num_req >= num_threads ? num_threads : num_req;
}


static void
xt_request_msgs_ebuf_wait_omp(Xt_request request)
{
  size_t num_threads_limit = xt_request_msgs_ebuf_get_num_omp_threads(request);
#pragma omp parallel num_threads(num_threads_limit)
  {
    Xt_request_msgs_ebuf request_msgs = (Xt_request_msgs_ebuf)request;
    int num_threads = omp_get_num_threads(),
      tid = omp_get_thread_num();
    int nreq = request_msgs->n_requests,
      start_req = (nreq * tid) / num_threads,
      nreq_ = (nreq * (tid+1)) / num_threads - start_req;

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
    xt_mpi_call(MPI_Waitall(nreq_, request_msgs->requests+start_req,
                            MPI_STATUSES_IGNORE), request_msgs->comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic pop
#endif
#pragma omp barrier
    finish_requests_mt(request_msgs);
  }
}

static int
xt_request_msgs_ebuf_test_omp(Xt_request request)
{
  bool flag = true;
  size_t num_threads_limit = xt_request_msgs_ebuf_get_num_omp_threads(request);
  /* cache_line_size_scale_factor * sizeof (int)
   * should equal the cache-line size and
   * cache_line_size_scale_factor >= 2 must be true */
  enum { cache_line_size_scale_factor = 16, auto_limit = 32 };
  int ops_comp_count_auto[auto_limit][cache_line_size_scale_factor],
    (*ops_comp_count)[cache_line_size_scale_factor];
  if (num_threads_limit <= auto_limit)
    ops_comp_count = ops_comp_count_auto;
  else
    ops_comp_count = xmalloc(sizeof (int) * num_threads_limit
                             * cache_line_size_scale_factor);
#pragma omp parallel num_threads(num_threads_limit) firstprivate(ops_comp_count)
  {
    Xt_request_msgs_ebuf request_msgs = (Xt_request_msgs_ebuf)request;

    size_t num_threads = (size_t)omp_get_num_threads(),
      tid = (size_t)omp_get_thread_num(),
      num_req = (size_t)request_msgs->n_requests,
      start_req = (num_req * tid) / num_threads,
      nreq = (num_req * (tid+1)) / num_threads - start_req;
    MPI_Request *reqs_copy, reqs_copy_[auto_limit],
      *restrict reqs_orig = request_msgs->requests;
    reqs_copy = nreq <= auto_limit
      ? reqs_copy_
      : xmalloc(sizeof (*reqs_copy) * nreq);
    int *ops_completed_buffer
      = (int *)(reqs_orig + num_req) + start_req;
    for (size_t i = 0; i < nreq; ++i)
      reqs_copy[i] = reqs_orig[start_req+i];
    struct test_some_partial_result test_res
      = xt_mpi_test_some_mt((int)nreq, reqs_copy,
                            ops_completed_buffer, request_msgs->comm);
    size_t done_count = (size_t)test_res.num_req_completed;
    for (size_t i = 0; i < done_count; ++i)
      ops_completed_buffer[i] += (int)start_req;
    ops_comp_count[tid][0] = test_res.num_req_completed;
    ops_comp_count[tid][1] = (int)nreq;
#pragma omp barrier
    size_t accum = 0;
    for (size_t i = 0; i < tid; ++i)
      accum += (size_t)ops_comp_count[i][0];
    size_t req_ofs = start_req - accum;
    for (size_t i = tid; i < num_threads; ++i)
      accum += (size_t)ops_comp_count[i][0];
    if (accum == num_req)
      finish_requests_mt(request_msgs);
    else {
      for (size_t i = 0, j = 0; i < nreq; ++i) {
        if (reqs_copy[i] != MPI_REQUEST_NULL)
          reqs_orig[req_ofs+j] = reqs_copy[i];
        j += reqs_copy[i] != MPI_REQUEST_NULL;
      }
#pragma omp master
      request_msgs->n_requests -= (int)accum;
    }
    if (nreq <= auto_limit) ; else
      free(reqs_copy);
#pragma omp master
    flag = accum == num_req;
  }
  if (num_threads_limit > auto_limit)
    free(ops_comp_count);
  return flag;
}
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
