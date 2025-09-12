/**
 * @file xt_request_msgs.c
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
#include "xt/xt_config.h"
#include "xt/xt_request_msgs.h"
#include "xt_request_msgs_internal.h"
#include "xt_mpi_internal.h"
#include "xt_config_internal.h"
#include "xt_request_internal.h"

static void xt_request_msgs_wait(Xt_request request);
static int xt_request_msgs_test(Xt_request request);

#ifdef _OPENMP
static void xt_request_msgs_wait_omp(Xt_request request);
static int xt_request_msgs_test_omp(Xt_request request);
#endif

static const struct Xt_request_vtable request_msgs_vtable = {

  .wait = xt_request_msgs_wait,
  .test = xt_request_msgs_test,
}
#ifdef _OPENMP
  , request_msgs_auto_omp_vtable = {
  .wait = xt_request_msgs_wait_omp,
  .test = xt_request_msgs_test_omp,
};
#endif
  ;

typedef struct Xt_request_msgs_ *Xt_request_msgs;

struct Xt_request_msgs_ {

  const struct Xt_request_vtable *vtable;
  int n, auto_free;
  MPI_Comm comm;
  MPI_Request requests[];
};

Xt_request
xt_request_msgs_alloc(int n, MPI_Comm comm, Xt_config config)
{
  Xt_request_msgs request
    = xmalloc(sizeof(*request) + (size_t)n * sizeof(MPI_Request) +
              (size_t)n * sizeof(int));
#ifndef _OPENMP
  (void)config;
#else
  int mthread_mode = xt_config_get_redist_mthread_mode(config);
  if (mthread_mode == XT_MT_OPENMP)
    request->vtable = &request_msgs_auto_omp_vtable;
  else
#endif
    request->vtable = &request_msgs_vtable;
  request->n = n;
  request->auto_free = 1;
  request->comm = comm;
  return (Xt_request)request;
}


Xt_request
xt_request_msgs_new(int n,
                    const MPI_Request requests[n],
                    MPI_Comm comm)
{
  return xt_request_msgs_custom_new(n, requests, comm, &xt_default_config);
}

Xt_request
xt_request_msgs_custom_new(int n,
                           const MPI_Request requests[n],
                           MPI_Comm comm,
                           Xt_config config)
{
  assert(n >= 0);
  Xt_request_msgs request
    = (Xt_request_msgs)xt_request_msgs_alloc(n, comm, config);

  memcpy(request->requests, requests, (size_t)n * sizeof(*request->requests));

  return (Xt_request)request;
}

MPI_Request *
xt_request_msgs_get_req_ptr(Xt_request request)
{
  Xt_request_msgs request_msgs = (Xt_request_msgs)request;
#ifdef _OPENMP
  assert(request->vtable == &request_msgs_vtable
         || request->vtable == &request_msgs_auto_omp_vtable);
#else
  assert(request->vtable == &request_msgs_vtable);
#endif
  return request_msgs->requests;
}

void
xt_request_msgs_set_auto_free(Xt_request request, int auto_free)
{
  Xt_request_msgs request_msgs = (Xt_request_msgs)request;
  request_msgs->auto_free = auto_free;
}


static void xt_request_msgs_wait(Xt_request request) {

  Xt_request_msgs request_msgs = (Xt_request_msgs)request;

#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
  /* GCC 11 has no means to specify that the special value pointer
   * MPI_STATUSES_IGNORE does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
  xt_mpi_call(MPI_Waitall(request_msgs->n, request_msgs->requests,
                          MPI_STATUSES_IGNORE), request_msgs->comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic pop
#endif
  if (request_msgs->auto_free)
    free(request_msgs);
}

static int xt_request_msgs_test(Xt_request request) {

  Xt_request_msgs request_msgs = (Xt_request_msgs)request;

  size_t n = (size_t)request_msgs->n;
  int *ops_completed_buffer = (int *)(request_msgs->requests + n);
  int flag = xt_mpi_test_some(&(request_msgs->n), request_msgs->requests,
                              ops_completed_buffer, request_msgs->comm);

  if (flag && request_msgs->auto_free) free(request_msgs);

  return flag;
}

#ifdef _OPENMP
static void
xt_request_msgs_wait_omp(Xt_request request)
{
#pragma omp parallel
  {
    Xt_request_msgs request_msgs = (Xt_request_msgs)request;

#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
  /* GCC 11 has no means to specify that the special value pointer
   * MPI_STATUSES_IGNORE does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif
    int num_threads = omp_get_num_threads(),
      tid = omp_get_thread_num();
    int nreq = request_msgs->n,
      start_req = (nreq * tid) / num_threads,
      nreq_ = (nreq * (tid+1)) / num_threads - start_req;
    xt_mpi_call(MPI_Waitall(nreq_, request_msgs->requests+start_req,
                            MPI_STATUSES_IGNORE), request_msgs->comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic pop
#endif
  }

  if (((Xt_request_msgs)request)->auto_free)
    free(request);
}

static int
xt_request_msgs_test_omp(Xt_request request) {

  bool flag = true;
  size_t num_threads = (size_t)omp_get_max_threads();
  /* cache_line_size_scale_factor * sizeof (int)
   * should equal the cache-line size and
   * cache_line_size_scale_factor >= 2 must be true */
  enum { cache_line_size_scale_factor = 16 };
  int ops_comp_count[num_threads][cache_line_size_scale_factor];
#pragma omp parallel
  {
    Xt_request_msgs request_msgs = (Xt_request_msgs)request;
    size_t num_req = (size_t)request_msgs->n,
      tid = (size_t)omp_get_thread_num(),
      start_req = (num_req * tid) / num_threads,
      nreq = (num_req * (tid+1)) / num_threads - start_req;
    MPI_Request reqs_copy[num_req],
      *restrict reqs_orig = request_msgs->requests;
    int *restrict ops_completed_buffer
      = (int *)(reqs_orig + num_req) + start_req;
    for (size_t i = 0; i < nreq; ++i)
      reqs_copy[i] = reqs_orig[start_req+i];
    struct test_some_partial_result test_res
      = xt_mpi_test_some_mt((int)nreq, reqs_copy,
                            ops_completed_buffer, request_msgs->comm);
    size_t done_count = (size_t)test_res.num_req_completed;
    ops_comp_count[tid][0] = test_res.num_req_completed;
    ops_comp_count[tid][1] = (int)nreq;
    for (size_t i = 0; i < done_count; ++i)
      ops_completed_buffer[i] += (int)start_req;
#pragma omp barrier
    size_t accum = 0;
    for (size_t i = 0; i < tid; ++i)
      accum += (size_t)ops_comp_count[i][0];
    size_t req_ofs = start_req - accum;
    for (size_t i = tid; i < num_threads; ++i)
      accum += (size_t)ops_comp_count[i][0];
    if (accum != num_req) {
      for (size_t i = 0, j = 0; i < nreq; ++i) {
        if (reqs_copy[i] != MPI_REQUEST_NULL)
          reqs_orig[req_ofs+j] = reqs_copy[i];
        j += reqs_copy[i] != MPI_REQUEST_NULL;
      }
    }
#pragma omp master
    {
      flag = accum == num_req;
      request_msgs->n -= (int)accum;
    }
  }
  if (flag && ((Xt_request_msgs)request)->auto_free) free(request);

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
