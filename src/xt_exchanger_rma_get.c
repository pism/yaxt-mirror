/**
 * @file xt_exchanger_rma_get.c
 *
 * @copyright Copyright  (C)  2026 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Nils-Arne Dreier <dreier@dkrz.de>
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_config_internal.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt_exchanger.h"
#include "xt_request_internal.h"
#include "xt_exchanger_rma_get.h"

static const char filename[] = "xt_exchanger_rma_get.c";

/* ------------------------------------------------------------------ */
/* Internal data structures                                            */
/* ------------------------------------------------------------------ */

typedef struct Xt_exchanger_rma_get_ *Xt_exchanger_rma_get;

struct Xt_exchanger_rma_get_ {
  const struct xt_exchanger_vtable *vtable;

  int nsend, nrecv;
  MPI_Comm comm;
  MPI_Aint win_size;       /* total byte size of our send window */
  size_t total_recv_size;     /* total staging buffer size (bytes) */
  MPI_Win win;             /* pre-allocated RMA window (persistent) */
  void *win_mem;           /* pointer to window memory */
  MPI_Group post_group;    /* group of send partners (who GET from us) */
  MPI_Group start_group;   /* group of recv partners (who we GET from) */
  bool exposure_epoch_open; /* true between MPI_Win_post and MPI_Win_wait */
  int ready_tag;            /* tag for point-to-point readiness messages */
  /*
   * send_win_displs[i]     = byte displacement in our window for send partner i.
   * recv_win_displs[i]     = byte displacement in recv partner i's window where
   *                          our data resides (communicated during construction).
   * recv_win_sizes[i]      = packed byte size from recv partner i.
   * recv_staging_displs[i] = byte offset into staging buffer for recv partner i.
   */
  MPI_Aint *send_win_displs;     /* [nsend] */
  MPI_Aint *recv_win_displs;     /* [nrecv] */
  int      *recv_win_sizes;      /* [nrecv] */
  int      *recv_staging_displs; /* [nrecv] */
  struct Xt_redist_msg msgs[];   /* nsend send msgs, then nrecv recv msgs */
};

/* ------------------------------------------------------------------ */
/* Async request object                                                */
/* ------------------------------------------------------------------ */

/*
 * After calling MPI_Win_complete / MPI_Win_wait the RMA epoch is closed and
 * all data is in the staging buffers.  The request object keeps the state
 * needed to finish (unpack + free) on wait/test.
 */
struct Xt_request_rma_get_ {
  const struct Xt_request_vtable *vtable;

  void *dst_data;
  void *staging_buf;
  Xt_exchanger_rma_get exchanger; /* back-pointer to clear exposure_epoch_open */
  int n_ready_reqs;           /* nsend + nrecv readiness notifications */
  MPI_Request *ready_reqs;    /* array of readiness notification requests */
  bool access_epoch_started;  /* true once Win_start has been called */
};

typedef struct Xt_request_rma_get_ *Xt_request_rma_get;

/* ------------------------------------------------------------------ */
/* Forward declarations                                                */
/* ------------------------------------------------------------------ */

static Xt_exchanger
xt_exchanger_rma_get_copy(Xt_exchanger exchanger,
                          MPI_Comm new_comm, int new_tag_offset);
static void
xt_exchanger_rma_get_delete(Xt_exchanger exchanger);
static void
xt_exchanger_rma_get_s_exchange(Xt_exchanger exchanger,
                                const void *src_data, void *dst_data);
static void
xt_exchanger_rma_get_a_exchange(Xt_exchanger exchanger,
                                const void *src_data, void *dst_data,
                                Xt_request *request);
static int
xt_exchanger_rma_get_get_msg_ranks(Xt_exchanger exchanger,
                                   enum xt_msg_direction direction,
                                   int *restrict *ranks);
static MPI_Datatype
xt_exchanger_rma_get_get_MPI_Datatype(Xt_exchanger exchanger, int rank,
                                      enum xt_msg_direction direction,
                                      bool do_dup);
static Xt_exchanger_omp_share
xt_exchanger_rma_get_create_omp_share(Xt_exchanger exchanger);

static void xt_request_rma_get_wait(Xt_request request);
static int  xt_request_rma_get_test(Xt_request request);

static void xt_request_rma_get_access_epoch(Xt_request_rma_get r){

  Xt_exchanger_rma_get e = r->exchanger;

  if(e->nrecv == 0)
    return;

  assert(e->total_recv_size > 0);

   /* --- Step 3: Open our access epoch on the send partners' windows ---
      We can use MPI_MODE_NOCHECK because we did it manually
    */
  xt_mpi_call(MPI_Win_start(e->start_group, MPI_MODE_NOCHECK, e->win), e->comm);

  /* --- Step 4: Allocate staging buffer and issue MPI_Get for each recv partner --- */
  r->staging_buf =
    xmalloc(e->total_recv_size);

  for (int i = 0; i < e->nrecv; i++) {
    xt_mpi_call(MPI_Get((char *)r->staging_buf + e->recv_staging_displs[i],
                        e->recv_win_sizes[i], MPI_BYTE,
                        e->msgs[e->nsend + i].rank,
                        e->recv_win_displs[i],
                        e->recv_win_sizes[i], MPI_BYTE,
                        e->win),
                e->comm);
  }

  /* --- Step 5: Complete our access epoch (all GETs issued) --- */
  xt_mpi_call(MPI_Win_complete(e->win), e->comm);

}

/* ------------------------------------------------------------------ */
/* vtable                                                               */
/* ------------------------------------------------------------------ */

static const struct Xt_request_vtable request_rma_get_vtable = {
  .wait = xt_request_rma_get_wait,
  .test = xt_request_rma_get_test,
};

const struct xt_exchanger_vtable xt_exchanger_rma_get_vtable = {
  .copy             = xt_exchanger_rma_get_copy,
  .delete           = xt_exchanger_rma_get_delete,
  .s_exchange       = xt_exchanger_rma_get_s_exchange,
  .a_exchange       = xt_exchanger_rma_get_a_exchange,
  .get_msg_ranks    = xt_exchanger_rma_get_get_msg_ranks,
  .get_MPI_Datatype = xt_exchanger_rma_get_get_MPI_Datatype,
  .create_omp_share = xt_exchanger_rma_get_create_omp_share,
};


static void
xt_request_rma_get_wait(Xt_request request)
{
  Xt_request_rma_get r = (Xt_request_rma_get)request;
  Xt_exchanger_rma_get e = r->exchanger;

  if (!r->access_epoch_started) {
    if (r->n_ready_reqs > 0)
      xt_mpi_call(MPI_Waitall(r->n_ready_reqs, r->ready_reqs,
                               MPI_STATUSES_IGNORE), e->comm);
    free(r->ready_reqs);
    xt_request_rma_get_access_epoch(r);
    r->access_epoch_started = true;
  }

  xt_mpi_call(MPI_Win_wait(e->win), e->comm);
  e->exposure_epoch_open = false;

  for (int i = 0; i < e->nrecv; i++) {
    int position = 0;
    xt_mpi_call(MPI_Unpack((char *)r->staging_buf + e->recv_staging_displs[i],
                           e->recv_win_sizes[i], &position,
                           r->dst_data,
                           1, e->msgs[e->nsend + i].datatype,
                           e->comm),
                e->comm);
  }

  free(r->staging_buf);

  free(r);
}

static int
xt_request_rma_get_test(Xt_request request)
{
  Xt_request_rma_get r = (Xt_request_rma_get)request;
  Xt_exchanger_rma_get e = r->exchanger;

  if (!r->access_epoch_started) {
    int flag = 0;
    if (r->n_ready_reqs > 0)
      xt_mpi_call(MPI_Testall(r->n_ready_reqs, r->ready_reqs, &flag,
                               MPI_STATUSES_IGNORE), e->comm);
    else
      flag = 1;
    if (!flag)
      return 0;
    free(r->ready_reqs);
    xt_request_rma_get_access_epoch(r);
    r->access_epoch_started = true;
  }

  int flag = 0;
  xt_mpi_call(MPI_Win_test(e->win, &flag), e->comm);
  if (!flag)
    return 0;

  /* Exposure epoch closed — clear flag, unpack and finish */
  r->exchanger->exposure_epoch_open = false;

  for (int i = 0; i < e->nrecv; i++) {
    int position = 0;
    xt_mpi_call(MPI_Unpack((char *)r->staging_buf + e->recv_staging_displs[i],
                           e->recv_win_sizes[i], &position,
                           r->dst_data,
                           1, e->msgs[e->nsend + i].datatype,
                           e->comm),
                e->comm);
  }

  free(r->staging_buf);

  free(r);
  return 1;
}

/* ------------------------------------------------------------------ */
/* Exchanger vtable implementations                                    */
/* ------------------------------------------------------------------ */

Xt_exchanger
xt_exchanger_rma_get_new(int nsend, int nrecv,
                         const struct Xt_redist_msg *send_msgs,
                         const struct Xt_redist_msg *recv_msgs,
                         MPI_Comm comm, int tag_offset,
                         Xt_config config)
{
  assert((nsend >= 0) & (nrecv >= 0));

  Xt_exchanger_rma_get exchanger = xmalloc(
    sizeof(*exchanger)
    + ((size_t)nsend + (size_t)nrecv) * sizeof(exchanger->msgs[0]));

  exchanger->vtable = &xt_exchanger_rma_get_vtable;
  exchanger->nsend  = nsend;
  exchanger->nrecv  = nrecv;
  exchanger->comm = comm;
  int tag = tag_offset + xt_mpi_tag_exchange_msg;

  bool dt_dup = !(config->flags & exch_no_dt_dup);
  xt_redist_msgs_strided_copy((size_t)nsend, send_msgs, sizeof(send_msgs[0]),
                              exchanger->msgs, sizeof(exchanger->msgs[0]),
                              exchanger->comm, dt_dup);
  xt_redist_msgs_strided_copy((size_t)nrecv, recv_msgs, sizeof(recv_msgs[0]),
                              exchanger->msgs + nsend,
                              sizeof(exchanger->msgs[0]),
                              exchanger->comm, dt_dup);

  /* Compute our send window layout: displacement per send partner */
  MPI_Aint *send_win_displs =
    xmalloc((size_t)(nsend > 0 ? nsend : 1) * sizeof(*send_win_displs));
  MPI_Aint total_win_size = 0;
  for (int i = 0; i < nsend; i++) {
    int packed_size;
    xt_mpi_call(MPI_Pack_size(1, exchanger->msgs[i].datatype,
                              exchanger->comm, &packed_size),
                exchanger->comm);
    send_win_displs[i] = total_win_size;
    total_win_size    += (MPI_Aint)packed_size;
  }
  exchanger->win_size        = total_win_size;
  exchanger->send_win_displs = send_win_displs;

  /*
   * Exchange displacement information:
   * - Tell each recv partner at what displacement in our window their data is.
   *   For recv partner j, we find the send_msg entry for rank j and send
   *   send_win_displs[i] to rank j.
   * - Receive from each send partner the displacement in their window where
   *   our data lives (they send us their send_win_displs entry for our rank).
   *
   * We use MPI_Send/Recv here (construction, not the exchange itself).
   */
  MPI_Aint *recv_win_displs =
    xmalloc((size_t)(nrecv > 0 ? nrecv : 1) * sizeof(*recv_win_displs));
  int *recv_win_sizes =
    xmalloc((size_t)(nrecv > 0 ? nrecv : 1) * sizeof(*recv_win_sizes));

  /* Compute recv sizes locally — MPI_Pack_size is deterministic */
  for (int i = 0; i < nrecv; i++) {
    int packed_size;
    xt_mpi_call(MPI_Pack_size(1, exchanger->msgs[nsend + i].datatype,
                              exchanger->comm, &packed_size),
                exchanger->comm);
    recv_win_sizes[i] = packed_size;
  }
  exchanger->recv_win_sizes = recv_win_sizes;

  /*
   * Post Irecvs for displacements from each send partner (the send partners
   * will tell us where in their window our data is).
   */
  MPI_Request *reqs =
    xmalloc((size_t)(nsend + nrecv) * sizeof(*reqs));
  for (int i = 0; i < nrecv; i++) {
    int rank = exchanger->msgs[nsend + i].rank;
    xt_mpi_call(MPI_Irecv(&recv_win_displs[i], 1, MPI_AINT,
                          rank, tag, exchanger->comm, &reqs[i]),
                exchanger->comm);
  }

  /*
   * For each recv partner, find the matching send_msg entry and send its
   * displacement.  The recv partner (who is our send target) needs to know
   * at which offset in our window we placed their message.
   */
  for (int i = 0; i < nsend; i++) {
    int rank = exchanger->msgs[i].rank;
    xt_mpi_call(MPI_Isend(&send_win_displs[i], 1, MPI_AINT,
                          rank, tag, exchanger->comm, &reqs[nrecv + i]),
                exchanger->comm);
  }

  xt_mpi_call(MPI_Waitall(nsend + nrecv, reqs, MPI_STATUSES_IGNORE),
              exchanger->comm);
  free(reqs);
  /* Compute staging buffer geometry */
  int *recv_staging_displs =
    xmalloc((size_t)(nrecv > 0 ? nrecv : 1) * sizeof(*recv_staging_displs));
  size_t total_recv_size = 0;
  for (int i = 0; i < nrecv; i++) {
    recv_staging_displs[i] = total_recv_size;
    total_recv_size       += recv_win_sizes[i];
  }
  exchanger->recv_staging_displs = recv_staging_displs;
  exchanger->total_recv_size     = total_recv_size;

  exchanger->recv_win_displs = recv_win_displs;

  /* Build PSCW groups */
  MPI_Group comm_group;
  xt_mpi_call(MPI_Comm_group(exchanger->comm, &comm_group), exchanger->comm);

  /* post_group: send partners (processes that GET from our window) */
  int *send_ranks = xmalloc((size_t)(nsend > 0 ? nsend : 1) * sizeof(*send_ranks));
  for (int i = 0; i < nsend; i++)
    send_ranks[i] = exchanger->msgs[i].rank;
  xt_mpi_call(MPI_Group_incl(comm_group, nsend, send_ranks,
                             &exchanger->post_group), exchanger->comm);
  free(send_ranks);

  /* start_group: recv partners (processes whose windows we GET from) */
  int *recv_ranks = xmalloc((size_t)(nrecv > 0 ? nrecv : 1) * sizeof(*recv_ranks));
  for (int i = 0; i < nrecv; i++)
    recv_ranks[i] = exchanger->msgs[nsend + i].rank;
  xt_mpi_call(MPI_Group_incl(comm_group, nrecv, recv_ranks,
                             &exchanger->start_group), exchanger->comm);
  free(recv_ranks);

  xt_mpi_call(MPI_Group_free(&comm_group), exchanger->comm);

  /* Allocate the persistent RMA window for our send data */
  xt_mpi_call(MPI_Win_allocate(total_win_size > 0 ? total_win_size : 1,
                               1 /* disp_unit = byte */,
                               MPI_INFO_NULL, exchanger->comm,
                               &exchanger->win_mem, &exchanger->win),
              exchanger->comm);

  exchanger->exposure_epoch_open = false;
  exchanger->ready_tag = tag;
  return (Xt_exchanger)exchanger;
}

static Xt_exchanger
xt_exchanger_rma_get_copy(Xt_exchanger exchanger,
                          MPI_Comm new_comm, int new_tag_offset)
{
  Xt_exchanger_rma_get e = (Xt_exchanger_rma_get)exchanger;
  return xt_exchanger_rma_get_new(e->nsend, e->nrecv,
                                  e->msgs, e->msgs + e->nsend,
                                  new_comm, new_tag_offset,
                                  &xt_default_config);
}

static void
xt_exchanger_rma_get_delete(Xt_exchanger exchanger)
{
  Xt_exchanger_rma_get e = (Xt_exchanger_rma_get)exchanger;
  size_t nmsg = (size_t)e->nsend + (size_t)e->nrecv;
  xt_redist_msgs_strided_destruct(nmsg, e->msgs, e->comm,
                                  sizeof(e->msgs[0]));
  xt_mpi_call(MPI_Group_free(&e->post_group), e->comm);
  xt_mpi_call(MPI_Group_free(&e->start_group), e->comm);
  free(e->send_win_displs);
  free(e->recv_win_displs);
  free(e->recv_win_sizes);
  free(e->recv_staging_displs);

  //MPI_Win_free is collective. This does not work with YAC currently
  //xt_mpi_call(MPI_Win_free(&e->win), e->comm);

  free(e);
}

static void
xt_exchanger_rma_get_s_exchange(Xt_exchanger exchanger,
                                const void *src_data, void *dst_data)
{
  Xt_request request;
  xt_exchanger_rma_get_a_exchange(exchanger, src_data, dst_data, &request);
  xt_request_rma_get_wait(request);
}

static void
xt_exchanger_rma_get_a_exchange(Xt_exchanger exchanger,
                                const void *src_data, void *dst_data,
                                Xt_request *request)
{
  Xt_exchanger_rma_get e = (Xt_exchanger_rma_get)exchanger;

  /* --- Step 1: Pack source data into the pre-allocated window memory --- */
  for (int i = 0; i < e->nsend; i++) {
    int position = 0;
    xt_mpi_call(MPI_Pack(CAST_MPI_SEND_BUF(src_data),
                         1, e->msgs[i].datatype,
                         (char *)e->win_mem + e->send_win_displs[i],
                         /* outcount */ (int)(e->win_size - e->send_win_displs[i]),
                         &position,
                         e->comm),
                e->comm);
  }

  /* --- Step 2: PSCW — expose our window to all processes that will GET from us --- */
  int post_grp_size;
  MPI_Group_size(e->post_group, &post_grp_size);

  /* Guard: abort if a previous exposure epoch is still open.
   * The caller must call xt_request_wait / xt_request_test before
   * starting a new exchange on the same exchanger. */
  if (e->exposure_epoch_open)
    Xt_abort(e->comm,
             "ERROR(xt_exchanger_rma_get_a_exchange): a_exchange called "
             "while previous exchange is still pending (wait not called)",
             filename, __LINE__);

  // Use MPI_MODE_NOCHECK, because we do manual (a)synchronization
  xt_mpi_call(MPI_Win_post(e->post_group, MPI_MODE_NOCHECK, e->win), e->comm);
  e->exposure_epoch_open = true;
  /* --- Step 6: Build and return the async request --- */
  Xt_request_rma_get req = xmalloc(sizeof(*req));
  req->vtable      = &request_rma_get_vtable;
  req->dst_data    = dst_data;
  req->staging_buf = NULL;
  req->exchanger   = e;

  /* Point-to-point readiness notifications: each partner confirms its
   * window is posted before we open our access epoch. */
  int n_ready = e->nsend + e->nrecv;
  MPI_Request *ready_reqs =
    xmalloc((size_t)(n_ready > 0 ? n_ready : 1) * sizeof(*ready_reqs));

  /* Receive readiness from each recv partner (they notify us when their
   * window is posted so we can safely MPI_Win_start). */
  for (int i = 0; i < e->nrecv; i++) {
    xt_mpi_call(MPI_Irecv(NULL, 0, MPI_BYTE,
                          e->msgs[e->nsend + i].rank, e->ready_tag,
                          e->comm, &ready_reqs[i]),
                e->comm);
  }

  /* Notify each send partner that our window is posted and ready. */
  for (int i = 0; i < e->nsend; i++) {
    xt_mpi_call(MPI_Isend(NULL, 0, MPI_BYTE,
                          e->msgs[i].rank, e->ready_tag,
                          e->comm, &ready_reqs[e->nrecv + i]),
                e->comm);
  }

  req->n_ready_reqs = n_ready;
  req->ready_reqs = ready_reqs;
  req->access_epoch_started = false;

  *request = (Xt_request)req;
}

static int
xt_exchanger_rma_get_get_msg_ranks(Xt_exchanger exchanger,
                                   enum xt_msg_direction direction,
                                   int *restrict *ranks)
{
  Xt_exchanger_rma_get e = (Xt_exchanger_rma_get)exchanger;
  int nmsg = (direction == SEND) ? e->nsend : e->nrecv;
  struct Xt_redist_msg *msgs =
    (direction == SEND) ? e->msgs : e->msgs + e->nsend;

  int *ranks_ = *ranks;
  if (!ranks_)
    ranks_ = *ranks = xmalloc((size_t)nmsg * sizeof(*ranks_));
  for (int i = 0; i < nmsg; i++)
    ranks_[i] = msgs[i].rank;
  return nmsg;
}

static MPI_Datatype
xt_exchanger_rma_get_get_MPI_Datatype(Xt_exchanger exchanger, int rank,
                                      enum xt_msg_direction direction,
                                      bool do_dup)
{
  Xt_exchanger_rma_get e = (Xt_exchanger_rma_get)exchanger;
  int nmsg = (direction == SEND) ? e->nsend : e->nrecv;
  struct Xt_redist_msg *msgs =
    (direction == SEND) ? e->msgs : e->msgs + e->nsend;

  for (int i = 0; i < nmsg; i++) {
    if (msgs[i].rank != rank)
      continue;
    if (!do_dup)
      return msgs[i].datatype;
    MPI_Datatype dt_copy;
    xt_mpi_call(MPI_Type_dup(msgs[i].datatype, &dt_copy), e->comm);
    return dt_copy;
  }
  return MPI_DATATYPE_NULL;
}

static Xt_exchanger_omp_share
xt_exchanger_rma_get_create_omp_share(Xt_exchanger exchanger)
{
  (void)exchanger;
  return NULL;
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
