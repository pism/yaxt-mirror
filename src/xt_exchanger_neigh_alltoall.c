/**
 * @file xt_exchanger_neigh_alltoall.c
 *
 * @copyright Copyright  (C)  2018 Jörg Behrens <behrens@dkrz.de>
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
#include <string.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt_config_internal.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/xt_xmap.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_request.h"
#include "xt/xt_request_msgs.h"
#include "xt_request_msgs_internal.h"
#include "xt_request_msgs_ebuf_internal.h"
#include "xt_exchanger.h"
#include "xt_exchanger_internal.h"
#include "xt_exchanger_neigh_alltoall.h"

static const char filename[] = "xt_exchanger_neigh_alltoall.c";

static Xt_exchanger
xt_exchanger_neigh_alltoall_copy(Xt_exchanger exchanger,
                                 MPI_Comm newComm, int new_tag_offset);
static void xt_exchanger_neigh_alltoall_delete(Xt_exchanger exchanger);
static void
xt_exchanger_neigh_alltoall_s_exchange(Xt_exchanger exchanger,
                                       const void *src_data,
                                       void *dst_data);
static void
xt_exchanger_neigh_alltoall_a_exchange(Xt_exchanger exchanger,
                                       const void *src_data, void *dst_data,
                                       Xt_request *request);
static int
xt_exchanger_neigh_alltoall_get_num_msg(Xt_exchanger exchanger,
                                        int direction);
static int
xt_exchanger_neigh_alltoall_get_msg_ranks(Xt_exchanger exchanger,
                                          enum xt_msg_direction direction,
                                          int *restrict *ranks);

static MPI_Datatype
xt_exchanger_neigh_alltoall_get_MPI_Datatype(Xt_exchanger exchanger,
                                             int rank,
                                             enum xt_msg_direction direction,
                                             bool do_dup);

static Xt_exchanger_iter
xt_exchanger_na2a_get_iterator(Xt_exchanger exchanger, Xt_config config);

struct xt_exchanger_neigh_alltoall_team_share
{
  int *ranks;
  MPI_Comm nb_comm;
  int nmsg[2];
};

static void xt_exchanger_neigh_alltoall_team_share_default_init(void *share)
{
  struct xt_exchanger_neigh_alltoall_team_share *team_share = share;
  team_share->ranks = NULL;
  team_share->nb_comm = MPI_COMM_NULL;
}

static void xt_exchanger_neigh_alltoall_team_share_destroy(void *share)
{
  struct xt_exchanger_neigh_alltoall_team_share *team_share = share;
  if (team_share->nb_comm != MPI_COMM_NULL)
    xt_mpi_call(MPI_Comm_free(&team_share->nb_comm), Xt_default_comm);
  if ((void *)team_share->ranks != (void *)(team_share+1))
    free(team_share->ranks);
}

const struct xt_exchanger_vtable xt_exchanger_neigh_alltoall_vtable = {
  .copy = xt_exchanger_neigh_alltoall_copy,
  .delete = xt_exchanger_neigh_alltoall_delete,
  .s_exchange = xt_exchanger_neigh_alltoall_s_exchange,
  .a_exchange = xt_exchanger_neigh_alltoall_a_exchange,
  .get_num_msg = xt_exchanger_neigh_alltoall_get_num_msg,
  .get_msg_ranks = xt_exchanger_neigh_alltoall_get_msg_ranks,
  .get_MPI_Datatype = xt_exchanger_neigh_alltoall_get_MPI_Datatype,
  .team_share_size = sizeof (struct xt_exchanger_neigh_alltoall_team_share),
  .team_share_default_init
  = xt_exchanger_neigh_alltoall_team_share_default_init,
  .team_share_destroy = xt_exchanger_neigh_alltoall_team_share_destroy,
  .get_iterator = xt_exchanger_na2a_get_iterator,
};

typedef struct Xt_exchanger_neigh_alltoall_ * Xt_exchanger_neigh_alltoall;

struct Xt_exchanger_neigh_alltoall_ {

  const struct xt_exchanger_vtable * vtable;

  MPI_Aint *displs;
  MPI_Datatype *datatypes;
  int *counts;
  struct xt_exchanger_neigh_alltoall_team_share *team_share,
    team_share_[];
};

static Xt_exchanger_neigh_alltoall
xt_exchanger_neigh_alltoall_alloc(size_t nsend, size_t nrecv,
                                  void *exchanger_team_share)
{
  size_t nmsg = nsend + nrecv;
  bool need_team_share_alloc = exchanger_team_share == NULL;
  size_t ranks_buf_size = sizeof (int) * (nmsg+(nmsg&1));
  Xt_exchanger_neigh_alltoall exchanger
    = xmalloc(sizeof(*exchanger)
              + (need_team_share_alloc
                 ? sizeof (*exchanger->team_share) + ranks_buf_size : 0)
              + sizeof (MPI_Aint) * nmsg
              + nmsg * sizeof(*exchanger->datatypes)
              + sizeof (int) * nmsg * ((unsigned)need_team_share_alloc + 1U));
  exchanger->displs
    = (void *)(
      (unsigned char *)(exchanger->team_share_ + need_team_share_alloc)
      + (need_team_share_alloc ? ranks_buf_size : 0));
  exchanger->datatypes = (void *)(exchanger->displs + nmsg);
  exchanger->counts = (void *)(exchanger->datatypes + nmsg);
  exchanger->vtable = &xt_exchanger_neigh_alltoall_vtable;
  if (need_team_share_alloc) {
    xt_exchanger_neigh_alltoall_team_share_default_init(exchanger->team_share_);
    exchanger->team_share = exchanger->team_share_;
    exchanger->team_share->ranks = (void *)(exchanger->team_share_ + 1);
  } else {
    exchanger->team_share = exchanger_team_share;
    if (!exchanger->team_share->ranks)
      exchanger->team_share->ranks = xmalloc(nmsg * sizeof(int));
  }
  return exchanger;
}

static void copy_dt(size_t n,
                    const struct Xt_msg_param *restrict msgs,
                    MPI_Datatype *restrict datatypes,
                    MPI_Comm comm, bool dt_dup) {

  for (size_t i = 0; i < n; ++i)
    if (dt_dup)
      xt_mpi_call(MPI_Type_dup(msgs[i].datatype, datatypes + i), comm);
    else
      datatypes[i] = msgs[i].datatype;
}

static void copy_ranks(size_t n,
                       const struct Xt_msg_param *restrict msgs,
                       int *restrict ranks)
{
  for (size_t i = 0; i < n; ++i)
    ranks[i] = msgs[i].rank;
}

static void copy_displs(size_t n,
                        const struct Xt_msg_param *restrict msgs,
                        MPI_Aint *restrict displs)
{
  for (size_t i = 0; i < n; ++i)
    displs[i] = msgs[i].displacement;
}

static void copy_counts(size_t n,
                        const struct Xt_msg_param *restrict msgs,
                        int *restrict count)
{
  for (size_t i = 0; i < n; ++i)
    count[i] = msgs[i].count;
}

Xt_exchanger
xt_exchanger_neigh_alltoall_new(int nsend, int nrecv,
                                const struct Xt_msg_param *send_msgs,
                                const struct Xt_msg_param *recv_msgs,
                                MPI_Comm comm, int tag_offset, Xt_config config)
{
  /** note: tag_offset + xt_mpi_tag_exchange_msg must not
   *        be used on @a comm by any other part of the program during the
   *        lifetime of the created exchanger object
   */
  (void)tag_offset;
  int flag;
  xt_mpi_call(MPI_Comm_test_inter(comm, &flag), comm);
  if (flag)
    Xt_abort(comm, "ERROR(xt_exchanger_neigh_alltoall_new): "
             "inter-communicator's are not defined for virtual topologies",
             filename, __LINE__);

  assert((nsend >= 0) & (nrecv >= 0));
  Xt_exchanger_neigh_alltoall exchanger
    = xt_exchanger_neigh_alltoall_alloc((size_t)nsend, (size_t)nrecv,
                                        config->exchanger_team_share);
  bool dt_dup = !(config->flags & exch_no_dt_dup);
  copy_displs((size_t)nsend, send_msgs, exchanger->displs);
  copy_displs((size_t)nrecv, recv_msgs, exchanger->displs + nsend);
  copy_dt((size_t)nsend, send_msgs, exchanger->datatypes, comm, dt_dup);
  copy_dt((size_t)nrecv, recv_msgs, exchanger->datatypes + nsend, comm, dt_dup);
  copy_counts((size_t)nsend, send_msgs, exchanger->counts);
  copy_counts((size_t)nrecv, recv_msgs, exchanger->counts + nsend);
  struct xt_exchanger_neigh_alltoall_team_share *team_share
    = exchanger->team_share;
  size_t nmsg = (size_t)nsend + (size_t)nrecv;
  if (!team_share->ranks)
    team_share->ranks = xmalloc(nmsg * sizeof(*team_share->ranks));
  if (team_share->nb_comm == MPI_COMM_NULL) {
    team_share->nmsg[RECV] = nrecv;
    team_share->nmsg[SEND] = nsend;
    copy_ranks((size_t)nsend, send_msgs, team_share->ranks);
    copy_ranks((size_t)nrecv, recv_msgs, team_share->ranks+nsend);
    enum { no_reorder = 0 }; // no reordering of ranks in new comm
#if __GNUC__ >= 11 && __GNUC__ <= 13
    /* GCC 11 has no means to specify that the special value pointer
     * MPI_UNWEIGHTED does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
    xt_mpi_call(
      MPI_Dist_graph_create_adjacent(
        comm, nrecv, team_share->ranks + nsend, MPI_UNWEIGHTED, nsend,
        team_share->ranks, MPI_UNWEIGHTED, MPI_INFO_NULL, no_reorder,
        &team_share->nb_comm), comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13
#pragma GCC diagnostic pop
#endif
  }
  return (Xt_exchanger)exchanger;
}

static Xt_exchanger
xt_exchanger_neigh_alltoall_copy(Xt_exchanger exchanger,
                                 MPI_Comm new_comm, int new_tag_offset)
{
  (void)new_tag_offset;
  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  const struct xt_exchanger_neigh_alltoall_team_share
    *team_share = exchanger_na->team_share;
  size_t nsend = (size_t)team_share->nmsg[SEND],
    nrecv = (size_t)team_share->nmsg[RECV];

  Xt_exchanger_neigh_alltoall
    exchanger_copy = xt_exchanger_neigh_alltoall_alloc(nsend, nrecv, NULL);

  struct xt_exchanger_neigh_alltoall_team_share
    *team_share_copy = exchanger_copy->team_share;
  size_t nmsg = nsend + nrecv;
  memcpy(exchanger_copy->displs, exchanger_na->displs,
         sizeof (MPI_Aint) * nmsg);
  memcpy(exchanger_copy->counts, exchanger_na->counts,
         sizeof (int) * nmsg);
  team_share_copy->nmsg[SEND] = (int)nsend;
  team_share_copy->nmsg[RECV] = (int)nrecv;
  memcpy(team_share_copy->ranks, team_share->ranks,
         nmsg * sizeof(*team_share_copy->ranks));
#ifndef NEED_MPICH_UNWEIGHTED_COMM_DUP_WORKAROUND
  xt_mpi_call(MPI_Comm_dup(team_share->nb_comm, &team_share_copy->nb_comm),
              new_comm);
#else
  /* MPICH up to version 3.4.2 at least cannot do MPI_Comm_dup
   * for topology communicators */
  enum { no_reorder = 0 }; // no reordering of ranks in new comm
#if __GNUC__ >= 11 && __GNUC__ <= 13
  /* GCC 11 has no means to specify that the special value pointer
   * MPI_UNWEIGHTED does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
  xt_mpi_call(
    MPI_Dist_graph_create_adjacent(
      new_comm, (int)nrecv, team_share_copy->ranks + nsend, MPI_UNWEIGHTED,
      (int)nsend, team_share_copy->ranks, MPI_UNWEIGHTED, MPI_INFO_NULL,
      no_reorder, &team_share_copy->nb_comm), new_comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13
#pragma GCC diagnostic pop
#endif
#endif
  for (size_t i = 0; i < nmsg; ++i)
    xt_mpi_call(MPI_Type_dup(exchanger_na->datatypes[i],
                             exchanger_copy->datatypes + i), new_comm);

  return (Xt_exchanger)exchanger_copy;
}

static void xt_exchanger_neigh_alltoall_delete(Xt_exchanger exchanger) {

  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  struct xt_exchanger_neigh_alltoall_team_share
    *team_share_na = exchanger_na->team_share;

  size_t nmsg = (size_t)team_share_na->nmsg[SEND]
    + (size_t)team_share_na->nmsg[RECV];

  for (size_t i = 0; i < nmsg; ++i)
    xt_mpi_call(MPI_Type_free(exchanger_na->datatypes + i),
                team_share_na->nb_comm);
  if (exchanger_na->team_share == exchanger_na->team_share_)
    xt_exchanger_neigh_alltoall_team_share_destroy(exchanger_na->team_share_);
  free(exchanger_na);
}

static void
xt_exchanger_neigh_alltoall_s_exchange(Xt_exchanger exchanger,
                                       const void *src_data, void *dst_data) {

  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  struct xt_exchanger_neigh_alltoall_team_share
    *team_share = exchanger_na->team_share;

  size_t nsend = (size_t)team_share->nmsg[SEND];
  int *count = exchanger_na->counts;
  xt_mpi_call(
    MPI_Neighbor_alltoallw(src_data, count,
                           exchanger_na->displs, exchanger_na->datatypes,
                           dst_data, count+nsend,
                           exchanger_na->displs, exchanger_na->datatypes +
                           nsend, team_share->nb_comm),
    team_share->nb_comm);
}

static void xt_exchanger_neigh_alltoall_a_exchange(Xt_exchanger exchanger,
                                                   const void * src_data,
                                                   void *dst_data,
                                                   Xt_request *request) {

  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  struct xt_exchanger_neigh_alltoall_team_share
    *team_share = exchanger_na->team_share;

  size_t nsend = (size_t)team_share->nmsg[SEND];
#if defined OMPI_MAJOR_VERSION && OMPI_MAJOR_VERSION == 4 \
  && OMPI_MINOR_VERSION == 0 && OMPI_RELEASE_VERSION == 2
  /* ugly work-around: Open MPI retains pointers to arguments
   * of MPI_Ineighbor_alltoallw, but the exchanger might have
   * been destroyed by the time Open MPI starts to use it.
   * Therefore, this work-around creates an extended request object and
   * stores copies of some arguments in the request object tail.
   * Not pretty, I know, but Open MPI 4.0.x is too recent to ignore. */
  size_t nmsg = nsend+(size_t)team_share->nmsg[RECV],
    displs_size = nmsg * sizeof (MPI_Aint),
    dt_size = nmsg * sizeof (MPI_Datatype),
    count_size = sizeof (int) * nmsg,
    tail_size = displs_size + dt_size + count_size;
  Xt_request request_ = *request
    = xt_request_msgs_ebuf_alloc(
      1, team_share->nb_comm, tail_size, &xt_default_config);
  MPI_Request *request_mpi = xt_request_msgs_ebuf_get_req_ptr(request_);
  MPI_Aint *displs = xt_request_msgs_ebuf_get_extra_buf(request_);
  MPI_Datatype *dt = (void *)(displs+nmsg);
  int *counts = (void *)(dt + nmsg);
  memcpy(displs, exchanger_na->displs, displs_size);
  memcpy(dt, exchanger_na->datatypes, dt_size);
  memcpy(counts, exchanger_na->counts, count_size);
#else
  Xt_request request_ = *request
    = xt_request_msgs_alloc(1, team_share->nb_comm, &xt_default_config);
  MPI_Request *request_mpi = xt_request_msgs_get_req_ptr(request_);
  MPI_Aint *displs = exchanger_na->displs;
  MPI_Datatype *dt = exchanger_na->datatypes;
  int *counts = exchanger_na->counts;
#endif
  xt_mpi_call(
    MPI_Ineighbor_alltoallw(src_data, counts, displs, dt,
                            dst_data, counts+nsend, displs, dt+nsend,
                            team_share->nb_comm, request_mpi),
    team_share->nb_comm);
}

static MPI_Datatype
xt_exchanger_neigh_alltoall_get_MPI_Datatype(Xt_exchanger exchanger,
                                             int rank,
                                             enum xt_msg_direction direction,
                                             bool do_dup)
{
  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  struct xt_exchanger_neigh_alltoall_team_share
    *team_share = exchanger_na->team_share;
  size_t nsend = (size_t)team_share->nmsg[SEND],
    nmsg = (size_t)team_share->nmsg[direction],
    ofs = direction == xt_direction_send ? 0 : nsend;
  int *restrict ranks = team_share->ranks + ofs;
  MPI_Datatype datatype_copy = MPI_DATATYPE_NULL;
  for (size_t i = 0; i < nmsg; ++i) {
    if (ranks[i] == rank) {
      if (do_dup)
        xt_mpi_call(MPI_Type_dup(exchanger_na->datatypes[i+ofs],
                                 &datatype_copy), team_share->nb_comm);
      else
        datatype_copy = exchanger_na->datatypes[i+ofs];
      break;
    }
  }
  return datatype_copy;
}

static int
xt_exchanger_neigh_alltoall_get_num_msg(Xt_exchanger exchanger,
                                        int direction)
{
  struct xt_exchanger_neigh_alltoall_team_share
    *team_share = ((Xt_exchanger_neigh_alltoall)exchanger)->team_share;
    int nmsg;
  assert(direction >= 0 && direction <= 2);
  switch (direction) {
  case xt_direction_send:
  case xt_direction_recv:
    nmsg = team_share->nmsg[direction];
    break;
  case 2:
    nmsg = team_share->nmsg[0] + team_share->nmsg[1];
    break;
#  if HAVE_DECL___BUILTIN_UNREACHABLE
  default:
    __builtin_unreachable();
#  endif
  }
  return nmsg;
}

static int
xt_exchanger_neigh_alltoall_get_msg_ranks(Xt_exchanger exchanger,
                                          enum xt_msg_direction direction,
                                          int *restrict *ranks)
{
  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  struct xt_exchanger_neigh_alltoall_team_share
    *team_share = exchanger_na->team_share;
  size_t nsend = (size_t)team_share->nmsg[SEND],
    nmsg = (size_t)team_share->nmsg[direction],
    ofs = direction == xt_direction_send ? 0 : nsend;
  int *ranks_ = *ranks;
  if (!ranks_)
    ranks_ = *ranks = xmalloc(nmsg * sizeof(**ranks));
  memcpy(ranks_, team_share->ranks + ofs, nmsg * sizeof(**ranks));
  return (int)nmsg;
}

typedef struct Xt_exchanger_na2a_iter_ *Xt_exchanger_na2a_iter;

struct Xt_exchanger_na2a_iter_ {
  const struct xt_exchanger_iter_vtable *vtable;
  Xt_exchanger_neigh_alltoall exchanger;
  unsigned msg_ofs;
  uint32_t config_flags;
};

static int
xt_exchanger_na2a_iter_next(Xt_exchanger_iter);
static void
xt_exchanger_na2a_iter_delete(Xt_exchanger_iter);
static MPI_Datatype
xt_exchanger_na2a_iter_get_MPI_Datatype(Xt_exchanger_iter iter);
static int
xt_exchanger_na2a_iter_get_rank(Xt_exchanger_iter iter);
static int
xt_exchanger_na2a_iter_get_count(Xt_exchanger_iter iter);
static int
xt_exchanger_na2a_iter_get_direction(Xt_exchanger_iter iter);
static MPI_Aint
xt_exchanger_na2a_iter_get_displacement(Xt_exchanger_iter iter);

static const struct xt_exchanger_iter_vtable
xt_exchanger_na2a_iterator_vtable = {
  .next = xt_exchanger_na2a_iter_next,
  .get_MPI_Datatype = xt_exchanger_na2a_iter_get_MPI_Datatype,
  .get_rank = xt_exchanger_na2a_iter_get_rank,
  .get_count = xt_exchanger_na2a_iter_get_count,
  .get_direction = xt_exchanger_na2a_iter_get_direction,
  .get_displacement = xt_exchanger_na2a_iter_get_displacement,
  .delete = xt_exchanger_na2a_iter_delete,
};


static Xt_exchanger_iter
xt_exchanger_na2a_get_iterator(Xt_exchanger exchanger, Xt_config config)
{
  Xt_exchanger_neigh_alltoall exchanger_na2a
    = (Xt_exchanger_neigh_alltoall)exchanger;
  Xt_exchanger_na2a_iter iter_na2a = xmalloc(sizeof (*iter_na2a));
  iter_na2a->vtable = &xt_exchanger_na2a_iterator_vtable;
  iter_na2a->exchanger = exchanger_na2a;
  iter_na2a->msg_ofs = 0;
  iter_na2a->config_flags = config->flags;
  return (Xt_exchanger_iter)iter_na2a;
}

static void
xt_exchanger_na2a_iter_delete(Xt_exchanger_iter iter)
{
  free(iter);
}

static int
xt_exchanger_na2a_iter_next(Xt_exchanger_iter iter)
{
  Xt_exchanger_na2a_iter iter_na2a = (Xt_exchanger_na2a_iter)iter;
  Xt_exchanger_neigh_alltoall exchanger_na2a = iter_na2a->exchanger;
  unsigned msg_ofs = iter_na2a->msg_ofs,
    nmsg = (unsigned)exchanger_na2a->team_share->nmsg[0]
    + (unsigned)exchanger_na2a->team_share->nmsg[1];
  if (msg_ofs+1 <= nmsg) {
    iter_na2a->msg_ofs = ++msg_ofs;
    return msg_ofs < nmsg;
  } else
    return 0;
}

static MPI_Datatype
xt_exchanger_na2a_iter_get_MPI_Datatype(Xt_exchanger_iter iter)
{
  Xt_exchanger_na2a_iter iter_na2a = (Xt_exchanger_na2a_iter)iter;
  Xt_exchanger_neigh_alltoall exchanger_na2a = iter_na2a->exchanger;
  unsigned msg_ofs = iter_na2a->msg_ofs,
    nmsg = (unsigned)exchanger_na2a->team_share->nmsg[0]
    + (unsigned)exchanger_na2a->team_share->nmsg[1];
  MPI_Datatype datatype_copy = MPI_DATATYPE_NULL;
  if (msg_ofs < nmsg) {
    if (iter_na2a->config_flags & exch_no_dt_dup)
      datatype_copy = exchanger_na2a->datatypes[msg_ofs];
    else
      xt_mpi_call(MPI_Type_dup(exchanger_na2a->datatypes[msg_ofs],
                               &datatype_copy),
                  exchanger_na2a->team_share->nb_comm);
  }
  return datatype_copy;
}

static int
xt_exchanger_na2a_iter_get_rank(Xt_exchanger_iter iter)
{
  Xt_exchanger_na2a_iter iter_na2a = (Xt_exchanger_na2a_iter)iter;
  struct xt_exchanger_neigh_alltoall_team_share *team_share
    = iter_na2a->exchanger->team_share;
  unsigned msg_ofs = iter_na2a->msg_ofs,
    nmsg = (unsigned)team_share->nmsg[0]
    + (unsigned)team_share->nmsg[1];
  return msg_ofs < nmsg ? team_share->ranks[msg_ofs] : -1;
}

static int
xt_exchanger_na2a_iter_get_count(Xt_exchanger_iter iter)
{
  Xt_exchanger_na2a_iter iter_na2a = (Xt_exchanger_na2a_iter)iter;
  Xt_exchanger_neigh_alltoall exchanger_na2a = iter_na2a->exchanger;
  unsigned msg_ofs = iter_na2a->msg_ofs,
    nmsg = (unsigned)exchanger_na2a->team_share->nmsg[0]
    + (unsigned)exchanger_na2a->team_share->nmsg[1];
  return msg_ofs < nmsg ? exchanger_na2a->counts[msg_ofs] : -1;
}

static int
xt_exchanger_na2a_iter_get_direction(Xt_exchanger_iter iter)
{
  Xt_exchanger_na2a_iter iter_na2a = (Xt_exchanger_na2a_iter)iter;
  struct xt_exchanger_neigh_alltoall_team_share *team_share
    = iter_na2a->exchanger->team_share;
  unsigned msg_ofs = iter_na2a->msg_ofs,
    nmsg = (unsigned)team_share->nmsg[0]
    + (unsigned)team_share->nmsg[1];
  if (msg_ofs < nmsg)
    return msg_ofs >= (unsigned)team_share->nmsg[0];
  else
    return -1;
}

static MPI_Aint
xt_exchanger_na2a_iter_get_displacement(Xt_exchanger_iter iter)
{
  Xt_exchanger_na2a_iter iter_na2a = (Xt_exchanger_na2a_iter)iter;
  Xt_exchanger_neigh_alltoall exchanger_na2a = iter_na2a->exchanger;
  unsigned msg_ofs = iter_na2a->msg_ofs,
    nmsg = (unsigned)exchanger_na2a->team_share->nmsg[0]
    + (unsigned)exchanger_na2a->team_share->nmsg[1];
  return msg_ofs < nmsg ? exchanger_na2a->displs[msg_ofs] : -1;
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
