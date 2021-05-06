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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <string.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt_redist_internal.h"
#include "xt/xt_xmap.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_request.h"
#include "xt/xt_request_msgs.h"
#include "xt_exchanger.h"
#include "xt_exchanger_neigh_alltoall.h"

#define MAX(a,b) ((a) >= (b) ? (a) : (b))

static Xt_exchanger
xt_exchanger_neigh_alltoall_copy(Xt_exchanger exchanger,
                                 MPI_Comm newComm, int new_tag_offset);
static void xt_exchanger_neigh_alltoall_delete(Xt_exchanger exchanger);
static void xt_exchanger_neigh_alltoall_s_exchange(Xt_exchanger exchanger,
                                                   const void * src_data,
                                                   void * dst_data);
static void xt_exchanger_neigh_alltoall_a_exchange(Xt_exchanger exchanger,
                                                   const void * src_data,
                                                   void * dst_data,
                                                   Xt_request *request);
static int
xt_exchanger_neigh_alltoall_get_msg_ranks(Xt_exchanger exchanger,
                                          enum xt_msg_direction direction,
                                          int *restrict *ranks);

static MPI_Datatype
xt_exchanger_neigh_alltoall_get_MPI_Datatype(Xt_exchanger exchanger,
                                             int rank,
                                             enum xt_msg_direction direction);


static const struct xt_exchanger_vtable exchanger_neigh_alltoall_vtable = {
  .copy = xt_exchanger_neigh_alltoall_copy,
  .delete = xt_exchanger_neigh_alltoall_delete,
  .s_exchange = xt_exchanger_neigh_alltoall_s_exchange,
  .a_exchange = xt_exchanger_neigh_alltoall_a_exchange,
  .get_msg_ranks = xt_exchanger_neigh_alltoall_get_msg_ranks,
  .get_MPI_Datatype = xt_exchanger_neigh_alltoall_get_MPI_Datatype,
};

typedef struct Xt_exchanger_neigh_alltoall_ * Xt_exchanger_neigh_alltoall;

struct Xt_exchanger_neigh_alltoall_ {

  const struct xt_exchanger_vtable * vtable;

  int nmsg[2];
  int tag_offset;
  MPI_Comm comm;
  int * ranks;
  int * one_counts;
  MPI_Aint * displs;
  MPI_Datatype * datatypes;
};

static Xt_exchanger_neigh_alltoall
xt_exchanger_neigh_alltoall_alloc(size_t nsend, size_t nrecv)
{
  size_t nmsg = nsend + nrecv;
  size_t max_msgs = MAX(nsend, nrecv);
  Xt_exchanger_neigh_alltoall exchanger = xmalloc(1 * sizeof(*exchanger));
  exchanger->ranks = xmalloc(nmsg * sizeof(*(exchanger->ranks)));
  exchanger->datatypes = xmalloc(nmsg * sizeof(*(exchanger->datatypes)));
  exchanger->one_counts = xmalloc(max_msgs * sizeof(*(exchanger->one_counts)));
  exchanger->displs = xmalloc(max_msgs * sizeof(*(exchanger->displs)));
  exchanger->vtable = &exchanger_neigh_alltoall_vtable;
  for (size_t i = 0; i < max_msgs; ++i) {
    exchanger->one_counts[i] = 1;
    exchanger->displs[i] = 0;
  }
  return exchanger;
}

static void copy_from_redist_msgs(size_t n,
                                  const struct Xt_redist_msg *restrict msgs,
                                  int *restrict ranks,
                                  MPI_Datatype *restrict datatypes,
                                  MPI_Comm comm) {

  for (size_t i = 0; i < n; ++i) {
    ranks[i] = msgs[i].rank;
    xt_mpi_call(MPI_Type_dup(msgs[i].datatype, datatypes + i), comm);
  }
}

Xt_exchanger
xt_exchanger_neigh_alltoall_new(int nsend, int nrecv,
                                const struct Xt_redist_msg *send_msgs,
                                const struct Xt_redist_msg *recv_msgs,
                                MPI_Comm comm, int tag_offset) {

  /** note: tag_offset + xt_mpi_tag_exchange_msg must not
   *        be used on @a comm by any other part of the program during the
   *        lifetime of the created exchanger object
   */

  int flag;
  xt_mpi_call(MPI_Comm_test_inter(comm, &flag), comm);
  if (flag)
    Xt_abort(comm, "ERROR(xt_exchanger_neigh_alltoall_new): "
             "inter-communicator's are not defined for virtual topologies",
             __FILE__, __LINE__);

  assert((nsend >= 0) & (nrecv >= 0));
  Xt_exchanger_neigh_alltoall exchanger
    = xt_exchanger_neigh_alltoall_alloc((size_t)nsend, (size_t)nrecv);
  exchanger->tag_offset = tag_offset;
  exchanger->nmsg[SEND] = nsend;
  copy_from_redist_msgs((size_t)nsend, send_msgs, exchanger->ranks,
                        exchanger->datatypes, comm);
  exchanger->nmsg[RECV] = nrecv;
  copy_from_redist_msgs((size_t)nrecv, recv_msgs, exchanger->ranks + nsend,
                        exchanger->datatypes + nsend, comm);

  int reorder = 0; // no reordering of ranks in new comm
  xt_mpi_call(
    MPI_Dist_graph_create_adjacent(
      comm, nrecv, exchanger->ranks + nsend, MPI_UNWEIGHTED, nsend,
      exchanger->ranks, MPI_UNWEIGHTED, MPI_INFO_NULL, reorder,
      &(exchanger->comm)), comm);

  return (Xt_exchanger)exchanger;
}

static Xt_exchanger
xt_exchanger_neigh_alltoall_copy(Xt_exchanger exchanger,
                                 MPI_Comm new_comm, int new_tag_offset)
{
  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  size_t nsend = (size_t)(exchanger_na->nmsg[SEND]),
    nrecv = (size_t)(exchanger_na->nmsg[RECV]),
    nmsg = nsend + nrecv;

  Xt_exchanger_neigh_alltoall
    exchanger_copy = xt_exchanger_neigh_alltoall_alloc(nsend, nrecv);

  exchanger_copy->nmsg[SEND] = (int)nsend;
  exchanger_copy->nmsg[RECV] = (int)nrecv;
  exchanger_copy->tag_offset = new_tag_offset;
  exchanger_copy->comm = new_comm;
  memcpy(exchanger_copy->ranks, exchanger_na->ranks,
         nmsg * sizeof(*(exchanger_copy->ranks)));
  for (size_t i = 0; i < nmsg; ++i)
    xt_mpi_call(MPI_Type_dup(exchanger_na->datatypes[i],
                             exchanger_copy->datatypes + i), new_comm);

  return (Xt_exchanger)exchanger_copy;
}

static void xt_exchanger_neigh_alltoall_delete(Xt_exchanger exchanger) {

  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;

  size_t nmsg = (size_t)exchanger_na->nmsg[SEND]
    + (size_t)exchanger_na->nmsg[RECV];
  MPI_Comm comm = exchanger_na->comm;

  free(exchanger_na->ranks);
  free(exchanger_na->one_counts);
  free(exchanger_na->displs);
  for (size_t i = 0; i < nmsg; ++i) {
    MPI_Datatype *dt = exchanger_na->datatypes + i;
    if (*dt != MPI_DATATYPE_NULL)
      xt_mpi_call(MPI_Type_free(dt), comm);
  }
  free(exchanger_na->datatypes);
  xt_mpi_call(MPI_Comm_free(&(exchanger_na->comm)), Xt_default_comm);
  free(exchanger_na);
}

static void xt_exchanger_neigh_alltoall_s_exchange(Xt_exchanger exchanger,
                                                   const void * src_data,
                                                   void * dst_data) {

  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;

  xt_mpi_call(
    MPI_Neighbor_alltoallw(src_data, exchanger_na->one_counts,
                           exchanger_na->displs, exchanger_na->datatypes,
                           dst_data, exchanger_na->one_counts,
                           exchanger_na->displs, exchanger_na->datatypes +
                           (size_t)(exchanger_na->nmsg[SEND]),
                           exchanger_na->comm),
    exchanger_na->comm);
}

static void xt_exchanger_neigh_alltoall_a_exchange(Xt_exchanger exchanger,
                                                   const void * src_data,
                                                   void * dst_data,
                                                   Xt_request *request) {

  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;

  MPI_Request tmp_request;

  xt_mpi_call(
    MPI_Ineighbor_alltoallw(src_data, exchanger_na->one_counts,
                            exchanger_na->displs, exchanger_na->datatypes,
                            dst_data, exchanger_na->one_counts,
                            exchanger_na->displs, exchanger_na->datatypes +
                            (size_t)(exchanger_na->nmsg[SEND]),
                            exchanger_na->comm, &tmp_request),
    exchanger_na->comm);

  *request = xt_request_msgs_new(1, &tmp_request, exchanger_na->comm);
}

static MPI_Datatype
xt_exchanger_neigh_alltoall_get_MPI_Datatype(Xt_exchanger exchanger,
                                             int rank,
                                             enum xt_msg_direction direction)
{
  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  size_t nsend = (size_t)exchanger_na->nmsg[SEND],
    nmsg = (size_t)exchanger_na->nmsg[direction],
    ofs = direction == SEND ? 0 : nsend;
  int *restrict ranks = exchanger_na->ranks + ofs;
  MPI_Datatype datatype_copy = MPI_DATATYPE_NULL;
  for (size_t i = 0; i < nmsg; ++i) {
    if (ranks[i] == rank) {
      xt_mpi_call(MPI_Type_dup(exchanger_na->datatypes[i+ofs], &datatype_copy),
                  exchanger_na->comm);
      break;
    }
  }
  return datatype_copy;
}

static int
xt_exchanger_neigh_alltoall_get_msg_ranks(Xt_exchanger exchanger,
                                          enum xt_msg_direction direction,
                                          int *restrict *ranks)
{
  Xt_exchanger_neigh_alltoall exchanger_na =
    (Xt_exchanger_neigh_alltoall)exchanger;
  size_t nsend = (size_t)exchanger_na->nmsg[SEND],
    nmsg = (size_t)exchanger_na->nmsg[direction],
    ofs = direction == SEND ? 0 : nsend;
  int *ranks_ = *ranks;
  if (!ranks_)
    ranks_ = *ranks = xmalloc(nmsg * sizeof(**ranks));
  memcpy(ranks_, exchanger_na->ranks + ofs, nmsg * sizeof(**ranks));
  return (int)nmsg;
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
