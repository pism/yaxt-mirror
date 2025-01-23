/**
 * @file xt_exchanger_simple_base.c
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
#include <mpi.h>
#include <stdbool.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_config_internal.h"
#include "xt_mpi_internal.h"
#include "xt_exchanger.h"
#include "xt_exchanger_internal.h"
#include "xt_exchanger_simple_base.h"

static const char filename[] = "xt_exchanger_simple_base.c";

static Xt_exchanger
xt_exchanger_simple_base_copy(Xt_exchanger exchanger,
                              MPI_Comm newComm, int new_tag_offset);
static void xt_exchanger_simple_base_delete(Xt_exchanger exchanger);
static void
xt_exchanger_simple_base_s_exchange(Xt_exchanger exchanger,
                                    const void *src_data,
                                    void *dst_data);
static void
xt_exchanger_simple_base_a_exchange(Xt_exchanger exchanger,
                                    const void *src_data, void *dst_data,
                                    Xt_request *request);
static int
xt_exchanger_simple_base_get_msg_ranks(Xt_exchanger exchanger,
                                       enum xt_msg_direction direction,
                                       int *restrict *ranks);

static int
xt_exchanger_simple_base_get_num_msg(Xt_exchanger exchanger,
                                     int direction);

static MPI_Datatype
xt_exchanger_simple_base_get_MPI_Datatype(Xt_exchanger exchanger,
                                          int rank,
                                          enum xt_msg_direction direction,
                                          bool do_dup);

static Xt_exchanger_thread_share
xt_exchanger_simple_base_create_thread_share(Xt_exchanger exchanger,
                                             Xt_config config);

static void
xt_exchanger_simple_base_destroy_thread_share(
  Xt_exchanger exchanger, Xt_exchanger_thread_share thread_share);

static Xt_exchanger_iter
xt_exchanger_simple_base_get_iterator(Xt_exchanger exchanger, Xt_config config);

const struct xt_exchanger_vtable xt_exchanger_simple_base_vtable = {
  .copy = xt_exchanger_simple_base_copy,
  .delete = xt_exchanger_simple_base_delete,
  .s_exchange = xt_exchanger_simple_base_s_exchange,
  .a_exchange = xt_exchanger_simple_base_a_exchange,
  .get_num_msg = xt_exchanger_simple_base_get_num_msg,
  .get_msg_ranks = xt_exchanger_simple_base_get_msg_ranks,
  .get_MPI_Datatype = xt_exchanger_simple_base_get_MPI_Datatype,
  .create_thread_share = xt_exchanger_simple_base_create_thread_share,
  .destroy_thread_share = xt_exchanger_simple_base_destroy_thread_share,
  .get_iterator = xt_exchanger_simple_base_get_iterator,
};

typedef struct Xt_exchanger_simple_base_ * Xt_exchanger_simple_base;

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

static Xt_exchanger_simple_base
xt_exchanger_simple_base_alloc(size_t nmsg)
{
  Xt_exchanger_simple_base exchanger;
  size_t header_size = sizeof(*exchanger),
    body_size = nmsg * sizeof (exchanger->msgs[0]);
  exchanger = xmalloc(header_size + body_size);
  exchanger->vtable = &xt_exchanger_simple_base_vtable;
  return exchanger;
}

static inline int
adjusted_rank(int r, int comm_rank, int comm_size)
{
  return r + (r <= comm_rank ? comm_size : 0);
}

#define XT_SORTFUNC_DECL static
#define SORT_TYPE struct Xt_msg_param
#define SORT_TYPE_SUFFIX msg_param
#define SORT_TYPE_CMP_LT(u,v,i,j)                       \
  (adjusted_rank((u).rank, comm_rank, comm_size)        \
   < adjusted_rank((v).rank, comm_rank, comm_size))
#define SORT_TYPE_CMP_LE(u,v,i,j)                       \
  (adjusted_rank((u).rank, comm_rank, comm_size)        \
   <= adjusted_rank((v).rank, comm_rank, comm_size))
#define SORT_TYPE_CMP_EQ(u,v,i,j) (((u).rank) == ((v).rank))
#define XT_SORT_EXTRA_ARGS_DECL , int comm_rank, int comm_size
#define XT_SORT_EXTRA_ARGS_PASS , comm_rank, comm_size
#define XT_SORT_VECSWAP_EXTRA_ARGS_DECL
#define XT_SORT_VECSWAP_EXTRA_ARGS_PASS

#include "xt_quicksort_base.h"

Xt_exchanger
xt_exchanger_simple_base_new(
  int nsend, int nrecv,
  const struct Xt_msg_param *send_msgs,
  const struct Xt_msg_param *recv_msgs,
  MPI_Comm comm, int tag_offset,
  xt_simple_s_exchange_func s_func,
  xt_simple_a_exchange_func a_func,
  xt_simple_create_thread_share_func create_thread_share_func,
  xt_simple_destroy_thread_share_func destroy_thread_share_func,
  Xt_config config)
{
  /** note: tag_offset + xt_mpi_tag_exchange_msg must not
   *        be used on @a comm by any other part of the program during the
   *        lifetime of the created exchanger object
   */

  if (s_func == 0)
    Xt_abort(comm, "ERROR(xt_exchanger_simple_base_new): invalid synchronous "
             "exchange function pointer", filename, __LINE__);

  assert((nsend >= 0) & (nrecv >= 0));
  size_t nmsg = (size_t)nsend + (size_t)nrecv;
  Xt_exchanger_simple_base exchanger
    = xt_exchanger_simple_base_alloc(nmsg);
  exchanger->comm = comm;
  exchanger->tag_offset = tag_offset;
  exchanger->nmsg[SEND] = nsend;
  exchanger->nmsg[RECV] = nrecv;
  bool dt_dup = !(config->flags & exch_no_dt_dup);
  xt_msg_params_copy((size_t)nsend, send_msgs, exchanger->msgs,
                    comm, dt_dup);
  xt_msg_params_copy((size_t)nrecv, recv_msgs, exchanger->msgs + nsend,
                     comm, dt_dup);
  exchanger->s_func = s_func;
  exchanger->a_func = a_func;
  exchanger->create_thread_share_func = create_thread_share_func;
  exchanger->destroy_thread_share_func = destroy_thread_share_func;

  {
    int comm_size, comm_rank, is_inter;
    xt_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
    xt_mpi_call(MPI_Comm_test_inter(comm, &is_inter), comm);
    int (*get_comm_size)(MPI_Comm, int *)
      = is_inter ? MPI_Comm_remote_size : MPI_Comm_size;
    xt_mpi_call(get_comm_size(comm, &comm_size), comm);
    /* In order to avoid congestion of messages, the order of send and
     * receive messages is changed. This is done by sorting the
     * messages according to the rank of the respective message
     * partner. Before the sorting to ranks that are smaller or equal
     * to the local rank the size of the communicator is added.
     *
     * example: process 5 is supposed to communicate with
     * processes: 9, 5, 2, 6, 1
     * 1. add comm_size(10): 9, 15, 12, 6, 11
     * 2. sort: 6, 9, 11, 12, 15
     * 3. subtrace comm_size again -> final order: 6, 9, 1, 2, 5
     */
    xt_quicksort_msg_param(exchanger->msgs, (size_t)nsend,
                           comm_rank, comm_size);
    xt_quicksort_msg_param(exchanger->msgs + (size_t)nsend, (size_t)nrecv,
                           comm_rank, comm_size);
  }

  return (Xt_exchanger)exchanger;
}

static Xt_exchanger
xt_exchanger_simple_base_copy(Xt_exchanger exchanger,
                              MPI_Comm new_comm, int new_tag_offset)
{
  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;
  int nsend = exchanger_sb->nmsg[SEND],
    nrecv = exchanger_sb->nmsg[RECV];
  size_t nmsg = (size_t)nsend + (size_t)nrecv;
  Xt_exchanger_simple_base exchanger_copy
    = xt_exchanger_simple_base_alloc(nmsg);
  exchanger_copy->nmsg[SEND] = nsend;
  exchanger_copy->nmsg[RECV] = nrecv;
  exchanger_copy->s_func = exchanger_sb->s_func;
  exchanger_copy->a_func = exchanger_sb->a_func;
  exchanger_copy->create_thread_share_func = exchanger_sb->create_thread_share_func;
  struct Xt_msg_param *restrict new_msgs = exchanger_copy->msgs,
    *restrict orig_msgs = exchanger_sb->msgs;
  xt_msg_params_copy(nmsg, orig_msgs, new_msgs, new_comm, true);
  exchanger_copy->comm = new_comm;
  exchanger_copy->tag_offset = new_tag_offset;
  return (Xt_exchanger)exchanger_copy;
}


static void xt_exchanger_simple_base_delete(Xt_exchanger exchanger) {

  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;

  size_t nmsg = (size_t)exchanger_sb->nmsg[SEND] + (size_t)exchanger_sb->nmsg[RECV];
  struct Xt_msg_param *restrict msgs = exchanger_sb->msgs;
  xt_msg_params_destruct(nmsg, msgs, exchanger_sb->comm);
  free(exchanger_sb);
}

static void
xt_exchanger_simple_base_s_exchange(Xt_exchanger exchanger,
                                    const void *src_data, void *dst_data)
{
  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;

  int nsend = exchanger_sb->nmsg[SEND];
  exchanger_sb->s_func(src_data, dst_data, nsend,
                       exchanger_sb->nmsg[RECV], exchanger_sb->msgs,
                       exchanger_sb->msgs + (size_t)nsend,
                       exchanger_sb->tag_offset, exchanger_sb->comm);
}

static void
xt_exchanger_simple_base_a_exchange(Xt_exchanger exchanger,
                                    const void *src_data, void *dst_data,
                                    Xt_request *request)
{
  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;

  if (exchanger_sb->a_func == NULL)
    Xt_abort(exchanger_sb->comm, "ERROR(xt_exchanger_simple_base_a_exchange): "
             "asynchronous exchange function not defined for current exchanger",
             filename, __LINE__);

  int nsend = exchanger_sb->nmsg[SEND];
  exchanger_sb->a_func(src_data, dst_data, nsend,
                       exchanger_sb->nmsg[RECV], exchanger_sb->msgs,
                       exchanger_sb->msgs + (size_t)nsend,
                       exchanger_sb->tag_offset, exchanger_sb->comm, request);
}

static MPI_Datatype
xt_exchanger_simple_base_get_MPI_Datatype(Xt_exchanger exchanger,
                                          int rank,
                                          enum xt_msg_direction direction,
                                          bool do_dup)
{
  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;
  size_t nsend = (size_t)exchanger_sb->nmsg[SEND],
    nmsg = (size_t)exchanger_sb->nmsg[direction],
    ofs = direction == xt_direction_send ? 0 : nsend;
  struct Xt_msg_param *restrict msgs = exchanger_sb->msgs + ofs;
  MPI_Datatype datatype_copy = MPI_DATATYPE_NULL;
  for (size_t i = 0; i < nmsg; ++i)
    if (msgs[i].rank == rank) {
      if (do_dup)
        xt_mpi_call(MPI_Type_dup(msgs[i].datatype, &datatype_copy),
                    exchanger_sb->comm);
      else
        datatype_copy = msgs[i].datatype;
      break;
    }
  return datatype_copy;
}

static int
xt_exchanger_simple_base_get_num_msg(Xt_exchanger exchanger,
                                     int direction)
{
  Xt_exchanger_simple_base exchanger_sb = (Xt_exchanger_simple_base)exchanger;
  int nmsg;
  assert(direction >= 0 && direction <= 2);
  switch (direction) {
  case xt_direction_send:
  case xt_direction_recv:
    nmsg = exchanger_sb->nmsg[direction];
    break;
  case 2:
    nmsg = exchanger_sb->nmsg[0] + exchanger_sb->nmsg[1];
    break;
#  if HAVE_DECL___BUILTIN_UNREACHABLE
  default:
    __builtin_unreachable();
#  endif
  }
  return nmsg;
}

static int
xt_exchanger_simple_base_get_msg_ranks(Xt_exchanger exchanger,
                                       enum xt_msg_direction direction,
                                       int *restrict *ranks)
{
  Xt_exchanger_simple_base exchanger_sb = (Xt_exchanger_simple_base)exchanger;
  size_t nmsg = (size_t)exchanger_sb->nmsg[direction];
  struct Xt_msg_param *restrict msgs = exchanger_sb->msgs
    + (direction == xt_direction_recv ? (size_t)exchanger_sb->nmsg[SEND] : 0);
  int *restrict ranks_ = *ranks;
  if (!ranks_)
    ranks_ = *ranks = xmalloc(nmsg * sizeof (*ranks_));
  for (size_t i = 0; i < nmsg; ++i)
    ranks_[i] = msgs[i].rank;
  return (int)nmsg;
}

static Xt_exchanger_thread_share
xt_exchanger_simple_base_create_thread_share(Xt_exchanger exchanger,
                                             Xt_config config)
{
  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;
  int nsend = exchanger_sb->nmsg[SEND];
  if (exchanger_sb->create_thread_share_func)
    return exchanger_sb->create_thread_share_func(
      nsend, exchanger_sb->nmsg[RECV],
      exchanger_sb->msgs, exchanger_sb->msgs + nsend,
      exchanger_sb->comm, config);
  else
    return NULL;
}

static void
xt_exchanger_simple_base_destroy_thread_share(
  Xt_exchanger exchanger, Xt_exchanger_thread_share thread_share)
{
  Xt_exchanger_simple_base exchanger_sb =
    (Xt_exchanger_simple_base)exchanger;
  /* assert that thread_share contains no in-flight requests */
  if (exchanger_sb->destroy_thread_share_func)
    exchanger_sb->destroy_thread_share_func(exchanger, thread_share);
}


typedef struct Xt_exchanger_simple_base_iter_ *Xt_exchanger_simple_base_iter;

struct Xt_exchanger_simple_base_iter_ {
  const struct xt_exchanger_iter_vtable *vtable;
  Xt_exchanger_simple_base exchanger;
  unsigned msg_ofs;
  uint32_t config_flags;
};

static int
xt_exchanger_simple_base_iter_next(Xt_exchanger_iter iter);
static MPI_Datatype
xt_exchanger_simple_base_iter_get_MPI_Datatype(Xt_exchanger_iter iter);
static int
xt_exchanger_simple_base_iter_get_rank(Xt_exchanger_iter iter);
static int
xt_exchanger_simple_base_iter_get_count(Xt_exchanger_iter iter);
static int
xt_exchanger_simple_base_iter_get_direction(Xt_exchanger_iter iter);
static MPI_Aint
xt_exchanger_simple_base_iter_get_displacement(Xt_exchanger_iter iter);
static void
xt_exchanger_simple_base_iter_delete(Xt_exchanger_iter iter);

static const struct xt_exchanger_iter_vtable
xt_exchanger_simple_base_iterator_vtable = {
  .next = xt_exchanger_simple_base_iter_next,
  .get_MPI_Datatype = xt_exchanger_simple_base_iter_get_MPI_Datatype,
  .get_rank = xt_exchanger_simple_base_iter_get_rank,
  .get_count = xt_exchanger_simple_base_iter_get_count,
  .get_direction = xt_exchanger_simple_base_iter_get_direction,
  .get_displacement = xt_exchanger_simple_base_iter_get_displacement,
  .delete = xt_exchanger_simple_base_iter_delete,
};


static Xt_exchanger_iter
xt_exchanger_simple_base_get_iterator(Xt_exchanger exchanger, Xt_config config)
{
  Xt_exchanger_simple_base exchanger_sb = (Xt_exchanger_simple_base)exchanger;
  Xt_exchanger_simple_base_iter iter_sb = xmalloc(sizeof (*iter_sb));
  iter_sb->vtable = &xt_exchanger_simple_base_iterator_vtable;
  iter_sb->exchanger = exchanger_sb;
  iter_sb->msg_ofs = 0;
  iter_sb->config_flags = config->flags;
  return (Xt_exchanger_iter)iter_sb;
}


static int
xt_exchanger_simple_base_iter_next(Xt_exchanger_iter iter)
{
  Xt_exchanger_simple_base_iter iter_sb = (Xt_exchanger_simple_base_iter)iter;
  unsigned msg_ofs = iter_sb->msg_ofs,
    nmsg = (unsigned)iter_sb->exchanger->nmsg[0]
    + (unsigned)iter_sb->exchanger->nmsg[1];
  if (msg_ofs+1 <= nmsg) {
    iter_sb->msg_ofs = ++msg_ofs;
    return msg_ofs < nmsg;
  } else
    return 0;
}

static MPI_Datatype
xt_exchanger_simple_base_iter_get_MPI_Datatype(Xt_exchanger_iter iter)
{
  Xt_exchanger_simple_base_iter iter_sb = (Xt_exchanger_simple_base_iter)iter;
  const struct Xt_msg_param *msgs = iter_sb->exchanger->msgs;
  unsigned msg_ofs = iter_sb->msg_ofs,
    nmsg = (unsigned)iter_sb->exchanger->nmsg[0]
    + (unsigned)iter_sb->exchanger->nmsg[1];
  MPI_Datatype datatype_copy = MPI_DATATYPE_NULL;
  if (msg_ofs < nmsg) {
    if (iter_sb->config_flags & exch_no_dt_dup)
      datatype_copy = msgs[msg_ofs].datatype;
    else
      xt_mpi_call(MPI_Type_dup(msgs[msg_ofs].datatype, &datatype_copy),
                  iter_sb->exchanger->comm);
  }
  return datatype_copy;
}

static int
xt_exchanger_simple_base_iter_get_rank(Xt_exchanger_iter iter)
{
  Xt_exchanger_simple_base_iter iter_sb = (Xt_exchanger_simple_base_iter)iter;
  const struct Xt_msg_param *msgs = iter_sb->exchanger->msgs;
  unsigned msg_ofs = iter_sb->msg_ofs,
    nmsg = (unsigned)iter_sb->exchanger->nmsg[0]
    + (unsigned)iter_sb->exchanger->nmsg[1];
  return msg_ofs < nmsg ? msgs[msg_ofs].rank : -1;
}

static int
xt_exchanger_simple_base_iter_get_count(Xt_exchanger_iter iter)
{
  Xt_exchanger_simple_base_iter iter_sb = (Xt_exchanger_simple_base_iter)iter;
  const struct Xt_msg_param *msgs = iter_sb->exchanger->msgs;
  unsigned msg_ofs = iter_sb->msg_ofs,
    nmsg = (unsigned)iter_sb->exchanger->nmsg[0]
    + (unsigned)iter_sb->exchanger->nmsg[1];
  return msg_ofs < nmsg ? msgs[msg_ofs].count : -1;
}

static int
xt_exchanger_simple_base_iter_get_direction(Xt_exchanger_iter iter)
{
  Xt_exchanger_simple_base_iter iter_sb = (Xt_exchanger_simple_base_iter)iter;
  unsigned msg_ofs = iter_sb->msg_ofs,
    nmsg = (unsigned)iter_sb->exchanger->nmsg[0]
    + (unsigned)iter_sb->exchanger->nmsg[1];
  if (msg_ofs < nmsg)
    return msg_ofs >= (unsigned)iter_sb->exchanger->nmsg[0];
  else
    return -1;
}

static MPI_Aint
xt_exchanger_simple_base_iter_get_displacement(Xt_exchanger_iter iter)
{
  Xt_exchanger_simple_base_iter iter_sb = (Xt_exchanger_simple_base_iter)iter;
  const struct Xt_msg_param *msgs = iter_sb->exchanger->msgs;
  unsigned msg_ofs = iter_sb->msg_ofs,
    nmsg = (unsigned)iter_sb->exchanger->nmsg[0]
    + (unsigned)iter_sb->exchanger->nmsg[1];
  return msg_ofs < nmsg ? msgs[msg_ofs].displacement : -1;
}

static void
xt_exchanger_simple_base_iter_delete(Xt_exchanger_iter iter)
{
  free(iter);
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
