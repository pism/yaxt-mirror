/**
 * @file xt_dmap_offset.c
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

#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "core/core.h"
#include "xt/xt_core.h"
#include "xt/xt_dmap_offset.h"
#include "core/ppm_xfuncs.h"
#include "xt_dmap_internal.h"
#include "xt/xt_xmap.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt/quicksort.h"
#include "xt/xt_redist.h"
#include "xt/xt_redist_single_array_base.h"

static void dmap_offset_repeat(Xt_dmap dmap, int num_repetitions,
                               const int src_displacements[num_repetitions],
                               const int dst_displacements[num_repetitions]);
static void dmap_offset_reorder(
  Xt_dmap dmap, enum xt_reorder_type type);
static Xt_redist dmap_offset_generate_redist(
  Xt_dmap dmap, MPI_Datatype base_datatype);
static void dmap_offset_delete(Xt_dmap dmap);

static const struct xt_dmap_vtable dmap_offset_vtable = {
  .repeat = dmap_offset_repeat,
  .reorder = dmap_offset_reorder,
  .generate_redist = dmap_offset_generate_redist,
  .delete = dmap_offset_delete,
};

struct Xt_dmap_msg {
  int rank;
  int * offsets;
  int count;
};

struct dmap_offset {
  const struct xt_dmap_vtable * vtable;

  int nsend, nrecv;
  MPI_Comm comm;
  struct Xt_dmap_msg msgs[];
};

static struct dmap_offset * dmap_offset_alloc(int nmsg) {

  struct dmap_offset * dmap;
  size_t header_size = sizeof(*dmap),
         body_size = (size_t)nmsg * sizeof (dmap->msgs[0]);
  dmap = xmalloc(header_size + body_size);
  dmap->vtable = &dmap_offset_vtable;
  return dmap;
}

static void generate_msg_infos(
  Xt_xmap_iter iter, const int *offsets, struct Xt_dmap_msg * msgs) {

  do {

    // get the points that need to be transfered
    const int *curr_transfer_pos = xt_xmap_iterator_get_transfer_pos(iter);
    int count = xt_xmap_iterator_get_num_transfer_pos(iter);
    int * curr_offsets = xmalloc(count * sizeof(*curr_offsets));

    // get the actual offsets of the points that are to be sent
    for (int i = 0; i < count; ++i)
      curr_offsets[i] = offsets[curr_transfer_pos[i]];

    msgs->offsets = curr_offsets;
    msgs->count = count;
    msgs->rank = xt_xmap_iterator_get_rank(iter);

    msgs++;

  // while there are messages to be processed
  } while (xt_xmap_iterator_next(iter));

  xt_xmap_iterator_delete(iter);
}

Xt_dmap xt_dmap_offset_new(
  Xt_xmap xmap, const int *src_offsets, const int *dst_offsets) {

  int nrecv = xt_xmap_get_num_sources(xmap);
  int nsend = xt_xmap_get_num_destinations(xmap);

  struct dmap_offset * dmap_off = dmap_offset_alloc(nsend + nrecv);
  dmap_off->nsend = nsend;
  dmap_off->nrecv = nrecv;
  dmap_off->comm = xt_xmap_get_communicator(xmap);

  struct Xt_dmap_msg * send_msgs = &(dmap_off->msgs[0]);
  struct Xt_dmap_msg * recv_msgs = send_msgs + (size_t)nsend;

  if (nrecv)
    generate_msg_infos(xt_xmap_get_in_iterator(xmap), dst_offsets, recv_msgs);

  if (nsend)
    generate_msg_infos(xt_xmap_get_out_iterator(xmap), src_offsets, send_msgs);

  return (Xt_dmap)dmap_off;
}

static void dmap_msg_offset_repeat(int msg_count, struct Xt_dmap_msg * msgs,
                                   int num_repetitions,
                                   const int displacements[num_repetitions]) {

  // for all messages
  for (int i = 0; i < msg_count; ++i) {

    struct Xt_dmap_msg * curr_msg = msgs + i;
    int count = curr_msg->count;
    int * new_offsets = xmalloc(num_repetitions * count * sizeof(*new_offsets));
    int * old_offsets = curr_msg->offsets;

    // add displacements to original offsets
    for (int j = 0; j < num_repetitions; ++j)
      for (int k = 0; k < count; ++k)
        new_offsets[j * count + k] = old_offsets[k] + displacements[j];

    free((void*)old_offsets);
    curr_msg->offsets = new_offsets;
    curr_msg->count *= num_repetitions;
  }
}

static void dmap_offset_repeat(Xt_dmap dmap, int num_repetitions,
                               const int src_displacements[num_repetitions],
                               const int dst_displacements[num_repetitions]) {

  struct dmap_offset * dmap_off = (struct dmap_offset*)dmap;

  int nsend = dmap_off->nsend;
  int nrecv = dmap_off->nrecv;
  struct Xt_dmap_msg * send_msgs = &(dmap_off->msgs[0]);
  struct Xt_dmap_msg * recv_msgs = send_msgs + (size_t)nsend;

  // apply offsets to send and recv messages
  dmap_msg_offset_repeat(nsend, send_msgs, num_repetitions, src_displacements);
  dmap_msg_offset_repeat(nrecv, recv_msgs, num_repetitions, dst_displacements);
}

static inline void
setup_buffer(int comm_size, int n, const struct Xt_dmap_msg * msgs,
             int ** buffer, int ** counts, int ** displs) {

  int buffer_size = 0;
  *counts = xcalloc((size_t)comm_size, sizeof(*counts));
  *displs = xmalloc((size_t)comm_size * sizeof(*displs));

  for (int i = 0; i < n; ++i) {
    buffer_size += msgs[i].count;
    (*counts)[msgs[i].rank] = msgs[i].count;
  }

  for (int i = 0, accu = 0; i < comm_size; accu += (*counts)[i++])
    (*displs)[i] = accu;

  *buffer = xmalloc((size_t)buffer_size * sizeof(*buffer));
}

static void reorder_offsets(
  int nsend, int nrecv, const struct Xt_dmap_msg * send_msgs,
  const struct Xt_dmap_msg *recv_msgs, MPI_Comm comm) {

  int comm_size;
  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int * send_buffer, * send_counts, * send_displs;
  int * recv_buffer, * recv_counts, * recv_displs;
  setup_buffer(comm_size, nsend, send_msgs,
               &send_buffer, &send_counts, &send_displs);
  setup_buffer(comm_size, nrecv, recv_msgs,
               &recv_buffer, &recv_counts, &recv_displs);

  for (int i = 0; i < nsend; ++i) {
    const struct Xt_dmap_msg * curr_msg = send_msgs + i;
    int count = curr_msg->count;
    if (count > 0) {
      // pack local send offsets into send buffer
      memcpy(send_buffer + send_displs[curr_msg->rank],
             curr_msg->offsets, count * sizeof(*send_buffer));
      // sort local send offsets
      xt_quicksort_int(curr_msg->offsets, count);
    }
  }

  // exchange local send offsets
  xt_mpi_call(MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_INT,
                            recv_buffer, recv_counts, recv_displs, MPI_INT,
                            comm), comm);

  for (int i = 0; i < nrecv; ++i) {
    const struct Xt_dmap_msg * curr_msg = recv_msgs + i;
    int count = curr_msg->count;
    if (count > 0)
      // reorder local receive offsets (curr_msg->offsets) based on remote send
      // offsets (recv_buffer + recv_displs[curr_msg->rank])
      xt_quicksort_index_int(
        recv_buffer + recv_displs[curr_msg->rank], count, curr_msg->offsets, 0);
  }

  free(send_buffer);
  free(send_counts);
  free(send_displs);
  free(recv_buffer);
  free(recv_counts);
  free(recv_displs);
}

static void dmap_offset_reorder(
  Xt_dmap dmap, enum xt_reorder_type type) {

  struct dmap_offset * dmap_off = (struct dmap_offset*)dmap;

  int nsend = dmap_off->nsend;
  int nrecv = dmap_off->nrecv;
  struct Xt_dmap_msg * send_msgs = &(dmap_off->msgs[0]);
  struct Xt_dmap_msg * recv_msgs = send_msgs + (size_t)nsend;
  MPI_Comm comm = dmap_off->comm;

  switch (type) {
    case (XT_REORDER_NONE):
      return;
    case (XT_REORDER_SEND_UP):
      reorder_offsets(nsend, nrecv, send_msgs, recv_msgs, comm);
      break;
    case (XT_REORDER_RECV_UP):
      reorder_offsets(nrecv, nsend, recv_msgs, send_msgs, comm);
      break;
    default:
      Xt_abort(comm, "ERROR(dmap_offset_reorder):invalid reorder type",
               __FILE__, __LINE__);
  };
}

static struct Xt_redist_msg * generate_redist_msgs(
  int count, struct Xt_dmap_msg * dmap_msgs, MPI_Datatype base_datatype,
  MPI_Comm comm) {

  struct Xt_redist_msg * redist_msgs = xmalloc(count * sizeof(*redist_msgs));

  for (int i = 0; i < count; ++i) {

    struct Xt_dmap_msg * curr_dmap_msg = dmap_msgs + i;

    redist_msgs[i].datatype =
      xt_mpi_generate_datatype(
        curr_dmap_msg->offsets, curr_dmap_msg->count, base_datatype, comm);
    redist_msgs[i].rank = curr_dmap_msg->rank;
  }

  return redist_msgs;
}

static Xt_redist dmap_offset_generate_redist(
  Xt_dmap dmap, MPI_Datatype base_datatype) {

  struct dmap_offset * dmap_off = (struct dmap_offset*)dmap;

  int nsend = dmap_off->nsend;
  int nrecv = dmap_off->nrecv;

  struct Xt_redist_msg * redist_msgs =
    generate_redist_msgs(
      nsend + nrecv, dmap_off->msgs, base_datatype, dmap_off->comm);

  Xt_redist redist =
    xt_redist_single_array_base_new(
      nsend, nrecv, redist_msgs, redist_msgs + nsend, dmap_off->comm);

  for (int i = 0; i < nsend + nrecv; ++i) {
    MPI_Datatype dt = redist_msgs[i].datatype;
    if (dt != MPI_DATATYPE_NULL)
      xt_mpi_call(MPI_Type_free(&dt), dmap_off->comm);
  }

  free(redist_msgs);

  return redist;
}

static void dmap_offset_delete(Xt_dmap dmap) {

  struct dmap_offset * dmap_off = (struct dmap_offset*)dmap;
  int nmsg = dmap_off->nsend + dmap_off->nrecv;
  for (int i = 0; i < nmsg; ++i) free(dmap_off->msgs[i].offsets);
  free(dmap_off);
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
