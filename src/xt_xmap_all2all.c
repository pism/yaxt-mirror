/**
 * @file xt_xmap_all2all.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
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
 * URL: https://redmine.dkrz.de/doc/yaxt/html/index.html
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <mpi.h>

#include "xt/xt_idxlist.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_xmap.h"
#include "xt/xt_mpi.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_xmap_intersection.h"

static void exchange_idxlists(Xt_idxlist ** src_intersections,
                              int ** src_ranks,
                              Xt_count * num_src_intersections,
                              Xt_idxlist ** dst_intersections,
                              int ** dst_ranks,
                              Xt_count * num_dst_intersections,
                              Xt_idxlist src_idxlist_local,
                              Xt_idxlist dst_idxlist_local,
                              MPI_Comm comm)
{

  /*
    Note: The meaning of source (src) and destination (dst) points can already be understood by
    looking at the serial case, where it is just a transformation of sequences of integers
    (called indices). The starting state (source sequence) is transformed into an end state
    (dst sequence). The transformation does not have to be bijective. For each position dpos of
    the dst sequence we have at least one position spos in the src sequence with:
    dst(dpos) = src(spos)
   */

  Xt_idxlist *src_idxlists_all, *dst_idxlists_all;

  unsigned i;
  int comm_size;

  xt_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // compute size of local index lists
  long long src_pack_size = xt_idxlist_get_pack_size(src_idxlist_local, comm);
  if (src_pack_size<0) die("src index list is too large");
  long long dst_pack_size = xt_idxlist_get_pack_size(dst_idxlist_local, comm);
  if (dst_pack_size<0) die("dst index list is too large");

  if (src_pack_size + dst_pack_size >= INT_MAX) die("local src+dst index lists are too large");

  int send_buffer_size;
  send_buffer_size  = src_pack_size + dst_pack_size;

  void * send_buffer;
  send_buffer = xmalloc(send_buffer_size);

  // pack local index lists
  int position = 0;
  xt_idxlist_pack(src_idxlist_local, send_buffer, send_buffer_size, &position, comm);
  xt_idxlist_pack(dst_idxlist_local, send_buffer, send_buffer_size, &position, comm);

  // exchange buffer sizes
  int send_buffer_header = send_buffer_size;
  int * recv_buffer_header;

  recv_buffer_header = xmalloc(comm_size * sizeof(*recv_buffer_header));

  xt_mpi_call(MPI_Allgather(&send_buffer_header, 1, MPI_INT,
                                recv_buffer_header, 1, MPI_INT, comm), comm);

  void * recv_buffer;
  int * displ;

  displ = xmalloc(comm_size * sizeof(*displ));

  displ[0] = 0;
  int recv_buffer_size = recv_buffer_header[0];
  for (i = 1; i < comm_size; ++i) {

     displ[i] = displ[i-1] + recv_buffer_header[i-1];
     recv_buffer_size += recv_buffer_header[i];
     if (recv_buffer_size<0) die("accumulated buffer sizes too big");
  }

  recv_buffer = xmalloc(recv_buffer_size);

  // exchange buffers
  xt_mpi_call(MPI_Allgatherv(send_buffer, send_buffer_size, MPI_PACKED,
                                 recv_buffer, recv_buffer_header, displ,
                                 MPI_PACKED, comm), comm);

  src_idxlists_all = xmalloc(comm_size * sizeof(*src_idxlists_all));

  dst_idxlists_all = xmalloc(comm_size * sizeof(*dst_idxlists_all));

  // unpack buffers
  for (i = 0, position = 0; i < comm_size; ++i) {
    src_idxlists_all[i] = xt_idxlist_unpack(recv_buffer, recv_buffer_size,
                                            &position, comm);
    dst_idxlists_all[i] = xt_idxlist_unpack(recv_buffer, recv_buffer_size,
                                            &position, comm);
  }

  // allocate memory for intersections
  (*dst_intersections) = xmalloc(comm_size * sizeof (**dst_intersections));
  (*dst_ranks) = xmalloc(comm_size * sizeof(**dst_ranks));
  (*num_dst_intersections) = 0;
  (*src_intersections) = xmalloc(comm_size * sizeof (**src_intersections));
  (*src_ranks) = xmalloc(comm_size * sizeof(**src_ranks));
  (*num_src_intersections) = 0;

  // compute intersections
  for (i = 0; i < comm_size; ++i) {

    Xt_idxlist intersect;

    intersect = xt_idxlist_get_intersection(src_idxlists_all[i],
                                            dst_idxlist_local);
    if (xt_idxlist_get_num_indices(intersect) > 0) {

      (*dst_intersections)[(*num_dst_intersections)] = intersect;
      (*dst_ranks)[(*num_dst_intersections)] = i;
      (*num_dst_intersections)++;
    }
    else
      xt_idxlist_delete(intersect);

    intersect = xt_idxlist_get_intersection(src_idxlist_local,
                                            dst_idxlists_all[i]);
    if (xt_idxlist_get_num_indices(intersect) > 0) {

      (*src_intersections)[(*num_src_intersections)] = intersect;
      (*src_ranks)[(*num_src_intersections)] = i;
      (*num_src_intersections)++;
    }
    else
      xt_idxlist_delete(intersect);
  }

  for (i = 0; i < comm_size; ++i) {

    xt_idxlist_delete(src_idxlists_all[i]);
    xt_idxlist_delete(dst_idxlists_all[i]);
  }

  free(src_idxlists_all);
  free(dst_idxlists_all);

  // clean up
  free(send_buffer);
  free(recv_buffer_header);
  free(recv_buffer);
  free(displ);
}

Xt_xmap xt_xmap_all2all_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist, MPI_Comm comm) {

  MPI_Comm curr_comm;

  xt_mpi_call(MPI_Comm_dup(comm, &curr_comm), comm);

  Xt_idxlist * src_intersections = NULL, * dst_intersections = NULL;
  int * src_ranks = NULL, * dst_ranks = NULL;
  Xt_count num_src_intersections = 0, num_dst_intersections = 0;

  // exchange index lists between all processes in comm
  exchange_idxlists(&src_intersections, &src_ranks, &num_src_intersections,
                    &dst_intersections, &dst_ranks, &num_dst_intersections,
                    src_idxlist, dst_idxlist, curr_comm);

  MPI_Comm_free(&curr_comm);

  Xt_xmap xmap;

  xmap = xt_xmap_intersection_new(src_intersections, src_ranks, num_src_intersections,
                                  dst_intersections, dst_ranks, num_dst_intersections,
                                  src_idxlist, dst_idxlist, comm);

  for (Xt_count i = 0; i < num_src_intersections; ++i)
    if (src_intersections[i] != NULL)
      xt_idxlist_delete(src_intersections[i]);
  for (Xt_count i = 0; i < num_dst_intersections; ++i)
    if (dst_intersections[i] != NULL)
      xt_idxlist_delete(dst_intersections[i]);
  free(src_intersections), free(dst_intersections);
  free(src_ranks), free(dst_ranks);

  return xmap;
}
