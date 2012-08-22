#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include "mpi.h"

#include "xt_idxlist.h"
#include "xt_idxvec.h"
#include "xt_xmap.h"
#include "xt_mpi.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"

static MPI_Comm     xmap_all2all_get_communicator(Xt_xmap xmap);
static int          xmap_all2all_get_num_destinations(Xt_xmap xmap);
static int          xmap_all2all_get_num_sources(Xt_xmap xmap);
static void         xmap_all2all_get_destination_ranks(Xt_xmap xmap, int * ranks);
static void         xmap_all2all_get_source_ranks(Xt_xmap xmap, int * ranks);
static Xt_xmap_iter xmap_all2all_get_destination_iterator(Xt_xmap xmap);
static Xt_xmap_iter xmap_all2all_get_source_iterator(Xt_xmap xmap);
static void         xmap_all2all_delete(Xt_xmap xmap);
static int          xmap_iterator_all2all_next(Xt_xmap_iter iter);
static int          xmap_all2all_iterator_get_rank(Xt_xmap_iter iter);
static Xt_idxlist   xmap_all2all_iterator_get_xlist(Xt_xmap_iter iter);
static void         xmap_all2all_iterator_delete(Xt_xmap_iter iter);
static int          xmap_all2all_get_max_src_pos(Xt_xmap xmap);
static int          xmap_all2all_get_max_dst_pos(Xt_xmap xmap);


static const struct Xt_xmap_iterator_vtable xmap_iterator_all2all_vtable = {
        .next      = xmap_iterator_all2all_next,
        .get_rank  = xmap_all2all_iterator_get_rank,
        .get_xlist = xmap_all2all_iterator_get_xlist,
        .delete    = xmap_all2all_iterator_delete};


struct Xt_xmap_iterator_all2all {

  const struct Xt_xmap_iterator_vtable * vtable;

  struct exchange_data * msg;
  int msgs_left;
};

static const struct Xt_xmap_vtable xmap_all2all_vtable = {
        .get_communicator         = xmap_all2all_get_communicator,
        .get_num_destinations     = xmap_all2all_get_num_destinations,
        .get_num_sources          = xmap_all2all_get_num_sources,
        .get_destination_ranks    = xmap_all2all_get_destination_ranks,
        .get_source_ranks         = xmap_all2all_get_source_ranks,
        .get_destination_iterator = xmap_all2all_get_destination_iterator,
        .get_source_iterator      = xmap_all2all_get_source_iterator,
        .delete                   = xmap_all2all_delete,
        .get_max_src_pos          = xmap_all2all_get_max_src_pos,
        .get_max_dst_pos          = xmap_all2all_get_max_dst_pos};

struct exchange_data {
  // list of relative positions in memory to send or receive
  Xt_idxlist xlist;
  int rank;
};

struct Xt_xmap_all2all {

  const struct Xt_xmap_vtable * vtable;

  struct exchange_data *dst_msg, *src_msg;
  int ndst, nsrc;

  // we need the max position in order to enable quick range-checks
  // for xmap-users like redist
  int max_src_pos; // max possible pos over all src xlists (always >= 0)
  int max_dst_pos; // same for dst

  MPI_Comm comm;
};

static MPI_Comm xmap_all2all_get_communicator(Xt_xmap xmap) {

  struct Xt_xmap_all2all * xmap_all2all;

  xmap_all2all = (struct Xt_xmap_all2all*)xmap;

  return xmap_all2all->comm;
}

static int xmap_all2all_get_num_destinations(Xt_xmap xmap) {

  struct Xt_xmap_all2all * xmap_all2all;

  xmap_all2all = (struct Xt_xmap_all2all*)xmap;

  return xmap_all2all->ndst;
}

static int xmap_all2all_get_num_sources(Xt_xmap xmap) {

  struct Xt_xmap_all2all * xmap_all2all;

  xmap_all2all = (struct Xt_xmap_all2all*)xmap;

  return xmap_all2all->nsrc;
}

static void xmap_all2all_get_destination_ranks(Xt_xmap xmap, int * ranks) {

  struct Xt_xmap_all2all * xmap_all2all;

  xmap_all2all = (struct Xt_xmap_all2all*)xmap;

  int i;

  for (i = 0; i < xmap_all2all->ndst; ++i)
    ranks[i] = xmap_all2all->dst_msg[i].rank;
}

static void xmap_all2all_get_source_ranks(Xt_xmap xmap, int * ranks) {

  struct Xt_xmap_all2all * xmap_all2all;

  xmap_all2all = (struct Xt_xmap_all2all*)xmap;

  int i;

  for (i = 0; i < xmap_all2all->nsrc; ++i)
    ranks[i] = xmap_all2all->src_msg[i].rank;
}

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

  xt_idxlist *src_idxlists_all, *dst_idxlists_all;

  unsigned i;
  int comm_size;

  mpi_err_handler(MPI_Comm_size(comm, &comm_size), comm);

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

  mpi_err_handler(MPI_Allgather(&send_buffer_header, 1, MPI_INT,
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
  mpi_err_handler(MPI_Allgatherv(send_buffer, send_buffer_size, MPI_PACKED,
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

static void generate_src_xlist(struct Xt_xmap_all2all * xmap,
                               Xt_idxlist * intersections,
                               int * ranks,
                               Xt_count num_intersections,
                               Xt_idxlist mypart_idxlist) {

  Xt_count intersection_size[num_intersections];
  int max_intersection_size = 0;
  const int single_match_only = 0;

  xmap->nsrc = num_intersections;
  xmap->src_msg = xmalloc(num_intersections * sizeof(*(xmap->src_msg)));

  // find max size of intersections:
  for (int i = 0; i < xmap->nsrc; ++i) {
    intersection_size[i]  =  xt_idxlist_get_num_indices(intersections[i]);
    if (intersection_size[i] > max_intersection_size)
      max_intersection_size = intersection_size[i];
  }

  Xt_idx *intersection_idxvec = xmalloc(max_intersection_size * sizeof(Xt_idx));
  int *intersection_pos = xmalloc(max_intersection_size * sizeof(int));

  int retval;
  for (int i = 0; i < xmap->nsrc; ++i) {
    xt_idxlist_get_indices(intersections[i], intersection_idxvec);
    retval = xt_idxlist_get_position_of_indices(mypart_idxlist, intersection_idxvec,
                                                intersection_size[i], intersection_pos,
                                                single_match_only);
    assert(retval != 1);

    xmap->src_msg[i].xlist = xt_idxvec_new(intersection_pos, intersection_size[i]);
    xmap->src_msg[i].rank = ranks[i];
  }

  free(intersection_pos);
  free(intersection_idxvec);
}

static void generate_dst_xlist(struct Xt_xmap_all2all * xmap,
                               Xt_idxlist * intersections,
                               int * ranks,
                               Xt_count num_intersections,
                               Xt_idxlist mypart_idxlist) {

  Xt_count intersection_size[num_intersections];
  int max_intersection_size = 0;
  const int single_match_only = 1;

  xmap->ndst = num_intersections;
  xmap->dst_msg = xmalloc(num_intersections * sizeof(*(xmap->dst_msg)));

  // find max size of intersections:
  for (int i = 0; i < xmap->ndst; ++i) {
    intersection_size[i]  =  xt_idxlist_get_num_indices(intersections[i]);
    if (intersection_size[i] > max_intersection_size)
      max_intersection_size = intersection_size[i];
  }

  Xt_idx *intersection_idxvec = xmalloc(max_intersection_size * sizeof(Xt_idx));
  int *intersection_pos = xmalloc(max_intersection_size * sizeof(int));

  int retval;
  for (int i = 0; i < xmap->ndst; ++i) {
    xt_idxlist_get_indices(intersections[i], intersection_idxvec);
    retval = xt_idxlist_get_position_of_indices(mypart_idxlist, intersection_idxvec,
                                                intersection_size[i], intersection_pos,
                                                single_match_only);
    assert(retval != 1);

    xmap->dst_msg[i].xlist = xt_idxvec_new(intersection_pos, intersection_size[i]);
    xmap->dst_msg[i].rank = ranks[i];
  }

  free(intersection_pos);
  free(intersection_idxvec);
}

static void generate_xlist(struct Xt_xmap_all2all * xmap,
                           Xt_idxlist * src_intersections,
                           int * src_ranks,
                           Xt_count num_src_intersections,
                           Xt_idxlist * dst_intersections,
                           int * dst_ranks,
                           Xt_count num_dst_intersections,
                           Xt_idxlist src_idxlist_local,
                           Xt_idxlist dst_idxlist_local) {

  generate_src_xlist(xmap, src_intersections, src_ranks,
                     num_src_intersections, src_idxlist_local);
  generate_dst_xlist(xmap, dst_intersections, dst_ranks,
                     num_dst_intersections, dst_idxlist_local);
}

Xt_xmap xt_xmap_all2all_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist, MPI_Comm comm) {

  struct Xt_xmap_all2all * xmap;

  xmap = xmalloc(1 * sizeof(*xmap));

  xmap->vtable = &xmap_all2all_vtable;

  mpi_err_handler(MPI_Comm_dup(comm, &(xmap->comm)), comm);

  Xt_idxlist * src_intersections = NULL, * dst_intersections = NULL;
  int * src_ranks = NULL, * dst_ranks = NULL;
  Xt_count num_src_intersections = 0, num_dst_intersections = 0;

  // exchange index lists between all processes in comm
  exchange_idxlists(&src_intersections, &src_ranks, &num_src_intersections,
                    &dst_intersections, &dst_ranks, &num_dst_intersections,
                    src_idxlist, dst_idxlist, xmap->comm);

  // generate exchange lists
  generate_xlist(xmap, src_intersections, src_ranks, num_src_intersections,
                 dst_intersections, dst_ranks, num_dst_intersections, src_idxlist,
                 dst_idxlist);

  for (Xt_count i = 0; i < num_src_intersections; ++i)
    xt_idxlist_delete(src_intersections[i]);
  for (Xt_count i = 0; i < num_dst_intersections; ++i)
    xt_idxlist_delete(dst_intersections[i]);
  free(src_intersections), free(dst_intersections);
  free(src_ranks), free(dst_ranks);

  // we could also calculate the (more precise) max pos using only xmap data
  // but using this simple estimate we are still okay for usage checks
  xmap->max_src_pos = xt_idxlist_get_num_indices(src_idxlist);
  xmap->max_dst_pos = xt_idxlist_get_num_indices(dst_idxlist);

  return (Xt_xmap)xmap;
}

static int xmap_all2all_get_max_src_pos(Xt_xmap xmap) {
  struct Xt_xmap_all2all *xmap_all2all = (struct Xt_xmap_all2all *)xmap;
  return xmap_all2all->max_src_pos;
}

static int xmap_all2all_get_max_dst_pos(Xt_xmap xmap) {
  struct Xt_xmap_all2all *xmap_all2all = (struct Xt_xmap_all2all *)xmap;
  return xmap_all2all->max_dst_pos;
}


static void xmap_all2all_delete(Xt_xmap xmap) {

  struct Xt_xmap_all2all * xmap_all2all;

  xmap_all2all = (struct Xt_xmap_all2all*)xmap;

  unsigned i;

  for (i = 0; i < xmap_all2all->ndst; ++i)
    xt_idxlist_delete(xmap_all2all->dst_msg[i].xlist);

  for (i = 0; i < xmap_all2all->nsrc; ++i)
    xt_idxlist_delete(xmap_all2all->src_msg[i].xlist);

  free(xmap_all2all->dst_msg);
  free(xmap_all2all->src_msg);

  mpi_err_handler(MPI_Comm_free(&(xmap_all2all->comm)), MPI_COMM_WORLD);
  free(xmap_all2all);
}

static Xt_xmap_iter xmap_all2all_get_destination_iterator(Xt_xmap xmap) {

  struct Xt_xmap_all2all * xmap_all2all;

  xmap_all2all = (struct Xt_xmap_all2all*)xmap;

  if (xmap_all2all->ndst == 0)
    return NULL;

  struct Xt_xmap_iterator_all2all * iter;

  iter = xmalloc(1 * sizeof(*iter));

  iter->vtable = &xmap_iterator_all2all_vtable;
  iter->msg = xmap_all2all->dst_msg;
  iter->msgs_left = xmap_all2all->ndst - 1;

  return (Xt_xmap_iter)iter;
}

static Xt_xmap_iter xmap_all2all_get_source_iterator(Xt_xmap xmap) {

  struct Xt_xmap_all2all * xmap_all2all;

  xmap_all2all = (struct Xt_xmap_all2all*)xmap;

  if (xmap_all2all->nsrc == 0)
    return NULL;

  struct Xt_xmap_iterator_all2all * iter;

  iter = xmalloc(1 * sizeof(*iter));

  iter->vtable = &xmap_iterator_all2all_vtable;
  iter->msg = xmap_all2all->src_msg;
  iter->msgs_left = xmap_all2all->nsrc - 1;

  return (Xt_xmap_iter)iter;
}

static int xmap_iterator_all2all_next(Xt_xmap_iter iter) {

  struct Xt_xmap_iterator_all2all * iter_all2all;

  iter_all2all = (struct Xt_xmap_iterator_all2all *)iter;

  if (iter_all2all == NULL || iter_all2all->msgs_left == 0)
    return 0;

  iter_all2all->msg++;
  iter_all2all->msgs_left--;

  return 1;
}

static int xmap_all2all_iterator_get_rank(Xt_xmap_iter iter) {

  struct Xt_xmap_iterator_all2all * iter_all2all;

  iter_all2all = (struct Xt_xmap_iterator_all2all *)iter;

  assert(iter_all2all != NULL);

  return iter_all2all->msg->rank;
}

static Xt_idxlist xmap_all2all_iterator_get_xlist(Xt_xmap_iter iter) {

  struct Xt_xmap_iterator_all2all * iter_all2all;

  iter_all2all = (struct Xt_xmap_iterator_all2all *)iter;

  assert(iter_all2all != NULL);

  return  xt_idxlist_copy(iter_all2all->msg->xlist);
}

static void xmap_all2all_iterator_delete(Xt_xmap_iter iter) {

  struct Xt_xmap_iterator_all2all * iter_all2all;

  iter_all2all = (struct Xt_xmap_iterator_all2all *)iter;

  free(iter_all2all);
}
