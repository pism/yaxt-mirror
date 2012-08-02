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
        .delete                   = xmap_all2all_delete};

struct exchange_data {
  // list of global indices to send or receive from partner rank
  Xt_idxlist intersect;
  // list of relative positions in memory to send or receive
  Xt_idxlist xlist;
  int rank;
};

struct Xt_xmap_all2all {

  const struct Xt_xmap_vtable * vtable;

  struct exchange_data *dst_msg, *src_msg;
  int ndst, nsrc;
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

static void exchange_idxlists(struct Xt_xmap_all2all * xmap,
                              Xt_idxlist src_idxlist_local,
                              Xt_idxlist dst_idxlist_local)
{

  xt_idxlist *src_idxlists_all, *dst_idxlists_all;

  unsigned i;
  int comm_size;

  mpi_err_handler(MPI_Comm_size(xmap->comm, &comm_size), xmap->comm);

  // compute size of local index lists
  long long src_pack_size = xt_idxlist_get_pack_size(src_idxlist_local, xmap->comm);
  if (src_pack_size<0) die("src index list is too large");
  long long dst_pack_size = xt_idxlist_get_pack_size(dst_idxlist_local, xmap->comm);
  if (dst_pack_size<0) die("dst index list is too large");

  if (src_pack_size + dst_pack_size >= INT_MAX) die("local src+dst index lists are too large");

  int send_buffer_size;
  send_buffer_size  = src_pack_size + dst_pack_size;

  void * send_buffer;
  send_buffer = xmalloc(send_buffer_size);

  // pack local index lists
  int position = 0;
  xt_idxlist_pack(src_idxlist_local, send_buffer, send_buffer_size, &position,
                  xmap->comm);
  xt_idxlist_pack(dst_idxlist_local, send_buffer, send_buffer_size, &position,
                  xmap->comm);

  // exchange buffer sizes
  int send_buffer_header = send_buffer_size;
  int * recv_buffer_header;

  recv_buffer_header = xmalloc(comm_size * sizeof(*recv_buffer_header));

  mpi_err_handler(MPI_Allgather(&send_buffer_header, 1, MPI_INT,
                                recv_buffer_header,
                                1, MPI_INT, xmap->comm), xmap->comm);

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
                                 MPI_PACKED,
                                 xmap->comm), xmap->comm);

  src_idxlists_all = xmalloc(comm_size * sizeof(*src_idxlists_all));

  dst_idxlists_all = xmalloc(comm_size * sizeof(*dst_idxlists_all));

  // unpack buffers
  for (i = 0, position = 0; i < comm_size; ++i) {
    src_idxlists_all[i] = xt_idxlist_unpack(recv_buffer, recv_buffer_size,
                                            &position, xmap->comm);
    dst_idxlists_all[i] = xt_idxlist_unpack(recv_buffer, recv_buffer_size,
                                            &position, xmap->comm);
  }

  // allocate memory for intersections
  xmap->dst_msg = xmalloc(comm_size * sizeof (*xmap->dst_msg));
  xmap->ndst = 0;
  xmap->src_msg = xmalloc(comm_size * sizeof (*xmap->src_msg));
  xmap->nsrc = 0;
  if (!xmap->dst_msg || !xmap->src_msg)
    abort();

  // compute intersections
  for (i = 0; i < comm_size; ++i)
  {
    Xt_idxlist intersect
      = xt_idxlist_get_intersection(src_idxlists_all[i], dst_idxlist_local);
    if (xt_idxlist_get_num_indices(intersect) > 0)
    {
      xmap->dst_msg[xmap->ndst].intersect = intersect;
      xmap->dst_msg[xmap->ndst].rank = i;
      xmap->ndst++;
    }
    else
      xt_idxlist_delete(intersect);
    intersect = xt_idxlist_get_intersection(src_idxlist_local,
                                            dst_idxlists_all[i]);
    if (xt_idxlist_get_num_indices(intersect) > 0)
    {
      xmap->src_msg[xmap->nsrc].rank = i;
      xmap->src_msg[xmap->nsrc].intersect = intersect;
      xmap->nsrc++;
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
                               Xt_idxlist mypart_idxlist) {

  Xt_count intersection_size[xmap->nsrc];
  int max_intersection_size = 0;
  const int single_match_only = 0;

  // find max size of intersections:
  for (int i = 0; i < xmap->nsrc; ++i) {
    intersection_size[i]  =  xt_idxlist_get_num_indices(xmap->src_msg[i].intersect);
    if (intersection_size[i] > max_intersection_size)  max_intersection_size = intersection_size[i];
  }

  Xt_idx *intersection_idxvec = xmalloc(max_intersection_size * sizeof(Xt_idx));
  int *intersection_pos = xmalloc(max_intersection_size * sizeof(int));

  int retval;
  for (int i = 0; i < xmap->nsrc; ++i) {
    xt_idxlist_get_indices(xmap->src_msg[i].intersect, intersection_idxvec);
    retval = xt_idxlist_get_position_of_indices(mypart_idxlist, intersection_idxvec,
                                                intersection_size[i], intersection_pos, single_match_only);
    assert(retval != 1);

    xmap->src_msg[i].xlist = xt_idxvec_new(intersection_pos, intersection_size[i]);

  }

  free(intersection_pos);
  free(intersection_idxvec);

}

static void generate_dst_xlist(struct Xt_xmap_all2all * xmap,
                               Xt_idxlist mypart_idxlist) {

  Xt_count intersection_size[xmap->ndst];
  int max_intersection_size = 0;
  const int single_match_only = 1;


  // find max size of intersections:
  for (int i = 0; i < xmap->ndst; ++i) {
    intersection_size[i]  =  xt_idxlist_get_num_indices(xmap->dst_msg[i].intersect);
    if (intersection_size[i] > max_intersection_size)  max_intersection_size = intersection_size[i];
  }

  Xt_idx *intersection_idxvec = xmalloc(max_intersection_size * sizeof(Xt_idx));
  int *intersection_pos = xmalloc(max_intersection_size * sizeof(int));

  int retval;
  for (int i = 0; i < xmap->ndst; ++i) {
    xt_idxlist_get_indices(xmap->dst_msg[i].intersect, intersection_idxvec);
    retval = xt_idxlist_get_position_of_indices(mypart_idxlist, intersection_idxvec,
                                                intersection_size[i], intersection_pos, single_match_only);
    assert(retval != 1);

    xmap->dst_msg[i].xlist = xt_idxvec_new(intersection_pos, intersection_size[i]);

  }

  free(intersection_pos);
  free(intersection_idxvec);

}

static void generate_xlist(struct Xt_xmap_all2all * xmap,
                           Xt_idxlist src_idxlist_local,
                           Xt_idxlist dst_idxlist_local) {

  generate_src_xlist(xmap, src_idxlist_local);
  generate_dst_xlist(xmap, dst_idxlist_local);
}

Xt_xmap xt_xmap_all2all_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist, MPI_Comm comm) {

  struct Xt_xmap_all2all * xmap;

  xmap = xmalloc(1 * sizeof(*xmap));

  xmap->vtable = &xmap_all2all_vtable;

  mpi_err_handler(MPI_Comm_dup(comm, &(xmap->comm)), comm);

  // exchange index lists between all processes in comm
  exchange_idxlists(xmap, src_idxlist, dst_idxlist);

  // generate exchange lists
  generate_xlist(xmap, src_idxlist, dst_idxlist);

  return (Xt_xmap)xmap;
}

static void xmap_all2all_delete(Xt_xmap xmap) {

  struct Xt_xmap_all2all * xmap_all2all;

  xmap_all2all = (struct Xt_xmap_all2all*)xmap;

  unsigned i;

  for (i = 0; i < xmap_all2all->ndst; ++i) {
    xt_idxlist_delete(xmap_all2all->dst_msg[i].intersect);
    xt_idxlist_delete(xmap_all2all->dst_msg[i].xlist);
  }

  for (i = 0; i < xmap_all2all->nsrc; ++i) {
    xt_idxlist_delete(xmap_all2all->src_msg[i].intersect);
    xt_idxlist_delete(xmap_all2all->src_msg[i].xlist);
  }

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
