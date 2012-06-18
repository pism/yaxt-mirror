#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "mpi.h"

#include "xt_idxlist.h"
#include "xt_idxvec.h"
#include "xt_xmap.h"
#include "xt_mpi.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"

struct exchange_data {
  // list of global indices to send or receive from partner rank
  Xt_idxlist intersect;
  // list of relative positions in memory to send or receive
  Xt_idxlist xlist;
  int rank;
};

struct xt_xmap {
  struct exchange_data *dst_msg, *src_msg;
  int ndst, nsrc;
  MPI_Comm comm;
};

struct Xt_xmap_iterator {

  struct exchange_data * msg;
  int msgs_left;
};

MPI_Comm xt_xmap_get_communicator(Xt_xmap xmap) {

  return xmap->comm;
}

int xt_xmap_get_num_destinations(Xt_xmap xmap) {

  return xmap->ndst;
}

int xt_xmap_get_num_sources(Xt_xmap xmap) {

  return xmap->nsrc;
}

void xt_xmap_get_destination_ranks(Xt_xmap xmap, int * ranks) {

  int i;

  for (i = 0; i < xmap->ndst; ++i)
    ranks[i] = xmap->dst_msg[i].rank;
}

void xt_xmap_get_source_ranks(Xt_xmap xmap, int * ranks) {

  int i;

  for (i = 0; i < xmap->nsrc; ++i)
    ranks[i] = xmap->src_msg[i].rank;
}

static void exchange_idxlists(xt_xmap xmap,
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


static void generate_src_xlist(xt_xmap xmap,
                               Xt_idxlist mypart_idxlist) {

  int intersection_size[xmap->nsrc];
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


static void generate_dst_xlist(xt_xmap xmap,
                               Xt_idxlist mypart_idxlist) {

  int intersection_size[xmap->ndst];
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

static void generate_xlist(xt_xmap xmap,
                           Xt_idxlist src_idxlist_local,
                           Xt_idxlist dst_idxlist_local) {

  generate_src_xlist(xmap, src_idxlist_local);
  generate_dst_xlist(xmap, dst_idxlist_local);
}

struct xt_xmap*
xt_xmap_new(xt_idxlist src_idxlist, xt_idxlist dst_idxlist, MPI_Comm comm) {

  struct xt_xmap * xmap;

  xmap = xmalloc(1 * sizeof(*xmap));

  mpi_err_handler(MPI_Comm_dup(comm, &(xmap->comm)), comm);

  // exchange index lists between all processes in comm
  exchange_idxlists(xmap, src_idxlist, dst_idxlist);

  // generate exchange lists
  generate_xlist(xmap, src_idxlist, dst_idxlist);

  return xmap;
}

void xt_xmap_delete(xt_xmap xmap) {

  unsigned i;

  for (i = 0; i < xmap->ndst; ++i) {
    xt_idxlist_delete(xmap->dst_msg[i].intersect);
    xt_idxlist_delete(xmap->dst_msg[i].xlist);
  }

  for (i = 0; i < xmap->nsrc; ++i) {
    xt_idxlist_delete(xmap->src_msg[i].intersect);
    xt_idxlist_delete(xmap->src_msg[i].xlist);
  }

  free(xmap->dst_msg);
  free(xmap->src_msg);

  mpi_err_handler(MPI_Comm_free(&(xmap->comm)), MPI_COMM_WORLD);
  free(xmap);
}

static Xt_idxlist get_xlist(struct exchange_data * msgs, int num_msgs, int rank) {

  int i;

  for (i = 0; i < num_msgs; ++i)
    if (rank == msgs[i].rank)
      return xt_idxlist_copy(msgs[i].xlist);

  return NULL;
}

Xt_idxlist xt_xmap_get_dst_xlist(Xt_xmap xmap, int rank) {

  return get_xlist(xmap->dst_msg, xmap->ndst, rank);
}

Xt_idxlist xt_xmap_get_src_xlist(Xt_xmap xmap, int rank) {

  return get_xlist(xmap->src_msg, xmap->nsrc, rank);
}

Xt_xmap_iter xt_xmap_get_destination_iterator(Xt_xmap xmap) {

  if (xmap->ndst == 0)
    return NULL;

  Xt_xmap_iter iter;

  iter = xmalloc(1 * sizeof(*iter));

  iter->msg = xmap->dst_msg;
  iter->msgs_left = xmap->ndst - 1;

  return iter;
}

Xt_xmap_iter xt_xmap_get_source_iterator(Xt_xmap xmap) {

  if (xmap->nsrc == 0)
    return NULL;

  Xt_xmap_iter iter;

  iter = xmalloc(1 * sizeof(*iter));

  iter->msg = xmap->src_msg;
  iter->msgs_left = xmap->nsrc - 1;

  return iter;
}

int xt_xmap_iterator_next(Xt_xmap_iter iter) {

  if (iter == NULL || iter->msgs_left == 0)
    return 0;

  iter->msg++;
  iter->msgs_left--;

  return 1;
}

int xt_xmap_iterator_get_rank(Xt_xmap_iter iter) {

  assert(iter != NULL);

  return iter->msg->rank;
}

Xt_idxlist xt_xmap_iterator_get_xlist(Xt_xmap_iter iter) {

  assert(iter != NULL);

  return  xt_idxlist_copy(iter->msg->xlist);
}

void xt_xmap_iterator_delete(Xt_xmap_iter iter) {

  free(iter);
}
