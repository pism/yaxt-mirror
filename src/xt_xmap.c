#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "mpi.h"

#include "xt_idxvec.h"
#include "xt_idxlist.h"
#include "xt_xmap.h"
#include "xt_mpi.h"

struct exchange_data {
  // list of global indices to send or receive from partner rank
  Xt_idxvec intersect;
  // list of relative positions in memory to send or receive
  Xt_idxvec xlist;
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

static void exchange_idxvecs(xt_xmap xmap,
                             Xt_idxvec src_idxvec_local,
                             Xt_idxvec dst_idxvec_local)
{

  xt_idxvec *src_idxvecs_all, *dst_idxvecs_all;

  unsigned i, j;
  int comm_size, comm_rank;

  mpi_err_handler(MPI_Comm_size(xmap->comm, &comm_size), xmap->comm);

  // compute size of local index vectors
  size_t send_buffer_size;
  send_buffer_size  = xt_idxvec_get_pack_size(src_idxvec_local, xmap->comm);
  send_buffer_size += xt_idxvec_get_pack_size(dst_idxvec_local, xmap->comm);

  void * send_buffer;

  send_buffer = malloc(send_buffer_size);
  if (send_buffer == NULL) {
    fputs("malloc failed for send buffer of index vector exchange\n", stderr);
     exit(EXIT_FAILURE);
  }

  // pack local index vectors
  int position = 0;
  xt_idxvec_pack(src_idxvec_local, send_buffer, send_buffer_size, &position,
                 xmap->comm);
  xt_idxvec_pack(dst_idxvec_local, send_buffer, send_buffer_size, &position,
                 xmap->comm);

  // exchange buffer sizes
  int send_buffer_header = send_buffer_size;
  int * recv_buffer_header;

  recv_buffer_header = malloc(comm_size * sizeof(*recv_buffer_header));
  if (recv_buffer_header == NULL) {
    fputs("malloc failed for receive buffer header"
          " of index vector exchange\n", stderr);
    exit(EXIT_FAILURE);
  }

  mpi_err_handler(MPI_Allgather(&send_buffer_header, 1, MPI_INT,
                                recv_buffer_header,
                                1, MPI_INT, xmap->comm), xmap->comm);

  void * recv_buffer;
  int * displ;

  displ = malloc(comm_size * sizeof(*displ));
  if (displ == NULL) {
    fputs("malloc failed for displacment of index vector exchange\n", stderr);
    exit(EXIT_FAILURE);
  }

  displ[0] = 0;
  size_t recv_buffer_size = recv_buffer_header[0];
  for (i = 1; i < comm_size; ++i) {

     displ[i] = displ[i-1] + recv_buffer_header[i-1];
     recv_buffer_size += recv_buffer_header[i];
  }

  recv_buffer = malloc(recv_buffer_size);
  if (recv_buffer == NULL) {
    fputs("malloc failed for receive buffer of index vector exchange\n",
          stderr);
    exit(EXIT_FAILURE);
  }

  // exchange buffers
  mpi_err_handler(MPI_Allgatherv(send_buffer, send_buffer_size, MPI_PACKED,
                                 recv_buffer, recv_buffer_header, displ,
                                 MPI_PACKED,
                                 xmap->comm), xmap->comm);

  src_idxvecs_all = malloc(comm_size * sizeof(*src_idxvecs_all));
  if (src_idxvecs_all == NULL) {
    fputs("malloc failed for source index vectors"
          " all of index vector exchange\n", stderr);
    exit(EXIT_FAILURE);
  }

  dst_idxvecs_all = malloc(comm_size * sizeof(*dst_idxvecs_all));
  if (dst_idxvecs_all == NULL) {
    fputs("malloc failed for destination index vectors"
          " all of index vector exchange\n", stderr);
    exit(EXIT_FAILURE);
  }

  // unpack buffers
  for (i = 0, position = 0; i < comm_size; ++i) {

    src_idxvecs_all[i] = xt_idxvec_unpack(recv_buffer, recv_buffer_size,
                                          &position, xmap->comm);
    dst_idxvecs_all[i] = xt_idxvec_unpack(recv_buffer, recv_buffer_size,
                                          &position, xmap->comm);
  }

  // allocate memory for intersections
  xmap->dst_msg = malloc(comm_size * sizeof (*xmap->dst_msg));
  xmap->ndst = 0;
  xmap->src_msg = malloc(comm_size * sizeof (*xmap->src_msg));
  xmap->nsrc = 0;
  if (!xmap->dst_msg || !xmap->src_msg)
    abort();

  // compute intersections
  for (i = 0; i < comm_size; ++i)
  {
    Xt_idxvec intersect
      = xt_idxvec_get_intersection(src_idxvecs_all[i], dst_idxvec_local);
    if (xt_idxvec_get_num_indices(intersect) > 0)
    {
      xmap->dst_msg[xmap->ndst].intersect = intersect;
      xmap->dst_msg[xmap->ndst].rank = i;
      xmap->ndst++;
    }
    else
      xt_idxvec_delete(intersect);
    intersect = xt_idxvec_get_intersection(src_idxvec_local,
                                           dst_idxvecs_all[i]);
    if (xt_idxvec_get_num_indices(intersect) > 0)
    {
      xmap->src_msg[xmap->nsrc].rank = i;
      xmap->src_msg[xmap->nsrc].intersect = intersect;
      xmap->nsrc++;
    }
    else
      xt_idxvec_delete(intersect);
  }

  for (i = 0; i < comm_size; ++i) {

    xt_idxvec_delete(src_idxvecs_all[i]);
    xt_idxvec_delete(dst_idxvecs_all[i]);
  }

  free(src_idxvecs_all);
  free(dst_idxvecs_all);

  // clean up
  free(send_buffer);
  free(recv_buffer_header);
  free(recv_buffer);
  free(displ);
}

static void generate_src_xlist(xt_xmap xmap,
                               Xt_idxvec idxvec) {

  unsigned i;
  int j;

  Xt_idx * index_list = NULL;
  int index_list_size = 0;

  // for all destination messages
  for (i = 0; i < xmap->nsrc; ++i) {

    int num_intersection_indices;

    num_intersection_indices =
      xt_idxvec_get_num_indices(xmap->src_msg[i].intersect);

    if (num_intersection_indices > index_list_size) {

      index_list_size = num_intersection_indices;

      index_list = realloc(index_list, index_list_size * sizeof(*index_list));
      if (index_list == NULL) {
        fputs("malloc failed for index list"
              " of the exchange list\n", stderr);
        exit(EXIT_FAILURE);
      }
    }

    // for all indices in the intersection
    for (j = 0; j < num_intersection_indices; ++j) {

      // look up current intersection index in the destination index vector
      int global_index;
      int local_index;

      assert(xt_idxvec_get_index(xmap->src_msg[i].intersect, j, &global_index) != 1);
      assert(xt_idxvec_get_position_of_index(idxvec, global_index, index_list+j) != 1);
    }

    xmap->src_msg[i].xlist = xt_idxlist_new(index_list, num_intersection_indices);
  }

  free(index_list);
}

static void generate_dst_xlist(xt_xmap xmap,
                               Xt_idxvec idxvec) {

  unsigned i;
  int j;

  Xt_idx * index_list = NULL;
  int index_list_size = 0;

  // for all destination messages
  for (i = 0; i < xmap->ndst; ++i) {

    int num_intersection_indices;

    num_intersection_indices =
      xt_idxvec_get_num_indices(xmap->dst_msg[i].intersect);

    if (num_intersection_indices > index_list_size) {

      index_list_size = num_intersection_indices;

      index_list = realloc(index_list, index_list_size * sizeof(*index_list));
      if (index_list == NULL) {
        fputs("malloc failed for index list"
              " of the exchange list\n", stderr);
        exit(EXIT_FAILURE);
      }
    }

    int prev_global_index;
    assert(xt_idxvec_get_index(xmap->dst_msg[i].intersect, 0, &prev_global_index) != 1);
    prev_global_index++;

    int repetition_count = 0;
    int repetition_offset = 0;

    // for all indices in the intersection
    for (j = 0; j < num_intersection_indices; ++j) {

      // look up current intersection index in the destination index vector
      int global_index;

      assert(xt_idxvec_get_index(xmap->dst_msg[i].intersect, j, &global_index) != 1);

      if (global_index != prev_global_index) {
        repetition_offset = 0;
        repetition_count = 0;
      }

      assert(xt_idxvec_get_position_of_index_off(idxvec, global_index,
                                                 index_list+j, repetition_offset) != 1);

      prev_global_index = global_index;
      repetition_offset = index_list[j] + 1;
      ++repetition_count;
    }

    xmap->dst_msg[i].xlist = xt_idxlist_new(index_list, num_intersection_indices);
  }

  free(index_list);
}

static void generate_xlist(xt_xmap xmap,
                           Xt_idxvec src_idxvec_local,
                           Xt_idxvec dst_idxvec_local) {

  generate_src_xlist(xmap, src_idxvec_local);
  generate_dst_xlist(xmap, dst_idxvec_local);
}

struct xt_xmap*
xt_xmap_new(xt_idxvec src_idxvec, xt_idxvec dst_idxvec, MPI_Comm comm) {

  struct xt_xmap * xmap;

  xmap = malloc(1 * sizeof(*xmap));
  if (xmap == NULL) {
    fprintf(stderr, "malloc failed for exchange map\n");
    exit(EXIT_FAILURE);
  }

  mpi_err_handler(MPI_Comm_dup(comm, &(xmap->comm)), comm);

  // exchange index vectors between all processes in comm
  exchange_idxvecs(xmap, src_idxvec, dst_idxvec);

  // generate exchange lists
  generate_xlist(xmap, src_idxvec, dst_idxvec);

  return xmap;
}

void xt_xmap_delete(xt_xmap xmap) {

  unsigned i;

  for (i = 0; i < xmap->ndst; ++i) {
    xt_idxvec_delete(xmap->dst_msg[i].intersect);
    xt_idxvec_delete(xmap->dst_msg[i].xlist);
  }

  for (i = 0; i < xmap->nsrc; ++i) {
    xt_idxvec_delete(xmap->src_msg[i].intersect);
    xt_idxvec_delete(xmap->src_msg[i].xlist);
  }

  free(xmap->dst_msg);
  free(xmap->src_msg);

  mpi_err_handler(MPI_Comm_free(&(xmap->comm)), MPI_COMM_WORLD);
  free(xmap);
}

static Xt_idxvec get_xlist(struct exchange_data * msgs, int num_msgs, int rank) {

  int i;

  for (i = 0; i < num_msgs; ++i)
    if (rank == msgs[i].rank)
      return xt_idxvec_copy(msgs[i].xlist);

  return NULL;
}

Xt_idxvec xt_xmap_get_dst_xlist(Xt_xmap xmap, int rank) {

  return get_xlist(xmap->dst_msg, xmap->ndst, rank);
}

Xt_idxvec xt_xmap_get_src_xlist(Xt_xmap xmap, int rank) {

  return get_xlist(xmap->src_msg, xmap->nsrc, rank);
}

Xt_xmap_iter xt_xmap_get_destination_iterator(Xt_xmap xmap) {

  if (xmap->ndst == 0)
    return NULL;

  Xt_xmap_iter iter;

  iter = malloc(1 * sizeof(*iter));
  if (iter == NULL) {
    fputs("malloc failed for exchange map destination iterator\n", stderr);
     exit(EXIT_FAILURE);
  }

  iter->msg = xmap->dst_msg;
  iter->msgs_left = xmap->ndst - 1;
}

Xt_xmap_iter xt_xmap_get_source_iterator(Xt_xmap xmap) {

  if (xmap->nsrc == 0)
    return NULL;

  Xt_xmap_iter iter;

  iter = malloc(1 * sizeof(*iter));
  if (iter == NULL) {
    fputs("malloc failed for exchange map source iterator\n", stderr);
     exit(EXIT_FAILURE);
  }

  iter->msg = xmap->src_msg;
  iter->msgs_left = xmap->nsrc - 1;
}

int xt_xmap_iterator_next(Xt_xmap_iter iter) {

  if (iter == NULL || iter->msgs_left == 0)
    return 0;

  iter->msg++;
  iter->msgs_left--;
}

int xt_xmap_iterator_get_rank(Xt_xmap_iter iter) {

  assert(iter != NULL);

  return iter->msg->rank;
}

Xt_idxvec xt_xmap_iterator_get_xlist(Xt_xmap_iter iter) {

  assert(iter != NULL);

  return  xt_idxvec_copy(iter->msg->xlist);
}

void xt_xmap_iterator_delete(Xt_xmap_iter iter) {

  free(iter);
}
