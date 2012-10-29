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
#include "xt_xmap_intersection.h"

static MPI_Comm     xmap_intersection_get_communicator(Xt_xmap xmap);
static int          xmap_intersection_get_num_destinations(Xt_xmap xmap);
static int          xmap_intersection_get_num_sources(Xt_xmap xmap);
static void         xmap_intersection_get_destination_ranks(Xt_xmap xmap, int * ranks);
static void         xmap_intersection_get_source_ranks(Xt_xmap xmap, int * ranks);
static Xt_xmap_iter xmap_intersection_get_destination_iterator(Xt_xmap xmap);
static Xt_xmap_iter xmap_intersection_get_source_iterator(Xt_xmap xmap);
static void         xmap_intersection_delete(Xt_xmap xmap);
static int          xmap_iterator_intersection_next(Xt_xmap_iter iter);
static int          xmap_intersection_iterator_get_rank(Xt_xmap_iter iter);
int const *         xmap_intersection_iterator_get_transfer_pos(Xt_xmap_iter iter);
Xt_count            xmap_intersection_iterator_get_num_transfer_pos(Xt_xmap_iter iter);
static void         xmap_intersection_iterator_delete(Xt_xmap_iter iter);
static int          xmap_intersection_get_max_src_pos(Xt_xmap xmap);
static int          xmap_intersection_get_max_dst_pos(Xt_xmap xmap);


static const struct Xt_xmap_iterator_vtable xmap_iterator_intersection_vtable = {
        .next                 = xmap_iterator_intersection_next,
        .get_rank             = xmap_intersection_iterator_get_rank,
        .get_transfer_pos     = xmap_intersection_iterator_get_transfer_pos,
        .get_num_transfer_pos = xmap_intersection_iterator_get_num_transfer_pos,
        .delete               = xmap_intersection_iterator_delete};


struct Xt_xmap_iterator_intersection {

  const struct Xt_xmap_iterator_vtable * vtable;

  struct exchange_data * msg;
  int msgs_left;
};

static const struct Xt_xmap_vtable xmap_intersection_vtable = {
        .get_communicator         = xmap_intersection_get_communicator,
        .get_num_destinations     = xmap_intersection_get_num_destinations,
        .get_num_sources          = xmap_intersection_get_num_sources,
        .get_destination_ranks    = xmap_intersection_get_destination_ranks,
        .get_source_ranks         = xmap_intersection_get_source_ranks,
        .get_destination_iterator = xmap_intersection_get_destination_iterator,
        .get_source_iterator      = xmap_intersection_get_source_iterator,
        .delete                   = xmap_intersection_delete,
        .get_max_src_pos          = xmap_intersection_get_max_src_pos,
        .get_max_dst_pos          = xmap_intersection_get_max_dst_pos};

struct exchange_data {
  // list of relative positions in memory to send or receive
  int * transfer_pos;
  Xt_count num_transfer_pos;
  int rank;
};

struct Xt_xmap_intersection {

  const struct Xt_xmap_vtable * vtable;

  struct exchange_data *dst_msg, *src_msg;
  int ndst, nsrc;

  // we need the max position in order to enable quick range-checks
  // for xmap-users like redist
  int max_src_pos; // max possible pos over all src transfer_pos (always >= 0)
  int max_dst_pos; // same for dst

  MPI_Comm comm;
};

static MPI_Comm xmap_intersection_get_communicator(Xt_xmap xmap) {

  struct Xt_xmap_intersection * xmap_intersection;

  xmap_intersection = (struct Xt_xmap_intersection*)xmap;

  return xmap_intersection->comm;
}

static int xmap_intersection_get_num_destinations(Xt_xmap xmap) {

  struct Xt_xmap_intersection * xmap_intersection;

  xmap_intersection = (struct Xt_xmap_intersection*)xmap;

  return xmap_intersection->ndst;
}

static int xmap_intersection_get_num_sources(Xt_xmap xmap) {

  struct Xt_xmap_intersection * xmap_intersection;

  xmap_intersection = (struct Xt_xmap_intersection*)xmap;

  return xmap_intersection->nsrc;
}

static void xmap_intersection_get_destination_ranks(Xt_xmap xmap, int * ranks) {

  struct Xt_xmap_intersection * xmap_intersection;

  xmap_intersection = (struct Xt_xmap_intersection*)xmap;

  int i;

  for (i = 0; i < xmap_intersection->ndst; ++i)
    ranks[i] = xmap_intersection->dst_msg[i].rank;
}

static void xmap_intersection_get_source_ranks(Xt_xmap xmap, int * ranks) {

  struct Xt_xmap_intersection * xmap_intersection;

  xmap_intersection = (struct Xt_xmap_intersection*)xmap;

  int i;

  for (i = 0; i < xmap_intersection->nsrc; ++i)
    ranks[i] = xmap_intersection->src_msg[i].rank;
}

/* compute list positions for one direction */
static void generate_dir_transfer_pos(Xt_idxlist * intersections,
                                      int * ranks,
                                      Xt_count num_intersections,
                                      Xt_idxlist mypart_idxlist,
                                      const int single_match_only,
                                      int *resCount,
                                      struct exchange_data **resSets)
{

  *resCount = num_intersections;
  *resSets = xmalloc(num_intersections * sizeof(**resSets));

  for (int i = 0; i < num_intersections; ++i) {

    const Xt_idx *intersection_idxvec
      = xt_idxlist_get_indices_const(intersections[i]);
    Xt_count intersection_size = xt_idxlist_get_num_indices(intersections[i]);
    int *intersection_pos = xmalloc(intersection_size
                                    * sizeof(*intersection_pos));

    int retval;
    retval = xt_idxlist_get_positions_of_indices(
      mypart_idxlist, intersection_idxvec, intersection_size, intersection_pos,
      single_match_only);
    assert(retval != 1);

    (*resSets)[i].transfer_pos = intersection_pos;
    (*resSets)[i].num_transfer_pos = intersection_size;
    (*resSets)[i].rank = ranks[i];
  }
}


static void generate_transfer_pos(struct Xt_xmap_intersection * xmap,
                                  Xt_idxlist * src_intersections,
                                  int * src_ranks,
                                  Xt_count num_src_intersections,
                                  Xt_idxlist * dst_intersections,
                                  int * dst_ranks,
                                  Xt_count num_dst_intersections,
                                  Xt_idxlist src_idxlist_local,
                                  Xt_idxlist dst_idxlist_local) {

  generate_dir_transfer_pos(src_intersections, src_ranks,
                            num_src_intersections, src_idxlist_local, 0,
                            &(xmap->nsrc), &(xmap->src_msg));
  generate_dir_transfer_pos(dst_intersections, dst_ranks,
                            num_dst_intersections, dst_idxlist_local, 1,
                            &(xmap->ndst), &(xmap->dst_msg));
}

Xt_xmap xt_xmap_intersection_new(Xt_idxlist * src_intersections, int * src_ranks,
                                 Xt_count num_src_intersections,
                                 Xt_idxlist * dst_intersections, int * dst_ranks,
                                 Xt_count num_dst_intersections,
                                 Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                                 MPI_Comm comm) {

  struct Xt_xmap_intersection * xmap;

  xmap = xmalloc(1 * sizeof(*xmap));

  xmap->vtable = &xmap_intersection_vtable;

  mpi_err_handler(MPI_Comm_dup(comm, &(xmap->comm)), comm);

  // generate exchange lists
  generate_transfer_pos(xmap, src_intersections, src_ranks, num_src_intersections,
                        dst_intersections, dst_ranks, num_dst_intersections, src_idxlist,
                        dst_idxlist);

  // we could also calculate the (more precise) max pos using only xmap data
  // but using this simple estimate we are still okay for usage checks
  xmap->max_src_pos = xt_idxlist_get_num_indices(src_idxlist);
  xmap->max_dst_pos = xt_idxlist_get_num_indices(dst_idxlist);

  return (Xt_xmap)xmap;
}

static int xmap_intersection_get_max_src_pos(Xt_xmap xmap) {
  struct Xt_xmap_intersection *xmap_intersection = (struct Xt_xmap_intersection *)xmap;
  return xmap_intersection->max_src_pos;
}

static int xmap_intersection_get_max_dst_pos(Xt_xmap xmap) {
  struct Xt_xmap_intersection *xmap_intersection = (struct Xt_xmap_intersection *)xmap;
  return xmap_intersection->max_dst_pos;
}


static void xmap_intersection_delete(Xt_xmap xmap) {

  struct Xt_xmap_intersection * xmap_intersection;

  xmap_intersection = (struct Xt_xmap_intersection*)xmap;

  unsigned i;

  for (i = 0; i < xmap_intersection->ndst; ++i)
    free(xmap_intersection->dst_msg[i].transfer_pos);

  for (i = 0; i < xmap_intersection->nsrc; ++i)
    free(xmap_intersection->src_msg[i].transfer_pos);

  free(xmap_intersection->dst_msg);
  free(xmap_intersection->src_msg);

  mpi_err_handler(MPI_Comm_free(&(xmap_intersection->comm)), MPI_COMM_WORLD);
  free(xmap_intersection);
}

static Xt_xmap_iter xmap_intersection_get_destination_iterator(Xt_xmap xmap) {

  struct Xt_xmap_intersection * xmap_intersection;

  xmap_intersection = (struct Xt_xmap_intersection*)xmap;

  if (xmap_intersection->ndst == 0)
    return NULL;

  struct Xt_xmap_iterator_intersection * iter;

  iter = xmalloc(1 * sizeof(*iter));

  iter->vtable = &xmap_iterator_intersection_vtable;
  iter->msg = xmap_intersection->dst_msg;
  iter->msgs_left = xmap_intersection->ndst - 1;

  return (Xt_xmap_iter)iter;
}

static Xt_xmap_iter xmap_intersection_get_source_iterator(Xt_xmap xmap) {

  struct Xt_xmap_intersection * xmap_intersection;

  xmap_intersection = (struct Xt_xmap_intersection*)xmap;

  if (xmap_intersection->nsrc == 0)
    return NULL;

  struct Xt_xmap_iterator_intersection * iter;

  iter = xmalloc(1 * sizeof(*iter));

  iter->vtable = &xmap_iterator_intersection_vtable;
  iter->msg = xmap_intersection->src_msg;
  iter->msgs_left = xmap_intersection->nsrc - 1;

  return (Xt_xmap_iter)iter;
}

static int xmap_iterator_intersection_next(Xt_xmap_iter iter) {

  struct Xt_xmap_iterator_intersection * iter_intersection;

  iter_intersection = (struct Xt_xmap_iterator_intersection *)iter;

  if (iter_intersection == NULL || iter_intersection->msgs_left == 0)
    return 0;

  iter_intersection->msg++;
  iter_intersection->msgs_left--;

  return 1;
}

static int xmap_intersection_iterator_get_rank(Xt_xmap_iter iter) {

  struct Xt_xmap_iterator_intersection * iter_intersection;

  iter_intersection = (struct Xt_xmap_iterator_intersection *)iter;

  assert(iter_intersection != NULL);

  return iter_intersection->msg->rank;
}

int const * xmap_intersection_iterator_get_transfer_pos(Xt_xmap_iter iter) {

  struct Xt_xmap_iterator_intersection * iter_intersection;

  iter_intersection = (struct Xt_xmap_iterator_intersection *)iter;

  assert(iter_intersection != NULL);

  return  iter_intersection->msg->transfer_pos;
}

Xt_count xmap_intersection_iterator_get_num_transfer_pos(Xt_xmap_iter iter) {

  struct Xt_xmap_iterator_intersection * iter_intersection;

  iter_intersection = (struct Xt_xmap_iterator_intersection *)iter;

  assert(iter_intersection != NULL);

  return iter_intersection->msg->num_transfer_pos;
}

static void xmap_intersection_iterator_delete(Xt_xmap_iter iter) {

  struct Xt_xmap_iterator_intersection * iter_intersection;

  iter_intersection = (struct Xt_xmap_iterator_intersection *)iter;

  free(iter_intersection);
}
