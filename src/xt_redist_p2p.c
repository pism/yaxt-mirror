#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <mpi.h>

#include "xt_mpi.h"
#include "xt_redist_p2p.h"
#include "xt_xmap.h"
#include "xt_idxlist.h"
#include "core/ppm_xfuncs.h"

static void
redist_p2p_delete(Xt_redist redist);

static void
redist_p2p_s_exchange(Xt_redist redist, void **src_data, unsigned num_src_arrays,
                                        void **dst_data, unsigned num_dst_arrays);

static const struct xt_redist_vtable redist_p2p_vtable = {redist_p2p_delete,
                                                          redist_p2p_s_exchange};

struct Xt_redist_p2p_msg {

  int rank;
  MPI_Datatype datatype;
};

struct Xt_redist_p2p {

  const struct xt_redist_vtable *vtable;

  int ndst, nsrc;

  struct Xt_redist_p2p_msg * send_msgs;
  struct Xt_redist_p2p_msg * recv_msgs;

  MPI_Comm comm;
};

static MPI_Datatype
generate_datatype(Xt_idxlist xlist, int *offsets,
                  MPI_Datatype base_datatype, MPI_Comm comm) {

  MPI_Datatype type;

  Xt_count num_indices;
  int * displ;

  num_indices = xt_idxlist_get_num_indices(xlist);

  displ = xmalloc(num_indices * sizeof(*displ));

  assert(sizeof(Xt_idx) == sizeof(int));

  xt_idxlist_get_indices(xlist, displ);

  if (offsets != NULL) {
    Xt_count i;
    for (i = 0; i < num_indices; ++i)
      displ[i] = offsets[displ[i]];
  }

  mpi_err_handler(MPI_Type_create_indexed_block(num_indices, 1, displ,
                                                base_datatype, &type), comm);
  mpi_err_handler(MPI_Type_commit(&type), comm);

  free(displ);

  return type;
}

static void
generate_msg_infos(int num_msgs, Xt_xmap_iter iter, int *offsets,
                   MPI_Datatype base_datatype, struct Xt_redist_p2p_msg ** msgs,
                   MPI_Comm comm) {

  if (num_msgs <= 0) {
    *msgs = NULL;
    return;
  }

  *msgs = xmalloc(num_msgs * sizeof(**msgs));

  struct Xt_redist_p2p_msg * curr_msg;

  curr_msg = *msgs;

  do {

    Xt_idxlist curr_xlist;

    curr_xlist = xt_xmap_iterator_get_xlist(iter);

    curr_msg->datatype = generate_datatype(curr_xlist, offsets, base_datatype, comm);
    curr_msg->rank = xt_xmap_iterator_get_rank(iter);

    curr_msg++;

    xt_idxlist_delete(curr_xlist);

  } while (xt_xmap_iterator_next(iter));

  xt_xmap_iterator_delete(iter);
}

Xt_redist xt_redist_p2p_off_new(Xt_xmap xmap, int *src_offsets,
                                int *dst_offsets, MPI_Datatype datatype) {

  struct Xt_redist_p2p *redist;

  redist = xmalloc(1 * sizeof(*redist));

  redist->vtable = &redist_p2p_vtable;

  MPI_Comm xmap_comm;

  xmap_comm = xt_xmap_get_communicator(xmap);

  mpi_err_handler(MPI_Comm_dup(xmap_comm, &(redist->comm)), xmap_comm);

  redist->ndst = xt_xmap_get_num_destinations(xmap);
  redist->nsrc = xt_xmap_get_num_sources(xmap);

  generate_msg_infos(redist->ndst, xt_xmap_get_destination_iterator(xmap),
                     dst_offsets, datatype, &(redist->recv_msgs), redist->comm);

  generate_msg_infos(redist->nsrc, xt_xmap_get_source_iterator(xmap),
                     src_offsets, datatype, &(redist->send_msgs), redist->comm);

  return (Xt_redist)redist;
}

Xt_redist xt_redist_p2p_new(Xt_xmap xmap, MPI_Datatype datatype) {

  return xt_redist_p2p_off_new(xmap, NULL, NULL, datatype);
}

static void
redist_p2p_delete(Xt_redist redist) {

  struct Xt_redist_p2p *redist_p2p;

  redist_p2p = (struct Xt_redist_p2p *)redist;

  int i;

  for (i = 0; i < redist_p2p->nsrc; ++i)
    MPI_Type_free(&(redist_p2p->send_msgs[i].datatype));
  free(redist_p2p->send_msgs);

  for (i = 0; i < redist_p2p->ndst; ++i)
    MPI_Type_free(&(redist_p2p->recv_msgs[i].datatype));
  free(redist_p2p->recv_msgs);

  MPI_Comm_free(&(redist_p2p->comm));

  free(redist_p2p);
}

static void
redist_p2p_s_exchange(Xt_redist redist, void **src_data, unsigned num_src_arrays,
                                        void **dst_data, unsigned num_dst_arrays) {

  struct Xt_redist_p2p *redist_p2p;

  redist_p2p = (struct Xt_redist_p2p *)redist;

  assert(num_src_arrays == 1 && num_dst_arrays == 1);

  MPI_Request * recv_request;

  recv_request = xmalloc(redist_p2p->ndst * sizeof(*recv_request));

  int i;

  for (i = 0; i < redist_p2p->ndst; ++i)
    mpi_err_handler(MPI_Irecv(dst_data[0], 1, redist_p2p->recv_msgs[i].datatype,
                              redist_p2p->recv_msgs[i].rank, 0, redist_p2p->comm,
                              recv_request+i), redist_p2p->comm);

  for (i = 0; i < redist_p2p->nsrc; ++i)
    mpi_err_handler(MPI_Send(src_data[0], 1, redist_p2p->send_msgs[i].datatype,
                             redist_p2p->send_msgs[i].rank, 0, redist_p2p->comm),
                    redist_p2p->comm);

  mpi_err_handler(MPI_Waitall(redist_p2p->ndst, recv_request,
                              MPI_STATUSES_IGNORE), redist_p2p->comm);

  free(recv_request);
}
