#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt_mpi.h"
#include "xt_redist_collection.h"
#include "ensure_array_size.h"
#include "xt_redist.h"

#define DATATYPE_CACHE_SIZE (16)

static void
redist_collection_delete(Xt_redist redist);

static void
redist_collection_s_exchange(Xt_redist redist, void **src_data,
                             unsigned num_src_arrays, void **dst_data,
                             unsigned num_dst_arrays);

static void
redist_collection_s_exchange1(Xt_redist redist, void *src_data, void *dst_data);

static MPI_Datatype
redist_collection_get_send_MPI_Datatype(Xt_redist redist, int rank);

static MPI_Datatype
redist_collection_get_recv_MPI_Datatype(Xt_redist redist, int rank);

static const struct xt_redist_vtable redist_collection_vtable = {
  .delete                = redist_collection_delete,
  .s_exchange            = redist_collection_s_exchange,
  .s_exchange1           = redist_collection_s_exchange1,
  .get_send_MPI_Datatype = redist_collection_get_send_MPI_Datatype,
  .get_recv_MPI_Datatype = redist_collection_get_recv_MPI_Datatype
};

struct redist_collection_msg {

  int rank;
  MPI_Datatype * datatypes; // datatypes of the redists (size == num_redists)
  MPI_Datatype datatype_cache[DATATYPE_CACHE_SIZE];
};

struct Xt_redist_collection {

  const struct xt_redist_vtable *vtable;

  unsigned num_redists;

  MPI_Aint * cached_src_displacements[DATATYPE_CACHE_SIZE];
  MPI_Aint * cached_dst_displacements[DATATYPE_CACHE_SIZE];
  int cache_src_token;
  int cache_dst_token;

  int ndst, nsrc;
  struct redist_collection_msg * send_msgs;
  struct redist_collection_msg * recv_msgs;

  MPI_Comm comm;
};

static void get_MPI_datatypes(struct redist_collection_msg ** msgs, int * nmsgs,
                              Xt_redist * redists, unsigned num_redists, MPI_Comm comm,
                              MPI_Datatype (*get_MPI_datatype)(Xt_redist,int)) {

  int msgs_array_size = 0;

  int comm_size;
  mpi_err_handler(MPI_Comm_size(comm, &comm_size), comm);

  MPI_Datatype datatypes[num_redists];

  for (int i = 0; i < comm_size; ++i) {

    int flag = 0;

    for (unsigned j = 0; j < num_redists; ++j)
      flag |= (datatypes[j] = get_MPI_datatype(redists[j], i)) != MPI_DATATYPE_NULL;

    if (flag) {

        ENSURE_ARRAY_SIZE(*msgs, msgs_array_size, *nmsgs+1);

        (*msgs)[*nmsgs].rank = i;
        (*msgs)[*nmsgs].datatypes = xmalloc(num_redists *
          sizeof(*((*msgs)[*nmsgs].datatypes)));
        memcpy((*msgs)[*nmsgs].datatypes, datatypes, num_redists * sizeof(*datatypes));

        for (unsigned k = 0; k < DATATYPE_CACHE_SIZE; ++k)
          (*msgs)[*nmsgs].datatype_cache[k] = MPI_DATATYPE_NULL;

        ++*nmsgs;
    }
  }

  if (*nmsgs > 0)
    *msgs = xrealloc(*msgs, *nmsgs * sizeof(**msgs));
}

Xt_redist xt_redist_collection_new(Xt_redist * redists, unsigned num_redists,
                                   MPI_Comm comm) {

  struct Xt_redist_collection * redist_coll;

  redist_coll = xmalloc(1 * sizeof(*redist_coll));

  redist_coll->vtable = &redist_collection_vtable;
  redist_coll->num_redists = num_redists;
  redist_coll->ndst = 0;
  redist_coll->nsrc = 0;
  redist_coll->send_msgs = NULL;
  redist_coll->recv_msgs = NULL;
  for (int i = 0; i < DATATYPE_CACHE_SIZE; ++i) {
    redist_coll->cached_src_displacements[i] = NULL;
    redist_coll->cached_dst_displacements[i] = NULL;
  }
  redist_coll->cache_src_token = 0;
  redist_coll->cache_dst_token = 0;
  mpi_err_handler(MPI_Comm_dup(comm, &(redist_coll->comm)), comm);

  get_MPI_datatypes(&(redist_coll->send_msgs), &(redist_coll->nsrc), redists,
                    num_redists, redist_coll->comm, xt_redist_get_send_MPI_Datatype);

  get_MPI_datatypes(&(redist_coll->recv_msgs), &(redist_coll->ndst), redists,
                    num_redists, redist_coll->comm, xt_redist_get_recv_MPI_Datatype);

  return (Xt_redist)redist_coll;
}

static void generate_datatype(MPI_Aint * displacements, int * block_lengths, 
                              struct redist_collection_msg * msg,
                              unsigned num_redists, int cache_index,
                              MPI_Comm comm) {

  MPI_Datatype * datatype;

  datatype = msg->datatype_cache+cache_index;

  if (*datatype != MPI_DATATYPE_NULL) mpi_err_handler(MPI_Type_free(datatype), comm);

  int num_datatypes = 0;

  for (int i = 0; i < num_redists; ++i)
    if (msg->datatypes[i] != MPI_DATATYPE_NULL)
      ++num_datatypes;

  MPI_Datatype * datatypes;
  MPI_Aint * displacements_;

  if (num_datatypes != num_redists) {

    datatypes = xmalloc(num_datatypes * sizeof(*datatypes));
    displacements_ = xmalloc(num_datatypes * sizeof(*displacements));

    num_datatypes = 0;

    for (int i = 0; i < num_redists; ++i) {
      if (msg->datatypes[i] != MPI_DATATYPE_NULL) {

        datatypes[num_datatypes] = msg->datatypes[i];
        displacements_[num_datatypes] = displacements[i];
        ++num_datatypes;
      }
    }
  } else {

    datatypes = msg->datatypes;
    displacements_ = displacements;
  }

  mpi_err_handler(MPI_Type_create_struct(num_datatypes, block_lengths, displacements_,
                                         datatypes, datatype), comm);

  mpi_err_handler(MPI_Type_commit(datatype), comm);

  if (num_datatypes != num_redists) {
    free(datatypes);
    free(displacements_);
  }
}

static int generate_datatypes(void ** data, struct redist_collection_msg * msgs, 
                              int num_messages, unsigned num_redists, int * cache_token,
                              MPI_Aint ** cached_displacements, MPI_Comm comm) {

  if (*cached_displacements == NULL)
    *cached_displacements = xmalloc(num_redists * sizeof(**cached_displacements));

  MPI_Aint addresses[num_redists];
  int block_lengths[num_redists];

  for (unsigned i = 0; i < num_redists; ++i)
    mpi_err_handler(MPI_Get_address(data[i], addresses+i), comm);

  for (unsigned i = 0; i < num_redists; ++i)
    (*cached_displacements)[i] = addresses[i] - addresses[0];

  for (unsigned i = 0; i < num_redists; ++i)
    block_lengths[i] = 1;

  for (int i = 0; i < num_messages; ++i)
    generate_datatype(*cached_displacements, block_lengths, msgs+i, num_redists,
                      *cache_token, comm);

  int cache_index = *cache_token;

  *cache_token = (*cache_token+1)%DATATYPE_CACHE_SIZE;

  return cache_index;
}

int lookup_cache_index(void ** data, MPI_Aint ** cached_displacements,
                       unsigned num_redists, MPI_Comm comm) {

  MPI_Aint addresses[num_redists];
  MPI_Aint displacements[num_redists];

  for (unsigned i = 0; i < num_redists; ++i)
    mpi_err_handler(MPI_Get_address(data[i], addresses+i), comm);

  for (unsigned i = 0; i < num_redists; ++i)
    displacements[i] = addresses[i] - addresses[0];

  for (int i = 0; i < DATATYPE_CACHE_SIZE; ++i) {
    if (cached_displacements[i] != NULL) {
      unsigned j;
      for (j = 0; j < num_redists; ++j)
        if (displacements[j] != cached_displacements[i][j]) break;
      if (j == num_redists) return i;
    } else {
      break;
    }
  }
  return -1;
}

static int get_cache_index(void ** data, struct redist_collection_msg * msgs, 
                           int num_messages, unsigned num_redists, int * cache_token,
                           MPI_Aint ** cached_displacements, MPI_Comm comm) {

  int cache_index = lookup_cache_index(data, cached_displacements, num_redists, comm);

  if (cache_index == -1)
    cache_index = generate_datatypes(data, msgs, num_messages, num_redists, cache_token,
                                     cached_displacements+*cache_token, comm);

  return cache_index;
}

static void
redist_collection_s_exchange(Xt_redist redist, void **src_data,
                             unsigned num_src_arrays, void **dst_data,
                             unsigned num_dst_arrays) {

  struct Xt_redist_collection * redist_coll;

  redist_coll = (struct Xt_redist_collection *) redist;

  if (num_src_arrays != redist_coll->num_redists ||
      num_dst_arrays != redist_coll->num_redists)
    Xt_abort(redist_coll->comm, "ERROR: wrong number of array in "
             "redist_collection_s_exchange", __FILE__, __LINE__);

  MPI_Request * recv_requests;

  recv_requests = xmalloc(redist_coll->ndst * sizeof(*recv_requests));

  int dst_cache_index = get_cache_index(dst_data, redist_coll->recv_msgs,
                                        redist_coll->ndst, redist_coll->num_redists,
                                        &(redist_coll->cache_dst_token),
                                        redist_coll->cached_dst_displacements,
                                        redist_coll->comm);

  for (int i = 0; i < redist_coll->ndst; ++i)
    mpi_err_handler(MPI_Irecv(dst_data[0], 1,
                              redist_coll->recv_msgs[i].datatype_cache[dst_cache_index],
                              redist_coll->recv_msgs[i].rank, 0, redist_coll->comm,
                              recv_requests+i), redist_coll->comm);

  int src_cache_index = get_cache_index(src_data, redist_coll->send_msgs,
                                        redist_coll->nsrc, redist_coll->num_redists,
                                        &(redist_coll->cache_src_token),
                                        redist_coll->cached_src_displacements,
                                        redist_coll->comm);

  for (int i = 0; i < redist_coll->nsrc; ++i)
    mpi_err_handler(MPI_Send(src_data[0], 1,
                             redist_coll->send_msgs[i].datatype_cache[src_cache_index],
                             redist_coll->send_msgs[i].rank, 0, redist_coll->comm),
                    redist_coll->comm);

  mpi_err_handler(MPI_Waitall(redist_coll->ndst, recv_requests,
                              MPI_STATUSES_IGNORE), redist_coll->comm);

  free(recv_requests);
}

static void
free_redist_collection_msgs(struct redist_collection_msg * msgs,
                            int nmsgs, unsigned num_redists, MPI_Comm comm) {

  for (int i = 0; i < nmsgs; ++i) {

    for (unsigned j = 0; j < num_redists; ++j)
      if (msgs[i].datatypes[j] != MPI_DATATYPE_NULL)
        mpi_err_handler(MPI_Type_free(msgs[i].datatypes+j), comm);
    free(msgs[i].datatypes);

    for (unsigned j = 0; j < DATATYPE_CACHE_SIZE; ++j) {

      if (msgs[i].datatype_cache[j] != MPI_DATATYPE_NULL)
        mpi_err_handler(MPI_Type_free(msgs[i].datatype_cache+j), comm);
      else
        break;
    }
  }
}

static void
redist_collection_delete(Xt_redist redist) {

  struct Xt_redist_collection * redist_coll;

  redist_coll = (struct Xt_redist_collection *)redist;

  free_redist_collection_msgs(redist_coll->send_msgs, redist_coll->nsrc,
                              redist_coll->num_redists, redist_coll->comm);
  free(redist_coll->send_msgs);

  free_redist_collection_msgs(redist_coll->recv_msgs, redist_coll->ndst,
                              redist_coll->num_redists, redist_coll->comm);
  free(redist_coll->recv_msgs);

  for (int i = 0; i < DATATYPE_CACHE_SIZE; ++i) {
    free(redist_coll->cached_src_displacements[i]);
    free(redist_coll->cached_dst_displacements[i]);
  }

  mpi_err_handler(MPI_Comm_free(&(redist_coll->comm)), MPI_COMM_WORLD);

  free(redist_coll);
}

static MPI_Datatype
redist_collection_get_send_MPI_Datatype(Xt_redist redist, int rank) {

  struct Xt_redist_collection * redist_coll;

  redist_coll = (struct Xt_redist_collection *)redist;

  Xt_abort(redist_coll->comm, "ERROR: get_send_MPI_Datatype is not"
           " supported for this xt_redist type (xt_redist_collection)",
           __FILE__, __LINE__);

  return MPI_DATATYPE_NULL;
}

static MPI_Datatype
redist_collection_get_recv_MPI_Datatype(Xt_redist redist, int rank) {

  struct Xt_redist_collection * redist_coll;

  redist_coll = (struct Xt_redist_collection *)redist;

  Xt_abort(redist_coll->comm, "ERROR: get_recv_MPI_Datatype is not"
           " supported for this xt_redist type (xt_redist_collection)",
           __FILE__, __LINE__);

  return MPI_DATATYPE_NULL;
}

static void
redist_collection_s_exchange1(Xt_redist redist, void *src_data, void *dst_data) {

  struct Xt_redist_collection * redist_coll;

  redist_coll = (struct Xt_redist_collection *)redist;

  Xt_abort(redist_coll->comm, "ERROR: s_exchange1 is not implemented for"
           " this xt_redist type (xt_redist_collection)", __FILE__, __LINE__);
}
