
#include "core/core.h"

#include "xt_idxlist.h"
#include "xt_idxvec.h"
#include "xt_xmap.h"
#include "xt_xmap_all2all.h"
#include "xt_idxstripes.h"
#include "xt_redist_p2p.h"

struct xt_idxlist_f {
   struct xt_idxlist *cptr;
};

struct xt_xmap_f {
   struct Xt_xmap *cptr;
};

struct xt_redist_f {
  struct xt_redist *cptr;
};

void xt_initialize_f(MPI_Fint *comm_f) {
  MPI_Comm comm_c;
  comm_c = MPI_Comm_f2c(*comm_f);
  Xt_initialize(comm_c);
}

void xt_abort_cmsl_f(MPI_Fint *comm_f, const char *msg, const char *source, int line) {
  MPI_Comm comm_c = MPI_Comm_f2c(*comm_f);
  Xt_abort(comm_c, msg, source, line);
}

void xt_abort_msl_f(const char *msg, const char *source, int line) {
  MPI_Comm comm_c = SymPrefix(default_comm);
  Xt_abort(comm_c, msg, source, line);
}

struct xt_idxlist_f xt_idxvec_new_f(Xt_idx *idxvec, Xt_count num_indices) {
  struct xt_idxlist_f idxlist_f;
  idxlist_f.cptr = xt_idxvec_new(idxvec, num_indices);
  return idxlist_f;
}

void xt_idxlist_delete_f(struct xt_idxlist_f *idxlist_f) {
  xt_idxlist_delete(idxlist_f->cptr);
}

struct xt_xmap_f xt_xmap_all2all_new_f(struct xt_idxlist_f *src_idxlist_f, struct xt_idxlist_f *dst_idxlist_f, MPI_Fint *comm_f) {
  MPI_Comm comm_c = MPI_Comm_f2c(*comm_f);
  struct xt_xmap_f xmap_f;
  xmap_f.cptr = xt_xmap_all2all_new(src_idxlist_f->cptr, dst_idxlist_f->cptr, comm_c);
  return xmap_f;
}

void xt_xmap_delete_f(struct xt_xmap_f *xmap_f) {
  xt_xmap_delete(xmap_f->cptr);
}

struct xt_redist_f xt_redist_p2p_blocks_new_f(struct xt_xmap_f *xmap_f,
                                              int *src_block_offsets, int *src_block_sizes,
                                              int *dst_block_offsets, int *dst_block_sizes,
                                              MPI_Fint *datatype_f) {
  MPI_Datatype datatype_c = MPI_Type_f2c(*datatype_f);
  struct xt_redist_f redist_f;

  redist_f.cptr =  xt_redist_p2p_blocks_new(xmap_f->cptr,
                                            src_block_offsets, src_block_sizes,
                                            dst_block_offsets, dst_block_sizes, datatype_c);

  return redist_f;
}

struct xt_redist_f xt_redist_p2p_off_new_f(struct xt_xmap_f *xmap_f, int *src_offsets,  int *dst_offsets, MPI_Fint *datatype_f) {
  MPI_Datatype datatype_c = MPI_Type_f2c(*datatype_f);
  struct xt_redist_f redist_f;

  redist_f.cptr =  xt_redist_p2p_off_new(xmap_f->cptr, src_offsets, dst_offsets, datatype_c);

  return redist_f;
}

struct xt_redist_f xt_redist_p2p_new_f(struct xt_xmap_f *xmap_f, MPI_Fint *datatype_f) {
  MPI_Datatype datatype_c = MPI_Type_f2c(*datatype_f);
  struct xt_redist_f redist_f;

  redist_f.cptr =  xt_redist_p2p_new(xmap_f->cptr, datatype_c);

  return redist_f;
}

void xt_redist_delete_f(struct xt_redist_f *redist_f) {
  xt_redist_delete(redist_f->cptr);
}

void xt_redist_s_exchange1_f(struct xt_redist_f *redist_f, void **src_data, void **dst_data) {
  xt_redist_s_exchange(redist_f->cptr, src_data, 1, dst_data, 1);
}

struct xt_idxlist_f xt_idxstripes_new_f(struct Xt_stripe const *stripes,  Xt_count num_stripes) {
  struct xt_idxlist_f idxlist_f;
  idxlist_f.cptr = xt_idxstripes_new(stripes, num_stripes);
  return idxlist_f;
}
