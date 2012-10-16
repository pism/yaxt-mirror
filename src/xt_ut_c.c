
//
// implementation of unitrans using yaxt
//

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "xt_mpi.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt_idxvec.h"
#include "xt_ut_c.h"
#include "xt_xmap.h"
#include "xt_xmap_all2all.h"
#include "xt_redist_p2p.h"

#include "xt_handles.h"

typedef struct {
  int decomp_size;
  int comm_tmpl_size;
  int comm_trans_size;
  int debug_lvl;
  int mode;
  int debug_unit;
} xt_ut_init_type;

static xt_ut_init_type init_val = {0,0,0,0,0,0};


typedef struct {
  Xt_idxlist idxvec;
  int domain_size;
} xt_ut_deco_type;

static const xt_ut_deco_type zero_deco = {.idxvec=NULL, .domain_size=0};

static struct xt_handle_set_type *deco_set;

typedef struct {
  Xt_xmap xmap;
} xt_ut_template_type;

static const xt_ut_template_type zero_template = {.xmap=NULL};

static struct xt_handle_set_type *template_set;


typedef struct {
  Xt_redist redist;
} xt_ut_trans_type;

static const xt_ut_trans_type zero_trans = {.redist = NULL};


static struct xt_handle_set_type *trans_set;

void xt_ut_abort(char *msg, char *source, int *line) {
  Xt_abort(MPI_COMM_WORLD, msg, source, *line);
}

void xt_ut_init(int decomp_size, int comm_tmpl_size, int comm_trans_size, int debug_lvl, int mode, int debug_unit) {

  init_val.decomp_size     = decomp_size;
  init_val.comm_tmpl_size  = comm_tmpl_size;
  init_val.comm_trans_size = comm_trans_size;
  init_val.debug_lvl       = debug_lvl;
  init_val.mode            = mode;
  init_val.debug_unit      = debug_unit;

  deco_set = xt_handle_set_new(init_val.decomp_size);

  template_set = xt_handle_set_new(init_val.comm_tmpl_size);

  trans_set = xt_handle_set_new(init_val.comm_trans_size);

}

static int new_ideco() {
  xt_ut_deco_type *data = xmalloc(sizeof(*data));
  return xt_handle_new(deco_set, data);
}

void xt_ut_destroy_decomposition(int id) {
  xt_ut_deco_type *deco = xt_handle2pointer(deco_set, id);
  if (!deco) die("xt_ut_destroy_decomposition: invalid handle");
  xt_idxlist_delete(deco->idxvec);
  xt_handle_delete(deco_set, id);
}


static int new_itemplate() {
  xt_ut_template_type *data = xmalloc(sizeof(*data));
  return xt_handle_new(template_set, data);
}

void xt_ut_destroy_transposition_template(int id) {
  xt_ut_template_type *template = xt_handle2pointer(template_set, id);
  if (!template) die("xt_ut_destroy_transposition_template: invalid handle");
  xt_xmap_delete(template->xmap);
  xt_handle_delete(template_set, id);
}


static int new_itrans() {
  xt_ut_trans_type *data = xmalloc(sizeof(*data));
  //if (!data) xt_die("malloc_faild", __FILE__, __LINE__);
  return xt_handle_new(trans_set, data);
}

void xt_ut_destroy_transposition(int id) {
  xt_ut_trans_type *trans = xt_handle2pointer(trans_set, id);
  if (!trans) die("xt_ut_destroy_transposition: invalid handle");
  xt_redist_delete(trans->redist);
  xt_handle_delete(trans_set, id);
}


MPI_Fint xt_ut_init_decomposition_1d(int *iv, int iv_n, int domain_size) {
  MPI_Fint id = new_ideco();
  xt_ut_deco_type *deco = xt_handle2pointer(deco_set, id);
  if (!deco) die("xt_ut_init_decomposition_1d: internal error");
  deco->idxvec = xt_idxvec_new(iv, iv_n);
  return id;
}


MPI_Fint xt_ut_init_oneway_transposition_template(int id_in, int id_out, int mpi_world, int icheck_unique) {
  xt_ut_deco_type *deco_in = xt_handle2pointer(deco_set, id_in);
  if(!deco_in) die("xt_ut_init_oneway_transposition_template: invalid handle (first argument)");

  xt_ut_deco_type *deco_out = xt_handle2pointer(deco_set, id_out);
  if(!deco_out) die("xt_ut_init_oneway_transposition_template: invalid handle (second argument)");

  int id = new_itemplate();
  xt_ut_template_type *template = xt_handle2pointer(template_set, id);
  if (!template) die("xt_ut_init_oneway_transposition_template: internal error");
  template->xmap = xt_xmap_all2all_new(deco_in->idxvec, deco_out->idxvec, MPI_COMM_WORLD);

  return id;
}


MPI_Fint xt_ut_init_transposition_simple(MPI_Fint itemplate, MPI_Fint f_datatype) {
  xt_ut_template_type *template = xt_handle2pointer(template_set, itemplate);
  if (!template) die("xt_ut_init_transposition_simple: invalid handle");

  MPI_Datatype datatype = MPI_Type_f2c(f_datatype);
  MPI_Fint itrans = new_itrans();
  xt_ut_trans_type *trans = xt_handle2pointer(trans_set, itrans);
  if (!trans) die("xt_ut_init_transposition_simple: internal error");
  trans->redist = xt_redist_p2p_off_new(template->xmap, NULL, NULL, datatype);

  return itrans;
}

MPI_Fint xt_ut_init_transposition(MPI_Fint itemplate,
                                  MPI_Fint offset_in[], MPI_Fint offset_in_size,
                                  MPI_Fint offset_out[], MPI_Fint offset_out_size,
                                  MPI_Fint f_datatype) {

  xt_ut_template_type *template = xt_handle2pointer(template_set, itemplate);
  if (!template) die("xt_ut_init_transposition: invalid handle");

  MPI_Fint itrans = new_itrans();
  xt_ut_trans_type *trans = xt_handle2pointer(trans_set, itrans);
  if (!trans) die("xt_ut_init_transposition: internal error");
  MPI_Datatype datatype = MPI_Type_f2c(f_datatype);
  trans->redist = xt_redist_p2p_off_new(template->xmap, offset_in, offset_out, datatype);

  return itrans;
}


void xt_create_id_positions(int positions[], int num) {
  for (int i=0; i<num; i++) {
    positions[i] = i;
  }
}

void xt_ut_transpose(void **pt_in, int itrans, int direction, void **pt_out) {
  xt_ut_trans_type *trans = xt_handle2pointer(trans_set, itrans);
  if (!trans) die("xt_ut_transpose1: invalid handle");
  xt_redist_s_exchange(trans->redist, pt_in, 1, pt_out, 1);
}
