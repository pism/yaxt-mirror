#include "xt/xt_xmap.h"

#include "mpi.h"


MPI_Comm xt_xmap_get_communicator(Xt_xmap xmap) {

  return xmap->vtable->get_communicator(xmap);
}

int xt_xmap_get_num_destinations(Xt_xmap xmap) {

  return xmap->vtable->get_num_destinations(xmap);
}

int xt_xmap_get_num_sources(Xt_xmap xmap) {

  return xmap->vtable->get_num_sources(xmap);
}

void xt_xmap_get_destination_ranks(Xt_xmap xmap, int * ranks) {

  xmap->vtable->get_destination_ranks(xmap, ranks);
}

void xt_xmap_get_source_ranks(Xt_xmap xmap, int * ranks) {

  xmap->vtable->get_source_ranks(xmap, ranks);
}

void xt_xmap_delete(Xt_xmap xmap) {

  xmap->vtable->delete(xmap);
}

Xt_xmap_iter xt_xmap_get_destination_iterator(Xt_xmap xmap) {

  return xmap->vtable->get_destination_iterator(xmap);
}

Xt_xmap_iter xt_xmap_get_source_iterator(Xt_xmap xmap) {

  return xmap->vtable->get_source_iterator(xmap);
}

int xt_xmap_iterator_next(Xt_xmap_iter iter) {

  return iter->vtable->next(iter);
}

int xt_xmap_iterator_get_rank(Xt_xmap_iter iter) {

  return iter->vtable->get_rank(iter);
}

int const * xt_xmap_iterator_get_transfer_pos(Xt_xmap_iter iter) {

  return iter->vtable->get_transfer_pos(iter);
}

Xt_count xt_xmap_iterator_get_num_transfer_pos(Xt_xmap_iter iter) {

  return iter->vtable->get_num_transfer_pos(iter);
}

void xt_xmap_iterator_delete(Xt_xmap_iter iter) {

  iter->vtable->delete(iter);
}

int xt_xmap_get_max_src_pos(Xt_xmap xmap) {
  return xmap->vtable->get_max_src_pos(xmap);
}

int xt_xmap_get_max_dst_pos(Xt_xmap xmap) {
  return xmap->vtable->get_max_dst_pos(xmap);
}
