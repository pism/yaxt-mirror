#include "yaxt.h"
#include "xt_redist.h"

struct Xt_redist {

  struct Xt_redist_vtable *vtable;
};

void xt_redist_delete(Xt_redist redist) {

  redist->vtable->delete(redist);
}

void xt_redist_s_exchange(Xt_redist redist, void **src_data, unsigned num_src_arrays,
                                            void **dst_data, unsigned num_dst_arrays) {

  redist->vtable->s_exchange(redist, src_data, num_src_arrays, dst_data, num_dst_arrays);
}
