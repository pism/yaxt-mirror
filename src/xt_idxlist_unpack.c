#include "xt/xt_idxlist.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxlist_collection.h"
#include "xt/xt_idxsection.h"
#include "xt/xt_idxstripes.h"

struct xt_idxlist_unpack_routine_table xt_idxlist_unpack_routines[]
= {{VECTOR, xt_idxvec_unpack},
   {COLLECTION, xt_idxlist_collection_unpack},
   {SECTION, xt_idxsection_unpack},
   {STRIPES, xt_idxstripes_unpack}};

unsigned xt_num_unpack_routines = sizeof(xt_idxlist_unpack_routines) /
                                  sizeof(xt_idxlist_unpack_routines[0]);
