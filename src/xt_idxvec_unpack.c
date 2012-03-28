#include "xt_idxvec.h"
#include "xt_idxlist.h"
#include "xt_idxvec_collection.h"

struct xt_idxvec_unpack_routine_table xt_idxvec_unpack_routines[] =
                                 {{LIST, xt_idxlist_unpack},
                                  {COLLECTION, xt_idxvec_collection_unpack}};

unsigned xt_num_unpack_routines = sizeof(xt_idxvec_unpack_routines) /
                                  sizeof(xt_idxvec_unpack_routines[0]);
