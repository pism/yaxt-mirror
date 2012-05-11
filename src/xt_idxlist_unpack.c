#include "xt_idxlist.h"
#include "xt_idxlist_t.h"
#include "xt_idxlist_collection.h"

struct xt_idxlist_unpack_routine_table xt_idxlist_unpack_routines[] =
                                 {{LIST, xt_idxlist_t_unpack},
                                  {COLLECTION, xt_idxlist_collection_unpack}};

unsigned xt_num_unpack_routines = sizeof(xt_idxlist_unpack_routines) /
                                  sizeof(xt_idxlist_unpack_routines[0]);
