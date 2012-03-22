#include "xt_idxset.h"
#include "xt_idxlist.h"
#include "xt_idxset_multi.h"

struct xt_idxset_unpack_routine_table xt_idxset_unpack_routines[] =
                                 {{LIST, xt_idxlist_unpack},
                                 {MULTI, xt_idxset_multi_unpack}};

unsigned xt_num_unpack_routines = sizeof(xt_idxset_unpack_routines) /
                                  sizeof(xt_idxset_unpack_routines[0]);
