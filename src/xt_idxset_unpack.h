#ifndef XT_IDX_UNPACK_H
#define XT_IDX_UNPACK_H

#include "xt_idxset.h"
#include "xt_idxlist.h"

struct {
   enum xt_idxset_type type;
   xt_idxset (*xt_idxset_unpack)(void*,int,int*,MPI_Comm);
} xt_idxset_unpack_routines[] = {{LIST, xt_idxlist_unpack}};

#endif // XT_IDX_UNPACK_H
