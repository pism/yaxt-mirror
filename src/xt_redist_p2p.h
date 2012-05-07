#ifndef XT_REDIST_P2P_H
#define XT_REDIST_P2P_H

#include <mpi.h>

#include "xt_redist.h"
#include "xt_xmap.h"

Xt_redist xt_redist_p2p_new(Xt_xmap xmap, MPI_Datatype datatype);

#endif // XT_REDIST_P2P_H
