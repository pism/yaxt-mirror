#ifndef XT_CORE_H
#define XT_CORE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <mpi.h>

/**
 * Distributed elements are enumerated with numbers of this signed
 * integral type.
 */
typedef XT_IDX Xt_idx;
#define Xt_idx_dt (XT_IDX_MPIDT)
/**
 * Local elements are addressed with numbers of this signed
 * integral type.
 */
typedef XT_COUNT Xt_count;
#define Xt_count_dt (XT_COUNT_MPIDT)

typedef struct Xt_idxlist *Xt_idxlist;
typedef struct Xt_xmap *Xt_xmap;
typedef struct xt_redist *Xt_redist;

/**
 * \file xt_core.h
 * \brief base definitions header file
 *
 * contains types used throughout the library and general initialization
 */

/**
 * initialize library
 * @param[in] default_comm communicator to use for collective aborts
 */
void
xt_initialize(MPI_Comm default_comm);

/**
 * destructor for an index list
 *
 * @param[in] idxlist index list that is to be destroyed
 */
void xt_idxlist_delete(Xt_idxlist idxlist);

/**
 * inquiry function to check if the inverse is defined \n
 * this call is collective among of processes that are
 * associated to the exchange map
 *
 * @param[in] xmap exchange map the is to be checked
 * @return         returns != 0 if the exchange map is bijective
 */
int xt_xmap_is_bijective(Xt_xmap xmap);

#endif
