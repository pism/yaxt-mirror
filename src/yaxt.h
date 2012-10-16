#ifndef YAXT_H
#define YAXT_H

#include "mpi.h"

typedef int Xt_idx;
#define Xt_idx_dt (MPI_INT)
typedef int Xt_count;
#define Xt_count_dt (MPI_INT)

typedef struct xt_idxlist *xt_idxlist;
typedef struct xt_idxlist *Xt_idxlist;
typedef struct Xt_xmap *Xt_xmap;
typedef struct xt_redist *Xt_redist;

/**
 * \file yaxt.h
 * \brief main header file
 *
 * contains all external interfaces
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
void xt_idxlist_delete(xt_idxlist idxlist);

/**
 * inquiry function to check if the inverse is defined \n
 * this call is collective among of processes that are
 * associated to the exchange map
 *
 * @param[in] xmap exchange map the is to be checked
 * @return         returns != 0 if the exchange map is bijective
 */
int xt_xmap_is_bijective(Xt_xmap xmap);

#endif // YAXT_H
