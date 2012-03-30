#ifndef YAXT_H
#define YAXT_H

#include "mpi.h"

typedef int Xt_idx;
#define Xt_idx_dt (MPI_INT)

typedef struct xt_idxvec *xt_idxvec;
typedef struct xt_idxvec *Xt_idxvec;
typedef struct xt_xmap *xt_xmap;

/**
 * \file yaxt.h
 * \brief main header file
 *
 * contains all external interfaces
 */

/**
 * constructor for an index list
 * 
 * @param[in] idxlist     array containing the global indices of points
 *                        that are available locally
 * @param[in] num_indices number of indices in idxlist
 * @return                returns a xt_idxvec
 *
 * \remark idxlist is only allowed to contain unique indices.
 */
xt_idxvec xt_idxvec_new(int* idxlist, unsigned num_indices);

/**
 * destructor for an index list
 *
 * @param[in] idxvec index vector that is to be destroyed
 */
void xt_idxvec_delete(xt_idxvec idxvec);


/**
 * constructor for an exchange map \n
 * this operation is collective over all processes in comm
 *
 * @param[in] src_idxvec source index vector
 * @param[in] dst_idxvec destination index vector
 * @param[in] comm       MPI communicator that contains all processes
 *                       that part in the exchange (xt_xmap_new will
 *                       make its of copy of comm)
 */
xt_xmap xt_xmap_new(xt_idxvec src_idxvec, xt_idxvec dst_idxvec, MPI_Comm comm);

/**
 * destructor for an exchange map
 *
 * @param[in] xmap exchange map to be destructed
 */
void xt_xmap_delete(xt_xmap xmap);

/**
 * inquiry function to check if the inverse is defined \n
 * this call is collective among of processes that are
 * associated to the exchange map
 *
 * @param[in] xmap exchange map the is to be checked
 * @return         returns != 0 if the exchange map is bijective
 */
int xt_xmap_is_bijective(xt_xmap xmap);

#endif // YAXT_H
