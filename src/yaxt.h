#ifndef YAXT_H
#define YAXT_H

#include "mpi.h"

typedef struct xt_idxset* xt_idxset;
typedef struct xt_xmap* xt_xmap;

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
 * @return                returns a xt_idxset
 */
xt_idxset xt_idxset_new(int* idxlist, unsigned num_indices);

/**
 * destructor for an index list
 *
 * @param[in] idxset index set that is to be destructed
 */
void xt_idxset_delete(xt_idxset idxset);


/**
 * constructor for an exchange map \n
 * this operation is collective over all processes in comm
 *
 * @param[in] src_idxset source index set
 * @param[in] dst_idxset destination index set
 * @param[in] comm       MPI communicator that contains all processes
 *                       that part in the exchange (xt_xmap_new will
 *                       make its of copy of comm)
 */
xt_xmap xt_xmap_new(xt_idxset src_idxset, xt_idxset dst_idxset, MPI_Comm comm);

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
