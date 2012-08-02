#ifndef XT_XMAP_H
#define XT_XMAP_H

#include "yaxt.h"

/**
 * \file xt_xmap.h
 * \brief exchange map declaration
 *
 * contains declaration for the exchange map data structure
 */

typedef struct Xt_xmap_iterator *Xt_xmap_iter;

struct Xt_xmap_iterator_vtable {
  int (*next)(Xt_xmap_iter iter);
  int (*get_rank)(Xt_xmap_iter iter);
  Xt_idxlist (*get_xlist)(Xt_xmap_iter iter);
  void (*delete)(Xt_xmap_iter iter);
};

struct Xt_xmap_iterator {
  struct Xt_xmap_iterator_vtable * vtable;
};

struct Xt_xmap_vtable {

  MPI_Comm (*get_communicator)(Xt_xmap);
  int (*get_num_destinations)(Xt_xmap);
  int (*get_num_sources)(Xt_xmap);
  void (*get_destination_ranks)(Xt_xmap, int*);
  void (*get_source_ranks)(Xt_xmap, int*);
  Xt_idxlist (*get_dst_xlist)(Xt_xmap, int);
  Xt_idxlist (*get_src_xlist)(Xt_xmap, int);
  Xt_xmap_iter (*get_destination_iterator)(Xt_xmap);
  Xt_xmap_iter (*get_source_iterator)(Xt_xmap);
  void (*delete)(Xt_xmap);
};

struct Xt_xmap {
  const struct Xt_xmap_vtable * vtable;
};

/**
 * gets the MPI communicator of the exchange map
 * @param[in] xmap exchange map
 * @return         MPI communicator of the exchange map
 */
MPI_Comm xt_xmap_get_communicator(Xt_xmap xmap);

/**
 * gets the number of processes that will require data from
 * the local process in an exchange operation
 *
 * @param[in] xmap exchange map
 * @return         number of prosses that will require data from
 *                 the local process in an exchange operation
 */
int xt_xmap_get_num_destinations(Xt_xmap xmap);

/**
 * gets the number of processes from which the local process will
 * required data in an exchange operation
 *
 * @param[in] xmap exchange map
 * @return number of prosses from which the local process will
 *         required data in an exchange operation
 */
int xt_xmap_get_num_sources(Xt_xmap xmap);

/**
 * gets the ranks of processes that will require data from
 * the local process in an exchange operation
 *
 * @param[in]  xmap  exchange map
 * @param[out] ranks ranks of processes that will require data
 *                   from the local process in an exchange
 *                   operation
 *
 * \remarks the user needs to provide an array for ranks that
 *          is sufficient in size to store all ranks
 */
void xt_xmap_get_destination_ranks(Xt_xmap xmap, int * ranks);

/**
 * gets the ranks of processes from which the local process will
 * required data in an exchange operation
 *
 * @param[in]  xmap  exchange map
 * @param[out] ranks ranks of processes from which the local process
 *                   will required data in an exchange operation
 *
 * \remarks the user needs to provide an array for ranks that
 *          is sufficient in size to store all ranks
 */
void xt_xmap_get_source_ranks(Xt_xmap xmap, int * ranks);

/**
 * gets the destination exchange list for a given rank
 *
 * @param[in] xmap exchange map
 * @param[in] rank rank of process for which the exchange list is to be retrieved
 * @return         destination exchange list for given rank
 *                 NULL in case there is no exchange with the given rank
 */
Xt_idxlist xt_xmap_get_dst_xlist(Xt_xmap xmap, int rank);

/**
 * gets the source exchange list for a given rank
 *
 * @param[in] xmap exchange map
 * @param[in] rank rank of process for which the exchange list is to be retrieved
 * @return         source exchange list for given rank
 *                 NULL in case there is no exchange with the given rank
 */
Xt_idxlist xt_xmap_get_src_xlist(Xt_xmap xmap, int rank);

/**
 * gets an iterator for all processes that will require data from
 * the local process in an exchange operation
 *
 * @param[in] xmap exchange map
 * @return         iterator for all destination processes
 *                 is NULL if there are no destination processes
 */
Xt_xmap_iter xt_xmap_get_destination_iterator(Xt_xmap xmap);

/**
 * gets an iterator for all processes from which the local process will
 * required data in an exchange operation
 *
 * @param[in] xmap exchange map
 * @return         iterator for all source processes
 *                 is NULL if there are no source processes
 */
Xt_xmap_iter xt_xmap_get_source_iterator(Xt_xmap xmap);

/**
 * sets an exchange map iterator to the next item
 *
 * @param[in,out] iter exchange map iterator
 * @return             returns 0 if there was not further item
 */
int xt_xmap_iterator_next(Xt_xmap_iter iter);

/**
 * gets the rank of the current item of an exchange map iterator
 *
 * @param[in] iter exchange map iterator
 * @return         rank of the current item
 */
int xt_xmap_iterator_get_rank(Xt_xmap_iter iter);

/**
 * get the exchange list of the current item of an exchange map iterator
 *
 * @param[in] iter exchange map iterator
 * @return         exchange list of the current item
 */
Xt_idxlist xt_xmap_iterator_get_xlist(Xt_xmap_iter iter);

/**
 * destructor for an exchange map iterator
 *
 * @param[in] iter exchange map iterator
 */
void xt_xmap_iterator_delete(Xt_xmap_iter iter);

/**
 * destructor for an exchange map
 *
 * @param[in] xmap exchange map to be destructed
 */
void xt_xmap_delete(Xt_xmap xmap);

#endif // XT_XMAP_H
