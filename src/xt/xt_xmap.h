/**
 * @file xt_xmap.h
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://redmine.dkrz.de/doc/yaxt/html/index.html
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef XT_XMAP_H
#define XT_XMAP_H

#include "xt/xt_core.h"

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
  int const * (*get_transfer_pos)(Xt_xmap_iter iter);
  Xt_count (*get_num_transfer_pos)(Xt_xmap_iter iter);
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
  Xt_xmap_iter (*get_destination_iterator)(Xt_xmap);
  Xt_xmap_iter (*get_source_iterator)(Xt_xmap);
  void (*delete)(Xt_xmap);
  int (*get_max_src_pos)(Xt_xmap);
  int (*get_max_dst_pos)(Xt_xmap);
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
 * get the positions of all element of the current item of an exchange
 * map iterator
 *
 * @param[in] iter exchange map iterator
 * @return         exchange list of the current item
 */
int const * xt_xmap_iterator_get_transfer_pos(Xt_xmap_iter iter);

/**
 * get the number of positions returned by
 * \ref xt_xmap_iterator_get_transfer_pos of the current item of an exchange
 * map iterator
 * @param[in] iter exchange map iterator
 * @return         number of positions returned by
 *                 \ref xt_xmap_iterator_get_transfer_pos for the current item
 */
Xt_count xt_xmap_iterator_get_num_transfer_pos(Xt_xmap_iter iter);

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

/**
 * gives a upper limit to the range of positions in the src part of the transfer_pos,
 *
 * @param[in] xmap exchange map
 * @return         upper limit of range of positions in the source list
 */
int xt_xmap_get_max_src_pos(Xt_xmap xmap);

/**
 * gives a upper limit to the range of positions in the dst part of the transfer_pos,
 *
 * @param[in] xmap exchange map
 * @return         upper limit of range of positions in the destination list
 */
int xt_xmap_get_max_dst_pos(Xt_xmap xmap);


#endif // XT_XMAP_H
