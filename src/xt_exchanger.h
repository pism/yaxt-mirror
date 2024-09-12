/**
 * @file xt_exchanger.h
 * @brief exchanging of data based on information provided by redist's
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#ifndef XT_EXCHANGER_H
#define XT_EXCHANGER_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>
#include <mpi.h>

#include "xt/xt_config.h"
#include "xt/xt_core.h"
#include "xt/xt_request.h"
#include "xt_redist_internal.h"
#include "core/ppm_visibility.h"


/** \example test_exchanger_parallel.c
 */

typedef struct Xt_exchanger_ *Xt_exchanger;

/**
 * Executes a synchronous data exchange.
 * @param[in]  exchanger exchanger object
 * @param[in]  src_data  source data
 * @param[out] dst_data  destination data
 */
PPM_DSO_INTERNAL void
xt_exchanger_s_exchange(Xt_exchanger exchanger, const void *src_data,
                        void *dst_data);

/**
 * Executes an asynchronous data exchange.
 * @param[in]  exchanger exchanger object
 * @param[in]  src_data  source data
 * @param[out] dst_data  destination data
 * @param[out] request   pointer that will reference request object
 *                       created for duration of exchange
 */
PPM_DSO_INTERNAL void
xt_exchanger_a_exchange(Xt_exchanger exchanger,
                        const void *src_data, void *dst_data,
                        Xt_request *request);

/**
 * Copies an exchange object.
 * @param[in] orig           exchanger object to be copied
 * @param[in] new_comm       communicator to be used by the new exchanger
 * @param[in] new_tag_offset new tag_offset
 */
PPM_DSO_INTERNAL Xt_exchanger
xt_exchanger_copy(Xt_exchanger orig, MPI_Comm new_comm, int new_tag_offset);

struct Xt_msg_param {
  MPI_Aint displacement;
  MPI_Datatype datatype;
  int rank, count;
};

/**
 * Constructor function for an exchanger, obtainable from an
 * \ref Xt_config configuration object */
typedef Xt_exchanger
(*Xt_exchanger_new)(int nsend, int nrecv,
                    const struct Xt_msg_param *send_msgs,
                    const struct Xt_msg_param *recv_msgs,
                    MPI_Comm comm, int tag_offset,
                    Xt_config config);

/**
 * Destructor for an exchanger.
 */
PPM_DSO_INTERNAL void
xt_exchanger_delete(Xt_exchanger);

/**
 * Query number of ranks that this rank sends messages to
 * (direction==xt_direction_send) or receives from
 * (direction==xt_direction_recv) when using this exchanger.
 * @param[in] exchanger exchanger object to query
 * @return number of destination ranks
 */
int
xt_exchanger_get_num_msg(Xt_exchanger exchanger,
                         int direction);


typedef struct Xt_exchanger_iter_ *Xt_exchanger_iter;

/**
 * gets an iterator for all messages of an
 * exchange operation either originating at or terminating at the
 * local rank
 *
 * @param[in] exchanger message exchanger object
 * @return         iterator for all local messages\n
 *                 is NULL if there are no local messages
 */
Xt_exchanger_iter
xt_exchanger_get_iterator(Xt_exchanger exchanger, Xt_config config);

/**
 * discards an exchanger iterator and frees associated memory
 *
 * @param[in,out] iter exchanger iterator to delete
 */
void
xt_exchanger_iter_delete(Xt_exchanger_iter iter);

/**
 * sets an exchanger iterator to the next item
 *
 * @param[in,out] iter exchanger iterator
 * @return             returns 0 if there was no further item
 */
int
xt_exchanger_iter_next(Xt_exchanger_iter iter);

/**
 * gets a copy of the MPI_Datatype used for the currently selected message
 * @param[in] iter      exchanger iterator object
 * @return copy (via MPI_Type_dup) of datatype in use for the current
 * message or MPI_DATATYPE_NULL if xt_exchanger_iter_next has already
 * returned 0 for this iterator
 */
MPI_Datatype
xt_exchanger_iter_get_MPI_Datatype(Xt_exchanger_iter iter);

/**
 * gets remote rank for the currently selected exchanger message
 * @param[in] iter      exchanger iterator object
 * @return rank that is source or destination (depending on direction)
 *         for the current message or -1 if xt_exchanger_iter_next has already
 *         returned 0 for this iterator
 */
int
xt_exchanger_iter_get_rank(Xt_exchanger_iter iter);

/**
 * gets element count for the currently selected exchanger message
 * @param[in] iter      exchanger iterator object
 * @return number of instances of datatype forming the data part of
 * the message or -1 if xt_exchanger_iter_next has already
 * returned 0 for this iterator
 */
int
xt_exchanger_iter_get_count(Xt_exchanger_iter iter);

/**
 * query direction for the currently selected exchanger message
 * @param[in] iter      exchanger iterator object
 * @return 0 if this rank is the source,
 *         1 if it's the destination
 *        -1 if the iterator has already been advanced
 *              past the last message
 */
int
xt_exchanger_iter_get_direction(Xt_exchanger_iter iter);

/**
 * gets byte offset for the currently selected exchanger message
 * @param[in] iter      exchanger iterator object
 * @return byte displacement at which the transfer will be initiated
 *          or -1 if xt_exchanger_iter_next has already
 *          returned 0 for this iterator (but which can also be a
 *          valid displacement in some cases)
 */
MPI_Aint
xt_exchanger_iter_get_displacement(Xt_exchanger_iter iter);



#endif // XT_EXCHANGER_H
/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
