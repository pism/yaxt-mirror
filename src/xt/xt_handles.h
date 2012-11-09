/**
 * @file xt_handles.h
 * @brief maps integers to data pointers
 *
 * contains utility routines for handle sets
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *
 * @version 1.0
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
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

#ifndef XT_HANDLES_H
#define XT_HANDLES_H

/** \example test_handles.c
 */

typedef struct xt_handle_set_type **xt_handle_set_type;

/**
 * checks if a given handle is valid within a handle set
 *
 * @param[in] hset       handle set
 * @param[in] handle     handle within handle set
 */
int xt_handle_is_valid(struct xt_handle_set_type *hset, int handle);

/**
 * constructor for handle sets
 *
 * @param[in] cap        start capacity (size of handle space)
 * @return returns an empty handle set of capacity cap
 */
struct xt_handle_set_type *xt_handle_set_new(int cap);

/**
 * destructor
 *
 * @param[in,out] hset handle set
 */
void xt_handle_set_delete(struct xt_handle_set_type *hset);

/**
 * registers user pointer with a new handle
 *
 * @param[in,out] hset     handle set
 * @param[in] p            pointer to user data
 * @return returns integer handle
 */
int xt_handle_new(struct xt_handle_set_type *hset, void *p);

/**
 * unregisters handle
 *
 * @param[in,out] hset     handle set
 * @param[in] handle       handle
 */
void xt_handle_delete(struct xt_handle_set_type *hset, int handle);

/**
 * inquires about user pointer connected with handle
 *
 * @param[in] hset     handle set
 * @param[in] handle       handle
 */
void *xt_handle2pointer(struct xt_handle_set_type *hset, int handle);

#endif
