/**
 * @file xt_config.h
 * @brief opaque configuration object for settings where the default
 * needs to be overridden
 *
 * @copyright Copyright  (C)  2020 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 *
 * URL: https://doc.redmine.dkrz.de/yaxt/html/
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

#ifndef XT_CONFIG_H
#define XT_CONFIG_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct Xt_config_ *Xt_config;

/**
 * constructor for configuration object
 *
 * @return returns a configuration object where every setting is set
 * to the corresponding default.
 */
Xt_config xt_config_new(void);

/**
 * destructor of configuration objects
 *
 * @param[in,out] config configuration object to destroy
 */
void xt_config_delete(Xt_config config);

enum Xt_exchangers {
  xt_exchanger_irecv_send,
  xt_exchanger_irecv_isend,
  xt_exchanger_irecv_isend_packed,
  xt_exchanger_mix_isend_irecv,
  xt_exchanger_neigh_alltoall,
};

/**
 * set exchanger to use when the \a config object is passed to constructors
 * @param[in,out] config configuration object to modify
 * @param method an entry from enum Xt_exchangers to signify the
 * desired exchanger for data transfers
 */
void xt_config_set_exchange_method(Xt_config config, int method);

/**
 * get exchanger used when the \a config object is passed to constructors
 * @param[in] config configuration object to modify
 * @return an entry from \a Xt_exchangers representing the method of
 * data transfer used
 */
int xt_config_get_exchange_method(Xt_config config);

#endif

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
