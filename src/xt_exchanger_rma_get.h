/**
 * @file xt_exchanger_rma_get.h
 *
 * @copyright Copyright  (C)  2026 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Nils-Arne Dreier <dreier@dkrz.de>
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
#ifndef XT_EXCHANGER_RMA_GET_H
#define XT_EXCHANGER_RMA_GET_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "core/ppm_visibility.h"
#include "xt/xt_config.h"
#include "xt/xt_core.h"
#include "xt_exchanger.h"
#include "xt_redist_internal.h"

/**
 * Constructor for an exchanger using MPI RMA (one-sided) GET operations.
 *
 * Each process allocates an MPI window over its source buffer via
 * MPI_Win_allocate, then receivers pull data from sources using MPI_Get.
 * Synchronization follows the PSCW (Post-Start-Complete-Wait) active-target
 * protocol:
 *   - Sources call MPI_Win_post  (expose window to receivers)
 *   - Receivers call MPI_Win_start (begin access epoch for sources)
 *   - Receivers issue MPI_Get    (pull data from source windows)
 *   - Receivers call MPI_Win_complete (end access epoch)
 *   - Sources call MPI_Win_wait  (wait for all GET operations to finish)
 *
 * For the asynchronous variant (a_exchange), the epoch is started here and
 * a custom Xt_request object is returned; calling xt_request_wait on it
 * calls MPI_Win_complete / MPI_Win_wait to close the epoch.
 *
 * @param[in] nsend      number of send messages (processes that will GET
 *                       from our source window)
 * @param[in] nrecv      number of receive messages (sources we will GET from)
 * @param[in] send_msgs  array with send messages (partners that read us)
 * @param[in] recv_msgs  array with receive messages (partners we read from)
 * @param[in] comm       MPI communicator for the exchange
 * @param[in] tag_offset unused (kept for API compatibility; PSCW needs no tag)
 * @param[in] config     optional customization parameters
 *
 * @remark This exchanger is not suitable for inter-communicators because
 *         MPI_Win_allocate requires an intra-communicator.
 * @remark MPI_Win_allocate copies the source buffer into RMA-allocated memory
 *         at the start of each exchange, which adds a copy overhead compared
 *         to point-to-point exchangers.
 */
PPM_DSO_INTERNAL Xt_exchanger
xt_exchanger_rma_get_new(int nsend, int nrecv,
                         const struct Xt_redist_msg *send_msgs,
                         const struct Xt_redist_msg *recv_msgs,
                         MPI_Comm comm, int tag_offset,
                         Xt_config config);

PPM_DSO_INTERNAL extern const struct xt_exchanger_vtable
xt_exchanger_rma_get_vtable;

#endif /* XT_EXCHANGER_RMA_GET_H */

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
