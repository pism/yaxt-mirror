/**
 * @file test_redist_common_f2c.c
 *
 * @copyright Copyright  (C)  2026 Jörg Behrens <behrens@dkrz.de>
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "yaxt.h"
#include "tests.h"
#include "test_redist_common.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"

struct xt_redist_f {
  Xt_redist cptr;
};

void
check_redist_f2c(struct xt_redist_f *redist,
                 int sync_mode,
                 int num_redists,
                 const void *const src[num_redists],
                 const size_t dst_num_elems[num_redists],
                 void *const dst[num_redists],
                 const prepare_dst dst_prep[num_redists],
                 const void *dst_prep_info,
                 const void *const ref_dst_data[num_redists],
                 const MPI_Fint dst_data_dt[num_redists],
                 const MPI_Fint ref_dst_data_dt[num_redists],
                 const char *file, int line);

void
check_redist_f2c(struct xt_redist_f *redist,
                 int sync_mode,
                 int num_redists,
                 const void *const src[num_redists],
                 const size_t dst_num_elems[num_redists],
                 void *const dst[num_redists],
                 const prepare_dst dst_prep[num_redists],
                 const void *dst_prep_info,
                 const void *const ref_dst_data[num_redists],
                 const MPI_Fint dst_data_dt[num_redists],
                 const MPI_Fint ref_dst_data_dt[num_redists],
                 const char *file, int line)
{
  MPI_Datatype dst_data_dt_[num_redists], ref_dst_data_dt_[num_redists];
  for (int i = 0; i < num_redists; ++i) {
    dst_data_dt_[i] = MPI_Type_f2c(dst_data_dt[i]);
    ref_dst_data_dt_[i] = MPI_Type_f2c(ref_dst_data_dt[i]);
  }
  check_redist_(redist->cptr, sync_mode, num_redists,
                src, dst_num_elems, dst, dst_prep, dst_prep_info,
                ref_dst_data, dst_data_dt_, ref_dst_data_dt_,
                file, line);
}



/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
