/**
 * @file xt_mpi.c
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *                                 Joerg Behrens <behrens@dkrz.de>
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 *         Joerg Behrens <behrens@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"

//! COMPACT_DT enables the anlysis of displacements in order to give a
//! more compact description to the datatype generators of MPI. For strong
//! enough MPI implementations this not required. Then you can undefine
//! COMPACT_DT and save some prcessing time within yaxt without losing communication
//! performance.
#define COMPACT_DT

//! Stripe describes a stripe of displacements (like xt_idxstripes for indices)
typedef struct stripe_struct {
  // disp[i] = start + i*stride
  // for 0 <= i < nstrides;
  int start;
  int stride;
  int nstrides;
} Stripe;

static const Stripe zero_stripe = { .start=0, .stride=0, .nstrides=0 };

static MPI_Datatype
xt_mpi_generate_compact_datatype_block(const int *disp, const int *blocklengths,
                                       int count, MPI_Datatype old_type);
static MPI_Datatype
xt_mpi_generate_compact_datatype(int const *disp, int disp_len,
                                 MPI_Datatype old_type);


//taken from http://beige.ucs.indiana.edu/I590/node85.html
void xt_mpi_error(int error_code, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  char error_string[MPI_MAX_ERROR_STRING];
  int length_of_error_string, error_class;

  MPI_Error_class(error_code, &error_class);
  MPI_Error_string(error_class, error_string, &length_of_error_string);
  fprintf(stderr, "%3d: %s\n", rank, error_string);
  MPI_Error_string(error_code, error_string, &length_of_error_string);
  fprintf(stderr, "%3d: %s\n", rank, error_string);
  MPI_Abort(comm, error_code);
}

static MPI_Datatype copy_mpi_datatype(MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_dup(old_type, &datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_simple(int displacement, MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_create_indexed_block(1, 1, &displacement, old_type,
                                                &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_contiguous(int displacement, int blocklength,
                            MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  if (displacement == 0)
    xt_mpi_call(MPI_Type_contiguous(blocklength, old_type, &datatype),
                    comm);
  else
    xt_mpi_call(MPI_Type_create_indexed_block(1, blocklength,
                                                  &displacement, old_type,
                                                  &datatype), comm);

  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;

}

static MPI_Datatype
gen_mpi_datatype_vector(int count, int blocklength, int stride,
                        int offset, MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_vector(count, blocklength, stride, old_type,
                                  &datatype), comm);
  if (offset != 0) {

    MPI_Datatype datatype_;
    int blocklength = 1;
    MPI_Aint old_type_size, old_type_lb;

    xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                        &old_type_size), comm);

    MPI_Aint displacement = offset * old_type_size;

    xt_mpi_call(MPI_Type_create_hindexed(1, &blocklength, &displacement,
                                             datatype, &datatype_),
                    comm);
    xt_mpi_call(MPI_Type_free(&datatype), comm);
    datatype = datatype_;
  }
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_indexed_block(int const * displacements, int blocklength,
                               int count, MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_create_indexed_block(count, blocklength,
                                                (void *)displacements,
                                                old_type, &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_indexed(const int *displacements, const int *blocklengths,
                         int count, MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_indexed(count, (int*)blocklengths, (void*)displacements,
                                   old_type, &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static inline int
check_for_vector_type(const int *displacements, const int *blocklengths,
                      int count) {

  int blocklength = blocklengths[0];

  for (int i = 1; i < count; ++i)
    if (blocklengths[i] != blocklength)
      return 0;

  int stride = displacements[1] - displacements[0];

  for (int i = 1; i + 1 < count; ++i)
    if (displacements[i+1] - displacements[i] != stride)
      return 0;

  return 1;
}

static inline int check_for_indexed_block_type(const int *blocklengths,
                                               int count) {

  int blocklength = blocklengths[0];

  for (int i = 1; i < count; ++i)
    if (blocklengths[i] != blocklength)
      return 0;

  return 1;
}

MPI_Datatype
xt_mpi_generate_datatype_block(const int *displacements,
                               const int *blocklengths,
                               int count, MPI_Datatype old_type,
                               MPI_Comm comm) {
#ifdef COMPACT_DT
  return xt_mpi_generate_compact_datatype_block(displacements, blocklengths, count, old_type);
#endif

  MPI_Datatype datatype;

  if (count == 0)
    datatype = MPI_DATATYPE_NULL;
  else if (count == 1 && blocklengths[0] == 1 && displacements[0] == 0)
    datatype = copy_mpi_datatype(old_type, comm);
  else if (count == 1 && blocklengths[0] == 1)
    datatype = gen_mpi_datatype_simple(displacements[0], old_type, comm);
  else if (count == 1)
    datatype = gen_mpi_datatype_contiguous(displacements[0], blocklengths[0],
                                           old_type, comm);
  else if (check_for_vector_type(displacements, blocklengths, count))
    datatype = gen_mpi_datatype_vector(count, blocklengths[0],
                                       displacements[1] - displacements[0],
                                       displacements[0], old_type, comm);
  else if (check_for_indexed_block_type(blocklengths, count))
    datatype = gen_mpi_datatype_indexed_block(displacements, blocklengths[0],
                                              count, old_type, comm);
  else
    datatype = gen_mpi_datatype_indexed(displacements, blocklengths, count,
                                        old_type, comm);

  return datatype;
}

MPI_Datatype xt_mpi_generate_datatype(int const * displacements, int count,
                                      MPI_Datatype old_type, MPI_Comm comm) {

  if (count <= 0)
    return MPI_DATATYPE_NULL;

#ifdef COMPACT_DT
  return xt_mpi_generate_compact_datatype(displacements, count, old_type);
#endif

  int * blocklengths = xmalloc((size_t)count * sizeof(*blocklengths));
  int new_count = 0;
  {
    int i = 0;
    do {
      int j = 1;
      while (i + j < count && displacements[i] + j == displacements[i + j])
        ++j;
      blocklengths[new_count++] = j;
      i += j;
    } while (i < count);
  }

  int * tmp_displ = NULL;
  const int *displ;

  if (new_count != count) {

    tmp_displ = xmalloc((size_t)new_count * sizeof(*tmp_displ));

    int offset = 0;

    for (int i = 0; i < new_count; ++i) {

      tmp_displ[i] = displacements[offset];
      offset += blocklengths[i];
    }

    displ = tmp_displ;
  } else
    displ = displacements;

  MPI_Datatype datatype;

  datatype = xt_mpi_generate_datatype_block(displ, blocklengths, new_count,
                                            old_type, comm);

  free(blocklengths);

  free(tmp_displ);

  return datatype;
}


static int
scan_stripe(int const *disp, int disp_len, Stripe *v, int vsize) {

  if (disp_len<1) return 0;

  Stripe x = zero_stripe;
  int i = 0;
  int p = 0;
  while (p < disp_len) {

    if (!x.nstrides) {
      x.start   = disp[p];
      x.stride = 1;
      x.nstrides = 1;
      p++; continue;
    }

    if (x.nstrides == 1) {
      x.stride = disp[p] - disp[p-1];
      x.nstrides = 2;
      p++; continue;
    }

    // x.nstrides >= 2:
    if (disp[p] - disp[p-1] == x.stride) {
      x.nstrides++;
      p++; continue;
    }

    if (x.nstrides > 2 || (x.nstrides == 2 && x.stride == 1) ) {
      // we accept small contiguous vectors (nstrides==2, stride==1)
      if (i >= vsize) die("scan_stripe: vsize too small\n");
      v[i]= x;
      i++;
      x = zero_stripe;
      continue;
    }

    if (x.nstrides == 2) {
      // break up trivial vec:
      if (i >= vsize) die("scan_stripe: vsize too small\n");
      v[i] = x;
      v[i].nstrides = 1;
      v[i].stride = 1;
      i++;
      x.start += x.stride;
      x.nstrides = 1;
      x.stride = 1;
      continue;
    }

  }

  // tail cases:
  if (x.nstrides > 2 || (x.nstrides == 2 && x.stride == 1) ) {
    if (i >= vsize) die("scan_stripe: vsize too small\n");
      v[i]= x;
      i++;
  } else if (x.nstrides == 2) {
    if (i+1 >= vsize) die("scan_stripe: vsize too small\n");
      v[i] = x;
      v[i].nstrides = 1;
      v[i].stride = 1;
      i++;
      v[i] = x;
      v[i].start += x.stride;
      v[i].nstrides = 1;
      v[i].stride = 1;
      i++;
  } else if (x.nstrides == 1) {
    if (i >= vsize) die("scan_stripe: vsize too small\n");
    v[i] = x;
    v[i].nstrides = 1;
    v[i].stride = 1;
    i++;
  }

  unsigned vn = (unsigned)i;

  // check:
  p = 0;
  for (unsigned i = 0; i < vn; i++) {
    for (int k = 0; k < v[i].nstrides; k++) {
      int d = v[i].start + k*v[i].stride;
      if (p > disp_len) die("scan_stripe: internal error (1)");
      if (disp[p] != d) die("scan_stripe: internal error (2)");
      p++;
    }
  }
  if (p != disp_len) die("scan_stripe: internal error (3)");

  return (int)vn;
}

static int
match_simple_vec(int *pstart, Stripe const *v, int vlen, MPI_Datatype old_type, int *disp, MPI_Datatype *dt) {
  // we only accept non-trivial matches (nsteps>2) with stride /= 1
  // using only one vector from v
  int p = *pstart;
  if (p >= vlen) return 0;
  int nstrides = v[p].nstrides;
  int stride = v[p].stride;
  if (nstrides < 2 || stride == 1 ) return 0;

  (*pstart)++;

  if (vlen > 1) {
    *disp = v[p].start;
  } else {
    *disp = 0;
  }

  int bl = 1;
  MPI_Datatype dt1;
  xt_mpi_call(MPI_Type_vector(nstrides, bl, stride, old_type, &dt1),
              Xt_default_comm);

  int start = v[p].start - *disp;
  if (!start) {
    *dt = dt1;
    return nstrides;
  }

  // (start != 0) => add offset:
  MPI_Aint old_type_size, old_type_lb;
  xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                  &old_type_size), Xt_default_comm);

  MPI_Aint displacement = start * old_type_size;
  int bl2 = 1;
  MPI_Datatype dt2;
  xt_mpi_call(MPI_Type_create_hindexed(1, &bl2, &displacement, dt1, &dt2),
              Xt_default_comm);

  xt_mpi_call(MPI_Type_free(&dt1), Xt_default_comm);

  *dt = dt2;
  return nstrides;
}

static int
match_block_vec(int *pstart, Stripe const *v, int vlen, MPI_Datatype old_type, int *disp, MPI_Datatype *dt) {
  // using at least 3 vectors
  int p = *pstart;
  if (p+2 >= vlen) return 0;
  if (v[p].stride != 1 || v[p+1].stride != 1 ) return 0;
  int bl = v[p].nstrides;
  if (bl < 1 || v[p+1].nstrides != bl) return 0;

  int vstride = v[p+1].start - v[p].start;

  p += 2;
  while( p < vlen && v[p].stride == 1 && v[p].nstrides == bl &&
         v[p].start - v[p-1].start == vstride ) {
    p++;
  }
  int n = p - *pstart;
  if (n<3) return 0;

  if (n == vlen) {
    *disp = 0;
  } else {
    *disp = v[*pstart].start;
  }

  MPI_Datatype dt1;
  xt_mpi_call(MPI_Type_vector(n, bl, vstride, old_type, &dt1),
              Xt_default_comm);

  int start = v[*pstart].start - *disp;

  *pstart = p;
  if (!start) {
    *dt = dt1;
    return n;
  }

  // (start != 0) => add offset:
  MPI_Aint old_type_size, old_type_lb;
  xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                  &old_type_size), Xt_default_comm);

  MPI_Aint displacement = start * old_type_size;
  int bl2 = 1;
  MPI_Datatype dt2;
  xt_mpi_call(MPI_Type_create_hindexed(1, &bl2, &displacement, dt1, &dt2),
              Xt_default_comm);

  xt_mpi_call(MPI_Type_free(&dt1), Xt_default_comm);

  *dt = dt2;
  return n;
}

static int
match_contiguous(int *pstart, Stripe const *v, int vlen, MPI_Datatype old_type, int *disp, MPI_Datatype *dt) {
  int p = *pstart;
  if (p >= vlen) return 0;
  if (v[p].stride != 1) return 0;
  if (v[p].nstrides < 2) return 0;

  if (vlen > 1) {
    *disp = v[*pstart].start;
  } else {
    *disp = 0;
  }

  int d = v[p].start - *disp;

  if (!d) {
    xt_mpi_call(MPI_Type_contiguous(v[p].nstrides, old_type, dt), Xt_default_comm) ;
  } else {
    xt_mpi_call(MPI_Type_create_indexed_block(1, v[p].nstrides, &d, old_type, dt),
                Xt_default_comm);
  }

  (*pstart)++;
  return v[p].nstrides;
}

static int
match_indexed(int *pstart, Stripe const *v, int vlen, MPI_Datatype old_type, int *disp, MPI_Datatype *dt) {
  // we only accept non-trivial matches
  int p = *pstart;
  if (p >= vlen) return 0;
  if (v[p].stride != 1) return 0;
  if (v[p].nstrides < 2) return 0;

  while( p < vlen && v[p].stride == 1 ) {
    p++;
  }
  int n = p - *pstart;

  if (n < 2) return 0;

  int start;
  if (n == vlen) {
    start = 0;
  } else {
    start = v[*pstart].start;
  }

  *disp = start;
  int bl[n];
  int d[n];
  int hom_bl = 1;
  for (int i = 0; i < n; i++) {
    int iv = *pstart+i;
    d[i] = v[iv].start - start;
    bl[i] = v[iv].nstrides;
    if (bl[i] != bl[0]) hom_bl = 0;
  }

  if (hom_bl) {
    xt_mpi_call(MPI_Type_create_indexed_block(n, bl[0], d, old_type, dt), Xt_default_comm);
  } else {
    xt_mpi_call(MPI_Type_indexed(n, bl, d, old_type, dt),
                Xt_default_comm);
  }

  *pstart = p;

  return n;
}

static int
gen_fallback_type(int *istart, int *iend, Stripe const *v, int vlen, MPI_Datatype old_type, int *offset, MPI_Datatype *dt) {
  if (*istart<0 || *iend < *istart || *iend >= vlen) return 0;

  int ia = *istart;
  int ib = *iend;

  int n = 0;
  for (int i=ia; i<=ib; i++) {
    n += v[i].nstrides;
  }
  if (n<1) return 0;

  int start;
  if (ia == 0 && ib == vlen-1) {
    // generate absolute datatype
    start = 0;
  } else {
    // generate relative datatype that gets embedded by the caller
    start = v[ia].start;
  }

  *offset = start;

  int *d = xmalloc(sizeof (*d) * (size_t)n);
  int p=0;

  for (int i=ia; i<=ib; i++) {
    for (int k=0; k<v[i].nstrides; k++) {
      d[p] = v[i].start + k * v[i].stride - start;
      p++;
    }
  }

  if (n==1 && d[0] == 0) {
    // At the moment we disable the embed aspect that was used to avoid MPI_Type_dup.
    // It turned out this creates more complexity when freeing the intermediate datatypes in the callin function.
    // maybe we just forget the whole idea - then we can also remove the embed logic here.
    //if (embed)
    //  *dt = old_type;
    //else
    xt_mpi_call(MPI_Type_dup(old_type, dt), Xt_default_comm);
  } else {
    xt_mpi_call(MPI_Type_create_indexed_block(n, 1, d, old_type, dt), Xt_default_comm);
  }
  free(d);

  *istart = *iend + 1;

  return n;
}

static MPI_Datatype
parse_stripe(Stripe const *v, int vlen, MPI_Datatype old_type) {

  int set_start = -1;
  int set_end = -2;
  MPI_Datatype *wdt = xmalloc(sizeof(*wdt) * (size_t)vlen);
  int *wdisp = xmalloc(sizeof (*wdisp) * (size_t)vlen);
  int wlen = 0;
  MPI_Datatype match_dt;
  int match_disp;
  int p = 0;
  int m = 0;
  while (p<vlen) {
    if ( match_block_vec(&p, v, vlen, old_type, &match_disp, &match_dt) ||
         match_indexed(&p, v, vlen, old_type, &match_disp, &match_dt) ||
         match_simple_vec(&p, v, vlen, old_type, &match_disp, &match_dt) ||
         match_contiguous(&p, v, vlen, old_type, &match_disp, &match_dt) ) {
      if (set_start <= set_end) {
        gen_fallback_type(&set_start, &set_end, v, vlen, old_type, &wdisp[m], &wdt[m]);
        m++;
      }
      wdt[m] = match_dt;
      wdisp[m] = match_disp;
      m++;
    } else {
      if (set_start>set_end) set_start = p;
      set_end = p;
      p++;
    }

  }
  if (set_start <=  set_end) {
    gen_fallback_type(&set_start, &set_end, v, vlen, old_type, &wdisp[m], &wdt[m]);
    m++;
  }
  wlen = m;
  MPI_Datatype result_dt;
  if (wlen == 1 ) {
    if (wdisp[0] == 0) {
      xt_mpi_call(MPI_Type_dup(wdt[0], &result_dt), Xt_default_comm);
    } else {
      die("parse_stripe: internal error; wlen == 1 && match_disp != 0\n");
    }
  } else {
    MPI_Aint old_type_lb, old_type_extent;
    MPI_Aint wbdisp[wlen];
    int wblocklength[wlen];
    xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                  &old_type_extent), Xt_default_comm);
    for(int i=0; i<wlen; i++) {
      wbdisp[i] = wdisp[i] * old_type_extent;
      wblocklength[i] = 1;
    }
    xt_mpi_call(MPI_Type_create_struct(wlen, wblocklength, wbdisp,
                                       wdt, &result_dt), Xt_default_comm);
  }
  xt_mpi_call(MPI_Type_commit(&result_dt), Xt_default_comm);
  for(m = 0; m<wlen; m++) {
    xt_mpi_call(MPI_Type_free(&wdt[m]), Xt_default_comm);
  }
  free(wdt);
  free(wdisp);
  return result_dt;
}

static MPI_Datatype
xt_mpi_generate_compact_datatype_block(const int *disp, const int *blocklengths,
                                       int count, MPI_Datatype old_type) {

  if (count < 1) return MPI_DATATYPE_NULL;
  Stripe *v = xmalloc(sizeof(*v) * (size_t)count);
  for (int i=0; i<count; i++) {
    v[i].start = disp[i];
    v[i].stride = 1;
    v[i].nstrides = blocklengths[i];
  }
  MPI_Datatype dt = parse_stripe(v, count, old_type);
  free(v);
  return dt;
}

static MPI_Datatype
xt_mpi_generate_compact_datatype(int const *disp, int disp_len,
                                 MPI_Datatype old_type) {

  if (disp_len < 1) return MPI_DATATYPE_NULL;

  Stripe *v = xmalloc(sizeof(*v) * (size_t)disp_len);
  int vlen = scan_stripe(disp, disp_len, v, disp_len);
  MPI_Datatype dt = parse_stripe(v, vlen, old_type);
  free(v);
  return dt;
}

