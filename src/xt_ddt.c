/**
 * @file xt_ddt.c
 *
 * @copyright Copyright  (C)  2022 Jörg Behrens <behrens@dkrz.de>
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
#include "config.h"
#endif

#include <string.h>
#include <mpi.h>

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt/xt_ddt.h"

/* attribute handle for attaching xt_ddt objects to MPI datatypes */
static int xt_ddt_internal_keyval = MPI_KEYVAL_INVALID;

/* attribute attached MPI_COMM_SELF which cleans up
 * xt_mpi_datatype_ddt_internal_keyval when finalize is called */
static int xt_ddt_cleanup_internal_keyval = MPI_KEYVAL_INVALID;

struct xt_ddt_data {
  size_t element_size; // size in byte of base element
  size_t displ_count;  // number of elements in displs
  size_t * displs;     // displacements in byte of elements
                       // (if required this array has to be available
                       //  in GPU memory)
};

struct Xt_ddt_ {

  int ref_count;

  // if exchange data is on GPU, one kernel per data entry is launched
  size_t count; // number of elements in data
  struct xt_ddt_data data[];
};

typedef struct Xt_ddt_ *Xt_ddt;

struct xt_ddt_tree_elem {
  MPI_Aint extent; // not necessarly the true extent
  enum {
    DTYPE,      // element is a basic MPI datatype
    SINGLE_SUB, // element contains a single sub element
    MULTI_SUB,  // element contains a dedicated sub element
                // for each displacement
  } type;       // type of this element
  union {
    MPI_Datatype dtype;                 // basic MPI datatype (type == DTYPE)
    union {
      struct xt_ddt_tree_elem * single; // single sub element (type == SINGLE_SUB)
      struct xt_ddt_tree_elem ** multi; // multiple sub elements (type == MULTI_SUB)
    } sub_elems;
  };
  int displ_count;    // number of displacement
  union {
    MPI_Aint single;  // single displacement (displ_count == 1)
    MPI_Aint * multi; // multiple displacements (displ_count > 1)
  } displs;
};

// returns a pointer to the displacements of an element, if displ_count > 0
// they will be allocated first
static inline MPI_Aint * xt_ddt_tree_elem_alloc_displs(
  struct xt_ddt_tree_elem * elem, int displ_count) {

  elem->displ_count = displ_count;
  if (displ_count <= 1)
    return &(elem->displs.single);
  else
    return
     ((elem->displs.multi =
         xmalloc((size_t)displ_count * sizeof(*(elem->displs.multi)))));
}

static struct xt_ddt_tree_elem * 
  mpi_ddt_2_xt_ddt_tree(MPI_Datatype mpi_ddt, int free_mpi_ddt);

static struct xt_ddt_tree_elem * xt_ddt_tree_elem_named_new(
  MPI_Aint extent, MPI_Datatype mpi_ddt) {

  struct xt_ddt_tree_elem * elem = xmalloc(1 * sizeof(*elem));

  elem->extent = extent;
  elem->type = DTYPE;
  elem->dtype = mpi_ddt;
  elem->displ_count = 1;
  elem->displs.single = 0;
  return elem;
}

static void free_unnamed_mpi_ddt(MPI_Datatype mpi_ddt) {

  // get parameters for decoding of MPI datatype
  int num_ints, num_adds, num_dtypes, combiner;
  xt_mpi_call(
    MPI_Type_get_envelope(
      mpi_ddt, &num_ints, &num_adds, &num_dtypes, &combiner), Xt_default_comm);

  if (combiner != MPI_COMBINER_NAMED)
    xt_mpi_call(MPI_Type_free(&mpi_ddt), Xt_default_comm);
}

// generates a ddt tree element with count contiguous entries of mpi_ddt
static struct xt_ddt_tree_elem * xt_ddt_tree_elem_cont_new(
  int count, MPI_Datatype mpi_ddt) {

  // if the sub element count is zero
  if (count == 0) {
    free_unnamed_mpi_ddt(mpi_ddt);
    return NULL;
  }

  struct xt_ddt_tree_elem * sub_elem = mpi_ddt_2_xt_ddt_tree(mpi_ddt, 1);

  // if the sub element is empty
  if (sub_elem == NULL) return NULL;

  // if there is only a single sub element
  if (count == 1) return sub_elem;

  struct xt_ddt_tree_elem * elem;
  MPI_Aint sub_elem_extent = sub_elem->extent;

  // if the sub element is a basic MPI_Datatype
  if ((sub_elem->type == DTYPE) &&
      (sub_elem->displ_count == 1) &&
      (sub_elem->displs.single == 0)) {

    elem = sub_elem;

  } else {

    elem = xmalloc(1 * sizeof(*elem));
    elem->type = SINGLE_SUB;
    elem->sub_elems.single = sub_elem;
  }

  elem->extent = (MPI_Aint)count * sub_elem_extent;
  MPI_Aint * displs = xt_ddt_tree_elem_alloc_displs(elem, count);
  for (int i = 0; i < count; ++i) displs[i] = (MPI_Aint)i * sub_elem_extent;

  return elem;
}

// generates a ddt tree element based on a hvector style memory layout
static struct xt_ddt_tree_elem * xt_ddt_tree_elem_hvector_new(
  MPI_Aint extent, int count, int blocklength, MPI_Aint stride,
  MPI_Datatype mpi_ddt) {

  // if the sub element count is zero
  if (count == 0) {
    free_unnamed_mpi_ddt(mpi_ddt);
    return NULL;
  }

  struct xt_ddt_tree_elem * sub_elem =
    xt_ddt_tree_elem_cont_new(blocklength, mpi_ddt);

  // if the sub element is empty
  if (sub_elem == NULL) return NULL;

  // if there is only a single sub element
  if (count == 1) return sub_elem;

  struct xt_ddt_tree_elem * elem = xmalloc(1 * sizeof(*elem));

  elem->extent = extent;
  elem->type = SINGLE_SUB;
  elem->sub_elems.single = sub_elem;

  MPI_Aint * displs = xt_ddt_tree_elem_alloc_displs(elem, count);
  for (int i = 0; i < count; ++i) displs[i] = (MPI_Aint)i * stride;

  return elem;
}

// generates a ddt tree element based on a vector style memory layout
static struct xt_ddt_tree_elem * xt_ddt_tree_elem_vector_new(
  MPI_Aint extent, int count, int blocklength, int stride,
  MPI_Datatype mpi_ddt) {

  MPI_Aint lb, base_extent;
  xt_mpi_call(
    MPI_Type_get_extent(mpi_ddt, &lb, &base_extent), Xt_default_comm);

  MPI_Aint hstride = base_extent * (MPI_Aint)stride;

  return
    xt_ddt_tree_elem_hvector_new(
      extent, count, blocklength, hstride, mpi_ddt);
}

// generates a ddt tree element based on a hindexed style memory layout
static struct xt_ddt_tree_elem * xt_ddt_tree_elem_hindexed_new(
  MPI_Aint extent, int count, int * blocklengths,
  MPI_Aint * displacements, MPI_Datatype mpi_ddt) {

  int displ_count = 0;
  for (int i = 0; i < count; ++i) displ_count += blocklengths[i];

  // if there are no entries
  if (displ_count == 0) {
    free_unnamed_mpi_ddt(mpi_ddt);
    return NULL;
  }

  MPI_Aint lb, base_extent;
  MPI_Type_get_extent(mpi_ddt, &lb, &base_extent);
  struct xt_ddt_tree_elem * sub_elem = mpi_ddt_2_xt_ddt_tree(mpi_ddt, 1);

  // if the sub element is empty
  if (sub_elem == NULL) return NULL;

  struct xt_ddt_tree_elem * elem = xmalloc(1 * sizeof(*elem));

  elem->extent = extent;
  elem->type = SINGLE_SUB;
  elem->sub_elems.single = sub_elem;
  MPI_Aint * displs = xt_ddt_tree_elem_alloc_displs(elem, displ_count);

  size_t idx = 0;
  for (int i = 0; i < count; ++i) {
    MPI_Aint block_displacement = displacements[i];
    for (int j = 0; j < blocklengths[i]; ++j, ++idx) {
      displs[idx] = block_displacement + (MPI_Aint)(j * base_extent);
    }
  }

  return elem;
}

// generates a ddt tree element based on an indexed style memory layout
static struct xt_ddt_tree_elem * xt_ddt_tree_elem_indexed_new(
  MPI_Aint extent, int count, int * blocklengths,
  int * displacements, MPI_Datatype mpi_ddt) {

  MPI_Aint lb, base_extent;
  MPI_Type_get_extent(mpi_ddt, &lb, &base_extent);

  MPI_Aint * hdisplacements = xmalloc((size_t)count * sizeof(*hdisplacements));
  for (int i = 0; i < count; ++i)
    hdisplacements[i] = (MPI_Aint)displacements[i] * base_extent;

  struct xt_ddt_tree_elem * elem =
    xt_ddt_tree_elem_hindexed_new(
      extent, count, blocklengths, hdisplacements, mpi_ddt);

  free(hdisplacements);

  return elem;
}

// generates a ddt tree element based on an hindexed block style memory layout
static struct xt_ddt_tree_elem * xt_ddt_tree_elem_hindexed_block_new(
  MPI_Aint extent, int count, int blocklength, MPI_Aint * displacements,
  MPI_Datatype mpi_ddt) {

  // if there are no blocks or if the blocks are empty 
  if ((count == 0) || (blocklength == 0)) {
    free_unnamed_mpi_ddt(mpi_ddt);
    return NULL;
  }

  struct xt_ddt_tree_elem * sub_elem =
    xt_ddt_tree_elem_cont_new(blocklength, mpi_ddt);

  // if the sub element is empty
  if (sub_elem == NULL) return NULL;

  struct xt_ddt_tree_elem * elem = xmalloc(1 * sizeof(*elem));

  elem->extent = extent;
  elem->type = SINGLE_SUB;
  elem->sub_elems.single = sub_elem;
  MPI_Aint * displs = xt_ddt_tree_elem_alloc_displs(elem, count);

  for (int i = 0; i < count; ++i) displs[i] = displacements[i];

  return elem;
}

// generates a ddt tree element based on an indexed block style memory layout
static struct xt_ddt_tree_elem * xt_ddt_tree_elem_indexed_block_new(
  MPI_Aint extent, int count, int blocklength, int * displacements,
  MPI_Datatype mpi_ddt) {

  MPI_Aint lb, base_extent;
  xt_mpi_call(
    MPI_Type_get_extent(mpi_ddt, &lb, &base_extent), Xt_default_comm);

  MPI_Aint * hdisplacements = xmalloc((size_t)count * sizeof(*hdisplacements));
  for (int i = 0; i < count; ++i)
    hdisplacements[i] = (MPI_Aint)displacements[i] * base_extent;

  struct xt_ddt_tree_elem * elem =
    xt_ddt_tree_elem_hindexed_block_new(
      extent, count, blocklength, hdisplacements, mpi_ddt);

  free(hdisplacements);

  return elem;
}

// generates a ddt tree element based on a struct style memory layout
static struct xt_ddt_tree_elem * xt_ddt_tree_elem_struct_new(
  MPI_Aint extent, int count, int * blocklengths, MPI_Aint * displacements,
  MPI_Datatype * mpi_ddts) {

  // remove empty blocks
  int new_count = 0;
  for (int i = 0; i < count; ++i) {

    MPI_Aint true_lb, true_extent;
    xt_mpi_call(
      MPI_Type_get_true_extent(mpi_ddts[i], &true_lb, &true_extent),
      Xt_default_comm);

    if ((blocklengths[i] != 0) && (true_extent != 0)) {

      if (new_count != count) {

        blocklengths[new_count] = blocklengths[i];
        displacements[new_count] = displacements[i];
        mpi_ddts[new_count] = mpi_ddts[i];
      }

      ++new_count;
    } else {
      free_unnamed_mpi_ddt(mpi_ddts[i]);
    }
  }
  if (new_count != count) count = new_count;

  // if there are no blocks
  if (count == 0) return NULL;

  // if there is only a single block
  if (count == 1)
    return
      xt_ddt_tree_elem_hindexed_block_new(
        extent, 1, blocklengths[0], displacements, *mpi_ddts);

  struct xt_ddt_tree_elem ** sub_elems = 
    xmalloc((size_t)count * sizeof(*sub_elems));

  struct xt_ddt_tree_elem * elem = xmalloc(1 * sizeof(*elem));

  elem->extent = extent;
  elem->type = MULTI_SUB;
  elem->sub_elems.multi = sub_elems;
  MPI_Aint * displs = xt_ddt_tree_elem_alloc_displs(elem, count);

  for (int i = 0; i < count; ++i) {
    sub_elems[i] = xt_ddt_tree_elem_cont_new(blocklengths[i], mpi_ddts[i]);
    displs[i] = displacements[i];
  }

  return elem;
}

// generate the displacements for a subarray
// (this routine is recursive)
static void xt_ddt_tree_elem_subarray_get_displs(
  MPI_Aint * displs, MPI_Aint dim_displ,
  size_t ndim, int * sizes, int * sub_sizes, int * starts, int order,
  MPI_Aint base_extent) {

  // if the base element is reached
  if (ndim == 0) {
    *displs = dim_displ;
    return;
  }

  // number of entries in the current dimension and index of the first entry
  int dim_sub_size, dim_start;

  if (order == MPI_ORDER_FORTRAN) {
    dim_sub_size = sub_sizes[ndim - 1];
    dim_start = starts[ndim - 1];
  } else {
    dim_sub_size = sub_sizes[0];
    dim_start = starts[0];
    sizes += 1;
    sub_sizes += 1;
    starts += 1;
  }

  size_t sub_array_size = 1, sub_array_extent = base_extent;
  for (size_t i = 0; i < ndim-1; ++i) {
    sub_array_size *= (size_t)(sub_sizes[i]);
    sub_array_extent *= (MPI_Aint)(sizes[i]);
  }

  // set current offset to the first entry in the current dimension
  dim_displ += (MPI_Aint)dim_start * sub_array_extent;

  // for the number of entries in the current dimension
  for (int i = 0; i < dim_sub_size; ++i) {

    // recursive call for sub dimensions
    xt_ddt_tree_elem_subarray_get_displs(
      displs, dim_displ, ndim - 1,
      sizes, sub_sizes, starts, order, base_extent);

    // set to next entry in current dimension
    displs += sub_array_size;
    dim_displ += sub_array_extent;
  }
}

static struct xt_ddt_tree_elem * xt_ddt_tree_elem_subarray_new(
  MPI_Aint extent, int ndim, int * sizes, int * sub_sizes,
  int * starts, int order, MPI_Datatype mpi_ddt) {

  struct xt_ddt_tree_elem * sub_elem = mpi_ddt_2_xt_ddt_tree(mpi_ddt, 1);

  // if the sub element is empty
  if (sub_elem == NULL) return NULL;

  MPI_Aint base_extent = sub_elem->extent;

  struct xt_ddt_tree_elem * elem = xmalloc(1 * sizeof(*elem));
  elem->extent = extent;
  elem->type = SINGLE_SUB;
  elem->sub_elems.single = sub_elem;

  size_t displ_count = 1;
  for (int i = 0; i < ndim; ++i)
    displ_count *= (size_t)(sub_sizes[i]);

  MPI_Aint * displs = xt_ddt_tree_elem_alloc_displs(elem, displ_count);
  xt_ddt_tree_elem_subarray_get_displs(
    displs, 0, (size_t)ndim, sizes, sub_sizes, starts, order, base_extent);

  return elem;
}

static struct xt_ddt_tree_elem * xt_ddt_tree_elem_resized_new(
  MPI_Aint extent, MPI_Datatype mpi_ddt) {

  struct xt_ddt_tree_elem * elem = mpi_ddt_2_xt_ddt_tree(mpi_ddt, 1);
  elem->extent = extent;

  return elem;
}

// generates a ddt tree element from a MPI datatype
static struct xt_ddt_tree_elem * 
  mpi_ddt_2_xt_ddt_tree(MPI_Datatype mpi_ddt, int free_mpi_ddt) {

  // get parameters for decoding of MPI datatype
  int num_ints, num_adds, num_dtypes, combiner;
  xt_mpi_call(
    MPI_Type_get_envelope(
      mpi_ddt, &num_ints, &num_adds, &num_dtypes, &combiner), Xt_default_comm);
  MPI_Datatype array_of_dtypes_single[1];
  MPI_Datatype * array_of_dtypes = array_of_dtypes_single;

  // get the extent of the MPI datatype
  MPI_Aint lb, extent;
  xt_mpi_call(MPI_Type_get_extent(mpi_ddt, &lb, &extent), Xt_default_comm);

  struct xt_ddt_tree_elem * elem;

  // decode MPI datatype
  switch (combiner) {

    default: {
      Xt_abort(
        Xt_default_comm,
        "ERROR(mpi_ddt_2_xt_ddt_tree): unsupported combiner type",
        __FILE__, __LINE__);
      return NULL;
    }

    case MPI_COMBINER_NAMED:

      elem = xt_ddt_tree_elem_named_new(extent, mpi_ddt);
      break;

    case MPI_COMBINER_DUP: {

      int * array_of_ints = NULL;
      MPI_Aint * array_of_adds = NULL;
      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, array_of_dtypes), Xt_default_comm);

      elem = mpi_ddt_2_xt_ddt_tree(array_of_dtypes[0], 1);
      break;
    }

    case MPI_COMBINER_CONTIGUOUS: {

      int array_of_ints[1];
      MPI_Aint * array_of_adds = NULL;

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, array_of_dtypes), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_cont_new(
          array_of_ints[0], array_of_dtypes[0]);
      break;
    }

    case MPI_COMBINER_VECTOR: {

      int array_of_ints[3];
      MPI_Aint * array_of_adds = NULL;

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, array_of_dtypes), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_vector_new(
          extent, array_of_ints[0], array_of_ints[1], array_of_ints[2],
          array_of_dtypes[0]);
      break;
    }

    case MPI_COMBINER_HVECTOR: {

      int array_of_ints[2];
      MPI_Aint array_of_adds[1];

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, array_of_dtypes), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_hvector_new(
          extent, array_of_ints[0], array_of_ints[1], array_of_adds[0],
          array_of_dtypes[0]);
      break;
    }
    case MPI_COMBINER_INDEXED:{

      int * array_of_ints = xmalloc((size_t)num_ints * sizeof(*array_of_ints));
      MPI_Aint * array_of_adds = NULL;

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, array_of_dtypes), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_indexed_new(
          extent, array_of_ints[0],
          array_of_ints + 1, array_of_ints + 1 + 1 * array_of_ints[0],
          array_of_dtypes[0]);

      free(array_of_ints);

      break;
    }
    case MPI_COMBINER_HINDEXED:{

      int * array_of_ints =
        xmalloc((size_t)num_ints * sizeof(*array_of_ints));
      MPI_Aint * array_of_adds =
        xmalloc((size_t)num_adds * sizeof(*array_of_adds));

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, array_of_dtypes), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_hindexed_new(
          extent, array_of_ints[0], array_of_ints + 1,
          array_of_adds, array_of_dtypes[0]);

      free(array_of_adds);
      free(array_of_ints);

      break;
    }
    case MPI_COMBINER_INDEXED_BLOCK:{

      int * array_of_ints = xmalloc((size_t)num_ints * sizeof(*array_of_ints));
      MPI_Aint * array_of_adds = NULL;

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, array_of_dtypes), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_indexed_block_new(
          extent, array_of_ints[0], array_of_ints[1], array_of_ints + 2,
          array_of_dtypes[0]);

      free(array_of_ints);

      break;
    }
    case MPI_COMBINER_HINDEXED_BLOCK:{

      int array_of_ints[2];
      MPI_Aint * array_of_adds =
        xmalloc((size_t)num_adds * sizeof(*array_of_adds));

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, array_of_dtypes), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_hindexed_block_new(
          extent, array_of_ints[0], array_of_ints[1], array_of_adds,
          array_of_dtypes[0]);

      free(array_of_adds);

      break;
    }
    case MPI_COMBINER_STRUCT:{

      int * array_of_ints =
        xmalloc((size_t)num_ints * sizeof(*array_of_ints));
      MPI_Aint * array_of_adds =
        xmalloc((size_t)num_adds * sizeof(*array_of_adds));
      array_of_dtypes =
        xmalloc((size_t)num_dtypes * sizeof(*array_of_dtypes));

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, array_of_dtypes), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_struct_new(
          extent, array_of_ints[0], array_of_ints + 1,
          array_of_adds, array_of_dtypes);

      free(array_of_dtypes);
      free(array_of_adds);
      free(array_of_ints);

      break;
    }
    case MPI_COMBINER_SUBARRAY:{

      int * array_of_ints =
        xmalloc((size_t)num_ints * sizeof(*array_of_ints));
      MPI_Aint * array_of_adds = NULL;

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, array_of_dtypes), Xt_default_comm);

      elem =
        xt_ddt_tree_elem_subarray_new(
          extent, array_of_ints[0],
          array_of_ints + 1 + 0 * array_of_ints[0],
          array_of_ints + 1 + 1 * array_of_ints[0],
          array_of_ints + 1 + 2 * array_of_ints[0],
          array_of_ints[1 + 3 * array_of_ints[0]],
          array_of_dtypes[0]);

      free(array_of_ints);

      break;
    }
    case MPI_COMBINER_RESIZED:{

      int * array_of_ints = NULL;
      MPI_Aint array_of_adds[2];

      xt_mpi_call(
        MPI_Type_get_contents(
          mpi_ddt, num_ints, num_adds, num_dtypes,
          array_of_ints, array_of_adds, array_of_dtypes), Xt_default_comm);

      elem = xt_ddt_tree_elem_resized_new(extent, array_of_dtypes[0]);
      break;
    }
  };

  if ((combiner != MPI_COMBINER_NAMED) && (free_mpi_ddt))
    xt_mpi_call(MPI_Type_free(&mpi_ddt), Xt_default_comm);

  return elem;
}

static void xt_ddt_tree_elem_delete(struct xt_ddt_tree_elem * elem) {

  if (elem == NULL) return;

  if (elem->displ_count > 1) free(elem->displs.multi);

  switch (elem->type) {

    default:
      break;
    case (SINGLE_SUB):
      xt_ddt_tree_elem_delete(elem->sub_elems.single);
      break;
    case (MULTI_SUB): {
      struct xt_ddt_tree_elem ** sub_elems = elem->sub_elems.multi;
      for (int i = 0; i < elem->displ_count; ++i)
        xt_ddt_tree_elem_delete(sub_elems[i]);
      free(sub_elems);
      break;
    }
  };
  free(elem);
}

static void xt_ddt_tree_delete(struct xt_ddt_tree_elem * tree) {

  if (tree != NULL) xt_ddt_tree_elem_delete(tree);
}

static inline size_t mpi_ddt_to_data_idx(MPI_Datatype mpi_ddt) {

  MPI_Aint true_lb, true_extent;
  xt_mpi_call(
    MPI_Type_get_true_extent(mpi_ddt, &true_lb, &true_extent),
    Xt_default_comm);

  switch(true_extent) {
    case 1:  return 5;
    case 2:  return 4;
    case 4:  return 3;
    case 8:  return 2;
    case 16: return 1;
    case 32: return 0;
    default:
      Xt_abort(
        Xt_default_comm,
        "ERROR(mpi_ddt_to_data_idx): datatype size", __FILE__, __LINE__);
  };
  return SIZE_MAX;
}

static void xt_ddt_tree_elem_get_data_sizes(
  struct xt_ddt_tree_elem * elem, MPI_Datatype * prev_dtype,
  size_t * prev_data_idx, struct xt_ddt_data * data) {

  if (elem == NULL) return;

  switch (elem->type) {

    case(DTYPE): {

      if (elem->dtype != *prev_dtype) {

        *prev_dtype = elem->dtype;
        *prev_data_idx = mpi_ddt_to_data_idx(elem->dtype);
      }

      data[*prev_data_idx].displ_count += (size_t)(elem->displ_count);
      break;
    } 
    case (SINGLE_SUB): {

      struct xt_ddt_tree_elem * sub_elem = elem->sub_elems.single;

      for (int i = 0; i < elem->displ_count; ++i)
        xt_ddt_tree_elem_get_data_sizes(
          sub_elem, prev_dtype, prev_data_idx, data);
      break;
    }
    case (MULTI_SUB): {

      struct xt_ddt_tree_elem ** sub_elems = elem->sub_elems.multi;

      for (int i = 0; i < elem->displ_count; ++i)
        xt_ddt_tree_elem_get_data_sizes(
          sub_elems[i], prev_dtype, prev_data_idx, data);
      break;
    }
  };
}

static void xt_ddt_tree_get_data_sizes(
  struct xt_ddt_tree_elem * tree, struct xt_ddt_data * data) {

  MPI_Datatype prev_dtype = MPI_DATATYPE_NULL;
  size_t prev_data_idx = SIZE_MAX;

  if (tree != NULL)
    xt_ddt_tree_elem_get_data_sizes(
      tree, &prev_dtype, &prev_data_idx, data);
}

static void xt_ddt_tree_elem_to_data(
  struct xt_ddt_tree_elem * elem, struct xt_ddt_data * data, size_t displ,
  MPI_Datatype * prev_dtype, size_t * prev_data_idx) {

  if (elem == NULL) return;

  size_t elem_displ_count = (size_t)(elem->displ_count);
  MPI_Aint * elem_displs =
    (elem_displ_count <= 1)?(&(elem->displs.single)):(elem->displs.multi);

  switch (elem->type) {

    case(DTYPE): {

      if (elem->dtype != *prev_dtype) {

        *prev_dtype = elem->dtype;
        *prev_data_idx = mpi_ddt_to_data_idx(elem->dtype);
      }

      size_t displ_count = data[*prev_data_idx].displ_count;
      size_t * displs = data[*prev_data_idx].displs;
      for (size_t i = 0; i < elem_displ_count; ++i, ++displ_count)
        displs[displ_count] = elem_displs[i] + displ;
      data[*prev_data_idx].displ_count = displ_count;

      break;
    } 
    case (SINGLE_SUB): {

      struct xt_ddt_tree_elem * sub_elem = elem->sub_elems.single;

      for (size_t i = 0; i < elem_displ_count; ++i)
        xt_ddt_tree_elem_to_data(
          sub_elem, data, displ + elem_displs[i],
          prev_dtype, prev_data_idx);
      break;
    }
    case (MULTI_SUB): {

      struct xt_ddt_tree_elem ** sub_elems = elem->sub_elems.multi;

      for (size_t i = 0; i < elem_displ_count; ++i)
        xt_ddt_tree_elem_to_data(
          sub_elems[i], data, displ + elem_displs[i],
          prev_dtype, prev_data_idx);
      break;
    }
  };
}

static void xt_ddt_tree_to_data(
  struct xt_ddt_tree_elem * tree, struct xt_ddt_data * data) {

  size_t displ = 0;
  MPI_Datatype prev_dtype = MPI_DATATYPE_NULL;
  size_t prev_data_idx = SIZE_MAX;

  if (tree != NULL)
    xt_ddt_tree_elem_to_data(
      tree, data, displ, &prev_dtype, &prev_data_idx);
}

static Xt_ddt xt_ddt_new(MPI_Datatype mpi_ddt) {

  // parse the MPI datatype
  struct xt_ddt_tree_elem * tree = mpi_ddt_2_xt_ddt_tree(mpi_ddt, 0);

  struct xt_ddt_data data[6];
  for (int i = 0; i < 6; ++i) {
    data[i].element_size = 1 << (5 - i);
    data[i].displ_count = 0;
    data[i].displs = NULL;
  }

  // determine the number of displacements for each elemental data size
  xt_ddt_tree_get_data_sizes(tree, data);

  // compute the total number of elemental data
  size_t total_displs_size = 0;
  size_t data_count = 0;
  for (size_t i = 0; i < 6; ++i) {
    total_displs_size += data[i].displ_count;
    if (data[i].displ_count > 0) ++data_count;
  }

  Xt_ddt ddt = xmalloc(1 * sizeof(*ddt) + data_count * sizeof(ddt->data[0]));
  ddt->ref_count = 1;
  ddt->count = data_count;

  if (total_displs_size > 0) {

    // allocate memory for all displacements
    size_t * displs = xmalloc(total_displs_size * sizeof(*displs));

    for (size_t i = 0, offset = 0; i < 6; ++i) {
      data[i].displs = displs + offset;
      offset += data[i].displ_count;
      data[i].displ_count = 0;
    }

    // determine all displacements
    xt_ddt_tree_to_data(tree, data);

    for (int i = 0, j = 0; i < 6; ++i)
      if (data[i].displ_count > 0)
        ddt->data[j++] = data[i];
  }

  // free internal representation of MPI datatype
  xt_ddt_tree_delete(tree);

  return ddt;
}

void xt_ddt_inc_ref_count(Xt_ddt ddt) {

  ddt->ref_count++;
}

void xt_ddt_delete(Xt_ddt ddt) {

  if (ddt == NULL) return;

  ddt->ref_count--;

  if (ddt->ref_count) return;

  // only the displacements of the first entry needs to be freed, all other
  // are part of the same allocation
  if (ddt->count > 0) free(ddt->data[0].displs);
  free(ddt);
}

static int xt_ddt_internal_keyval_copy(
  MPI_Datatype XT_UNUSED(dtype), int XT_UNUSED(dtype_keyval),
  void *XT_UNUSED(extra_state), void *attribute_val_in,
  void *attribute_val_out, int *flag) {

  xt_ddt_inc_ref_count((Xt_ddt)attribute_val_in);
  *(Xt_ddt*)attribute_val_out = (Xt_ddt)attribute_val_in;
  *flag = 1;
  return MPI_SUCCESS;
}

static int xt_ddt_internal_keyval_delete(
  MPI_Datatype XT_UNUSED(dtype), int XT_UNUSED(dtype_keyval),
  void *attribute_val, void *XT_UNUSED(extra_state)) {

  xt_ddt_delete((Xt_ddt)attribute_val);
  return MPI_SUCCESS;
}

static int xt_ddt_cleanup_internal_keyval_delete(
  MPI_Comm comm, int XT_UNUSED(dtype_keyval),
  void *XT_UNUSED(attribute_val), void *XT_UNUSED(extra_state))
{

  if (xt_ddt_internal_keyval != MPI_KEYVAL_INVALID) {
    xt_mpi_call(
      MPI_Type_free_keyval(&xt_ddt_internal_keyval), comm);
    xt_ddt_internal_keyval = MPI_KEYVAL_INVALID;
  }

  xt_mpi_call(
    MPI_Comm_free_keyval(&xt_ddt_cleanup_internal_keyval), comm);
  xt_ddt_cleanup_internal_keyval = MPI_KEYVAL_INVALID;
  
  return MPI_SUCCESS;
}

static Xt_ddt get_xt_ddt_from_mpi_ddt(MPI_Datatype mpi_ddt) {

  // if xt_ddt_internal_keyval has not yet been created
  if (xt_ddt_internal_keyval == MPI_KEYVAL_INVALID) {

    // create keyval that will store xt_ddt's in MPI_Datatype's
    xt_mpi_call(
      MPI_Type_create_keyval(
        // MPI_TYPE_NULL_COPY_FN,
        xt_ddt_internal_keyval_copy,
        xt_ddt_internal_keyval_delete,
        &xt_ddt_internal_keyval, NULL),
      Xt_default_comm);

    // register a callback in MPI_Finalize (via an attribute to MPI_COMM_SELF)
    // to clean up xt_ddt_internal_keyval
    xt_mpi_call(
      MPI_Comm_create_keyval(
        MPI_COMM_NULL_COPY_FN,
        xt_ddt_cleanup_internal_keyval_delete,
        &xt_ddt_cleanup_internal_keyval, NULL),
      Xt_default_comm);
    xt_mpi_call(
      MPI_Comm_set_attr(
        MPI_COMM_SELF, xt_ddt_cleanup_internal_keyval, NULL),
      Xt_default_comm);
  }

  // get xt_ddt from MPI datatype, if available
  void * attr;
  int flag;
  xt_mpi_call(
    MPI_Type_get_attr(
      mpi_ddt, xt_ddt_internal_keyval, &attr, &flag), Xt_default_comm);

  // if the MPI datatype had already a xt_ddt attached to itself
  if (flag) return (Xt_ddt)attr;

  // generate a xt_ddt from the MPI datatype
  Xt_ddt ddt = xt_ddt_new(mpi_ddt);

  // attach this xt_ddt to the MPI datatype for later use
  xt_mpi_call(
    MPI_Type_set_attr(
      mpi_ddt, xt_ddt_internal_keyval, (void*)ddt), Xt_default_comm);

  return ddt;
}

size_t xt_ddt_get_pack_size(MPI_Datatype mpi_ddt) {

  Xt_ddt ddt = get_xt_ddt_from_mpi_ddt(mpi_ddt);

  if (ddt == NULL) return 0;

  size_t pack_size = 0;

  for (size_t i = 0; i < ddt->count; ++i)
    pack_size += ddt->data[i].element_size * ddt->data[i].displ_count;

  return pack_size;
}

void xt_ddt_pack(MPI_Datatype mpi_ddt, void const * src, void * dst) {

  Xt_ddt ddt = get_xt_ddt_from_mpi_ddt(mpi_ddt);

  size_t dst_offset = 0;

  // for all sections with the same elemental datatype extent
  for (size_t i = 0; i < ddt->count; ++i) {

    size_t element_size = ddt->data[i].element_size;
    size_t displ_count = ddt->data[i].displ_count;
    size_t * displs = ddt->data[i].displs;

    // for all displacements of the current section
    for (size_t j = 0; j < displ_count; ++j)
      memcpy(
        (unsigned char *)dst + dst_offset + j * element_size,
        (unsigned char *)src + displs[j], element_size);

    dst_offset += displ_count * element_size;
  }
}

void xt_ddt_unpack(MPI_Datatype mpi_ddt, void const * src, void * dst) {

  Xt_ddt ddt = get_xt_ddt_from_mpi_ddt(mpi_ddt);

  size_t src_offset = 0;

  // for all sections with the same elemental datatype extent
  for (size_t i = 0; i < ddt->count; ++i) {

    size_t element_size = ddt->data[i].element_size;
    size_t displ_count = ddt->data[i].displ_count;
    size_t * displs = ddt->data[i].displs;

    // for all displacements of the current section
    for (size_t j = 0; j < displ_count; ++j)
      memcpy(
        ((unsigned char *)dst) + displs[j],
        ((unsigned char *)src) + src_offset + j * element_size, element_size);

    src_offset += displ_count * element_size;
  }
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
