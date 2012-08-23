#ifndef XT_REDIST_H
#define XT_REDIST_H

#include <mpi.h>

#include "yaxt.h"

/**
 * \file xt_redist.h
 * \brief redistribution of data
 *
 * contains declaration the redistribution data structure, which
 * is generate from one or more xt_xmaps
 */

struct xt_redist_vtable {

  void (*delete)(Xt_redist);
  void (*s_exchange)(Xt_redist, void **, unsigned, void **, unsigned);
  void (*s_exchange1)(Xt_redist, void *, void *);
  MPI_Datatype (*get_send_MPI_Datatype)(Xt_redist, int);
  MPI_Datatype (*get_recv_MPI_Datatype)(Xt_redist, int);
};

struct xt_redist {
  struct xt_redist_vtable *vtable;
};

/**
 * destructor
 *
 * @param[in,out] redist redistribution structure
 */
void xt_redist_delete(Xt_redist redist);

/**
 * synchronouse redistribution of data
 *
 * @param[in]     redist         redistribution structure
 * @param[in]     src_data       array containing the base addresses of the input arrays
 * @param[in]     num_src_arrays number of base addresses in src_data
 * @param[in,out] dst_data       array containing the base addresses of the output arrays
 * @param[in]     num_dst_arrays number of base addresses in dst_data
 */
void xt_redist_s_exchange(Xt_redist redist, void **src_data, unsigned num_src_arrays,
                                            void **dst_data, unsigned num_dst_arrays);

/**
 * synchronouse redistribution of data - single array case
 *
 * @param[in]     redist         redistribution structure
 * @param[in]     src_data       input array
 * @param[in,out] dst_data       output array
 */
void xt_redist_s_exchange1(Xt_redist redist, void *src_data, void *dst_data);

/**
 * gets a copy of the MPI_Datatype used for the data of the send operation with
 * the given rank
 *
 * @param[in] redist redistribution structure
 * @param[in] rank   MPI rank
 * \return MPI_Datatype for the data of the send operation with the given rank
 *
 * \remarks returns MPI_DATATYPE_NULL if there is no send operation with the given rank
 */
MPI_Datatype xt_redist_get_send_MPI_Datatype(Xt_redist redist, int rank);

/**
 * gets a copy of the MPI_Datatype used for the data of the recv operation with
 * the given rank
 *
 * @param[in] redist redistribution structure
 * @param[in] rank   MPI rank
 * \return MPI_Datatype for the data of the recv operation with the given rank
 *
 * \remarks returns MPI_DATATYPE_NULL if there is no recv operation with the given rank
 */
MPI_Datatype xt_redist_get_recv_MPI_Datatype(Xt_redist redist, int rank);

#endif // XT_REDIST_H
