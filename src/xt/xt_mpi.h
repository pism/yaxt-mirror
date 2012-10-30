#ifndef XT_MPI_H
#define XT_MPI_H

#include <mpi.h>

/**
 * \file xt_mpi.h
 * \brief utility routines for MPI
 *
 * contains utility routines for handling MPI
 */

/**
 * checks the return value of an MPI call
 *
 * @param[in] error_code return code of an MPI call
 * @param[in] comm       communicator which was used for the respective MPI call
 */
void mpi_err_handler(int error_code, MPI_Comm comm);

/**
 * generates an MPI datatype
 *
 * @param[in] displacements array of displacements
 * @param[in] count         number of elements
 * @param[in] old_type      base MPI datatype of all elements
 * @param[in] comm          MPI communicator
 * @return                  MPI datatype for the given data layout
 *
 * \remarks the returned datatype needs to be freed by the user (MPI_Type_free)
 */
MPI_Datatype xt_mpi_generate_datatype(int const * displacements, int count,
                                      MPI_Datatype old_type, MPI_Comm comm);

/**
 * generates an MPI datatype
 *
 * @param[in] displacements array of displacements
 * @param[in] blocklengths  array of block sizes
 * @param[in] count         number of blocks
 * @param[in] old_type      base MPI datatype of all elements
 * @param[in] comm          MPI communicator
 * @return                  MPI datatype for the given data layout
 *
 * \remarks the returned datatype needs to be freed by the user (MPI_Type_free)
 */
MPI_Datatype xt_mpi_generate_datatype_block(int const * displacements, int * blocklengths,
                                            int count, MPI_Datatype old_type,
                                            MPI_Comm comm);

#endif // XT_MPI_H
