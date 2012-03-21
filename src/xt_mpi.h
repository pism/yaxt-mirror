#ifndef XT_MPI_H
#define XT_MPI_H

#include "mpi.h"

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

#endif // XT_MPI_H
