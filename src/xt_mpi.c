#include <stdlib.h>
#include <stdio.h>

#include "xt_mpi.h"

//taken from http://beige.ucs.indiana.edu/I590/node85.html
void mpi_err_handler(int error_code, MPI_Comm comm) {


   if (error_code != MPI_SUCCESS) {

      int rank;

      MPI_Comm_rank(comm, &rank);

      char error_string[1024];
      int length_of_error_string, error_class;

      MPI_Error_class(error_code, &error_class);
      MPI_Error_string(error_class, error_string, &length_of_error_string);
      fprintf(stderr, "%3d: %s\n", rank, error_string);
      MPI_Error_string(error_code, error_string, &length_of_error_string);
      fprintf(stderr, "%3d: %s\n", rank, error_string);
      MPI_Abort(comm, error_code);
   }
}
