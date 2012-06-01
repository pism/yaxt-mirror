dnl ######################################################################
dnl Basic autotools settings.
dnl ######################################################################
AC_PREREQ([2.64])
AC_INIT([yaxt],[0.0.1],[hanke@dkrz.de])
AC_CONFIG_SRCDIR([src/yaxt.h])
AC_CONFIG_AUX_DIR([config])
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE([1.10 foreign])
AM_MAINTAINER_MODE
dnl ######################################################################
dnl Checks for programs.
dnl ######################################################################
AC_PROG_CC
AC_PROG_CC_C99
AC_PROG_FC
AC_FC_SRCEXT([f90])
AC_FC_FREEFORM
AC_LANG_PUSH([Fortran])
AC_PROG_FPP
AC_LANG_POP([Fortran])
AC_PROG_RANLIB
dnl ######################################################################
dnl Checks for header files.
dnl ######################################################################
dnl
AC_HEADER_STDC

#  ----------------------------------------------------------------------
# Compile with MPI support
dnl
dnl ######################################################################
dnl Checks for MPI.
dnl ######################################################################
ACX_FORTRAN_PACKAGE([MPI],[mpif.h],,,
  [AC_MSG_ERROR([Required include mpif.h not found or not compilable.])],
  [mpi_waitall],[mpi mpi_f90 mpi_f77 mpich],[[-lmpi_f77 -lmpi],[-lmpi]],,
  [AC_MSG_ERROR([Cannot link Fortran MPI programs.])],[])
ACX_C_PACKAGE([MPI],[mpi.h],,,
  [AC_MSG_ERROR([Required header mpi.h not found or not compilable.])],
  [MPI_Waitall],[mpi mpich],,,
  [AC_MSG_ERROR([Cannot link C MPI programs.])])

# the following is needed for source parts from libraries with
# optional MPI dependencies
AC_DEFINE([USE_MPI],[1],[use MPI parallel process environment])
AC_CHECK_PROGS([MPI_LAUNCH],[mpirun mpiexec],[true])
AS_IF([test x$MPI_LAUNCH = xtrue],
  [AC_MSG_WARN([MPI launch command unavailable])])

AC_CONFIG_HEADER([include/config.h])

AC_CONFIG_FILES([tests/test_redist_p2p_parallel_run	\
	tests/test_xmap_parallel_run			\
	tests/test_idxvec_run				\
	tests/test_idxlist_collection_run		\
	tests/test_idxsection_run			\
	tests/test_redist_p2p_run			\
	tests/test_xmap_run				\
	tests/test_handles_run				\
	tests/test_ut_run				\
	util/serialrun],[chmod a+x "$ac_file"])
AC_OUTPUT([Makefile src/Makefile tests/Makefile])
