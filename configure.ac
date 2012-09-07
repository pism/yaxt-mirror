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
dnl Properties of the compiler
dnl ######################################################################
AM_PROG_CC_C_O
AC_PROG_FC_C_O
ACX_SL_FC_MOD_SUFFIX(,[AC_MSG_ERROR(dnl
[Fortran module file extension must be known to install F90 modules.])])
AM_CONDITIONAL([FC_MOD_UPPERCASE],
  [test x$ac_fc_mod_uppercase = xyes])
AC_LANG_PUSH([Fortran])
ACX_SL_FC_CHECK_MOD_PATH_FLAG
AC_LANG_POP([Fortran])
ACX_FC_VOLATILE([FCFLAGS="${FCFLAGS} ${FPP_DEFOPT}HAVE_VOLATILE"])
AC_FC_WRAPPERS

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
ACX_C_PACKAGE([MPI],[mpi.h],,,
  [AC_MSG_ERROR([Required header mpi.h not found or not compilable.])],
  [MPI_Waitall],[mpi mpich],,,
  [AC_MSG_ERROR([Cannot link C MPI programs.])])
ACX_FORTRAN_PACKAGE([MPI],[mpif.h],,,
  [AC_MSG_ERROR([Required include mpif.h not found or not compilable.])],
  [mpi_waitall],[mpi mpi_f90 mpi_f77 mpich],[[-lmpi_f77 -lmpi],[-lmpi]],,
  [AC_MSG_ERROR([Cannot link Fortran MPI programs.])],[])
ACX_F90_PACKAGE([MPI],[mpi],,,
  [AC_MSG_ERROR([Required Fortran 90 module not found or incompatible.])],
  [mpi_waitall],[mpi mpi_f90 mpi_f77 mpich],[[-lmpi_f77 -lmpi],[-lmpi]],,
  [AC_MSG_ERROR([Cannot link Fortran 90 MPI programs.])],,
[      INTEGER :: req(1), stat(mpi_status_size, 1), ierror
],[(1, req, stat, ierror)])

# the following is needed for source parts from libraries with
# optional MPI dependencies
AC_DEFINE([USE_MPI],[1],[use MPI parallel process environment])
AC_CHECK_PROGS([MPI_LAUNCH],[mpirun mpiexec],[true])
AS_IF([test x$MPI_LAUNCH = xtrue],
  [AC_MSG_WARN([MPI launch command unavailable])])

AC_CONFIG_HEADER([include/config.h])

AC_CONFIG_FILES([tests/test_redist_p2p_parallel_run		\
	tests/test_redist_collection_parallel_run		\
	tests/test_redist_collection_static_parallel_run	\
	tests/test_xmap_all2all_parallel_run			\
	tests/test_xmap_dist_dir_parallel_run			\
	tests/test_idxvec_run					\
	tests/test_idxlist_collection_run			\
	tests/test_idxsection_run				\
	tests/test_idxstripes_run				\
	tests/test_redist_collection_run			\
	tests/test_redist_collection_static_run			\
	tests/test_redist_p2p_run				\
	tests/test_xmap_dist_dir_run				\
	tests/test_xmap_all2all_run				\
	tests/test_handles_run					\
	tests/test_ut_run					\
        tests/test_perf_run              			\
        tests/test_perf_stripes_run              		\
        tests/test_sort_run              			\
        tests/test_yaxt_run              			\
	tests/test_mpi_generate_datatype_run			\
	util/serialrun],[chmod a+x "$ac_file"])
AC_OUTPUT([Makefile src/Makefile tests/Makefile])
