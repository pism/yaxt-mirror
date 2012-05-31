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
AC_PROG_FC
AC_FC_SRCEXT([f90])
AC_FC_FREEFORM
AC_PROG_RANLIB
dnl ######################################################################
dnl Checks for header files.
dnl ######################################################################
dnl
AC_HEADER_STDC

#  ----------------------------------------------------------------------
# Compile with MPI support
AC_ARG_ENABLE(mpi,AS_HELP_STRING([--enable-mpi],[Compile with MPI compiler [default=no]]),enable_mpi=${enableval}
,enable_mpi=no)
AS_IF([test x"${enable_mpi}" = x"yes"],
      [USE_MPI=yes])
AS_IF([test x"$USE_MPI" = xyes],
  [AC_DEFINE([USE_MPI],[1],[parallel I/O requested and available])
   AC_CHECK_PROGS([MPI_LAUNCH],[mpirun mpiexec],[])
   AS_IF([test x$MPI_LAUNCH = x],[])],
  [MPI_LAUNCH="`pwd`/util/serialrun"
   AC_SUBST([MPI_LAUNCH])])
AM_CONDITIONAL([USE_MPI],[test x"$USE_MPI" = xyes])
AC_SUBST([USE_MPI])

AC_CONFIG_HEADER([include/config.h])

AC_CONFIG_FILES([tests/test_redist_p2p_parallel_run	\
	tests/test_xmap_parallel_run			\
	tests/test_idxvec_run				\
	tests/test_idxlist_collection_run		\
	tests/test_idxsection_run			\
	tests/test_redist_p2p_run			\
	tests/test_xmap_run				\
	util/serialrun],[chmod a+x "$ac_file"])
AC_OUTPUT([Makefile src/Makefile tests/Makefile])

