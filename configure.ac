dnl ######################################################################
dnl Basic autotools settings.
dnl ######################################################################
AC_PREREQ([2.64])
dnl Update the corresponding entry in doc/Doxyfile when changing this
AC_INIT([yaxt],[0.1.3],[hanke@dkrz.de])
AC_CONFIG_SRCDIR([src/xt_idxlist.c])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE([1.10 foreign])
AM_MAINTAINER_MODE
dnl ######################################################################
dnl Checks for programs.
dnl ######################################################################
AC_PROG_CC
AC_PROG_CC_C99
AX_PTHREAD([AC_DEFINE([HAVE_PTHREAD],1,[Define if you have POSIX threads libraries and header files.])
   LIBS="$PTHREAD_LIBS $LIBS"
   CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
   CC="$PTHREAD_CC"])
AC_PROG_FC
AC_FC_SRCEXT([f90])
AC_FC_FREEFORM
AC_LANG_PUSH([Fortran])
AC_PROG_FPP
AC_LANG_POP([Fortran])
AC_PROG_MKDIR_P
AX_WITH_PERL([/bin/false], $PATH:$prefix/bin:/usr/bin:/usr/local/bin)
AX_PROG_PERL_VERSION([5.0005],[],
        AC_MSG_ERROR([Perl (version>=5.0005) is required]))
LT_INIT
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
AC_FC_DUMMY_MAIN(,[:])
AC_FC_WRAPPERS
dnl
dnl Check wether char is unsigned by default
AC_C_CHAR_UNSIGNED
AS_IF([test x$ac_cv_c_char_unsigned = xyes],
  [C_CHAR_IS_UNSIGNED=1],
  [C_CHAR_IS_UNSIGNED=0])
AC_SUBST([C_CHAR_IS_UNSIGNED])
ACX_CHECK_SIZEOF_RELATION([long], [int])
ACX_CHECK_SIZEOF_RELATION([long long], [long])

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
AC_PATH_PROGS([MPI_LAUNCH],[mpirun mpiexec],[true])
AS_IF([test x"$cross_compiling" = xno],
  [AC_MSG_CHECKING([if $MPI_LAUNCH works])
   saved_CFLAGS=$CFLAGS
   saved_LIBS=$LIBS
   CFLAGS="$CFLAGS $MPI_C_INCLUDE"
   LIBS="$LIBS $MPI_C_LIB"
   AC_LINK_IFELSE([AC_LANG_SOURCE([
@%:@include <stdio.h>
@%:@include <stdlib.h>

@%:@include <mpi.h>

@%:@define xmpi(ret)           \\
  do {                      \\
    if (ret != MPI_SUCCESS) \\
      exit(EXIT_FAILURE);   \\
  } while (0)

int
main(int argc, char **argv)
{
  xmpi(MPI_Init(&argc, &argv));
  char *numarg = argv@<:@1@:>@;
  int cmdnum = atoi(numarg);
  int procnum = 1;
  xmpi(MPI_Allreduce(MPI_IN_PLACE, &procnum, 1, MPI_INT, MPI_SUM,
                     MPI_COMM_WORLD));
  xmpi(MPI_Finalize());
  return (procnum == cmdnum)?EXIT_SUCCESS:EXIT_FAILURE;
}
])],
     [AS_IF([$MPI_LAUNCH -n 4 ./conftest$EXEEXT 4],
        [AC_MSG_RESULT([yes])],
        [AC_MSG_RESULT([no])
         AC_MSG_FAILURE([mpirun doesn't work])])],
     [AC_MSG_FAILURE([Cannot compile simple MPI program])])
   CFLAGS=$saved_CFLAGS
   LIBS=$saved_LIBS
  ])
AS_IF([test "x$MPI_LAUNCH" = xtrue],
  [AC_MSG_WARN([MPI launch command unavailable])])

# enable user to select an alternative enumeration type in case
# 2^31 elements isn't enough or too much
AC_ARG_WITH([idxtype],
  AS_HELP_STRING([--with-idxtype],[select type to use for global index values @<:@default=int@:>@]),
  [XT_INT=$withval
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[[char a[((($XT_INT)-1)<0)];
char b[($XT_INT)((($XT_INT)-1.5)+0.5)]]])],,
     [AC_MSG_ERROR([Selected index type not compilable or unsigned])])],
  [XT_INT=int])
m4_foreach([substtype],[[XT_INT]],
  [AS_CASE(["$substtype"],
     [int],
     substtype[_MPIDT=MPI_INT
     ACX_FORTRAN_TYPE_KIND([integer],[c_int],[int],
       [use iso_c_binding],,substtype[_FC_KIND=$acx_fortran_kind_subst])
     ],
     [long],
     substtype[_MPIDT=MPI_LONG
     ACX_FORTRAN_TYPE_KIND([integer],[c_long],[long],
       [use iso_c_binding],,substtype[_FC_KIND=$acx_fortran_kind_subst])
     ],
     [short],
     substtype[_MPIDT=MPI_SHORT
     ACX_FORTRAN_TYPE_KIND([integer],[c_short],[short],
       [use iso_c_binding],,substtype[_FC_KIND=$acx_fortran_kind_subst])
     ],
     ["long long"],
     substtype[_MPIDT=MPI_LONG_LONG
     ACX_FORTRAN_TYPE_KIND([integer],[c_long_long],[long long],
       [use iso_c_binding],,substtype[_FC_KIND=$acx_fortran_kind_subst])
     ],
     [AC_MSG_ERROR([Invalid type selected for substtype])])
   FCFLAGS="$FCFLAGS ${FPP_DEFOPT}substtype[]_FC_KIND=$substtype[]_FC_KIND"
  ])

AC_DEFINE_UNQUOTED([XT_INT],[$XT_INT],[Type to use for enumeration of distributed elements.])
AC_DEFINE_UNQUOTED([XT_INT_MPIDT],[$XT_INT_MPIDT],[Type to use for enumeration of distributed elements.])

AC_ARG_VAR([BUILD_CFLAGS],
  [append to CFLAGS during build but not in configure phase])
AC_CONFIG_COMMANDS_PRE([CFLAGS="$CFLAGS${BUILD_CFLAGS+ $BUILD_CFLAGS}"])
AC_ARG_VAR([BUILD_FCFLAGS],
  [append to FCFLAGS during build but not in configure phase])
AC_CONFIG_COMMANDS_PRE([FCFLAGS="$FCFLAGS${BUILD_FCFLAGS+ $BUILD_FCFLAGS}"])
AC_ARG_VAR([BUILD_LDFLAGS],
  [append to LDFLAGS during build but not in configure phase])
AC_CONFIG_COMMANDS_PRE([LDFLAGS="$LDFLAGS${BUILD_LDFLAGS+ $BUILD_LDFLAGS}"])

AC_CONFIG_HEADER([include/config.h])

AC_CONFIG_FILES([						\
	perf/perf_idxsection_get_positions_of_indices_run	\
	tests/test_redist_p2p_parallel_run			\
	tests/test_redist_collection_parallel_run		\
	tests/test_redist_collection_static_parallel_run	\
	tests/test_xmap_all2all_parallel_run			\
	tests/test_xmap_dist_dir_parallel_run			\
	tests/test_idxempty_run					\
	tests/test_idxvec_run					\
	tests/test_idxlist_collection_run			\
	tests/test_idxmod_run                                  	\
	tests/test_idxsection_run				\
	tests/test_idxstripes_run				\
	tests/test_redist_collection_run			\
	tests/test_redist_collection_static_run			\
	tests/test_redist_p2p_run				\
	tests/test_xmap_dist_dir_run				\
	tests/test_xmap_all2all_run				\
	tests/test_xmap_all2all_fail_run			\
	tests/test_handles_run					\
	tests/test_ut_run					\
        tests/test_perf_run              			\
        tests/test_perf_stripes_run              		\
        tests/test_sort_run              			\
        tests/test_uid_run              			\
        tests/test_yaxt_run              			\
	tests/test_mpi_generate_datatype_run			\
	examples/rrobin_run					\
	examples/tr_gp2fs_run					\
	util/serialrun],[chmod a+x "$ac_file"])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile examples/Makefile \
  inst_headers/Makefile inst_headers/f90/Makefile \
  src/pkgconfig/yaxt.pc perf/Makefile])
# force make to regenerate FC.deps file (because it depends on Makefile)
AC_CONFIG_COMMANDS([force_dep_regeneration],[sleep 1
touch tests/Makefile])

AC_OUTPUT
