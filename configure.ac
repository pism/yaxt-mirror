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
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT

