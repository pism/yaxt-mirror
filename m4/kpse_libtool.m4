# Public macros for the TeX Live (TL) tree.
# Copyright (C) 1995 - 2009 Karl Berry <tex-live@tug.org>
# Copyright (C) 2009, 2010 Peter Breitenlohner <tex-live@tug.org>
#
# This file is free software; the copyright holders
# give unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# serial 1

# _KPSE_USE_LIBTOOL()
# Switch link tests over to use libtool so as not to require dependent
# libraries to be listed explicitly.
# Extended for Fortran by Thomas Jahns <jahns@dkrz.de>, 2015
# -------------------
AC_DEFUN([_KPSE_USE_LIBTOOL],
[## $0: Generate a libtool script for use in configure tests
AC_PROVIDE_IFELSE([LT_INIT], ,
                  [m4_fatal([$0: requires libtool])])[]dnl
LT_OUTPUT
AC_CONFIG_COMMANDS_PRE([ac_objext=${acx_lt_saved_ac_objext}])
acx_lt_saved_ac_objext=$ac_objext
ac_objext=lo
m4_append([AC_LANG(C)],
[ac_link="$ac_pwd/libtool --mode=link --tag=CC $ac_link"
ac_compile="$ac_pwd/libtool --mode=compile --tag=CC $ac_compile"
])[]dnl
AC_PROVIDE_IFELSE([AC_PROG_CXX],
[m4_append([AC_LANG(C++)],
[ac_link="$ac_pwd/libtool --mode=link --tag=CXX $ac_link"
ac_compile="$ac_pwd/libtool --mode=compile --tag=CXX $ac_compile"
])])[]dnl
AC_PROVIDE_IFELSE([AC_PROG_FC],
[m4_append([AC_LANG(Fortran)],
[ac_link="$ac_pwd/libtool --mode=link --tag=FC $ac_link"
ac_compile="$ac_pwd/libtool --mode=compile --tag=FC $ac_compile"
])])[]dnl
AC_PROVIDE_IFELSE([AC_PROG_F77],
[m4_append([AC_LANG(Fortran 77)],
[ac_link="$ac_pwd/libtool --mode=link --tag=F77 $ac_link"
ac_compile="$ac_pwd/libtool --mode=compile --tag=F77 $ac_compile"
])])[]dnl
AC_LANG(_AC_LANG)[]dnl
dnl After compiling/linking checks, now also $top_builddir/$objdir
dnl needs to be cleaned.
ac_clean_files_save="$ac_clean_files_save $ac_pwd/$objdir"
ac_clean_files="$ac_clean_files $ac_pwd/$objdir"
]) # _KPSE_USE_LIBTOOL
