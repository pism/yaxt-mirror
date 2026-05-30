dnl acx_fc_c_sizeof.m4 --- check if Fortran compiler supports
dnl                        iso_c_bindings c_sizeof interface or the
dnl                        equivalent sizeof intrinsic of gfortran
dnl
dnl
dnl Copyright  (C)  2026  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords:
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://swprojects.dkrz.de/redmine/projects/scales-ppm
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modification, are  permitted provided that the following conditions are
dnl met:
dnl
dnl Redistributions of source code must retain the above copyright notice,
dnl this list of conditions and the following disclaimer.
dnl
dnl Redistributions in binary form must reproduce the above copyright
dnl notice, this list of conditions and the following disclaimer in the
dnl documentation and/or other materials provided with the distribution.
dnl
dnl Neither the name of the DKRZ GmbH nor the names of its contributors
dnl may be used to endorse or promote products derived from this software
dnl without specific prior written permission.
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
dnl IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
dnl TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
dnl PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
dnl OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
dnl EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
dnl PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
dnl PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
dnl LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
dnl NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
dnl SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl
dnl Commentary:
dnl
dnl
dnl
dnl Code:
dnl
dnl _ACX_FC_SIZE_QUERY([QUERY],[QUERY-DESCRIPTION],
dnl                    [ACTION-IF-WORKING],[ACTION-IF-NOT-WORKING],
dnl                    [PRE])
dnl -----------------
dnl Checks whether the compiler supports the QUERY function
dnl
AC_DEFUN([_ACX_FC_SIZE_QUERY],
  [AC_CHECK_SIZEOF([int])
   AS_VAR_PUSHDEF([acx_cv_fc_query],[acx_cv_fc_]$1)
   AC_CACHE_CHECK([$2],
     [acx_cv_fc_query],
     [AC_LANG_PUSH([Fortran])
      acx_cv_fc_query=no
      AC_COMPILE_IFELSE([      PROGRAM conftest
      USE iso_c_binding, ONLY: c_int, c_size_t
      $5
      IMPLICIT NONE
      INTERFACE
        SUBROUTINE posix_exit(code) BIND(c, name="exit")
          IMPORT :: c_int
          INTEGER(c_int), VALUE, INTENT(in) :: code
        END SUBROUTINE posix_exit
      END INTERFACE
      TYPE, BIND(C) :: pineapple
        INTEGER(c_int) :: banana
      END TYPE pineapple
      TYPE(pineapple) :: apple
      INTEGER(c_size_t) :: apple_size
      apple_size = $1(apple)
      IF (apple_size /= $ac_cv_sizeof_int) CALL posix_exit(1_c_int)
      END PROGRAM conftest],
        [AS_VAR_SET([acx_cv_fc_query],[yes])])
      AC_LANG_POP([Fortran])])
dnl
   AS_IF([test x"$acx_cv_fc_query" != xno],
     [$3],
     [m4_default([$4],
       [AC_MSG_WARN([Could not build Fortran program containing $1 query!])])])
   AS_VAR_POPDEF([acx_cv_fc_query])
dnl
])dnl ACX_FC_C_SIZEOF
#
# ACX_FC_C_SIZEOF([ACTION-IF-WORKING],[ACTION-IF-NOT-WORKING])
# -----------------
# Checks whether the compiler supports the C_SIZEOF function of
# Fortran 2008 ISO-C bindings.
#
AC_DEFUN([ACX_FC_C_SIZEOF],
  [_ACX_FC_SIZE_QUERY([c_sizeof],[for iso_c_binding c_sizeof],
     [$1],[$2],[USE iso_c_binding, ONLY: c_sizeof])])
dnl ACX_FC_C_SIZEOF
#
# ACX_FC_SIZEOF([ACTION-IF-WORKING],[ACTION-IF-NOT-WORKING])
# -----------------
# Checks whether the compiler supports the GNU gfortran SIZEOF intrinsic
#
AC_DEFUN([ACX_GFORTRAN_SIZEOF],
  [_ACX_FC_SIZE_QUERY([sizeof],[for GNU gfortran sizeof],
     [$1],[$2])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
