dnl acx_fc_c_funloc.m4 --- check if Fortran compiler supports
dnl                        iso_c_bindings c_funloc interface
dnl
dnl
dnl Copyright  (C)  2024  Thomas Jahns <jahns@dkrz.de>
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
# ACX_FC_C_FUNLOC([ACTION-IF-WORKING],[ACTION-IF-NOT-WORKING])
# -----------------
# Checks whether the compiler supports the C_FUNLOC and C_FUNPTR parts of
# Fortran 2003 ISO-C bindings.
#
AC_DEFUN([ACX_FC_C_FUNLOC],
  [AC_CACHE_CHECK([for iso_c_binding c_funloc/c_funptr],
     [acx_cv_fc_c_funloc],
     [AC_LANG_PUSH([Fortran])
      AC_COMPILE_IFELSE([      PROGRAM conftest
      USE iso_c_binding, ONLY: c_int, c_ptr, c_funloc, c_funptr, &
        c_associated
      INTERFACE
        SUBROUTINE free(p) BIND(c, name="free")
          IMPORT :: c_ptr
          TYPE(c_ptr), VALUE :: p
        END SUBROUTINE
        SUBROUTINE posix_exit(code) BIND(c, name="exit")
          IMPORT :: c_int
          INTEGER(c_int), VALUE, INTENT(in) :: code
        END SUBROUTINE posix_exit
      END INTERFACE
      TYPE(c_funptr) :: p_free_c, p_exit_c
      p_free_c = C_FUNLOC(free)
      p_exit_c = C_FUNLOC(posix_exit)
      IF (c_associated(p_free_c, p_exit_c)) CALL posix_exit(1_c_int)
      END PROGRAM conftest],
        [acx_cv_fc_c_funloc=yes],
        [acx_cv_fc_c_funloc=no])
      AC_LANG_POP([Fortran])])
dnl
   AS_IF([test x"$acx_cv_fc_c_funloc" = xyes],
     [$1],
     [m4_default([$2],
       [AC_MSG_WARN([Could not build Fortran program containing C function pointer!])])])
dnl
])dnl ACX_FC_C_FUNLOC
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
