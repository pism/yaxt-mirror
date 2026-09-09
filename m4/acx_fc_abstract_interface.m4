dnl acx_fc_abstract_interface.m4 --- check if Fortran compiler
dnl                                  permits ABSTRACT INTERFACE blocks
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
# ACX_FC_ABSTRACT_INTERFACE([ACTION-IF-WORKING],[ACTION-IF-NOT-WORKING])
# -----------------
# Checks whether the compiler supports abstract interface blocks
# to declare function signatures.
#
#
#
AC_DEFUN([ACX_FC_ABSTRACT_INTERFACE],
  [AC_CACHE_CHECK([for abstract interface support],
     [acx_cv_fc_abstract_interface],
     [AC_LANG_PUSH([Fortran])
      FCFLAGS_save=$FCFLAGS
      FCFLAGS="$FCFLAGS $FCFLAGS_f90"
      acx_cv_fc_abstract_interface=no
      AC_COMPILE_IFELSE(
[#ifdef __G95__
      ! g95 succeeds in the above test but fails when compiling the actual
      ! code for unknown reasons, this part of the test ensures a failure
      ! todo: find test that more closely matches the actual code
      MODULE conftest_module_o
      IMPLICIT NONE
      INTERFACE
        FUNCTION p(j)
        END FUNCTION p
        SUBROUTINE as
          IMPORT p
        END SUBROUTINE as
      END INTERFACE
      END MODULE conftest_module_o
#else
      MODULE conftest_module_a
      USE iso_c_binding, ONLY: c_ptr
      IMPLICIT NONE
      TYPE, BIND(C) :: eggplant
        TYPE(c_ptr) :: avocado
      END TYPE eggplant
      END MODULE conftest_module_a
      MODULE conftest_module_b
      USE conftest_module_a, ONLY: eggplant
      USE iso_c_binding, ONLY: c_ptr, c_associated, c_float, c_loc
      IMPLICIT NONE
      PRIVATE
      TYPE, BIND(C) :: peach
        TYPE(c_ptr) :: citrus
      END TYPE peach
      ABSTRACT INTERFACE
        SUBROUTINE apple(a, b, c, d)
          IMPORT :: eggplant, peach, c_float
          TYPE(eggplant), VALUE :: a
          TYPE(peach), INTENT(in) :: b
          REAL, INTENT(in) :: c
          REAL(c_float), TARGET, INTENT(out) :: d
        END SUBROUTINE apple
      END INTERFACE
      PUBLIC :: peach, cherry, apple
      CONTAINS
      SUBROUTINE banana(a, b, c, d)
        TYPE(eggplant), VALUE :: a
        TYPE(peach), INTENT(in) :: b
        REAL, INTENT(in) :: c
        REAL(c_float), TARGET, INTENT(out) :: d
        d = REAL(c * 0.5 + MERGE(1.0, 0.0, c_associated(c_loc(d), &
          a%avocado)), c_float)
      END SUBROUTINE banana
      FUNCTION cherry() RESULT(b_ptr)
        PROCEDURE(apple), POINTER :: b_ptr
        b_ptr => banana
      END FUNCTION cherry
      END MODULE conftest_module_b
#endif
],
        [AC_COMPILE_IFELSE(
[      MODULE conftest_module_c
      USE iso_c_binding, ONLY: c_int, c_float, c_null_ptr, c_loc, &
        c_associated
      USE conftest_module_a, ONLY: eggplant
      USE conftest_module_b, ONLY: peach, cherry, apple
      IMPLICIT NONE
      INTERFACE
        SUBROUTINE posix_exit(code) BIND(c, name="exit")
          IMPORT :: c_int
          INTEGER(c_int), VALUE, INTENT(in) :: code
        END SUBROUTINE posix_exit
      END INTERFACE
      TYPE(eggplant) :: berry
      TYPE(peach) :: melon
      CONTAINS
      SUBROUTINE test_main
      LOGICAL :: l
      berry%avocado = c_null_ptr
      l = orange()
      IF (.NOT. l) CALL posix_exit(1_c_int)
      END SUBROUTINE test_main
      FUNCTION orange()
        PROCEDURE(apple), POINTER :: strawberry
        LOGICAL :: orange
        REAL :: x
        strawberry => cherry()
        CALL strawberry(berry, melon, 2.0, x)
        orange = x <= 1.0 .AND. x >= 1.0
      END FUNCTION orange
      END MODULE conftest_module_c
      MODULE conftest_module_n
      IMPLICIT NONE
      ABSTRACT INTERFACE
        FUNCTION apricot(kiwi) RESULT(rhubarp)
          INTEGER :: rhubarp
          INTEGER, INTENT(in) :: kiwi
        END FUNCTION apricot
      END INTERFACE
      INTERFACE
        SUBROUTINE durian(init, destroy)
          IMPORT :: apricot
          PROCEDURE(apricot) :: init
          INTERFACE
            SUBROUTINE destroy(kiwi)
              INTEGER, INTENT(in) :: kiwi
            END SUBROUTINE destroy
          END INTERFACE
        END SUBROUTINE durian
      END INTERFACE
      PROCEDURE(apricot), POINTER :: papaya => NULL()
      CONTAINS
      FUNCTION mango(kiwi) RESULT(rhubarp)
      INTEGER :: rhubarp
      INTEGER, INTENT(in) :: kiwi
      rhubarp = 0
      END FUNCTION mango
      SUBROUTINE coconut(kiwi)
        INTEGER, INTENT(in) :: kiwi
      END SUBROUTINE coconut
      SUBROUTINE tomato()
        CALL durian(mango, coconut)
      END SUBROUTINE tomato
      END MODULE conftest_module_n
#ifdef __INTEL_COMPILER
      MODULE conftest_module_i
        IMPLICIT NONE
        ABSTRACT INTERFACE
          FUNCTION a(b)
            REAL :: b, a
          END FUNCTION a
        END INTERFACE
      CONTAINS
        FUNCTION c RESULT(d)
          PROCEDURE(a), POINTER :: d
          CALL e(d)
        END FUNCTION c
      END MODULE conftest_module_i
#endif
],
           [acx_cv_fc_abstract_interface=yes])])
   FCFLAGS=$FCFLAGS_save
   AC_LANG_POP([Fortran])])
dnl
   AS_IF([test x"$acx_cv_fc_abstract_interface" = xyes],
     [$1],
     [m4_default([$2],
       [AC_MSG_WARN([Could not build Fortran program containing abstract interface!])])])
dnl
])dnl ACX_FC_ABSTRACT_INTERFACE
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
