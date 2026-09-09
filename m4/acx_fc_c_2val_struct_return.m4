dnl acx_fc_c_2fld_struct_return.m4 --- check if Fortran compiler supports
dnl                                    returning of structs with two
dnl                                    integral fields for functions
dnl                                    using bind(c) interface
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
# ACX_FC_C_2FLD_STRUCT_RETURN([ACTION-IF-WORKING],
#                             [ACTION-IF-NOT-WORKING],
#                             [ACTION-IF-WORKING-WITH-WORKAROUND])
# -----------------
# Checks whether the compiler correctly returns structs containing two
# values for functions using Fortran 2003 ISO-C bindings.
# This is broken for some compiler versions on x86_64 that don't
# correctly use %rax, %rdx but instead pass a pointer to a result
# struct via %rdi.
#
AC_DEFUN([ACX_FC_C_2FLD_STRUCT_RETURN],
  [AC_CACHE_CHECK([for bind(c) return of struct containing two fields],
     [acx_cv_fc_c_2fld_struct_return],
     [acx_cv_fc_c_2fld_struct_return=failed
      AC_LANG_PUSH([Fortran])
      AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[      MODULE conftest_module_c
      USE iso_c_binding, ONLY: c_long
      IMPLICIT NONE
      TYPE, PUBLIC, BIND(C) :: cl_pair
        INTEGER(c_long) :: u, v
      END TYPE cl_pair
      TYPE(cl_pair), PUBLIC, PARAMETER :: clp_ref(2) = &
        (/ cl_pair(23_c_long,     42_c_long), &
        &  cl_pair(4711_c_long, 0815_c_long) /)
      PUBLIC :: test_clp_ret
      CONTAINS
      FUNCTION test_clp_ret(idx) BIND(c, name="test_clp_ret") &
        RESULT(pair)
        INTEGER(c_long), VALUE :: idx
        TYPE(cl_pair) :: pair
        pair = clp_ref(idx)
      END FUNCTION test_clp_ret
      END MODULE conftest_module_c])],
        [ACX_MV_OBJ([conftest],[conftest_f])
         AC_LANG_PUSH([C])
         save_CFLAGS=$CFLAGS
         save_LIBS=$LIBS
         while : ; do
         AC_COMPILE_IFELSE([AC_LANG_SOURCE([@%:@include <stdbool.h>
@%:@include <stdio.h>
@%:@include <stdlib.h>

struct cl_pair {
  long u, v;
};

@%:@ifndef NEED_FC_PAIR_STRUCT_WORKAROUND
extern struct cl_pair
test_clp_ret(long idx);
@%:@else
extern void
test_clp_ret(struct cl_pair *rv_out, long idx);
@%:@endif

bool
test_clp_check(void)
{
  struct cl_pair agg;
@%:@ifndef NEED_FC_PAIR_STRUCT_WORKAROUND
  agg = test_clp_ret(2L);
@%:@else
  test_clp_ret(&agg, 2L);
@%:@endif

  bool res = (agg.u == 4711 && agg.v == 815);
  if (!res) {
    fprintf(stderr, "encountered unexpected values: u=%ld != 4711, "
            "v=%ld != 815\n", agg.u, agg.v);
  }
  return res;
}
])],
           [ACX_MV_OBJ([conftest],[conftest_c])
            AS_CASE([$LIBS],
              [conftest_c.*conftest_f.*],[],
              [LIBS="conftest_c.$ac_objext conftest_f.$ac_objext $LIBS"])
            AC_LANG_PUSH([Fortran])
            AC_RUN_IFELSE([AC_LANG_SOURCE(
[      PROGRAM conftest
      USE iso_c_binding, ONLY: c_int, c_bool
      INTERFACE
        FUNCTION test_clp_check() BIND(c, name="test_clp_check") &
          RESULT(returned_valid_struct)
          ! this function is implemented by C code to check if the
          ! struct returned from test_clp_ret made it through the
          ! language barrier
          IMPORT :: c_bool
          LOGICAL(c_bool) :: returned_valid_struct
        END FUNCTION
        SUBROUTINE posix_exit(code) BIND(c, name="exit")
          IMPORT :: c_int
          INTEGER(c_int), VALUE, INTENT(in) :: code
        END SUBROUTINE posix_exit
      END INTERFACE
      IF (.NOT. test_clp_check()) CALL posix_exit(1_c_int)
      END PROGRAM conftest])],
             [AS_CASE([$CFLAGS],
                [*\ -DNEED_FC_PAIR_STRUCT_WORKAROUND],
                [acx_cv_fc_c_2fld_struct_return="with workaround"],
                [acx_cv_fc_c_2fld_struct_return=yes])],
             [acx_cv_fc_c_2fld_struct_return=failed],
             [AC_MSG_NOTICE([Skipping run test for 2-field-struct return in cross-compilation mode,])
	      AC_MSG_NOTICE([2-field-struct-return link-test succeeded for $FC.])
              acx_cv_fc_c_2fld_struct_return=yes])
            AC_LANG_POP([Fortran])])
         AS_CASE(["$CFLAGS:$acx_cv_fc_c_2fld_struct_return"],
           [*\ -DNEED_FC_PAIR_STRUCT_WORKAROUND:*|*:yes],[break])
         CFLAGS="$CFLAGS -DNEED_FC_PAIR_STRUCT_WORKAROUND"
         done
         CFLAGS=$save_CFLAGS
         LIBS=$save_LIBS
         AC_LANG_POP([C])])
      AC_LANG_POP([Fortran])])
dnl
   AS_CASE([$acx_cv_fc_c_2fld_struct_return],
     [yes],[$1],
m4_ifval(     [$3],[[with\ workaround],[$3],]
)dnl
     [*],
     [m4_default([$2],
       [AC_MSG_WARN([Could not return 2-field struct to Fortran program via ISO C bindings!])])])
dnl
])dnl ACX_FC_C_2FLD_STRUCT_RETURN
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
