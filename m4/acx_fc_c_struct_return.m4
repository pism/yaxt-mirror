dnl acx_fc_c_struct_return.m4 --- check if Fortran compiler supports
dnl                               returning of structs with a single
dnl                               pointer field for functions
dnl                               via bind(c) interface
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
# ACX_FC_C_STRUCT_RETURN_FROM_C([ACTION-IF-WORKING],
#                               [ACTION-IF-NOT-WORKING],
#                               [ACTION-IF-WORKING-WITH-WORKAROUND])
# -----------------
# Checks whether the compiler correctly calls functions returning
# structs containing one pointer for functions using Fortran 2003 ISO-C bindings.
# This is broken for some compiler versions on x86_64 that don't
# correctly use %rax, but instead pass a pointer to a result
# struct passed in as hidden first argument via %rdi.
#
AC_DEFUN([ACX_FC_C_STRUCT_RETURN_FROM_C],
  [AC_CACHE_CHECK([for bind(c) return of struct containing one field from C],
     [acx_cv_fc_c_struct_return_from_c],
     [acx_cv_fc_c_struct_return_from_c=failed
      AC_LANG_PUSH([C])
      save_CFLAGS=$CFLAGS
      save_LIBS=$LIBS
      LIBS="conftest_c.$ac_objext $LIBS"
      while : ; do
      AC_COMPILE_IFELSE([AC_LANG_SOURCE([@%:@include <stdbool.h>
@%:@include <stdio.h>
@%:@include <stdlib.h>


int banana@<:@2@:>@ = { 23, 42 };

struct ptr_cnt {
  void *p;
};

@%:@ifndef NEED_FC_STRUCT_RET_FROM_C_WORKAROUND
struct ptr_cnt
test_clp_ret(int idx, int ofs);
@%:@else
void
test_clp_ret(struct ptr_cnt *rv_out, int idx, int ofs);
@%:@endif

@%:@ifndef NEED_FC_STRUCT_RET_FROM_C_WORKAROUND
struct ptr_cnt
test_clp_ret(int idx, int ofs)
@%:@else
void
test_clp_ret(struct ptr_cnt *rv_out, int idx, int ofs)
@%:@endif
{
  struct ptr_cnt agg = { .p = banana + idx - ofs };
@%:@ifndef NEED_FC_STRUCT_RET_FROM_C_WORKAROUND
  return agg;
@%:@else
  *rv_out = agg;
@%:@endif
}
])],
        [ACX_MV_OBJ([conftest],[conftest_c])
         AC_LANG_PUSH([Fortran])
         AC_RUN_IFELSE([AC_LANG_SOURCE(
[      MODULE conftest_b
      USE iso_c_binding, ONLY: c_int
      IMPLICIT NONE
      PRIVATE
      INTEGER(c_int), PUBLIC, BIND(c,name='banana'), TARGET :: banana(2)
      END MODULE conftest_b
      PROGRAM conftest
      USE iso_c_binding, ONLY: c_int, c_ptr, c_null_ptr, c_loc, c_associated
      USE conftest_b, ONLY: banana
      IMPLICIT NONE
      TYPE, BIND(C) :: ptr_cnt
        TYPE(c_ptr) :: p = c_null_ptr
      END TYPE ptr_cnt
      INTEGER(c_int), PARAMETER :: ref_idx = 2_c_int
      INTERFACE
        FUNCTION test_clp_ret(idx, ofs) BIND(c, name="test_clp_ret") &
          RESULT(p_int)
          IMPORT :: ptr_cnt, c_int
          INTEGER(c_int), VALUE :: idx, ofs
          TYPE(ptr_cnt) :: p_int
        END FUNCTION test_clp_ret
        SUBROUTINE posix_exit(code) BIND(c, name="exit")
          IMPORT :: c_int
          INTEGER(c_int), VALUE, INTENT(in) :: code
        END SUBROUTINE posix_exit
      END INTERFACE
      TYPE(ptr_cnt) :: rv
      TYPE(c_ptr) :: ref_p
      rv = test_clp_ret(ref_idx, 1_c_int)
      ref_p = C_LOC(banana(ref_idx))
      IF (.NOT. c_associated(rv%p, ref_p)) CALL posix_exit(1_c_int)
      END PROGRAM conftest])],
           [AS_CASE([$CFLAGS],
              [*\ -DNEED_FC_STRUCT_RET_FROM_C_WORKAROUND],
              [acx_cv_fc_c_struct_return_from_c="with workaround"],
              [acx_cv_fc_c_struct_return_from_c=yes])],
           [acx_cv_fc_c_struct_return_from_c=failed],
           [AC_MSG_NOTICE([Skipping run test for struct return in cross-compilation mode,])
            AC_MSG_NOTICE([struct-return link-test succeeded for $FC.])
            acx_cv_fc_c_struct_return_from_c=yes])
         AC_LANG_POP([Fortran])])
      AS_CASE(["$CFLAGS:$acx_cv_fc_c_struct_return_from_c"],
        [*\ -DNEED_FC_STRUCT_RET_FROM_C_WORKAROUND:*|*:yes],[break])
      CFLAGS="$CFLAGS -DNEED_FC_STRUCT_RET_FROM_C_WORKAROUND"
      done
      CFLAGS=$save_CFLAGS
      LIBS=$save_LIBS
      AC_LANG_POP([C])])
dnl
   AS_CASE([$acx_cv_fc_c_struct_return_from_c],
     [yes],[$1],
m4_ifval(     [$3],[[with\ workaround],[$3],]
)dnl
     [*],
     [m4_default([$2],
       [AC_MSG_WARN([Could not return 2-field struct to Fortran program via ISO C bindings!])])])
dnl
])dnl ACX_FC_C_STRUCT_RETURN_FROM_C
#
# ACX_FC_C_STRUCT_RETURN_FROM_FORTRAN([ACTION-IF-WORKING],
#                                     [ACTION-IF-NOT-WORKING],
#                                     [ACTION-IF-WORKING-WITH-WORKAROUND])
# -----------------
# Checks whether the compiler correctly returns structs containing one
# pointer for functions using Fortran 2003 ISO-C bindings.
# This is broken for some compiler versions on x86_64 that don't
# correctly use %rax, but instead pass a pointer to a result
# struct passed in as hidden first argument via %rdi.
#
AC_DEFUN([ACX_FC_C_STRUCT_RETURN_FROM_FORTRAN],
  [AC_CACHE_CHECK([for bind(c) return of struct containing one field],
     [acx_cv_fc_c_struct_return_from_fortran],
     [acx_cv_fc_c_struct_return_from_fortran=failed
      AC_LANG_PUSH([Fortran])
      AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[      MODULE conftest_module_c
      USE iso_c_binding, ONLY: c_ptr, c_int, c_loc, &
        c_associated, c_null_ptr
      IMPLICIT NONE
      TYPE, PUBLIC, BIND(C) :: ptr_cnt
        TYPE(c_ptr) :: p = c_null_ptr
      END TYPE ptr_cnt
      INTEGER(c_int), BIND(c, name='banana'), TARGET :: banana(2)
      PUBLIC :: test_clp_ret
      CONTAINS
      FUNCTION test_clp_ret(table, idx, ofs) BIND(c, name="test_clp_ret") &
        RESULT(p_int)
        TYPE(ptr_cnt), INTENT(in) :: table
        INTEGER(c_int), VALUE :: idx, ofs
        TYPE(ptr_cnt) :: p_int
        p_int%p = C_LOC(banana(idx+ofs))
        IF (C_ASSOCIATED(table%p, p_int%p)) WRITE (0, '(a)') 'aliasing'
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


int banana@<:@2@:>@ = { 23, 42 };

struct ptr_cnt {
  void *p;
};

@%:@ifndef NEED_FC_STRUCT_RET_FROM_FORTRAN_WORKAROUND
extern struct ptr_cnt
test_clp_ret(struct ptr_cnt *alias_check, int idx, int ofs);
@%:@else
extern void
test_clp_ret(struct ptr_cnt *rv_out, struct ptr_cnt *alias_check,
             int idx, int ofs);
@%:@endif

@%:@define BANANA_OFS 1
int
test_clp_check(void)
{
  struct ptr_cnt agg = { .p = banana };
@%:@ifndef NEED_FC_STRUCT_RET_FROM_FORTRAN_WORKAROUND
  agg = test_clp_ret(&agg, BANANA_OFS, 1);
@%:@else
  test_clp_ret(&agg, &agg, BANANA_OFS, 1);
@%:@endif

  bool res = (agg.p == banana+BANANA_OFS);
  if (!res) {
    fprintf(stderr, "encountered unexpected values: p=%p != (banana+%ld=%p)\n",
            agg.p, BANANA_OFS-1, banana+BANANA_OFS-1);
  }
  return res ? 4711 : 0xbeef;
}
])],
           [ACX_MV_OBJ([conftest],[conftest_c])
            AS_CASE([$LIBS],
              [conftest_c.*conftest_f.*],[],
              [LIBS="conftest_c.$ac_objext conftest_f.$ac_objext $LIBS"])
            AC_LANG_PUSH([Fortran])
            AC_RUN_IFELSE([AC_LANG_SOURCE(
[      PROGRAM conftest
      USE iso_c_binding, ONLY: c_int
      INTERFACE
        FUNCTION test_clp_check() BIND(c, name="test_clp_check") &
          RESULT(retcode)
          ! this function is implemented by C code to check if the
          ! struct returned from test_clp_ret made it through the
          ! language barrier
          IMPORT :: c_int
          INTEGER(c_int) :: retcode
        END FUNCTION
        SUBROUTINE posix_exit(code) BIND(c, name="exit")
          IMPORT :: c_int
          INTEGER(c_int), VALUE, INTENT(in) :: code
        END SUBROUTINE posix_exit
      END INTERFACE
      IF (test_clp_check() /= 4711_c_int) CALL posix_exit(1_c_int)
      END PROGRAM conftest])],
             [AS_CASE([$CFLAGS],
                [*\ -DNEED_FC_STRUCT_RET_FROM_FORTRAN_WORKAROUND],
                [acx_cv_fc_c_struct_return_from_fortran="with workaround"],
                [acx_cv_fc_c_struct_return_from_fortran=yes])],
             [acx_cv_fc_c_struct_return_from_fortran=failed],
             [AC_MSG_NOTICE([Skipping run test for struct return in cross-compilation mode,])
	      AC_MSG_NOTICE([struct-return link-test succeeded for $FC.])
              acx_cv_fc_c_struct_return_from_fortran=yes])
            AC_LANG_POP([Fortran])])
         AS_CASE(["$CFLAGS:$acx_cv_fc_c_struct_return_from_fortran"],
           [*\ -DNEED_FC_STRUCT_RET_FROM_FORTRAN_WORKAROUND:*|*:yes],[break])
         CFLAGS="$CFLAGS -DNEED_FC_STRUCT_RET_FROM_FORTRAN_WORKAROUND"
         done
         CFLAGS=$save_CFLAGS
         LIBS=$save_LIBS
         AC_LANG_POP([C])])
      AC_LANG_POP([Fortran])])
dnl
   AS_CASE([$acx_cv_fc_c_struct_return_from_fortran],
     [yes],[$1],
m4_ifval(     [$3],[[with\ workaround],[$3],]
)dnl
     [*],
     [m4_default([$2],
       [AC_MSG_WARN([Could not return 2-field struct to Fortran program via ISO C bindings!])])])
dnl
])dnl ACX_FC_C_STRUCT_RETURN_FROM_FORTRAN
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
