dnl acx_use_libtool_configuration.m4 --- prevent problematic libtool build configurations
dnl
dnl Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords:
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://www.dkrz.de/redmine/projects/scales-ppm
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
dnl _ACX_LT_FLAGS_MANGLE
AC_DEFUN([_ACX_LT_FORT_FLAGS_MANGLE],
  [_AC_FORTRAN_ASSERT
   AC_LANG_CASE([Fortran],[acx_FC=$FC acx_FCFLAGS=$FCFLAGS],
     [Fortran 77],[acx_FC=$F77 acx_FCFLAGS=$FFLAGS])
   acx_temp=`$acx_FC -V 2>&1`
dnl fix problems from NAG compiler
   AS_IF([echo "$acx_temp" | grep '^NAG Fortran Compiler Release' >/dev/null],
     [acx_FCFLAGS=`echo "$acx_FCFLAGS" | tr ' ' '\n' | sed -e '/^-W@<:@lc@:>@/{' \
-e 's/^\(-Wl\)/-Xlinker \1/;s/^\(-Wc\)/-Xcompiler \1/' -e '}'     | tr '\n' ' '`])
dnl fix conflicting use of -module by libtool and ifort
   AS_CASE(["x${FC_MODOUT}x"],
     [x'-module 'x|x'-mod 'x],
     [FC_MODOUT="-Xcompiler ${FC_MODOUT}-Xcompiler "])
   AC_LANG_CASE([Fortran],[FC=$acx_FC FCFLAGS=$acx_FCFLAGS],
     [Fortran 77],[FC=$acx_F77 FCFLAGS=$acx_FFLAGS])])
dnl
dnl ACX_USE_LIBTOOL_CONFIGURATION
dnl Switch compiler to libtool wrapper and prevent occurrence of
dnl problematic setups
AC_DEFUN([ACX_USE_LIBTOOL_CONFIGURATION],
  [dnl before switching on libtool, identify compilers that prevent us from
   dnl certain build configurations
   ACX_LT_PROBLEMS
   AC_REQUIRE([LT_INIT])
   dnl _KPSE_USE_LIBTOOL ensures libtool is also used for configure-time tests,
   dnl which deduces dependent libraries automatically
   _KPSE_USE_LIBTOOL
   AC_PROVIDE_IFELSE([AC_PROG_FC],
     [AS_IF([test -n "$FC" -a X"$FC" != Xno],
        [AC_LANG_PUSH([Fortran])
         _ACX_LT_FORT_FLAGS_MANGLE
         AC_LANG_POP([Fortran])])])
   AC_PROVIDE_IFELSE([AC_PROG_F77],
     [AS_IF([test -n "$F77" -a X"$F77" != Xno],
        [AC_LANG_PUSH([Fortran 77])
         _ACX_LT_FORT_FLAGS_MANGLE
         AC_LANG_POP([Fortran 77])])])
dnl substitute -shared-intel if present
   AS_FOR([acx_flag_var],[acx_flag_var_],[CFLAGS CXXFLAGS FCFLAGS F77FLAGS LDFLAGS],
     [AS_IF([eval test x\$\{acx_flag_var+set\} = xset],
        [eval acx_temp="\" \$$acx_flag_var_ \""
         AS_CASE([$acx_temp],[*\ -shared-intel\ *],
           [eval acx_flag_var=\"`echo "$acx_temp" | sed -e 's/ -shared-intel / -Xcompiler -shared-intel /' -e 's/^ //;s/ $//'`\"])
dnl take care of ifort/icc/icpc two-part options
         eval acx_temp="\" \$$acx_flag_var_ \""
         eval acx_flag_var=\"`echo "$acx_temp" | sed -e 's/ -\(align\|allow\|assume\|ccdefault\|check\|convert\|debug\|debug-parameters\|diag-type\|diag-enable\|double-size\|dynamic-linker\|dyncom\|export-dir\|extend-source\|fp-model\|fpscomp\|gen-interfaces\|heap-arrays\|imacros\|integer-size\|iprefix\|iquote\|iwithprefixbefore\|module\|names\|opt-report\|opt-streaming-stores\|pch-dir\|pch-use\|prof-dir\|prof-file\|real-size\|reentrancy\|stand\|tcollect-filter\|tune\|warn\|watch\) \(@<:@^-@:>@@<:@^ @:>@*\) / -Xcompiler -\1 -Xcompiler \2 /g' -e 's/^ //;s/ $//'`\"])])
   AC_LANG_PUSH([C])
   _KPSE_CHECK_LIBTOOL
   AC_LANG_POP([C])
   AC_PROVIDE_IFELSE([AC_PROG_FC],
     [AS_IF([test -n "$FC" -a X"$FC" != Xno],
        [AC_LANG_PUSH([Fortran])
         _KPSE_CHECK_LIBTOOL
         AC_LANG_POP([Fortran])])])
   AC_PROVIDE_IFELSE([AC_PROG_F77],
     [AS_IF([test -n "$F77" -a X"$F77" != Xno],
        [AC_LANG_PUSH([Fortran 77])
         _KPSE_CHECK_LIBTOOL
         AC_LANG_POP([Fortran 77])])])
   AC_PROVIDE_IFELSE([AC_PROG_CXX],
     [AS_IF([test -n "$CXX" -a X"$CXX" != Xno],
        [AC_LANG_PUSH([C++])
         _KPSE_CHECK_LIBTOOL
         AC_LANG_POP([C++])])])])dnl
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
