!>
!> @file xt_idxmod_f.f90
!> @brief Fortran interface to yaxt xt_idxmod functions
!>
!> @copyright Copyright  (C)  2018 Jörg Behrens <behrens@dkrz.de>
!>                                 Moritz Hanke <hanke@dkrz.de>
!>                                 Thomas Jahns <jahns@dkrz.de>
!>
!> @author Jörg Behrens <behrens@dkrz.de>
!>         Moritz Hanke <hanke@dkrz.de>
!>         Thomas Jahns <jahns@dkrz.de>
!>

!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://doc.redmine.dkrz.de/yaxt/html/
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
MODULE xt_idxmod
  USE xt_core, ONLY: i2, i4, xt_int_kind, xt_abort
  USE xt_idxlist_abstract, ONLY: xt_idxlist, xt_idxlist_c2f, xt_idxlist_f2c
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_int, c_null_ptr, c_loc
  IMPLICIT NONE
  PRIVATE

  TYPE, BIND(C), PUBLIC :: xt_modifier
    TYPE(xt_idxlist) :: extract
    TYPE(xt_idxlist) :: subst
    INTEGER(c_int) :: mask
  END TYPE xt_modifier

  INTERFACE xt_idxmod_new
    MODULE PROCEDURE xt_idxmod_new_a1d
    MODULE PROCEDURE xt_idxmod_new_a1d_a1d
    MODULE PROCEDURE xt_idxmod_new_a1d_i4
    MODULE PROCEDURE xt_idxmod_new_a1d_i4_a1d
    MODULE PROCEDURE xt_idxmod_new_a1d_i4_a2d
  END INTERFACE xt_idxmod_new

  PUBLIC :: xt_idxmod_new

  INTERFACE
    FUNCTION xt_idxmod_new_c(patch, modifier, modifier_num, mstate_ptr) &
         BIND(C, name='xt_idxmod_new') RESULT(res)
      IMPORT :: xt_int_kind, xt_modifier, c_ptr, c_int
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in) :: patch
      INTEGER(c_int), VALUE, INTENT(in) :: modifier_num
      TYPE(xt_modifier), INTENT(in) :: modifier(modifier_num)
      TYPE(c_ptr), VALUE, INTENT(in) :: mstate_ptr
      TYPE(c_ptr) :: res
    END FUNCTION xt_idxmod_new_c
  END INTERFACE

CONTAINS

  FUNCTION xt_idxmod_new_a1d(patch, modifier) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: patch
    TYPE(xt_modifier), INTENT(in) :: modifier(:)
    TYPE(Xt_idxlist) :: res

    INTEGER :: num_modifier
    INTEGER(c_int) :: num_modifier_c
    num_modifier = SIZE(modifier)
    IF (num_modifier > HUGE(1_c_int)) &
         CALL xt_abort("number of modifiers too high", &
         __FILE__, &
         __LINE__)
    num_modifier_c = INT(num_modifier, c_int)
    res = xt_idxlist_c2f(xt_idxmod_new_c(xt_idxlist_f2c(patch), modifier, &
         num_modifier_c, c_null_ptr))
  END FUNCTION xt_idxmod_new_a1d

  FUNCTION xt_idxmod_new_a1d_a1d(patch, modifier, mstate) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: patch
    TYPE(xt_modifier), INTENT(in) :: modifier(:)
    INTEGER(c_int), TARGET, INTENT(inout) :: mstate(*)
    TYPE(Xt_idxlist) :: res

    INTEGER :: num_modifier
    INTEGER(c_int) :: num_modifier_c
    num_modifier = SIZE(modifier)
    IF (num_modifier > HUGE(1_c_int)) &
         CALL xt_abort("number of modifiers too high", &
         __FILE__, &
         __LINE__)
    num_modifier_c = INT(num_modifier, c_int)
    res = xt_idxlist_c2f(xt_idxmod_new_c(xt_idxlist_f2c(patch), modifier, &
         num_modifier_c, C_LOC(mstate)))
  END FUNCTION xt_idxmod_new_a1d_a1d

  FUNCTION xt_idxmod_new_a1d_i4(patch, modifier, num_modifier) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: patch
    TYPE(xt_modifier), INTENT(in) :: modifier(*)
    INTEGER(i4), INTENT(in) :: num_modifier
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_modifier_c

    num_modifier_c = INT(num_modifier, c_int)
    res = xt_idxlist_c2f(xt_idxmod_new_c(xt_idxlist_f2c(patch), modifier, &
         num_modifier_c, c_null_ptr))
  END FUNCTION xt_idxmod_new_a1d_i4

  FUNCTION xt_idxmod_new_a1d_i4_a1d(patch, modifier, num_modifier, mstate) &
       RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: patch
    TYPE(xt_modifier), INTENT(in) :: modifier(*)
    INTEGER(i4), INTENT(in) :: num_modifier
    INTEGER(c_int), TARGET, INTENT(inout) :: mstate(*)
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_modifier_c

    num_modifier_c = INT(num_modifier, c_int)
    res = xt_idxlist_c2f(xt_idxmod_new_c(xt_idxlist_f2c(patch), modifier, &
         num_modifier_c, C_LOC(mstate)))
  END FUNCTION xt_idxmod_new_a1d_i4_a1d

  FUNCTION xt_idxmod_new_a1d_i4_a2d(patch, modifier, num_modifier, mstate) &
       RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: patch
    TYPE(xt_modifier), INTENT(in) :: modifier(*)
    INTEGER(i4), INTENT(in) :: num_modifier
    INTEGER(c_int), TARGET, INTENT(inout) :: mstate(1,*)
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_modifier_c

    num_modifier_c = INT(num_modifier, c_int)
    res = xt_idxlist_c2f(xt_idxmod_new_c(xt_idxlist_f2c(patch), modifier, &
         num_modifier_c, C_LOC(mstate)))
  END FUNCTION xt_idxmod_new_a1d_i4_a2d

END MODULE xt_idxmod
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
