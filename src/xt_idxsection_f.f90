!> @file xt_idxsection_f.f90
!! @brief Fortran interface to yaxt implementation
!!
!! @copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!

!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://redmine.dkrz.de/doc/yaxt/html/index.html
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
MODULE xt_idxsection
  USE iso_c_binding, ONLY: c_int, c_ptr
  USE xt_core, ONLY: xt_int_kind, xt_abort, xt_get_default_comm, i2, i4, i8
  USE xt_idxlist_abstract, ONLY: xt_idxlist, xt_idxlist_c2f
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xt_idxsection_new

  INTERFACE xt_idxsection_new
    MODULE PROCEDURE xt_idxsection_new_a
    MODULE PROCEDURE xt_idxsection_new_i2
    MODULE PROCEDURE xt_idxsection_new_i4
    MODULE PROCEDURE xt_idxsection_new_i8
  END INTERFACE xt_idxsection_new

  INTERFACE
    FUNCTION xt_idxsection_new_c(start, num_dimensions, global_size, &
         local_size, local_start) BIND(c, name='xt_idxsection_new') &
         RESULT(idxsection)
      IMPORT :: c_int, c_ptr, xt_idxlist, xt_int_kind
      INTEGER(xt_int_kind), VALUE, INTENT(in) :: start
      INTEGER(c_int), VALUE, INTENT(in) :: num_dimensions
      INTEGER(xt_int_kind), INTENT(in) :: global_size(num_dimensions), &
           local_size(num_dimensions), local_start(num_dimensions)
      TYPE(c_ptr) :: idxsection
    END FUNCTION xt_idxsection_new_c
  END INTERFACE

CONTAINS

  FUNCTION xt_idxsection_new_a(start, global_size, local_size, local_start) &
       RESULT(idxsection)
    INTEGER(xt_int_kind), INTENT(in) :: start, global_size(:), local_size(:), &
         local_start(:)
    TYPE(xt_idxlist) :: idxsection
    INTEGER :: num_dimensions

    num_dimensions = SIZE(global_size)
    IF (SIZE(local_size) /= num_dimensions &
         .OR. SIZE(local_start) /= num_dimensions) &
         CALL xt_abort(xt_get_default_comm(), "non-matching array sizes", &
         __FILE__, &
         __LINE__)
    idxsection = xt_idxlist_c2f(&
         xt_idxsection_new_c(start, INT(num_dimensions, c_int), &
         &                   global_size, local_size, local_start))
  END FUNCTION xt_idxsection_new_a

  FUNCTION xt_idxsection_new_i2(start, num_dimensions, global_size, &
       local_size, local_start) RESULT(idxsection)
    INTEGER(i2), INTENT(in) :: num_dimensions
    INTEGER(xt_int_kind), INTENT(in) :: start, global_size(num_dimensions), &
         local_size(num_dimensions), local_start(num_dimensions)
    TYPE(xt_idxlist) :: idxsection
    INTEGER(c_int), PARAMETER :: dummy = 1

    idxsection = xt_idxlist_c2f(&
         xt_idxsection_new_c(start, INT(num_dimensions, c_int), &
         &                   global_size, local_size, local_start))
  END FUNCTION xt_idxsection_new_i2

  FUNCTION xt_idxsection_new_i4(start, num_dimensions, global_size, &
       local_size, local_start) RESULT(idxsection)
    INTEGER(i4), INTENT(in) :: num_dimensions
    INTEGER(xt_int_kind), INTENT(in) :: start, global_size(num_dimensions), &
         local_size(num_dimensions), local_start(num_dimensions)
    TYPE(xt_idxlist) :: idxsection
    INTEGER(c_int), PARAMETER :: dummy = 1

    IF (num_dimensions > HUGE(dummy)) &
         CALL xt_abort(xt_get_default_comm(), "num_dimensions too large", &
         __FILE__, &
         __LINE__)
    idxsection = xt_idxlist_c2f(&
         xt_idxsection_new_c(start, INT(num_dimensions, c_int), &
         &                   global_size, local_size, local_start))
  END FUNCTION xt_idxsection_new_i4

  FUNCTION xt_idxsection_new_i8(start, num_dimensions, global_size, &
       local_size, local_start) RESULT(idxsection)
    INTEGER(i8), INTENT(in) :: num_dimensions
    INTEGER(xt_int_kind), INTENT(in) :: start, global_size(num_dimensions), &
         local_size(num_dimensions), local_start(num_dimensions)
    TYPE(xt_idxlist) :: idxsection
    INTEGER(c_int), PARAMETER :: dummy = 1

    IF (num_dimensions > HUGE(dummy)) &
         CALL xt_abort(xt_get_default_comm(), "num_dimensions too large", &
         __FILE__, &
         __LINE__)
    idxsection = xt_idxlist_c2f(&
         xt_idxsection_new_c(start, INT(num_dimensions, c_int), &
         &                   global_size, local_size, local_start))
  END FUNCTION xt_idxsection_new_i8

END MODULE xt_idxsection
