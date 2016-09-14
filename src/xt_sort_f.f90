!>
!! @file xt_sort_f.f90
!! @brief Fortran interface to yaxt sort declarations
!!
!! @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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
MODULE xt_sort
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_size_t
  USE xt_core, ONLY: xt_abort, xt_int_kind
  IMPLICIT NONE
  PRIVATE

  INTERFACE xt_sort_int
    SUBROUTINE xt_sort_int_f2c(a, n) BIND(c, name='xt_sort_int_f2c')
      IMPORT :: c_int, c_size_t
      INTEGER(c_size_t), VALUE :: n
      INTEGER(c_int) :: a(n)
    END SUBROUTINE xt_sort_int_f2c
    MODULE PROCEDURE xt_sort_int_a
  END INTERFACE xt_sort_int
  PUBLIC :: xt_sort_int

  INTERFACE xt_sort_index
    SUBROUTINE xt_sort_index_f2c(a, n, positions, reset_positions) &
          BIND(c, name='xt_sort_index_f2c')
      IMPORT :: c_int
      INTEGER(c_int), VALUE :: n, reset_positions
      INTEGER(c_int) :: a(n), positions(n)
    END SUBROUTINE xt_sort_index_f2c
    MODULE PROCEDURE xt_sort_index_a_a_l
  END INTERFACE xt_sort_index
  PUBLIC :: xt_sort_index

  TYPE, BIND(c), PUBLIC :: xt_idxpos
    INTEGER(xt_int_kind) :: idx
    INTEGER(c_int) :: pos
  END TYPE xt_idxpos

  INTERFACE xt_sort_idxpos
    SUBROUTINE xt_sort_idxpos_f2c(a, n) BIND(c, name='xt_sort_idxpos_f2c')
      IMPORT :: c_size_t, xt_idxpos
      INTEGER(c_size_t), VALUE :: n
      TYPE(xt_idxpos) :: a(n)
    END SUBROUTINE xt_sort_idxpos_f2c
    MODULE PROCEDURE xt_sort_idxpos_a
  END INTERFACE xt_sort_idxpos
  PUBLIC :: xt_sort_idxpos

CONTAINS
  SUBROUTINE xt_sort_int_a(a)
    INTEGER(c_int), INTENT(inout) :: a(:)
    INTEGER :: a_size
    a_size = SIZE(a)
    IF (a_size > 1) CALL xt_sort_int_f2c(a, INT(a_size, c_size_t))
  END SUBROUTINE xt_sort_int_a

  SUBROUTINE xt_sort_index_a_a_l(a, positions, reset_positions)
    INTEGER(c_int), INTENT(inout) :: a(:), positions(:)
    LOGICAL, INTENT(in) :: reset_positions
    INTEGER :: a_size
    a_size = SIZE(a)
    IF (a_size > SIZE(positions)) &
         CALL xt_abort("positions array too small", &
         "xt_sort.f90", __LINE__)
    IF (a_size > 1) CALL xt_sort_index_f2c(a, INT(a_size, c_int), positions, &
         MERGE(1_c_int, 0_c_int, reset_positions))
  END SUBROUTINE xt_sort_index_a_a_l

  SUBROUTINE xt_sort_idxpos_a(a)
    TYPE(xt_idxpos), INTENT(inout) :: a(:)
    INTEGER :: a_size
    a_size = SIZE(a)
    IF (a_size > 1) CALL xt_sort_idxpos_f2c(a, INT(a_size, c_size_t))
  END SUBROUTINE xt_sort_idxpos_a
END MODULE xt_sort
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
