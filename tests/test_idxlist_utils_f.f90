!>
!> @file test_idxlist_utils_f.f90
!>
!> @copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
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
MODULE test_idxlist_utils
  USE yaxt, ONLY: xt_int_kind, xt_idxlist, xt_idxlist_f2c
  USE ftest_common, ONLY: test_abort
  USE iso_c_binding, ONLY: c_ptr
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: check_idxlist
CONTAINS
  SUBROUTINE check_idxlist(idxlist, ref_indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices(:)

    INTEGER :: num_ref_indices
    INTEGER(xt_int_kind) :: dummy(1)

    INTERFACE
      SUBROUTINE check_idxlist_c(idxlist, ref_indices, ref_num_indices) &
           BIND(c, name='check_idxlist')
        IMPORT :: xt_int_kind, c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(xt_int_kind), INTENT(in) :: ref_indices(*)
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: ref_num_indices
      END SUBROUTINE check_idxlist_c
    END INTERFACE

    num_ref_indices = SIZE(ref_indices)
    IF (num_ref_indices > 0) THEN
      CALL check_idxlist_c(xt_idxlist_f2c(idxlist), ref_indices, &
           INT(num_ref_indices, xt_int_kind))
    ELSE
      CALL check_idxlist_c(xt_idxlist_f2c(idxlist), dummy, &
           0_xt_int_kind)
    END IF
  END SUBROUTINE check_idxlist

END MODULE test_idxlist_utils

