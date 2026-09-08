!>
!! @file test_idxempty_f.f90
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
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#include "fc_feature_defs.inc"
#ifdef VERBOSE
#define PUT_ERR(msg) WRITE (0, '(a)') msg ; err_count__ = err_count__ + 1
#else
#define PUT_ERR(msg) err_count__ = err_count__ + 1
#endif
PROGRAM test_idxempty
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_idxlist, xt_idxempty_new, &
       xt_int_kind, xt_stripe, xt_idxlist_get_intersection, &
       xt_idxlist_get_index_stripes, xt_idxlist_delete, &
       xt_bounds, xt_idxlist_get_bounding_box, &
       xt_idxlist_get_index_at_position, &
       xt_idxlist_get_position_of_index, &
       xt_idxlist_get_indices_at_positions, &
       xt_idxlist_get_positions_of_indices
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort, err_count__
  USE test_idxlist_utils, ONLY: check_idxlist, test_err_count, &
       idxlist_pack_unpack_copy
  IMPLICIT NONE

  TYPE(xt_idxlist) :: idxempty, idxempty_copy
  INTEGER(xt_int_kind) :: no_idx(1)
  TYPE(xt_stripe), ALLOCATABLE :: stripes(:)
  CHARACTER(len=*), PARAMETER :: filename = 'test_idxempty_f.f90'


  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  idxempty = xt_idxempty_new()

  CALL check_idxlist(idxempty, no_idx(1:0))

  idxempty_copy = idxlist_pack_unpack_copy(idxempty)

  ! check the computed intersection, must be identical to original list
  CALL check_idxlist(idxempty_copy, no_idx(1:0))

  CALL check_intersection

  CALL xt_idxlist_get_index_stripes(idxempty, stripes)

  IF (ALLOCATED(stripes)) &
       CALL test_abort("unexpected non-zero amount of stripes for &
       &empty index set", &
       filename, __LINE__)

  CALL check_bounding_box

  CALL check_idxpos_queries

  CALL xt_idxlist_delete(idxempty)
  CALL xt_idxlist_delete(idxempty_copy)

  CALL xt_finalize
  CALL finish_mpi

  IF (test_err_count() /= 0) CALL test_abort("non-zero error count", &
       filename, __LINE__)

CONTAINS

  SUBROUTINE check_intersection
    TYPE(xt_idxlist) :: intersection

    intersection = xt_idxlist_get_intersection(idxempty, idxempty_copy)
    CALL check_idxlist(intersection, no_idx(1:0))
    CALL xt_idxlist_delete(intersection)

  END SUBROUTINE check_intersection

  SUBROUTINE check_bounding_box
    INTEGER, PARAMETER :: ndims = 3
    INTEGER(xt_int_kind), PARAMETER :: global_start_index = 0
    INTEGER(xt_int_kind) :: global_size(ndims)
    TYPE(xt_bounds) :: bounds(ndims)

    global_size = 10
    bounds = xt_idxlist_get_bounding_box(idxempty, global_size, &
         global_start_index)
    IF (ANY(bounds%size /= 0)) &
         CALL test_abort("ERROR: non-zero boundings box for xt_idxempty in &
         &xt_idxlist_get_bounding_box", &
         filename, __LINE__)
  END SUBROUTINE check_bounding_box

  ! check that position and index queries return the expected
  ! (missing) values
  SUBROUTINE check_idxpos_queries
    INTEGER, PARAMETER :: num_pos = 3, undef_pat = 57005
    INTEGER(xt_int_kind) :: idx(num_pos)
    INTEGER(xt_int_kind), PARAMETER :: undef_idx = -5, &
         undef_idx2 = 2989
    INTEGER :: position(num_pos), rc, i
    LOGICAL :: found
    INTEGER, PARAMETER :: undef_pos = 48879

    idx(1) = undef_idx
    found = .NOT. xt_idxlist_get_index_at_position(idxempty, 0, idx(1))
    IF (found .OR. idx(1) /= undef_idx) THEN
      PUT_ERR("retrieved index position from empty list")
    END IF
    position = (/ -1, 0, 2 /)
    idx = 7
    rc = xt_idxlist_get_indices_at_positions(&
         idxempty, position, num_pos, idx, undef_idx)
    IF (rc /= num_pos) THEN
      PUT_ERR("retrieved index position from empty list")
    END IF
    IF (ANY(idx /= undef_idx)) THEN
      PUT_ERR("incorrectly initialized index from empty list")
    END IF
    position(1) = undef_pos
    found = .NOT. xt_idxlist_get_position_of_index(&
         idxempty, 0_xt_int_kind, position(1))
    IF (found .OR. position(1) /= -1) THEN
      PUT_ERR("incorrect position query result")
    END IF
    DO i = 1, num_pos
      idx(i) = IEOR(INT(i-1, xt_int_kind), undef_idx2)
      position(i) = -IEOR(i+1, undef_pos)
    END DO
    rc = xt_idxlist_get_positions_of_indices(&
         idxempty, idx, position, .FALSE.)
    IF (rc /= num_pos) THEN
      PUT_ERR("empty index list does not reject all indices")
      DO i = 1, num_pos
        IF (idx(i) /= IEOR(INT(i, xt_int_kind), undef_idx2) &
             .OR.  position(i) /= -IEOR(i+1, undef_pos)) THEN
          PUT_ERR("encountered unexpected result")
        END IF
      END DO
    END IF
  END SUBROUTINE check_idxpos_queries

END PROGRAM test_idxempty
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
