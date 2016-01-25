!>
!! @file test_redist_collection_displace.f90
!! @brief Fortran cache displacement test of redist_collection class
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
MODULE redist_collection_displace
  USE mpi
  USE ftest_common, ONLY: test_abort
  USE yaxt, ONLY: xt_int_kind, xi => xt_int_kind, &
       xt_xmap, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_collection_new, &
       Xt_redist_copy, xt_redist_delete, xt_redist_s_exchange, &
       xt_slice_c_loc
  USE test_redist_common, ONLY: build_odd_selection_xmap
  USE iso_c_binding, ONLY: c_ptr
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: cache_size = 16, cache_overrun = 2
    INTEGER(xt_int_kind), PARAMETER :: num_slice = 3_xi, dst_step = 2_xi
    INTEGER, PARAMETER :: src_slice_len = 5
    INTEGER, PARAMETER :: dst_slice_len &
         = (src_slice_len + dst_step - 1)/dst_step
  PUBLIC :: test_displacement_variations
CONTAINS
  ! test with one redist used three times (with different input
  ! data displacements until the cache is full)
  ! set up data
  SUBROUTINE test_displacement_variations
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redists(num_slice), redist_coll, &
         redist_coll_copy

    xmap = build_odd_selection_xmap(src_slice_len)
    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! generate redist_collection
    redists = redist

    redist_coll = xt_redist_collection_new(redists, INT(num_slice), &
         cache_size, mpi_comm_world)

    CALL xt_redist_delete(redist)

    CALL run_displacement_check(redist_coll)
    redist_coll_copy = xt_redist_copy(redist_coll)
    CALL run_displacement_check(redist_coll_copy)

    ! clean up
    CALL xt_redist_delete(redist_coll)
    CALL xt_redist_delete(redist_coll_copy)
  END SUBROUTINE test_displacement_variations

  SUBROUTINE run_displacement_check(redist_coll)
    TYPE(xt_redist), INTENT(in) :: redist_coll
    INTEGER(xt_int_kind) :: i, j
    INTEGER :: k
    DOUBLE PRECISION, TARGET, SAVE :: src_data(src_slice_len, num_slice) &
         = RESHAPE((/ (DBLE(i), i = 1_xi, src_slice_len*num_slice) /), &
         (/ INT(src_slice_len), INT(num_slice) /))
    DOUBLE PRECISION, TARGET :: dst_data(dst_slice_len, num_slice)
    DOUBLE PRECISION, TARGET :: &
         src_data_(src_slice_len + cache_size + cache_overrun), &
         dst_data_(dst_slice_len + cache_size + cache_overrun)
    TYPE(c_ptr) :: src_data_p(num_slice), dst_data_p(num_slice)
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(dst_slice_len, num_slice) = &
         RESHAPE((/ ((DBLE(i + j * src_slice_len), &
         &            i = 1_xi, src_slice_len, dst_step), &
         &           j = 0_xi, num_slice - 1_xi) /), &
         &       (/ INT(dst_slice_len), INT(num_slice) /))

    dst_data = -1.0d0

    DO i = 1, num_slice - 1
      CALL xt_slice_c_loc(src_data(:, i), src_data_p(i))
      CALL xt_slice_c_loc(dst_data(:, i), dst_data_p(i))
    END DO

    ! test exchange
    DO k = 1, cache_size + cache_overrun
      src_data_(k:k+src_slice_len-1) = src_data(:,num_slice)
      dst_data_(k:k+dst_slice_len-1) = -1.0d0
      dst_data = -1.0d0

      CALL xt_slice_c_loc(src_data_(k:k+src_slice_len-1), src_data_p(3))
      CALL xt_slice_c_loc(dst_data_(k:k+dst_slice_len-1), dst_data_p(3))

      CALL xt_redist_s_exchange(redist_coll, INT(num_slice), src_data_p, &
           dst_data_p)

      IF (ANY(ref_dst_data(:, 1:num_slice-1) /= dst_data(:, 1:2)) &
           .OR. ANY(ref_dst_data(:,3) /= dst_data_(k:k+dst_slice_len-1))) &
           CALL test_abort("error in xt_redist_s_exchange", &
           __FILE__, &
           __LINE__)
    END DO
  END SUBROUTINE run_displacement_check

END MODULE redist_collection_displace
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
