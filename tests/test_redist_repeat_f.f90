!>
!! @file test_redist_repeatection_static_f.f90
!! @brief Fortran test of redist_repeatection_static class
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
PROGRAM test_redist_repeat
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_xmap, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_p2p_off_new, &
       xt_redist_repeat_new, &
       xt_redist_delete, xt_redist_s_exchange1
  USE test_redist_common, ONLY: build_odd_selection_xmap, check_redist
  USE iso_c_binding, ONLY: c_loc, c_int
  IMPLICIT NONE
  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  CALL simple_test
  CALL test_repeated_redist
  CALL test_repeated_redist_with_gap
  CALL test_repeated_overlapping_redist

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE simple_test
    ! general test with one redist
    ! set up data
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_repeat
    INTEGER, PARAMETER :: src_slice_len = 5, dst_slice_len = 3
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(dst_slice_len) &
         = (/ 1.0d0, 3.0d0, 5.0d0 /), &
         src_data(src_slice_len) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /)
    DOUBLE PRECISION :: dst_data(dst_slice_len)
    INTEGER(mpi_address_kind) :: src_extent, dst_extent
    INTEGER(mpi_address_kind) :: base_address, temp_address
    INTEGER(c_int) :: displacements(1) = 0
    INTEGER :: ierror

    xmap = build_odd_selection_xmap(src_slice_len)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    CALL mpi_get_address(src_data(1), base_address, ierror)
    CALL mpi_get_address(src_data(2), temp_address, ierror)
    src_extent = (temp_address - base_address) * src_slice_len
    CALL mpi_get_address(dst_data(1), base_address, ierror)
    CALL mpi_get_address(dst_data(2), temp_address, ierror)
    dst_extent = (temp_address - base_address) * dst_slice_len

    ! generate redist_repeat
    redist_repeat = xt_redist_repeat_new(redist, src_extent, dst_extent, 1, &
                                         displacements)

    CALL xt_redist_delete(redist)

    ! test exchange
    CALL check_redist(redist_repeat, src_data, dst_data, ref_dst_data)

    ! clean up
    CALL xt_redist_delete(redist_repeat)
  END SUBROUTINE simple_test

  SUBROUTINE test_repeated_redist_ds1(redist_repeat)
    TYPE(xt_redist), INTENT(in) :: redist_repeat
    INTEGER :: i, j
    DOUBLE PRECISION, PARAMETER :: src_data(5, 3) = RESHAPE((/&
         (DBLE(i), i = 1, 15)/), (/ 5, 3 /))
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(3, 3) &
         = RESHAPE((/ ((DBLE(i + j), i = 1,5,2), j = 0,10,5) /), (/ 3, 3 /))
    DOUBLE PRECISION :: dst_data(3, 3)

    CALL check_redist(redist_repeat, src_data, dst_data, ref_dst_data)
  END SUBROUTINE test_repeated_redist_ds1

  SUBROUTINE test_repeated_redist_ds1_with_gap(redist_repeat)
    TYPE(xt_redist), INTENT(in) :: redist_repeat
    INTEGER :: i, j
    DOUBLE PRECISION, PARAMETER :: src_data(5, 5) = RESHAPE((/&
         (DBLE(i), i = 1, 25)/), (/ 5, 5 /))
    DOUBLE PRECISION :: dst_data(3, 5)
#ifndef __PGI
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(3, 5) &
         = RESHAPE((/ ((DBLE((i + j)*MOD(j+1,2)-MOD(j,2)), i = 1,5,2), &
         j = 0,20,5) /), (/ 3, 5 /))
#else
    DOUBLE PRECISION :: ref_dst_data(3, 5)
    ref_dst_data &
         = RESHAPE((/ ((DBLE((i + j)*MOD(j+1,2)-MOD(j,2)), i = 1,5,2), &
         j = 0,20,5) /), (/ 3, 5 /))
#endif
    CALL check_redist(redist_repeat, src_data, dst_data, ref_dst_data)
  END SUBROUTINE test_repeated_redist_ds1_with_gap

  SUBROUTINE test_repeated_redist_ds2(redist_repeat)
    TYPE(xt_redist), INTENT(in) :: redist_repeat
    INTEGER :: i, j
    DOUBLE PRECISION, PARAMETER :: src_data(5, 3) = RESHAPE((/&
         (DBLE(i), i = 20, 34)/), (/ 5, 3 /))
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(3, 3) &
         = RESHAPE((/ ((DBLE(i + j), i = 1,5,2), j = 19,33,5) /), (/ 3, 3 /))
    DOUBLE PRECISION :: dst_data(3, 3)

    CALL check_redist(redist_repeat, src_data, dst_data, ref_dst_data)
  END SUBROUTINE test_repeated_redist_ds2

  SUBROUTINE test_repeated_redist
    ! test with one redist used three times (with two different input data
    ! displacements -> test of cache) (with default cache size)
    ! set up data
    INTEGER, PARAMETER :: num_slice = 3
    INTEGER, PARAMETER :: src_slice_len = 5
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_repeat
    INTEGER(mpi_address_kind) :: src_extent, dst_extent
    INTEGER(mpi_address_kind) :: base_address, temp_address
    INTEGER(c_int) :: displacements(3)
    DOUBLE PRECISION, TARGET :: src_template(5, 3), dst_template(3, 3)
    INTEGER :: ierror

    xmap = build_odd_selection_xmap(src_slice_len)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! generate redist_repeat
    CALL mpi_get_address(src_template(1,1), base_address, ierror)
    CALL mpi_get_address(src_template(1,2), temp_address, ierror)
    src_extent = temp_address - base_address
    CALL mpi_get_address(dst_template(1,1), base_address, ierror)
    CALL mpi_get_address(dst_template(1,2), temp_address, ierror)
    dst_extent = temp_address - base_address
    displacements = (/0,1,2/)

    redist_repeat = xt_redist_repeat_new(redist, src_extent, dst_extent, &
         num_slice, displacements)
    CALL xt_redist_delete(redist)

    ! test exchange
    CALL test_repeated_redist_ds1(redist_repeat)
    ! test exchange
    CALL test_repeated_redist_ds2(redist_repeat)
    ! clean up
    CALL xt_redist_delete(redist_repeat)
  END SUBROUTINE test_repeated_redist

  SUBROUTINE test_repeated_redist_with_gap
    ! test with one redist used three times (with two different input data
    ! displacements -> test of cache) (with default cache size)
    ! set up data
    INTEGER, PARAMETER :: num_slice = 3
    INTEGER, PARAMETER :: src_slice_len = 5
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_repeat
    INTEGER(mpi_address_kind) :: src_extent, dst_extent
    INTEGER(mpi_address_kind) :: base_address, temp_address
    INTEGER(c_int), PARAMETER :: displacements(3) = (/0,2,4/)
    DOUBLE PRECISION, TARGET :: src_template(5, 3), dst_template(3, 3)
    INTEGER :: ierror

    xmap = build_odd_selection_xmap(src_slice_len)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! generate redist_repeat
    CALL mpi_get_address(src_template(1,1), base_address, ierror)
    CALL mpi_get_address(src_template(1,2), temp_address, ierror)
    src_extent = temp_address - base_address
    CALL mpi_get_address(dst_template(1,1), base_address, ierror)
    CALL mpi_get_address(dst_template(1,2), temp_address, ierror)
    dst_extent = temp_address - base_address

    redist_repeat = xt_redist_repeat_new(redist, src_extent, dst_extent, &
         num_slice, displacements)
    CALL xt_redist_delete(redist)

    ! test exchange
    CALL test_repeated_redist_ds1_with_gap(redist_repeat)
    ! clean up
    CALL xt_redist_delete(redist_repeat)
  END SUBROUTINE test_repeated_redist_with_gap

  SUBROUTINE test_repeated_overlapping_redist
    ! test with one redist used three times (with two different input data
    ! displacements -> test of cache) (with default cache size)
    ! set up data
    INTEGER, PARAMETER :: npt = 9, selection_len = 6
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_repeat
    INTEGER(mpi_address_kind) :: src_extent, dst_extent
    INTEGER(mpi_address_kind) :: base_address, temp_address
    INTEGER(c_int), PARAMETER :: displacements(2) = (/ 0, 1 /)
    INTEGER :: i, j, ierror
    INTEGER, PARAMETER :: src_pos(npt) = (/ (i, i=1,npt) /), &
         dst_pos(npt) = (/ (2*i, i = 0, npt-1) /)
    DOUBLE PRECISION, TARGET :: src_data(npt), dst_data(npt)
#if __INTEL_COMPILER >= 1600 && __INTEL_COMPILER <= 1602 || defined __PGI
    DOUBLE PRECISION :: ref_dst_data(npt)
#else
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(npt) &
         = (/ ((DBLE(((2-j)*3+i+101)*((ABS(j)+j)/ABS(j+1)) &
         &           +(j-1-ABS(j-1))/2), &
         &     i=1,3 ),j=2,0,-1) /)
#endif
    DOUBLE PRECISION, TARGET :: src_template(2), dst_template(2)

    xmap = build_odd_selection_xmap(selection_len)

    redist = xt_redist_p2p_off_new(xmap, src_pos, dst_pos, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! init data
#if __INTEL_COMPILER >= 1600 && __INTEL_COMPILER <= 1602 || defined __PGI
    DO j = 2, 0, -1
      DO i = 1, 3
        ref_dst_data(i + (2-j)*3) = DBLE(((2-j)*3+i+101)*((ABS(j)+j)/ABS(j+1)) &
             &                           +(j-1-ABS(j-1))/2)
      END DO
    END DO
#endif
    DO i = 1, npt
      src_data(i) = 1.0d2 + DBLE(i)
    END DO
    dst_data = -1.0d0

    ! test individual redists
    CALL redist_dbl(redist, src_data, dst_data)
    CALL redist_dbl(redist, src_data(2:), dst_data(2:))
    ! check individual redists to have desired effect
    IF (ANY(dst_data /= ref_dst_data)) &
         CALL test_abort("error in xt_redist_s_exchange1", &
         __FILE__, &
         __LINE__)
    dst_data = -1.0d0
    ! generate redist_repeat
    CALL mpi_get_address(src_template(1), base_address, ierror)
    CALL mpi_get_address(src_template(2), temp_address, ierror)
    src_extent = temp_address - base_address
    CALL mpi_get_address(dst_template(1), base_address, ierror)
    CALL mpi_get_address(dst_template(2), temp_address, ierror)
    dst_extent = temp_address - base_address

    redist_repeat = xt_redist_repeat_new(redist, src_extent, dst_extent, &
         displacements)
    CALL xt_redist_delete(redist)

    ! test exchange
    CALL check_redist(redist_repeat, src_data, SIZE(dst_data), &
         dst_data, ref_dst_data)
    ! clean up
    CALL xt_redist_delete(redist_repeat)
  END SUBROUTINE test_repeated_overlapping_redist

  ! work-around for F2003 restriction on slice c_loc
  SUBROUTINE redist_dbl(redist, src_data, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    DOUBLE PRECISION, TARGET, INTENT(in) :: src_data(*)
    DOUBLE PRECISION, TARGET, INTENT(inout) :: dst_data(*)
    CALL xt_redist_s_exchange1(redist, C_LOC(src_data), C_LOC(dst_data))
  END SUBROUTINE redist_dbl

END PROGRAM test_redist_repeat
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
