!>
!! @file test_redist_common_f.f90
!! @brief common routines for Fortran test of redist classes
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
MODULE test_redist_common
  USE iso_c_binding, ONLY: c_loc
  USE mpi
  USE yaxt, ONLY: xt_idxlist, xt_int_kind, xt_idxvec_new, xt_idxlist_delete, &
       xt_xmap, xt_xmap_all2all_new, xt_redist, &
       xt_redist_s_exchange, xt_redist_s_exchange1, &
       xi => xt_int_kind
#ifdef __PGI
  ! PGI up to at least 15.4 has a bug that prevents proper import of
  ! multiply extended generics. This is a separate bug from the one exhibited
  ! in 12.7 and older (see test_xmap_intersection_parallel_f.f90 for that)
  USE xt_redist_real_dp, ONLY: xt_redist_s_exchange
  USE xt_redist_int_i2, ONLY: xt_redist_s_exchange
  USE xt_redist_int_i4, ONLY: xt_redist_s_exchange
  USE xt_redist_int_i8, ONLY: xt_redist_s_exchange
#endif
  USE ftest_common, ONLY: test_abort
  IMPLICIT NONE
  PRIVATE
  INTERFACE check_redist
    MODULE PROCEDURE check_redist_dp
    MODULE PROCEDURE check_redist_dp2
    MODULE PROCEDURE check_redist_dp2_2d
    MODULE PROCEDURE check_redist_xi
    MODULE PROCEDURE check_redist_xi2
  END INTERFACE check_redist
  PUBLIC :: build_odd_selection_xmap, check_redist, communicators_are_congruent
CONTAINS
  ! build xmap for destination list containing all odd elements of
  ! source list dimensioned 1 to src_slice_len
  FUNCTION build_odd_selection_xmap(src_slice_len) RESULT(xmap)
    INTEGER, INTENT(in) :: src_slice_len
    TYPE(xt_xmap) :: xmap
    INTEGER :: i, j, dst_slice_len
    INTEGER, PARAMETER :: dst_step = 2
    INTEGER(xt_int_kind), ALLOCATABLE :: index_list(:)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist

    dst_slice_len = (src_slice_len + dst_step - 1)/dst_step
    ALLOCATE(index_list(src_slice_len))
    DO i = 1, src_slice_len
      index_list(i) = INT(i, xt_int_kind)
    END DO
    src_idxlist = xt_idxvec_new(index_list)
    j = 1
    DO i = 1, src_slice_len, dst_step
      index_list(j) = INT(i, xt_int_kind)
      j = j + 1
    END DO
    dst_idxlist = xt_idxvec_new(index_list, dst_slice_len)
    DEALLOCATE(index_list)

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, mpi_comm_world)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)
  END FUNCTION build_odd_selection_xmap

  FUNCTION communicators_are_congruent(comm1, comm2) RESULT(congruent)
    INTEGER, INTENT(in) :: comm1, comm2
    LOGICAL :: congruent

    INTEGER :: ierror, rcode

    CALL mpi_comm_compare(comm1, comm2, rcode, ierror)
    congruent = ((rcode == mpi_ident) .OR. (rcode == mpi_congruent))
  END FUNCTION communicators_are_congruent

  SUBROUTINE check_redist_dp(redist, src, dst_size, dst, ref_dst)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, INTENT(in) :: dst_size
    DOUBLE PRECISION, TARGET, INTENT(in) :: src(*)
    DOUBLE PRECISION, INTENT(in) :: ref_dst(dst_size)
    DOUBLE PRECISION, TARGET, INTENT(inout) :: dst(dst_size)
    dst = -1.0d0
    CALL xt_redist_s_exchange1(redist, C_LOC(src), C_LOC(dst))
    IF (ANY(dst /= ref_dst)) &
         CALL test_abort("error in xt_redist_s_exchange1", &
         __FILE__, &
         __LINE__)
  END SUBROUTINE check_redist_dp

  SUBROUTINE check_redist_dp2(redist, src, dst, ref_dst)
    TYPE(xt_redist), INTENT(in) :: redist
    DOUBLE PRECISION, INTENT(in) :: src(:)
    DOUBLE PRECISION, INTENT(in) :: ref_dst(:)
    DOUBLE PRECISION, TARGET, INTENT(inout) :: dst(:)

    INTEGER :: dst_size, ref_dst_size
    dst_size = SIZE(dst)
    ref_dst_size = SIZE(ref_dst)
    IF (dst_size /= ref_dst_size) &
         CALL test_abort("error: ref_dst larger than dst", &
         __FILE__, &
         __LINE__)
    dst = -1.0d0
    CALL xt_redist_s_exchange(redist, src, dst)
    IF (ANY(dst /= ref_dst)) &
         CALL test_abort("error in xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)
  END SUBROUTINE check_redist_dp2

  SUBROUTINE check_redist_dp2_2d(redist, src, dst, ref_dst)
    TYPE(xt_redist), INTENT(in) :: redist
    DOUBLE PRECISION, INTENT(in) :: src(:,:)
    DOUBLE PRECISION, INTENT(in) :: ref_dst(:,:)
    DOUBLE PRECISION, TARGET, INTENT(inout) :: dst(:,:)

    INTEGER :: dst_size, ref_dst_size
    dst_size = SIZE(dst)
    ref_dst_size = SIZE(ref_dst)
    IF (dst_size /= ref_dst_size) &
         CALL test_abort("error: ref_dst larger than dst", &
         __FILE__, &
         __LINE__)
    dst = -1.0d0
    CALL xt_redist_s_exchange(redist, src, dst)
    IF (ANY(dst /= ref_dst)) &
         CALL test_abort("error in xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)
  END SUBROUTINE check_redist_dp2_2d

  SUBROUTINE check_redist_xi(redist, src, dst_size, dst, ref_dst)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, INTENT(in) :: dst_size
    INTEGER(xi), TARGET, INTENT(in) :: src(*)
    INTEGER(xi), INTENT(in) :: ref_dst(dst_size)
    INTEGER(xi), TARGET, INTENT(inout) :: dst(dst_size)
    dst = -1_xi
    CALL xt_redist_s_exchange1(redist, C_LOC(src), C_LOC(dst))
    IF (ANY(dst /= ref_dst)) &
         CALL test_abort("error in xt_redist_s_exchange1", &
         __FILE__, &
         __LINE__)
  END SUBROUTINE check_redist_xi

  SUBROUTINE check_redist_xi2(redist, src, dst, ref_dst)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(xi), INTENT(in) :: src(:)
    INTEGER(xi), INTENT(in) :: ref_dst(:)
    INTEGER(xi), TARGET, INTENT(inout) :: dst(:)

    INTEGER :: dst_size, ref_dst_size
    dst_size = SIZE(dst)
    ref_dst_size = SIZE(ref_dst)
    IF (dst_size /= ref_dst_size) &
         CALL test_abort("error: ref_dst larger than dst", &
         __FILE__, &
         __LINE__)
    CALL xt_redist_s_exchange(redist, src, dst)
    IF (ANY(dst /= ref_dst)) &
         CALL test_abort("error in xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)
  END SUBROUTINE check_redist_xi2

END MODULE test_redist_common
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
