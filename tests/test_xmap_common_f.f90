!>
!! @file test_xmap_common_f.f90
!! @brief generic Fortran xmap test procedures
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
MODULE test_xmap_common
  USE iso_c_binding, ONLY: c_ptr, c_loc
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, &
       xt_idxlist, xt_idxlist_delete, xt_idxvec_new, &
       xt_xmap, xt_xmap_copy, xt_xmap_delete, &
       xt_xmap_get_num_destinations, xt_xmap_get_num_sources, &
       xt_xmap_get_destination_ranks, xt_xmap_get_source_ranks, &
       xt_mpi_comm_mark_exclusive
  IMPLICIT NONE
  PRIVATE
  INTEGER :: my_rank
  INTEGER, PARAMETER :: xi = xt_int_kind
  PUBLIC :: xmap_self_test_main, xmap_self_test_main_dble
CONTAINS
  SUBROUTINE xmap_self_test_main(xmap_new)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER :: ierror, i
    INTEGER :: comms(2)

    CALL init_mpi
    CALL xt_initialize(mpi_comm_world)
    CALL mpi_comm_rank(mpi_comm_world, my_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", &
         __FILE__, &
         __LINE__)

    comms(1) = mpi_comm_world
    CALL mpi_comm_dup(mpi_comm_world, comms(2), ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", &
         __FILE__, &
         __LINE__)
    CALL xt_mpi_comm_mark_exclusive(comms(2))

    DO i = 1, SIZE(comms)
      CALL test_xmap1(xmap_new, comms(i))
      CALL test_xmap2(xmap_new, comms(i))
    END DO

    CALL mpi_comm_free(comms(2), ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", &
         __FILE__, &
         __LINE__)

    IF (test_err_count() /= 0) &
         CALL test_abort("non-zero error count!", &
         __FILE__, &
         __LINE__)
    CALL xt_finalize
    CALL finish_mpi
  END SUBROUTINE xmap_self_test_main

  SUBROUTINE xmap_self_test_main_dble(xmap_new)
    INTERFACE
      FUNCTION xmap_new(src_ids, num_src_ids, dst_ids, num_dst_ids, id_size, &
                        comm) RESULT(res)
        IMPORT :: xt_xmap, c_ptr
        IMPLICIT NONE
        TYPE(c_ptr), INTENT(in) :: src_ids, dst_ids
        INTEGER, VALUE, INTENT(in) :: num_src_ids, num_dst_ids, id_size
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER :: ierror, i
    INTEGER :: comms(2)

    CALL init_mpi
    CALL xt_initialize(mpi_comm_world)
    CALL mpi_comm_rank(mpi_comm_world, my_rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", &
         __FILE__, &
         __LINE__)

    comms(1) = mpi_comm_world
    CALL mpi_comm_dup(mpi_comm_world, comms(2), ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", &
         __FILE__, &
         __LINE__)
    CALL xt_mpi_comm_mark_exclusive(comms(2))

    DO i = 1, SIZE(comms)
      CALL test_xmap1_dble(xmap_new, comms(i))
      CALL test_xmap2_dble(xmap_new, comms(i))
    END DO

    CALL mpi_comm_free(comms(2), ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("MPI error!", &
         __FILE__, &
         __LINE__)

    IF (test_err_count() /= 0) &
         CALL test_abort("non-zero error count!", &
         __FILE__, &
         __LINE__)
    CALL xt_finalize
    CALL finish_mpi
  END SUBROUTINE xmap_self_test_main_dble

  SUBROUTINE shift_idx(idx, offset)
    INTEGER(xt_int_kind), INTENT(inout) :: idx(:)
    INTEGER(xt_int_kind), INTENT(in) :: offset
    INTEGER :: i
    DO i = 1, SIZE(idx)
      idx(i) = idx(i) + INT(my_rank, xi) * offset
    END DO
  END SUBROUTINE shift_idx

  SUBROUTINE shift_ids_dble(ids, offset)
    DOUBLE PRECISION, INTENT(inout) :: ids(:)
    DOUBLE PRECISION, INTENT(in) :: offset
    INTEGER :: i
    DO i = 1, SIZE(ids)
      ids(i) = ids(i) + DBLE(my_rank) * offset
    END DO
  END SUBROUTINE shift_ids_dble

  SUBROUTINE assert_xmap_is_to_self(xmap)
    TYPE(xt_xmap) :: xmap
    INTEGER :: rank(1)
    IF (xt_xmap_get_num_destinations(xmap) /= 1) &
         CALL test_abort("error in xmap construction", &
         __FILE__, &
         __LINE__)

    IF (xt_xmap_get_num_sources(xmap) /= 1) &
         CALL test_abort("error in xt_xmap_get_num_sources", &
         __FILE__, &
         __LINE__)
    CALL xt_xmap_get_destination_ranks(xmap, rank)
    IF (rank(1) /= my_rank) &
         CALL test_abort("error in xt_xmap_get_destination_ranks", &
         __FILE__, &
         __LINE__)

    CALL xt_xmap_get_source_ranks(xmap, rank)
    IF (rank(1) /= my_rank) &
         CALL test_abort("error in xt_xmap_get_source_ranks", &
         __FILE__, &
         __LINE__)

  END SUBROUTINE assert_xmap_is_to_self

  SUBROUTINE test_xmap(src_index_list, dst_index_list, xmap_new, comm)
    INTEGER(xt_int_kind), INTENT(in) :: src_index_list(:), dst_index_list(:)
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER, INTENT(inout) :: comm

    TYPE(xt_xmap) :: xmap, xmap_copy
    src_idxlist = xt_idxvec_new(src_index_list)
    dst_idxlist = xt_idxvec_new(dst_index_list)

    xmap = xmap_new(src_idxlist, dst_idxlist, comm)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

    CALL assert_xmap_is_to_self(xmap)
    xmap_copy = xt_xmap_copy(xmap)
    CALL assert_xmap_is_to_self(xmap_copy)

    CALL xt_xmap_delete(xmap)
    CALL xt_xmap_delete(xmap_copy)
  END SUBROUTINE test_xmap

  SUBROUTINE test_xmap_dble(src_ids, num_src_ids, dst_ids, num_dst_ids, &
                            id_size, xmap_new, comm)
    DOUBLE PRECISION, INTENT(in) :: src_ids(:), dst_ids(:)
    INTEGER, INTENT(in) :: num_src_ids, num_dst_ids, id_size
    INTERFACE
      FUNCTION xmap_new(src_ids, num_src_ids, dst_ids, num_dst_ids, id_size, &
                        comm) RESULT(res)
        IMPORT :: xt_xmap, c_ptr
        IMPLICIT NONE
        TYPE(c_ptr), INTENT(in) :: src_ids, dst_ids
        INTEGER, VALUE, INTENT(in) :: num_src_ids, num_dst_ids, id_size
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER, INTENT(inout) :: comm

    TYPE(xt_xmap) :: xmap, xmap_copy

    xmap = xmap_new(c_loc(src_ids(1)), num_src_ids, c_loc(dst_ids(1)), &
                    num_dst_ids, id_size, comm)

    CALL assert_xmap_is_to_self(xmap)
    xmap_copy = xt_xmap_copy(xmap)
    CALL assert_xmap_is_to_self(xmap_copy)

    CALL xt_xmap_delete(xmap)
    CALL xt_xmap_delete(xmap_copy)
  END SUBROUTINE test_xmap_dble

  SUBROUTINE test_xmap1(xmap_new, comm)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER, INTENT(inout) :: comm

    INTEGER(xt_int_kind) :: i
    INTEGER(xt_int_kind), PARAMETER :: num_src_idx = 7, num_dst_idx = 7
    INTEGER(xt_int_kind) :: src_index_list(num_src_idx), &
         dst_index_list(num_dst_idx)
    DO i = 1_xi, num_src_idx
      src_index_list(i) = i
    END DO
    CALL shift_idx(src_index_list, num_src_idx)
    DO i = 1_xi, num_dst_idx
      dst_index_list(i) = num_dst_idx - i + 1_xi
    END DO
    CALL shift_idx(dst_index_list, num_src_idx)
    CALL test_xmap(src_index_list, dst_index_list, xmap_new, comm)
  END SUBROUTINE test_xmap1

  SUBROUTINE test_xmap1_dble(xmap_new, comm)
    INTERFACE
      FUNCTION xmap_new(src_ids, num_src_ids, dst_ids, num_dst_ids, id_size, &
                        comm) RESULT(res)
        IMPORT :: xt_xmap, c_ptr
        IMPLICIT NONE
        TYPE(c_ptr), INTENT(in) :: src_ids, dst_ids
        INTEGER, VALUE, INTENT(in) :: num_src_ids, num_dst_ids, id_size
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER, INTENT(inout) :: comm

    INTEGER :: i
    INTEGER, PARAMETER :: num_src_ids = 7, num_dst_ids = 7
    DOUBLE PRECISION :: src_ids(num_src_ids), dst_ids(num_dst_ids)
    DO i = 1, num_src_ids
      src_ids(i) = DBLE(i)
    END DO
    CALL shift_ids_dble(src_ids, DBLE(num_src_ids))
    DO i = 1, num_dst_ids
      dst_ids(i) = DBLE(num_dst_ids - i + 1)
    END DO
    CALL shift_ids_dble(dst_ids, DBLE(num_src_ids))
    CALL test_xmap_dble(src_ids, num_src_ids, dst_ids, num_dst_ids, 8, &
                        xmap_new, comm)
  END SUBROUTINE test_xmap1_dble

  SUBROUTINE test_xmap2(xmap_new, comm)
    INTERFACE
      FUNCTION xmap_new(src_idxlist, dst_idxlist, comm) RESULT(res)
        IMPORT :: xt_idxlist, xt_xmap
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER, INTENT(inout) :: comm

    INTEGER(xt_int_kind) :: src_index_list(14), dst_index_list(13)
    src_index_list = &
         (/ 5_xi, 67_xi, 4_xi, 5_xi, 13_xi, &
         &  9_xi,  2_xi, 1_xi, 0_xi, 96_xi, &
         & 13_xi, 12_xi, 1_xi, 3_xi /)
    dst_index_list = &
         (/ 5_xi, 4_xi, 3_xi, 96_xi, 1_xi, &
         &  5_xi, 4_xi, 5_xi,  4_xi, 3_xi, &
         & 13_xi, 2_xi, 1_xi /)
    CALL test_xmap(src_index_list, dst_index_list, xmap_new, comm)
  END SUBROUTINE test_xmap2

  SUBROUTINE test_xmap2_dble(xmap_new, comm)
    INTERFACE
      FUNCTION xmap_new(src_ids, num_src_ids, dst_ids, num_dst_ids, id_size, &
                        comm) RESULT(res)
        IMPORT :: xt_xmap, c_ptr
        IMPLICIT NONE
        TYPE(c_ptr), INTENT(in) :: src_ids, dst_ids
        INTEGER, VALUE, INTENT(in) :: num_src_ids, num_dst_ids, id_size
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(xt_xmap) :: res
      END FUNCTION xmap_new
    END INTERFACE
    INTEGER, INTENT(inout) :: comm

    DOUBLE PRECISION :: src_ids(14), dst_ids(13)
    src_ids = &
         (/ 5d+0, 67d+0, 4d+0, 5d+0, 13d+0, &
         &  9d+0,  2d+0, 1d+0, 0d+0, 96d+0, &
         & 13d+0, 12d+0, 1d+0, 3d+0 /)
    dst_ids = &
         (/ 5d+0, 4d+0, 3d+0, 96d+0, 1d+0, &
         &  5d+0, 4d+0, 5d+0,  4d+0, 3d+0, &
         & 13d+0, 2d+0, 1d+0 /)
    CALL test_xmap_dble(src_ids, 14, dst_ids, 13, 8, xmap_new, comm)
  END SUBROUTINE test_xmap2_dble

END MODULE test_xmap_common
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
