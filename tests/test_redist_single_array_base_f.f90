!>
!! @file test_redist_single_array_base_f.f90
!!
!! @copyright Copyright  (C)  2017 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!
!
!  Keywords:
!  Maintainer: Jörg Behrens <behrens@dkrz.de>
!              Moritz Hanke <hanke@dkrz.de>
!              Thomas Jahns <jahns@dkrz.de>
!  URL: https://doc.redmine.dkrz.de/yaxt/html/
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are  permitted provided that the following conditions are
!  met:
!
!  Redistributions of source code must retain the above copyright notice,
!  this list of conditions and the following disclaimer.
!
!  Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer in the
!  documentation and/or other materials provided with the distribution.
!
!  Neither the name of the DKRZ GmbH nor the names of its contributors
!  may be used to endorse or promote products derived from this software
!  without specific prior written permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
PROGRAM test_redist_single_array_base_f
  USE mpi
  USE yaxt, ONLY: xt_redist, xt_offset_ext, xt_initialize, xt_finalize, &
       xt_redist_copy, xt_redist_delete,  xt_redist_get_mpi_comm, &
       xt_redist_msg, xt_redist_single_array_base_new
  USE xt_core, ONLY: i2, i4, i8

  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_redist_common, ONLY: check_redist, communicators_are_congruent
  USE test_idxlist_utils, ONLY: test_err_count
  IMPLICIT NONE

  ! init mpi
  CALL init_mpi

  CALL xt_initialize(mpi_comm_world)

  ! single double
  call test_single_double
  ! reverse order of some doubles
  call test_reverse_doubles

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)

  CALL xt_finalize
  CALL finish_mpi

CONTAINS

  SUBROUTINE test_single_double

    TYPE(xt_redist_msg) :: send_msgs(1)
    TYPE(xt_redist_msg) :: recv_msgs(1)

    INTEGER, PARAMETER :: num_elem = 1
    DOUBLE PRECISION, PARAMETER :: src_data(num_elem) &
         = (/ 0.0d0 /)
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(num_elem) &
         = (/ 0.0d0 /)

    send_msgs(1)%rank = 0
    send_msgs(1)%datatype = MPI_DOUBLE_PRECISION
    recv_msgs(1)%rank = 0
    recv_msgs(1)%datatype = MPI_DOUBLE_PRECISION

    CALL test_single_array_base(send_msgs, recv_msgs, src_data, ref_dst_data)

  END SUBROUTINE test_single_double

  SUBROUTINE test_reverse_doubles

    TYPE(xt_redist_msg) :: send_msgs(1)
    TYPE(xt_redist_msg) :: recv_msgs(1)

    INTEGER :: i, ierr
    INTEGER, PARAMETER :: num_elem = 10
    DOUBLE PRECISION, PARAMETER :: src_data(num_elem) &
         = (/ (DBLE(i), i = 1, num_elem) /)
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(num_elem) &
         = (/ (DBLE(i), i = num_elem, 1, -1) /)

    send_msgs(1)%rank = 0
    CALL MPI_Type_contiguous( &
      num_elem, MPI_DOUBLE_PRECISION, send_msgs(1)%datatype, ierr)
    CALL MPI_Type_commit(send_msgs(1)%datatype, ierr)
    recv_msgs(1)%rank = 0
    CALL MPI_Type_create_indexed_block( &
      num_elem, 1, (/ (i, i = num_elem - 1, 0, -1) /), MPI_DOUBLE_PRECISION, &
      recv_msgs(1)%datatype, ierr);
    CALL MPI_Type_commit(recv_msgs(1)%datatype, ierr)

    CALL test_single_array_base(send_msgs, recv_msgs, src_data, ref_dst_data)

    CALL MPI_Type_free(recv_msgs(1)%datatype, ierr)
    CALL MPI_Type_free(send_msgs(1)%datatype, ierr)

  END SUBROUTINE test_reverse_doubles

  SUBROUTINE test_single_array_base( &
      send_msgs, recv_msgs, src_data, ref_dst_data)
    TYPE(xt_redist_msg), INTENT(IN) :: send_msgs(:)
    TYPE(xt_redist_msg), INTENT(IN) :: recv_msgs(:)
    DOUBLE PRECISION, INTENT(IN) :: src_data(:)
    DOUBLE PRECISION, INTENT(IN) :: ref_dst_data(:)

    TYPE(xt_redist) :: redist

    redist = &
      xt_redist_single_array_base_new(send_msgs, recv_msgs, MPI_COMM_WORLD)
    call check_redist_extended(redist, src_data, ref_dst_data)

    redist = &
      xt_redist_single_array_base_new( &
        INT(size(send_msgs), i2), INT(size(recv_msgs), i2), &
        send_msgs, recv_msgs, MPI_COMM_WORLD)
    call check_redist_extended(redist, src_data, ref_dst_data)

    redist = &
      xt_redist_single_array_base_new( &
        INT(size(send_msgs), i4), INT(size(recv_msgs), i4), &
        send_msgs, recv_msgs, MPI_COMM_WORLD)
    call check_redist_extended(redist, src_data, ref_dst_data)

    redist = &
      xt_redist_single_array_base_new( &
        INT(size(send_msgs), i8), INT(size(recv_msgs), i8), &
        send_msgs, recv_msgs, MPI_COMM_WORLD)
    call check_redist_extended(redist, src_data, ref_dst_data)

  END SUBROUTINE test_single_array_base

  SUBROUTINE check_redist_extended(redist, src_data, ref_dst_data)
    TYPE(xt_redist), INTENT(INOUT) :: redist
    DOUBLE PRECISION, INTENT(IN) :: src_data(:)
    DOUBLE PRECISION, INTENT(IN) :: ref_dst_data(:)

    DOUBLE PRECISION :: dst_data(SIZE(ref_dst_data))

    TYPE(xt_redist) :: redist_copy

    ! test communicator of redist
    IF (.NOT. communicators_are_congruent(xt_redist_get_mpi_comm(redist), &
         MPI_COMM_WORLD)) &
         CALL test_abort("error in xt_redist_get_mpi_comm", &
         __FILE__, &
         __LINE__)

    ! test exchange
    CALL check_redist(redist, src_data, dst_data, ref_dst_data)
    redist_copy = xt_redist_copy(redist)
    CALL xt_redist_delete(redist)
    CALL check_redist(redist_copy, src_data, dst_data, ref_dst_data)
    CALL xt_redist_delete(redist_copy)

  END SUBROUTINE check_redist_extended

END PROGRAM test_redist_single_array_base_f
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
