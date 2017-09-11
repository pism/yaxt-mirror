!
! @file test_redist_single_array_base_parallel_f.f90
!
! @copyright Copyright  (C)  2017 Jörg Behrens <behrens@dkrz.de>
!                                 Moritz Hanke <hanke@dkrz.de>
!                                 Thomas Jahns <jahns@dkrz.de>
!
! @author Jörg Behrens <behrens@dkrz.de>
!         Moritz Hanke <hanke@dkrz.de>
!         Thomas Jahns <jahns@dkrz.de>
!
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
PROGRAM test_redist_single_array_base_parallel_f
  USE mpi
  USE yaxt, ONLY: xt_redist, xt_offset_ext, xt_initialize, xt_finalize, &
       xt_redist_copy, xt_redist_delete,  xt_redist_get_mpi_comm, &
       xt_redist_msg, xt_redist_single_array_base_new
  USE xt_core, ONLY: i2, i4, i8

  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_redist_common, ONLY: check_redist, communicators_are_congruent
  USE test_idxlist_utils, ONLY: test_err_count
  IMPLICIT NONE

  INTEGER :: comm_rank, comm_size, ierror

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  CALL mpi_comm_rank(mpi_comm_world, comm_rank, ierror)
  IF (ierror /= mpi_success) &
       CALL test_abort("MPI error!", &
       __FILE__, &
       __LINE__)

  CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
  IF (ierror /= mpi_success) &
       CALL test_abort("MPI error!", &
       __FILE__, &
       __LINE__)

  CALL test_round_robin
  CALL test_allgather
  CALL test_scatter

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)
  CALL xt_finalize
  CALL finish_mpi


CONTAINS

  SUBROUTINE test_round_robin

    TYPE(xt_redist_msg) :: send_msgs(1)
    TYPE(xt_redist_msg) :: recv_msgs(1)

    INTEGER, PARAMETER :: num_elem = 1
    DOUBLE PRECISION :: src_data(num_elem)
    DOUBLE PRECISION :: ref_dst_data(num_elem)

    send_msgs(1)%rank = MOD(comm_rank + 1, comm_size)
    send_msgs(1)%datatype = MPI_DOUBLE_PRECISION
    recv_msgs(1)%rank = MOD(comm_rank + comm_size - 1, comm_size)
    recv_msgs(1)%datatype = MPI_DOUBLE_PRECISION

    src_data(1) = comm_rank
    ref_dst_data(1) = MOD(comm_rank + comm_size - 1, comm_size)

    CALL test_single_array_base(send_msgs, recv_msgs, src_data, ref_dst_data)

  END SUBROUTINE test_round_robin

  SUBROUTINE test_allgather

    TYPE(xt_redist_msg) :: send_msgs(comm_size)
    TYPE(xt_redist_msg) :: recv_msgs(comm_size)

    DOUBLE PRECISION :: src_data(1)
    DOUBLE PRECISION :: ref_dst_data(comm_size)

    INTEGER :: i, ierr

    DO i = 1, comm_size
      send_msgs(i)%rank = i - 1
      send_msgs(i)%datatype = MPI_DOUBLE_PRECISION
      recv_msgs(i)%rank = i - 1
      CALL MPI_Type_create_indexed_block( &
        1, 1, (/i - 1/), MPI_DOUBLE_PRECISION, recv_msgs(i)%datatype, ierr)
      CALL MPI_Type_commit(recv_msgs(i)%datatype, ierr)
    END DO

    src_data(1) = comm_rank
    ref_dst_data(:) = (/ (i-1, i = 1, comm_size) /)

    CALL test_single_array_base(send_msgs, recv_msgs, src_data, ref_dst_data)

    DO i = 1, comm_size
      CALL MPI_Type_free(recv_msgs(i)%datatype, ierr)
    END DO

  END SUBROUTINE test_allgather

  SUBROUTINE test_scatter

    TYPE(xt_redist_msg), ALLOCATABLE :: send_msgs(:)
    TYPE(xt_redist_msg) :: recv_msgs(1)

    DOUBLE PRECISION, ALLOCATABLE :: src_data(:)
    DOUBLE PRECISION :: ref_dst_data(1)

    INTEGER :: i, ierr, nsend

    nsend = MERGE(comm_size, 0, comm_rank == 0)
    ALLOCATE(send_msgs(nsend))
    IF (comm_rank == 0) THEN
      DO i = 1, comm_size
        send_msgs(i)%rank = i - 1
        CALL MPI_Type_create_indexed_block( &
          1, 1, (/i - 1/), MPI_DOUBLE_PRECISION, send_msgs(i)%datatype, ierr)
        CALL MPI_Type_commit(send_msgs(i)%datatype, ierr)
      END DO
    END IF
    recv_msgs(1)%rank = 0
    recv_msgs(1)%datatype = MPI_DOUBLE_PRECISION

    ALLOCATE(src_data(nsend))
    DO i = 1, nsend
      src_data(i) = DBLE(i-1)
    END DO
    ref_dst_data(1) = comm_rank

    CALL test_single_array_base(send_msgs, recv_msgs, src_data, ref_dst_data)

    IF (comm_rank == 0) THEN
      DO i = 1, comm_size
        CALL MPI_Type_free(send_msgs(i)%datatype, ierr)
      END DO
    END IF

  END SUBROUTINE test_scatter

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

END PROGRAM test_redist_single_array_base_parallel_f
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
