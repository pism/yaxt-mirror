!>
!> @file test_idxempty.f90
!>
!> @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
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
PROGRAM test_idxempty
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_idxlist, xt_idxempty_new, &
       xt_int_kind, xt_idxlist_get_pack_size, xt_stripe, &
       xt_idxlist_pack, xt_idxlist_unpack, xt_idxlist_get_intersection, &
       xt_idxlist_get_index_stripes, xt_idxlist_delete, &
       xt_bounds, xt_idxlist_get_bounding_box
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: check_idxlist, test_err_count
  IMPLICIT NONE

  TYPE(xt_idxlist) :: idxempty, idxempty_copy
  INTEGER(xt_int_kind) :: no_idx(1)
  INTEGER :: buffer_size, position, rank, request, ierror
  CHARACTER(1), ALLOCATABLE :: send_buffer(:), recv_buffer(:)
  TYPE(xt_stripe), ALLOCATABLE :: stripes(:)


  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  idxempty = xt_idxempty_new()

  CALL check_idxlist(idxempty, no_idx(1:0))

  buffer_size = xt_idxlist_get_pack_size(idxempty, mpi_comm_world)

  ALLOCATE(send_buffer(buffer_size), recv_buffer(buffer_size))

  position = 0
  CALL xt_idxlist_pack(idxempty, send_buffer, buffer_size, position, &
       mpi_comm_world)

  CALL mpi_comm_rank(mpi_comm_world, rank, ierror)
  IF (ierror /= mpi_success) &
       CALL test_abort("mpi_comm_rank failed", &
       __FILE__, &
       __LINE__)
  ! send the buffer
  CALL mpi_isend(send_buffer, buffer_size, mpi_packed, rank, 0, &
       mpi_comm_world, request, ierror)
  IF (ierror /= mpi_success) &
       CALL test_abort("mpi_isend failed", &
       __FILE__, &
       __LINE__)

  CALL mpi_request_free(request, ierror)
  IF (ierror /= mpi_success) &
       CALL test_abort("mpi_request_free failed", &
       __FILE__, &
       __LINE__)

  ! receive buffer
  CALL mpi_recv(recv_buffer, buffer_size, mpi_packed, rank, 0, &
       mpi_comm_world, mpi_status_ignore, ierror)
  IF (ierror /= mpi_success) &
       CALL test_abort("mpi_recv failed", &
       __FILE__, &
       __LINE__)

  position = 0
  CALL xt_idxlist_unpack(idxempty_copy, recv_buffer, buffer_size, &
       position, mpi_comm_world)

  ! check the computed intersection, must be identical to original list
  CALL check_idxlist(idxempty_copy, no_idx(1:0))

  CALL check_intersection

  CALL xt_idxlist_get_index_stripes(idxempty, stripes)

  IF (ALLOCATED(stripes)) &
       CALL test_abort("unexpected non-zero amount of stripes for &
       &empty index set", &
       __FILE__, &
       __LINE__)

  CALL check_bounding_box

  CALL xt_idxlist_delete(idxempty)
  CALL xt_idxlist_delete(idxempty_copy)

  CALL finish_mpi

  IF (test_err_count() /= 0) STOP 1

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
         __FILE__, &
         __LINE__)
  END SUBROUTINE check_bounding_box

END PROGRAM test_idxempty
