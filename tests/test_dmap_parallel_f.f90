!>
!! @file test_dmap_parallel_f.f90
!!
!! @copyright Copyright  (C)  2018 Jörg Behrens <behrens@dkrz.de>
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
PROGRAM test_dmap_parallel_f
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, &
       xt_idxlist, xt_idxvec_new, xt_idxlist_delete, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_dmap, xt_dmap_offset_new, xt_dmap_delete, xt_dmap_repeat, &
       xt_dmap_reorder, xt_dmap_generate_redist, &
       XT_REORDER_NONE, XT_REORDER_SEND_UP, XT_REORDER_RECV_UP, &
       xt_redist, xt_redist_s_exchange, xt_redist_delete, &
       xi => xt_int_kind, &
       xt_sort_permutation, &
       xt_redist_get_send_mpi_datatype, xt_redist_get_recv_mpi_datatype
#if defined __PGI && __PGIC__ <= 15
  USE xt_redist_int_i4, ONLY: xt_redist_s_exchange
  USE xt_redist_int_i8, ONLY: xt_redist_s_exchange
#endif
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort, cmp_arrays
  IMPLICIT NONE
  INTEGER :: comm_size, comm_rank, ierror
  CHARACTER(len=*), PARAMETER :: filename = 'test_dmap_parallel_f.f90'
  INTEGER, PARAMETER :: num_indices_per_rank = 4
  INTEGER, PARAMETER :: reorder_type(3) &
       = (/xt_reorder_none, xt_reorder_send_up, xt_reorder_recv_up/)

  ! init mpi
  CALL init_mpi

  CALL xt_initialize(mpi_comm_world)
  CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
  IF (ierror /= MPI_SUCCESS) CALL test_abort('mpi_size: MPI_COMM_SIZE failed', &
       filename, __LINE__)
  CALL mpi_comm_rank(mpi_comm_world, comm_rank, ierror)
  IF (ierror /= MPI_SUCCESS) CALL test_abort('mpi_rank: MPI_COMM_RANK failed', &
       filename, __LINE__)

  CALL test1(mpi_comm_world, comm_rank, comm_size)
  CALL test2(mpi_comm_world, comm_rank, comm_size)

  CALL xt_finalize
  CALL finish_mpi

CONTAINS

  SUBROUTINE test1(comm, comm_rank, comm_size)
    INTEGER, INTENT(in) :: comm, comm_rank, comm_size
    INTEGER :: i
    TYPE(xt_idxlist) :: src_index_list, dst_index_list
    TYPE(xt_xmap) :: xmap

    INTEGER(xi) :: src_index_vector(num_indices_per_rank)
    INTEGER(xi), ALLOCATABLE :: dst_index_vector(:)

    ALLOCATE(dst_index_vector(comm_size * num_indices_per_rank))

    DO i = 1, num_indices_per_rank
      src_index_vector(i) = INT(num_indices_per_rank*comm_rank + i-1, xi)
    END DO
    DO i = 1, comm_size*num_indices_per_rank
      dst_index_vector(i) = INT(i-1, xi)
    END DO
    src_index_list = xt_idxvec_new(src_index_vector)
    dst_index_list = xt_idxvec_new(dst_index_vector)

    xmap = xt_xmap_all2all_new(src_index_list, dst_index_list, comm)

    DO i = 1, SIZE(reorder_type)
      CALL test_dmap(xmap, reorder_type(i), comm_size, comm_rank)
    END DO

    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(dst_index_list)
    CALL xt_idxlist_delete(src_index_list)

  END SUBROUTINE test1

  SUBROUTINE test_dmap(xmap, reorder_type, comm_size, comm_rank)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: reorder_type, comm_size, comm_rank
    INTEGER, PARAMETER :: num_ranks = 4

    INTEGER, PARAMETER :: src_offsets(num_indices_per_rank, 0:3) = &
         RESHAPE((/ 0,1,2,3, &
         3,2,1,0, &
         0,2,3,1, &
         3,1,2,0 /), (/ 4, 4 /)), &
         dst_offsets(num_ranks*num_indices_per_rank, 0:3) = &
         RESHAPE((/ 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15, &
         15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0, &
         10,15,13,4,3,1,8,14,12,9,0,2,6,7,5,11, &
         7,12,6,9,5,4,11,14,1,8,10,2,13,3,15,0 /), &
         (/ num_ranks*num_indices_per_rank, 4 /))
    TYPE(xt_dmap) :: dmap
    TYPE(xt_redist) :: redist
    INTEGER :: other_rank

    dmap = xt_dmap_offset_new(xmap, src_offsets(:, comm_rank), &
         dst_offsets(:, comm_rank))

    CALL xt_dmap_reorder(dmap, reorder_type)

    redist = xt_dmap_generate_redist(dmap, mpi_integer)

    ! test mpi datatypes in redist
    DO other_rank = 0, comm_size-1

      CALL test_send_dt(redist, reorder_type, &
           comm_size, comm_rank, other_rank, src_offsets, dst_offsets)
      CALL test_recv_dt(redist, reorder_type, &
           comm_size, comm_rank, other_rank, src_offsets, dst_offsets)
    END DO

    CALL xt_redist_delete(redist)
    CALL xt_dmap_delete(dmap)

  END SUBROUTINE test_dmap

  ! test send datatype (used to extract requested data from the local
  ! source data)
  SUBROUTINE test_send_dt(redist, reorder_type, &
       comm_size, comm_rank, other_rank, src_offsets, dst_offsets)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, INTENT(in) :: comm_size, comm_rank, other_rank
    INTEGER, INTENT(in) :: src_offsets(num_indices_per_rank, 0:comm_size-1), &
         dst_offsets(comm_size*num_indices_per_rank, 0:comm_size-1)
    INTEGER, INTENT(in) :: reorder_type
    !> offsets of the local points that need to be sent
    INTEGER :: local_send_offsets(num_indices_per_rank), &
         !> offsets on the remote process for receiving the data
         remote_recv_offsets(num_indices_per_rank), &
         !> local send data
         send_buffer(num_indices_per_rank), &
         !> contains data as it would arrive at the remote process
         recv_buffer(num_indices_per_rank), &
         !> expected contents of recv_buffer
         ref_recv_buffer(num_indices_per_rank)
    INTEGER :: i, ierror
    INTEGER :: send_datatype

    local_send_offsets = src_offsets(:, comm_rank)

    DO i = 1, num_indices_per_rank
      send_buffer(src_offsets(i, comm_rank)+1) &
           = comm_rank * num_indices_per_rank + i-1
    END DO
    DO i = 1, num_indices_per_rank
      remote_recv_offsets(i) &
           = dst_offsets(comm_rank * num_indices_per_rank + i, other_rank)
    END DO
    recv_buffer = -1
    DO i = 1, num_indices_per_rank
      ref_recv_buffer(i) = comm_rank * num_indices_per_rank + i-1
    END DO

    ! by default yaxt sends data in an order that is based on the sorted
    ! indices
    ! apply reordering to ref_recv_buffer (based on send or receive offsets)
    SELECT CASE (reorder_type)
    CASE (xt_reorder_send_up)
      CALL xt_sort_permutation(local_send_offsets, ref_recv_buffer)
    CASE (xt_reorder_recv_up)
      CALL xt_sort_permutation(remote_recv_offsets, ref_recv_buffer)
    END SELECT

    send_datatype = xt_redist_get_send_mpi_datatype(redist, other_rank)
    CALL mpi_type_commit(send_datatype, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("ERROR: committing send datatype failed", &
         filename, __LINE__)

    CALL mpi_sendrecv(&
         send_buffer, 1, send_datatype, 0, 0, &
         recv_buffer, num_indices_per_rank, mpi_integer, 0, 0, &
         mpi_comm_self, mpi_status_ignore, ierror)

    IF (ierror /= mpi_success .OR. ANY(ref_recv_buffer /= recv_buffer)) &
         CALL test_abort("ERROR: wrong send datatype", filename, __LINE__)

    CALL mpi_type_free(send_datatype, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("ERROR: freeing send datatype failed", &
         filename, __LINE__)

  END SUBROUTINE test_send_dt

  ! test recv datatype
  SUBROUTINE test_recv_dt(redist, reorder_type, &
       comm_size, comm_rank, other_rank, src_offsets, dst_offsets)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, INTENT(in) :: comm_size, comm_rank, other_rank
    INTEGER, INTENT(in) :: src_offsets(num_indices_per_rank, 0:comm_size-1), &
         dst_offsets(num_indices_per_rank*comm_size, 0:comm_size-1)
    INTEGER, INTENT(in) :: reorder_type

    !> offsets on the local process for receiving the data
    INTEGER :: local_recv_offsets(num_indices_per_rank), &
         !> offsets of the remote points that need to be sent
         remote_send_offsets(num_indices_per_rank), &
         !> contains data as it would arrive from the remote process
         send_buffer(num_indices_per_rank), &
         !> contains the destination array
         recv_buffer(comm_size*num_indices_per_rank), &
         !> expected contents of recv_buffer
         ref_recv_buffer(comm_size*num_indices_per_rank)
    INTEGER :: i
    INTEGER :: recv_datatype

    ref_recv_buffer = -1
    DO  i = 1, num_indices_per_rank
      local_recv_offsets(i) = &
           dst_offsets(other_rank * num_indices_per_rank + i, comm_rank)
      ref_recv_buffer(local_recv_offsets(i)+1) &
           = other_rank * num_indices_per_rank + i-1
    END DO
    remote_send_offsets = src_offsets(:, other_rank)
    DO i = 1, num_indices_per_rank
      send_buffer(i) = other_rank * num_indices_per_rank + i-1
    END DO
    recv_buffer = -1

    ! by default yaxt sends data in an order that is based on the sorted
    ! indices
    ! apply reordering to send_buffer (based on send or receive offsets)
    SELECT CASE (reorder_type)
    CASE (XT_REORDER_SEND_UP)
      CALL xt_sort_permutation(remote_send_offsets, send_buffer)
    CASE (xt_reorder_recv_up)
      CALL xt_sort_permutation(local_recv_offsets, send_buffer)
    END SELECT

    recv_datatype = xt_redist_get_recv_mpi_datatype(redist, other_rank)
    CALL mpi_type_commit(recv_datatype, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("ERROR: committing recv datatype failed", &
         filename, __LINE__)

    CALL mpi_sendrecv(&
         send_buffer, num_indices_per_rank, mpi_integer, 0, 0, &
         recv_buffer, 1, recv_datatype, 0, 0, &
         mpi_comm_self, mpi_status_ignore, ierror)

    IF (ierror /= mpi_success .OR. ANY(ref_recv_buffer /= recv_buffer)) &
         CALL test_abort("ERROR: wrong receive datatype", &
         filename, __LINE__)
    CALL mpi_type_free(recv_datatype, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("ERROR: freeing receive datatype failed", &
         filename, __LINE__)
  END SUBROUTINE test_recv_dt

  ! test setup:
  ! The source array of each process contains the indices in the range
  !   [rank*NPROMA*NBLK;(rank+1)*NPROMA*NBLK[
  ! The destination array contains all available indices, so that test generates
  ! a allgather communication pattern.
  ! To simulate the data layout of ICON, a level dimension is introduced.
  ! Processes always have the whole level column for each 2D point.
  ! The exchange map is generated for the 2D domain, the third dimension is
  ! through the use xt_dmap_repeat.
  SUBROUTINE test2(comm, comm_rank, comm_size)
    INTEGER, INTENT(in) :: comm, comm_rank, comm_size

    INTEGER, PARAMETER :: nproma = 4, nlev = 20, nblk = 6

    INTEGER(xi) :: src_index_vector(nproma*nblk)
    INTEGER :: src_offsets(NPROMA*NBLK)
    INTEGER(xi), ALLOCATABLE :: dst_index_vector(:)
    INTEGER, ALLOCATABLE :: dst_offsets(:)
    INTEGER :: lev_displacements(nlev)
    INTEGER :: src_data(nproma,nlev,nblk)
    INTEGER, ALLOCATABLE :: ref_dst_data(:,:,:), dst_data(:,:,:)
    TYPE(xt_dmap) :: dmap
    TYPE(xt_redist) :: redist
    INTEGER :: i, j, k
    TYPE(xt_idxlist) :: src_index_list, dst_index_list
    TYPE(xt_xmap) :: xmap

    ALLOCATE(dst_index_vector(nproma*nblk*comm_size), &
         dst_offsets(nproma*nblk*comm_size), &
         dst_data(nproma, nlev, nblk*comm_size), &
         ref_dst_data(nproma, nlev, nblk*comm_size))

    DO i = 1, nproma*nblk
      src_index_vector(i) = INT(nproma*nblk*comm_rank + i-1, xi)
    END DO
    DO i = 1, nproma*nblk*comm_size
      dst_index_vector(i) = INT(i-1, xi)
    END DO
    src_index_list = xt_idxvec_new(src_index_vector, nproma*nblk)
    dst_index_list = xt_idxvec_new(dst_index_vector, nproma*nblk*comm_size)

    xmap = xt_xmap_all2all_new(src_index_list, dst_index_list, comm)

    DO i = 1, nblk*nproma
      src_offsets(i) = MOD(i-1, nproma) + ((i-1)/nproma) * nproma * nlev
    END DO
    DO i = 1, nblk*comm_size*nproma
      dst_offsets(i) = MOD(i-1, nproma) + ((i-1)/nproma) * nproma * nlev
    END DO

    DO i = 1, nlev
      lev_displacements(i) = (i-1) * nproma
    END DO

    DO i = 1, nlev
      DO j = 1, nblk
        DO k = 1, nproma
          src_data(k,i,j) = comm_rank * nproma * nblk + k + j * nproma &
               + i * comm_size * nproma * nblk
        END DO
      END DO
      DO j = 1, nblk * comm_size
        DO k = 1, nproma
          ref_dst_data(k,i,j) = k + j * nproma + i * comm_size * nproma * nblk
        END DO
      END DO
    END DO


    DO i = 1, SIZE(reorder_type)
      dmap = xt_dmap_offset_new(xmap, src_offsets, dst_offsets)
      CALL xt_dmap_repeat(dmap, nlev, lev_displacements, lev_displacements)
      CALL xt_dmap_reorder(dmap, reorder_type(i))

      redist = xt_dmap_generate_redist(dmap, mpi_integer)
      dst_data = -1

      CALL xt_redist_s_exchange(redist, src_data, dst_data)

      IF (cmp_arrays(dst_data, ref_dst_data)) &
        CALL test_abort("ERROR: xt_dmap_repeat failed", filename, __LINE__)

      CALL xt_redist_delete(redist)
      CALL xt_dmap_delete(dmap)
    END DO

    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(dst_index_list)
    CALL xt_idxlist_delete(src_index_list)
  END SUBROUTINE test2


END PROGRAM test_dmap_parallel_f
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
