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
       xi => xt_int_kind
#if defined __PGI && __PGIC__ <= 15
  USE xt_redist_int_i4, ONLY: xt_redist_s_exchange
  USE xt_redist_int_i8, ONLY: xt_redist_s_exchange
#endif
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort, cmp_arrays
  IMPLICIT NONE
  INTEGER :: comm_size, comm_rank, ierror
  CHARACTER(len=*), PARAMETER :: filename = 'test_dmap_parallel_f.f90'
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

  CALL xt_finalize
  CALL finish_mpi

CONTAINS

  SUBROUTINE test1(comm, comm_rank, comm_size)
    INTEGER, INTENT(in) :: comm, comm_rank, comm_size
    INTEGER :: i, j, k
    TYPE(xt_idxlist) :: src_index_list, dst_index_list
    TYPE(xt_xmap) :: xmap
    INTEGER, PARAMETER :: nproma = 4, nlev = 20, nblk = 6
    INTEGER :: src_offsets(nproma*nblk)
    INTEGER :: lev_displacements(nlev)
    INTEGER :: src_data(nproma,nlev,nblk)

    TYPE(xt_dmap) :: dmap
    TYPE(xt_redist) :: redist
    INTEGER, ALLOCATABLE :: dst_data(:, :, :), dst_offsets(:), &
         ref_dst_data(:,:,:)


    ALLOCATE(dst_data(nproma, nlev, nblk*comm_size), &
         ref_dst_data(nproma, nlev, nblk*comm_size), &
         dst_offsets(nproma*nblk*comm_size))

    src_index_list = xt_idxvec_new((/ (INT(i + comm_rank * nproma * nblk, xi), &
         i = 1, nproma*nblk) /))
    dst_index_list = xt_idxvec_new((/ (INT(i, xi), i=1,nproma*nblk*comm_size)/))

    xmap = xt_xmap_all2all_new(src_index_list, dst_index_list, mpi_comm_world)

    DO i = 1, nblk
      DO j = 1, nproma
        src_offsets(j + (i-1) * nproma) = j - 1 + (i-1) * nproma * nlev
      END DO
    END DO
    DO i = 1, nblk*comm_size
      DO j = 1, nproma
        dst_offsets(j + (i-1)*nproma) = (j-1) + (i-1)*nproma*nlev
      END DO
    END DO

    DO i = 1, nlev
      lev_displacements(i) = (i-1)*nproma
    END DO

    DO i = 1, nlev
      DO j = 1, nblk
        DO k = 1, nproma
          src_data(k,i,j) = comm_rank * nproma * nblk &
               + (k-1) + (j-1)*nproma + (i-1)*comm_size*nproma*nblk
        END DO
      END DO
    END DO
    DO i = 1, nlev
      DO j = 1, nblk*comm_size
        DO k = 1, nproma
          ref_dst_data(k,i,j) = &
            (k-1) + (j-1)*nproma + (i-1)*comm_size*nproma*nblk
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
        CALL test_abort("error in xt_dmap", filename, __LINE__)

      CALL xt_redist_delete(redist)
      CALL xt_dmap_delete(dmap)

    END DO

    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(dst_index_list)
    CALL xt_idxlist_delete(src_index_list)

  END SUBROUTINE test1

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
