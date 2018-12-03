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
#define NUM_RANKS (4)

PROGRAM test_dmap_parallel_f
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, &
       xt_idxlist, xt_idxvec_new, xt_idxlist_delete, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_dmap, xt_dmap_offset_new, xt_dmap_delete, xt_dmap_repeat, &
       xt_dmap_reorder, xt_dmap_generate_redist, &
       XT_REORDER_NONE, XT_REORDER_SEND_UP, XT_REORDER_RECV_UP, &
       xt_redist, xt_redist_s_exchange, xt_redist_delete
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort, cmp_arrays
  USE test_idxlist_utils, ONLY: test_err_count
  IMPLICIT NONE
  INTEGER :: comm_size
  CHARACTER(len=*), PARAMETER :: filename = &
__FILE__

  ! init mpi
  CALL init_mpi

  CALL xt_initialize(mpi_comm_world)

  IF (mpi_size(mpi_comm_world) /= NUM_RANKS) &
       CALL test_abort("non-zero error count!", &
       filename, __LINE__)

  CALL test1()

  CALL xt_finalize
  CALL finish_mpi

CONTAINS

  INTEGER FUNCTION mpi_rank(comm)
    INTEGER, INTENT(in) :: comm
    INTEGER :: ierror

    CALL mpi_comm_rank(comm, mpi_rank, ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort('mpi_rank: MPI_COMM_RANK failed', &
         __FILE__, &
         __LINE__)

  END FUNCTION mpi_rank

  INTEGER FUNCTION mpi_size(comm)
    INTEGER, INTENT(in) :: comm
    INTEGER :: ierror

    CALL mpi_comm_size(comm, mpi_size, ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort('mpi_size: MPI_COMM_SIZE failed', &
         __FILE__, &
         __LINE__)

  END FUNCTION mpi_size

#define NPROMA (4)
#define NLEV (20)
#define NBLK (6)

  SUBROUTINE test1()
    INTEGER :: i, j, k, l, rank
    TYPE(xt_idxlist) :: src_index_list, dst_index_list
    TYPE(xt_xmap) :: xmap
    INTEGER :: src_offsets(NPROMA*NBLK)
    INTEGER :: dst_offsets(NPROMA*NBLK*NUM_RANKS)
    INTEGER :: lev_displacements(NLEV)
    INTEGER :: src_data(NPROMA,NLEV,NBLK)
    INTEGER :: ref_dst_data(NPROMA,NLEV,NBLK*NUM_RANKS)
    INTEGER :: reorder_type(3)
    TYPE(xt_dmap) :: dmap
    TYPE(xt_redist) :: redist
    INTEGER :: dst_data(NPROMA, NLEV, NBLK*NUM_RANKS)

    rank = mpi_rank(mpi_comm_world)

    src_index_list = xt_idxvec_new((/(i+rank*NPROMA*NBLK,i=1,NPROMA*NBLK)/))
    dst_index_list = xt_idxvec_new((/(i,i=1,NPROMA*NBLK*NUM_RANKS)/))

    xmap = xt_xmap_all2all_new(src_index_list, dst_index_list, mpi_comm_world)

    k = 1
    DO i = 1, NBLK
      DO j = 1, NPROMA
        src_offsets(k) = (j-1) + (i-1) * NPROMA * NLEV
        k = k + 1
      END DO
    END DO
    k = 1
    DO i = 1, NBLK*NUM_RANKS
      DO j = 1, NPROMA
        dst_offsets(k) = (j-1) + (i-1) * NPROMA * NLEV
        k = k + 1
      END DO
    END DO

    lev_displacements = (/((i-1)*NPROMA,i=1,NLEV)/)

    DO i = 0, NLEV-1
      DO j = 0, NBLK-1
        DO k = 0, NPROMA-1
          src_data(k+1,i+1,j+1) = rank * NPROMA * NBLK + k + j * NPROMA + &
                                  i * NUM_RANKS * NPROMA * NBLK
        END DO
      END DO
    END DO
    DO i = 0, NLEV-1
      DO j = 0, NBLK*NUM_RANKS-1
        DO k = 0, NPROMA-1
          ref_dst_data(k+1,i+1,j+1) = &
            k + j * NPROMA + i * NUM_RANKS * NPROMA * NBLK
        END DO
      END DO
    END DO

    reorder_type = (/XT_REORDER_NONE, XT_REORDER_SEND_UP, XT_REORDER_RECV_UP/)

    DO i = 1, SIZE(reorder_type)

      dmap = xt_dmap_offset_new(xmap, src_offsets, dst_offsets)
      CALL xt_dmap_repeat(dmap, NLEV, lev_displacements, lev_displacements)
      CALL xt_dmap_reorder(dmap, reorder_type(i))

      redist = xt_dmap_generate_redist(dmap, MPI_INTEGER)

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
