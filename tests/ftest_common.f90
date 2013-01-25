!>
!> @file ftest_common.f90
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
MODULE ftest_common
  USE mpi
  IMPLICIT NONE
  PRIVATE
  INTEGER, PUBLIC, PARAMETER :: dp = SELECTED_REAL_KIND(12, 307)

  TYPE timer
    CHARACTER(len=20) :: label = 'undef'
    INTEGER  :: istate  = -1
    REAL(dp) :: t0      = 0.0_dp
    REAL(dp) :: dt_work = 0.0_dp
  END TYPE timer

  INTERFACE test_abort
    MODULE PROCEDURE test_abort_cmsl_f
    MODULE PROCEDURE test_abort_msl_f
  END INTERFACE test_abort

  REAL(dp) :: sync_dt_sum = 0.0_dp
  LOGICAL, PARAMETER :: debug = .FALSE.
  LOGICAL :: verbose = .FALSE.

  PUBLIC :: init_mpi
  PUBLIC :: timer, treset, tstart, tstop, treport, sync
  PUBLIC :: id_map
  PUBLIC :: test_abort
CONTAINS
  SUBROUTINE init_mpi
    CHARACTER(len=*), PARAMETER :: context = 'init_mpi: '
    INTEGER :: ierror
    CALL mpi_init(ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_INIT failed', &
         __FILE__, __LINE__)
  END SUBROUTINE init_mpi

  PURE SUBROUTINE treset(t, label)
    TYPE(timer), INTENT(inout) :: t
    CHARACTER(len=*), INTENT(in) :: label
    t%label   = label
    t%istate  = 0
    t%t0      = 0.0_dp
    t%dt_work = 0.0_dp
  END SUBROUTINE treset

  SUBROUTINE tstart(t)
    TYPE(timer), INTENT(inout) :: t
    IF (debug) WRITE(0,*) 'tstart: ',t%label
    CALL sync
    t%istate = 1
    t%t0 = work_time()
  END SUBROUTINE tstart

  SUBROUTINE tstop(t)
    TYPE(timer), INTENT(inout) :: t
    REAL(dp) :: t1
    IF (debug) WRITE(0,*) 'tstop: ',t%label
    t1 = work_time()
    t%dt_work = t%dt_work + (t1 - t%t0)
    t%istate = 0
    CALL sync

  END SUBROUTINE tstop

  SUBROUTINE treport(t,extra_label,comm)
    TYPE(timer), INTENT(in) :: t
    CHARACTER(len=*), INTENT(in) :: extra_label
    INTEGER, INTENT(in) :: comm

    CHARACTER(len=*), PARAMETER :: context = 'treport: '
    REAL(dp) :: work_sum, work_max, work_avg, e
    REAL(dp) :: sbuf
    REAL(dp), ALLOCATABLE :: rbuf(:)
    INTEGER :: nprocs, rank, ierror

    sbuf = t%dt_work
    CALL mpi_comm_rank(comm, rank, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL test_abort(context//'MPI_COMM_RANK failed', &
         __FILE__, &
         __LINE__)
    CALL mpi_comm_size(comm, nprocs, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL test_abort(context//'MPI_COMM_RANK failed', &
         __FILE__, &
         __LINE__)
    ALLOCATE(rbuf(0:nprocs-1))
    rbuf = -1.0_dp
    CALL mpi_gather(sbuf, 1, MPI_DOUBLE_PRECISION, &
         &  rbuf, 1, MPI_DOUBLE_PRECISION, &
         &  0, comm, ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_GATHER failed', &
         __FILE__, &
         __LINE__)

    IF (rank == 0) THEN
      IF (rbuf(0) /= sbuf) CALL test_abort(context//'internal error (1)', &
           __FILE__, &
           __LINE__)
      IF (ANY(rbuf < 0.0_dp)) CALL test_abort(context//'internal error (2)', &
           __FILE__, &
           __LINE__)
      work_sum = SUM(rbuf)
      work_max = MAXVAL(rbuf)
      work_avg = work_sum / REAL(nprocs, dp)
      e = work_avg / (work_max + 1.e-20_dp)

      IF (verbose) WRITE(0,'(A,I4,2X,A16,3F18.8)') &
           'nprocs, label, wmax, wavg, e =', &
           nprocs, extra_label//':'//t%label, &
           work_max, work_avg, e
    ENDIF

  END SUBROUTINE treport

  SUBROUTINE sync
    CHARACTER(len=*), PARAMETER :: context = 'sync: '
    INTEGER :: ierror
    REAL(dp) :: t0, dt

    t0 = mpi_wtime()

    CALL mpi_barrier(MPI_COMM_WORLD, IERROR)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_BARRIER failed', &
         __FILE__, &
         __LINE__)

    dt = (MPI_WTIME() - t0)
    sync_dt_sum = sync_dt_sum + dt

  END SUBROUTINE sync

  REAL(dp) FUNCTION work_time()
    work_time = MPI_WTIME() - sync_dt_sum
    RETURN
  END FUNCTION work_time

  PURE SUBROUTINE id_map(map)
    INTEGER, INTENT(out) :: map(:,:)

    INTEGER :: i,j,p

    p = 0
    DO j = 1, SIZE(map,2)
      DO i = 1, SIZE(map,1)
        p = p + 1
        map(i,j) = p
      ENDDO
    ENDDO

  END SUBROUTINE id_map

  SUBROUTINE test_abort_msl_f(msg, source, line)
    CHARACTER(*), INTENT(in) :: msg
    CHARACTER(*), INTENT(in) :: source
    INTEGER, VALUE, INTENT(in) :: line
    CALL test_abort_cmsl_f(mpi_comm_world, msg, source, line)
  END SUBROUTINE test_abort_msl_f

  SUBROUTINE test_abort_cmsl_f(comm, msg, source, line)
    INTEGER, INTENT(in):: comm
    CHARACTER(*), INTENT(in) :: msg
    CHARACTER(*), INTENT(in) :: source
    INTEGER, VALUE, INTENT(in) :: line

    INTERFACE
      SUBROUTINE abort() BIND(c, name='abort')
      END SUBROUTINE abort
    END INTERFACE

    INTEGER :: ierror
    LOGICAL :: flag

    WRITE (0, '(3a,i0,2a)') 'Fatal error in ', source, ', line ', line, &
         ': ', msg
    FLUSH(0)
    CALL mpi_initialized(flag, ierror)
    IF (ierror == mpi_success .AND. flag) &
         CALL mpi_abort(comm, 1, ierror)
    CALL abort
  END SUBROUTINE test_abort_cmsl_f

END MODULE ftest_common
