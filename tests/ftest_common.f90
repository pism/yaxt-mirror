!>
!> @file ftest_common.f90
!>
!> @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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
MODULE ftest_common
  USE mpi
  USE xt_core, ONLY: i2, i4, i8
  USE iso_c_binding, ONLY: c_int, c_bool, c_double, &
       c_int16_t, c_int32_t, c_int64_t
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

  INTERFACE icmp
    MODULE PROCEDURE icmp_2d
    MODULE PROCEDURE icmp_3d
  END INTERFACE icmp

  INTERFACE cmp_arrays

    FUNCTION cmp_dbl_arrays_(asize, a, b) BIND(c, name='cmp_dbl_arrays')
      IMPORT :: c_bool, c_int, c_double
#ifdef HAVE_FC_LOGICAL_C_BOOL
      LOGICAL(c_bool) :: cmp_dbl_arrays_
#define cmp_dbl_arrays(asize,a,b) LOGICAL(cmp_dbl_arrays_(asize,a,b))
#else
      INTEGER(c_int) :: cmp_dbl_arrays_
#define cmp_dbl_arrays(asize,a,b) (cmp_dbl_arrays_(asize,a,b) /= 0)
#endif
      INTEGER(c_int), VALUE, INTENT(in) :: asize
      REAL(c_double), INTENT(in) :: a(asize), b(asize)
    END FUNCTION cmp_dbl_arrays_

    FUNCTION cmp_int16_arrays_(asize, a, b) BIND(c, name='cmp_int16_arrays')
      IMPORT :: c_bool, c_int, c_int16_t
#ifdef HAVE_FC_LOGICAL_C_BOOL
      LOGICAL(c_bool) :: cmp_int16_arrays_
#define cmp_int16_arrays(asize,a,b) LOGICAL(cmp_int16_arrays_(asize,a,b))
#else
      INTEGER(c_int) :: cmp_int16_arrays_
#define cmp_int16_arrays(asize,a,b) (cmp_int16_arrays_(asize,a,b) /= 0)
#endif
      INTEGER(c_int), VALUE, INTENT(in) :: asize
      INTEGER(c_int16_t), INTENT(in) :: a(asize), b(asize)
    END FUNCTION cmp_int16_arrays_

    FUNCTION cmp_int32_arrays_(asize, a, b) BIND(c, name='cmp_int32_arrays')
      IMPORT :: c_bool, c_int, c_int32_t
#ifdef HAVE_FC_LOGICAL_C_BOOL
      LOGICAL(c_bool) :: cmp_int32_arrays_
#define cmp_int32_arrays(asize,a,b) LOGICAL(cmp_int32_arrays_(asize,a,b))
#else
      INTEGER(c_int) :: cmp_int32_arrays_
#define cmp_int32_arrays(asize,a,b) (cmp_int32_arrays_(asize,a,b) /= 0)
#endif
      INTEGER(c_int), VALUE, INTENT(in) :: asize
      INTEGER(c_int32_t), INTENT(in) :: a(asize), b(asize)
    END FUNCTION cmp_int32_arrays_

    FUNCTION cmp_int64_arrays_(asize, a, b) BIND(c, name='cmp_int64_arrays')
      IMPORT :: c_bool, c_int, c_int64_t
#ifdef HAVE_FC_LOGICAL_C_BOOL
      LOGICAL(c_bool) :: cmp_int64_arrays_
#define cmp_int64_arrays(asize,a,b) LOGICAL(cmp_int64_arrays_(asize,a,b))
#else
      INTEGER(c_int) :: cmp_int64_arrays_
#define cmp_int64_arrays(asize,a,b) (cmp_int64_arrays_(asize,a,b) /= 0)
#endif
      INTEGER(c_int), VALUE, INTENT(in) :: asize
      INTEGER(c_int64_t), INTENT(in) :: a(asize), b(asize)
    END FUNCTION cmp_int64_arrays_

    MODULE PROCEDURE cmp_dbl_arrays_a1d_a1d
    MODULE PROCEDURE cmp_dbl_arrays_a2d_a2d
    MODULE PROCEDURE cmp_dbl_arrays_a3d_a3d
    MODULE PROCEDURE cmp_i2_arrays_a1d_a1d
    MODULE PROCEDURE cmp_i4_arrays_a1d_a1d
    MODULE PROCEDURE cmp_i8_arrays_a1d_a1d
    MODULE PROCEDURE cmp_i2_arrays_a2d_a2d
    MODULE PROCEDURE cmp_i4_arrays_a2d_a2d
    MODULE PROCEDURE cmp_i8_arrays_a2d_a2d
    MODULE PROCEDURE cmp_i2_arrays_a3d_a3d
    MODULE PROCEDURE cmp_i4_arrays_a3d_a3d
    MODULE PROCEDURE cmp_i8_arrays_a3d_a3d
  END INTERFACE cmp_arrays

  INTERFACE id_map
    MODULE PROCEDURE id_map_i2, id_map_i4, id_map_i8
  END INTERFACE id_map

  INTERFACE icbrt
    MODULE PROCEDURE icbrt_i2, icbrt_i4, icbrt_i8
  END INTERFACE icbrt

  REAL(dp) :: sync_dt_sum = 0.0_dp
  LOGICAL, PARAMETER :: debug = .FALSE.
  LOGICAL :: verbose = .FALSE.

  PUBLIC :: init_mpi, finish_mpi
  PUBLIC :: timer, treset, tstart, tstop, treport, mysync
  PUBLIC :: id_map, icbrt, icmp, factorize, regular_deco
  PUBLIC :: test_abort, set_verbose, get_verbose
  PUBLIC :: cmp_arrays
  PUBLIC :: run_randomized_tests, init_fortran_random
  CHARACTER(len=*), PARAMETER :: filename = 'ftest_common.f90'
CONTAINS

  SUBROUTINE init_mpi
    CHARACTER(len=*), PARAMETER :: context = 'init_mpi: '
    INTEGER :: ierror
    CALL mpi_init(ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_INIT failed', &
         filename, __LINE__)
  END SUBROUTINE init_mpi

  SUBROUTINE finish_mpi
    CHARACTER(len=*), PARAMETER :: context = 'finish_mpi: '
    INTEGER :: ierror
    CALL MPI_FINALIZE(ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_FINALIZE failed', &
         filename, &
         __LINE__)
  END SUBROUTINE finish_mpi

  SUBROUTINE set_verbose(verb)
    LOGICAL, INTENT(in) :: verb
    verbose = verb
  END SUBROUTINE set_verbose

  SUBROUTINE get_verbose(verb)
    LOGICAL, INTENT(out) :: verb
    verb = verbose
  END SUBROUTINE get_verbose

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
    CALL mysync
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
    CALL mysync

  END SUBROUTINE tstop

  SUBROUTINE treport(t,extra_label,comm)
    TYPE(timer), INTENT(in) :: t
    CHARACTER(len=*), INTENT(in) :: extra_label
    INTEGER, INTENT(in) :: comm

    CHARACTER(len=*), PARAMETER :: context = 'treport: '
    REAL(dp) :: work_sum, work_max, work_avg, e
    REAL(dp) :: sbuf(1)
    REAL(dp), ALLOCATABLE :: rbuf(:)
    INTEGER :: nprocs, rank, ierror

    sbuf(1) = t%dt_work
    CALL mpi_comm_rank(comm, rank, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL test_abort(context//'mpi_comm_rank failed', filename, __LINE__)
    CALL mpi_comm_size(comm, nprocs, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL test_abort(context//'mpi_comm_size failed', filename, __LINE__)
    ALLOCATE(rbuf(0:nprocs-1))
    rbuf = -1.0_dp
    CALL mpi_gather(sbuf, 1, MPI_DOUBLE_PRECISION, &
         &  rbuf, 1, MPI_DOUBLE_PRECISION, &
         &  0, comm, ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_GATHER failed', &
         filename, __LINE__)

    IF (rank == 0) THEN
      IF (cmp_dbl_arrays(1, rbuf, sbuf)) &
           CALL test_abort(context//'internal error (1)', &
           filename, __LINE__)
      IF (ANY(rbuf < 0.0_dp)) CALL test_abort(context//'internal error (2)', &
           filename, __LINE__)
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

  SUBROUTINE mysync
    CHARACTER(len=*), PARAMETER :: context = 'mysync: '
    INTEGER :: ierror
    REAL(dp) :: t0, dt

    t0 = mpi_wtime()

    CALL mpi_barrier(MPI_COMM_WORLD, IERROR)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_BARRIER failed', &
         filename, __LINE__)

    dt = (MPI_WTIME() - t0)
    sync_dt_sum = sync_dt_sum + dt

  END SUBROUTINE mysync

  REAL(dp) FUNCTION work_time()
    work_time = MPI_WTIME() - sync_dt_sum
    RETURN
  END FUNCTION work_time

  PURE SUBROUTINE id_map_i2(map)
    INTEGER(i2), INTENT(out) :: map(:,:)

    INTEGER :: i, j, m, n

    m = SIZE(map, 1)
    n = SIZE(map, 2)
    DO j = 1, n
      DO i = 1, m
        map(i,j) = INT((j - 1) * m + i, i2)
      ENDDO
    ENDDO

  END SUBROUTINE id_map_i2

  PURE SUBROUTINE id_map_i4(map)
    INTEGER(i4), INTENT(out) :: map(:,:)

    INTEGER :: i, j, m, n

    m = SIZE(map, 1)
    n = SIZE(map, 2)
    DO j = 1, n
      DO i = 1, m
        map(i,j) = INT((j - 1) * m + i, i4)
      ENDDO
    ENDDO

  END SUBROUTINE id_map_i4

  PURE SUBROUTINE id_map_i8(map)
    INTEGER(i8), INTENT(out) :: map(:,:)

    INTEGER :: i, j, m, n

    m = SIZE(map, 1)
    n = SIZE(map, 2)
    DO j = 1, n
      DO i = 1, m
        map(i,j) = INT((j - 1) * m + i, i8)
      ENDDO
    ENDDO

  END SUBROUTINE id_map_i8

  SUBROUTINE test_abort_msl_f(msg, source, line)
    CHARACTER(*), INTENT(in) :: msg
    CHARACTER(*), INTENT(in) :: source
    INTEGER, INTENT(in) :: line
    CALL test_abort_cmsl_f(mpi_comm_world, msg, source, line)
  END SUBROUTINE test_abort_msl_f

  SUBROUTINE test_abort_cmsl_f(comm, msg, source, line)
    INTEGER, INTENT(in):: comm
    CHARACTER(*), INTENT(in) :: msg
    CHARACTER(*), INTENT(in) :: source
    INTEGER, INTENT(in) :: line

    INTERFACE
      SUBROUTINE c_abort() BIND(c, name='abort')
      END SUBROUTINE c_abort
    END INTERFACE

    INTEGER :: ierror
    LOGICAL :: flag

    WRITE (0, '(3a,i0,2a)') 'Fatal error in ', source, ', line ', line, &
         ': ', TRIM(msg)
    FLUSH(0)
    CALL mpi_initialized(flag, ierror)
    IF (ierror == mpi_success .AND. flag) &
         CALL mpi_abort(comm, 1, ierror)
    CALL c_abort
  END SUBROUTINE test_abort_cmsl_f

  SUBROUTINE icmp_2d(label, f,g, rank)
    CHARACTER(len=*), PARAMETER :: context = 'ftest_common::icmp_2d: '
    CHARACTER(len=*), INTENT(in) :: label
    INTEGER, INTENT(in) :: f(:,:), g(:,:)
    INTEGER, INTENT(in) :: rank

    INTEGER :: i, j, n1, n2

    n1 = SIZE(f,1)
    n2 = SIZE(f,2)
    IF (SIZE(g,1) /= n1 .OR. SIZE(g,2) /= n2) &
         CALL test_abort(context//'shape mismatch error', filename, __LINE__)

    DO j = 1, n2
      DO i = 1, n1
        IF (f(i,j) /= g(i,j)) THEN
          WRITE(0,'(2a,4(a,i0))') context, label, ' test failed: i=', &
               i, ', j=', j, ', f(i,j)=', f(i,j), ', g(i,j)=', g(i,j)
          CALL test_abort(context//label//' test failed', filename, __LINE__)
        ENDIF
      ENDDO
    ENDDO
    IF (verbose) WRITE(0,*) rank,':',context//label//' passed'
  END SUBROUTINE icmp_2d

  SUBROUTINE icmp_3d(label, f,g, rank)
    CHARACTER(len=*), PARAMETER :: context = 'ftest_common::icmp_3d: '
    CHARACTER(len=*), INTENT(in)  :: label
    INTEGER, INTENT(in)  :: f(:,:,:), g(:,:,:)
    INTEGER, INTENT(in) :: rank

    INTEGER :: i1, i2, i3, n1, n2, n3

    n1 = SIZE(f,1)
    n2 = SIZE(f,2)
    n3 = SIZE(f,3)
    IF (SIZE(g,1) /= n1 .OR. SIZE(g,2) /= n2 .OR. SIZE(g,3) /= n3) &
      CALL test_abort(context//label//'shape mismatch', filename, __LINE__)

    DO i3 = 1, n3
      DO i2 = 1, n2
        DO i1 = 1, n1
          IF (f(i1,i2,i3) /= g(i1,i2,i3)) THEN
            WRITE(0,*) context,label,&
                 ' test failed: i1, i2, i3, f(i1,i2,i3), g(i1,i2,i3) =', &
                 i1, i2, i3, f(i1,i2,i3), g(i1,i2,i3)
            CALL test_abort(context//label//' test failed', filename, __LINE__)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    IF (verbose) WRITE(0,*) rank,':',context//label//' passed'
  END SUBROUTINE icmp_3d


  FUNCTION cmp_dbl_arrays_a1d_a1d(a, b) RESULT(differ)
    DOUBLE PRECISION, INTENT(in) :: a(:), b(:)
    LOGICAL :: differ
    INTEGER :: asize, bsize
    INTEGER(c_int) :: asize_c

    asize = SIZE(a)
    bsize = SIZE(b)
    IF (asize /= bsize) THEN
      WRITE (0, '(a)') 'warning: comparing arrays of different size'
      differ = .TRUE.
    ELSE IF (asize > 0) THEN
      asize_c = INT(asize, c_int)
      differ = cmp_dbl_arrays(asize_c, a, b)
    ELSE
      differ = .FALSE.
    END IF
  END FUNCTION cmp_dbl_arrays_a1d_a1d

  FUNCTION cmp_dbl_arrays_a2d_a2d(a, b) RESULT(differ)
    DOUBLE PRECISION, INTENT(in) :: a(:,:), b(:,:)
    LOGICAL :: differ
    INTEGER :: asize, bsize
    INTEGER(c_int) :: asize_c

    asize = SIZE(a)
    bsize = SIZE(b)
    IF (asize /= bsize) THEN
      WRITE (0, '(a)') 'warning: comparing arrays of different size'
      differ = .TRUE.
    ELSE IF (asize > 0) THEN
      asize_c = INT(asize, c_int)
      differ = cmp_dbl_arrays(asize_c, a, b)
    ELSE
      differ = .FALSE.
    END IF
  END FUNCTION cmp_dbl_arrays_a2d_a2d

  FUNCTION cmp_dbl_arrays_a3d_a3d(a, b) RESULT(differ)
    DOUBLE PRECISION, INTENT(in) :: a(:,:,:), b(:,:,:)
    LOGICAL :: differ
    INTEGER :: asize, bsize
    INTEGER(c_int) :: asize_c

    asize = SIZE(a)
    bsize = SIZE(b)
    IF (asize /= bsize) THEN
      WRITE (0, '(a)') 'warning: comparing arrays of different size'
      differ = .TRUE.
    ELSE IF (asize > 0) THEN
      asize_c = INT(asize, c_int)
      differ = cmp_dbl_arrays(asize_c, a, b)
    ELSE
      differ = .FALSE.
    END IF
  END FUNCTION cmp_dbl_arrays_a3d_a3d

  FUNCTION cmp_i2_arrays_a1d_a1d(a, b) RESULT(differ)
    INTEGER(i2), INTENT(in) :: a(:), b(:)
    LOGICAL :: differ
    INTEGER :: asize, bsize
    INTEGER(c_int) :: asize_c

    asize = SIZE(a)
    bsize = SIZE(b)
    IF (asize /= bsize) THEN
      WRITE (0, '(a)') 'warning: comparing arrays of different size'
      differ = .TRUE.
    ELSE IF (asize > 0) THEN
      asize_c = INT(asize, c_int)
      differ = cmp_int16_arrays(asize_c, a, b)
    ELSE
      differ = .FALSE.
    END IF
  END FUNCTION cmp_i2_arrays_a1d_a1d

  FUNCTION cmp_i4_arrays_a1d_a1d(a, b) RESULT(differ)
    INTEGER(i4), INTENT(in) :: a(:), b(:)
    LOGICAL :: differ
    INTEGER :: asize, bsize
    INTEGER(c_int) :: asize_c

    asize = SIZE(a)
    bsize = SIZE(b)
    IF (asize /= bsize) THEN
      WRITE (0, '(a)') 'warning: comparing arrays of different size'
      differ = .TRUE.
    ELSE IF (asize > 0) THEN
      asize_c = INT(asize, c_int)
      differ = cmp_int32_arrays(asize_c, a, b)
    ELSE
      differ = .FALSE.
    END IF
  END FUNCTION cmp_i4_arrays_a1d_a1d

  FUNCTION cmp_i8_arrays_a1d_a1d(a, b) RESULT(differ)
    INTEGER(i8), INTENT(in) :: a(:), b(:)
    LOGICAL :: differ
    INTEGER :: asize, bsize
    INTEGER(c_int) :: asize_c

    asize = SIZE(a)
    bsize = SIZE(b)
    IF (asize /= bsize) THEN
      WRITE (0, '(a)') 'warning: comparing arrays of different size'
      differ = .TRUE.
    ELSE IF (asize > 0) THEN
      asize_c = INT(asize, c_int)
      differ = cmp_int64_arrays(asize_c, a, b)
    ELSE
      differ = .FALSE.
    END IF
  END FUNCTION cmp_i8_arrays_a1d_a1d

  FUNCTION cmp_i2_arrays_a2d_a2d(a, b) RESULT(differ)
    INTEGER(i2), INTENT(in) :: a(:,:), b(:,:)
    LOGICAL :: differ
    INTEGER :: asize, bsize
    INTEGER(c_int) :: asize_c

    asize = SIZE(a)
    bsize = SIZE(b)
    IF (asize /= bsize) THEN
      WRITE (0, '(a)') 'warning: comparing arrays of different size'
      differ = .TRUE.
    ELSE IF (asize > 0) THEN
      asize_c = INT(asize, c_int)
      differ = cmp_int16_arrays(asize_c, a, b)
    ELSE
      differ = .FALSE.
    END IF
  END FUNCTION cmp_i2_arrays_a2d_a2d

  FUNCTION cmp_i4_arrays_a2d_a2d(a, b) RESULT(differ)
    INTEGER(i4), INTENT(in) :: a(:,:), b(:,:)
    LOGICAL :: differ
    INTEGER :: asize, bsize
    INTEGER(c_int) :: asize_c

    asize = SIZE(a)
    bsize = SIZE(b)
    IF (asize /= bsize) THEN
      WRITE (0, '(a)') 'warning: comparing arrays of different size'
      differ = .TRUE.
    ELSE IF (asize > 0) THEN
      asize_c = INT(asize, c_int)
      differ = cmp_int32_arrays(asize_c, a, b)
    ELSE
      differ = .FALSE.
    END IF
  END FUNCTION cmp_i4_arrays_a2d_a2d

  FUNCTION cmp_i8_arrays_a2d_a2d(a, b) RESULT(differ)
    INTEGER(i8), INTENT(in) :: a(:,:), b(:,:)
    LOGICAL :: differ
    INTEGER :: asize, bsize
    INTEGER(c_int) :: asize_c

    asize = SIZE(a)
    bsize = SIZE(b)
    IF (asize /= bsize) THEN
      WRITE (0, '(a)') 'warning: comparing arrays of different size'
      differ = .TRUE.
    ELSE IF (asize > 0) THEN
      asize_c = INT(asize, c_int)
      differ = cmp_int64_arrays(asize_c, a, b)
    ELSE
      differ = .FALSE.
    END IF
  END FUNCTION cmp_i8_arrays_a2d_a2d

  FUNCTION cmp_i2_arrays_a3d_a3d(a, b) RESULT(differ)
    INTEGER(i2), INTENT(in) :: a(:,:,:), b(:,:,:)
    LOGICAL :: differ
    INTEGER :: asize, bsize
    INTEGER(c_int) :: asize_c

    asize = SIZE(a)
    bsize = SIZE(b)
    IF (asize /= bsize) THEN
      WRITE (0, '(a)') 'warning: comparing arrays of different size'
      differ = .TRUE.
    ELSE IF (asize > 0) THEN
      asize_c = INT(asize, c_int)
      differ = cmp_int16_arrays(asize_c, a, b)
    ELSE
      differ = .FALSE.
    END IF
  END FUNCTION cmp_i2_arrays_a3d_a3d

  FUNCTION cmp_i4_arrays_a3d_a3d(a, b) RESULT(differ)
    INTEGER(i4), INTENT(in) :: a(:,:,:), b(:,:,:)
    LOGICAL :: differ
    INTEGER :: asize, bsize
    INTEGER(c_int) :: asize_c

    asize = SIZE(a)
    bsize = SIZE(b)
    IF (asize /= bsize) THEN
      WRITE (0, '(a)') 'warning: comparing arrays of different size'
      differ = .TRUE.
    ELSE IF (asize > 0) THEN
      asize_c = INT(asize, c_int)
      differ = cmp_int32_arrays(asize_c, a, b)
    ELSE
      differ = .FALSE.
    END IF
  END FUNCTION cmp_i4_arrays_a3d_a3d

  FUNCTION cmp_i8_arrays_a3d_a3d(a, b) RESULT(differ)
    INTEGER(i8), INTENT(in) :: a(:,:,:), b(:,:,:)
    LOGICAL :: differ
    INTEGER :: asize, bsize
    INTEGER(c_int) :: asize_c

    asize = SIZE(a)
    bsize = SIZE(b)
    IF (asize /= bsize) THEN
      WRITE (0, '(a)') 'warning: comparing arrays of different size'
      differ = .TRUE.
    ELSE IF (asize > 0) THEN
      asize_c = INT(asize, c_int)
      differ = cmp_int64_arrays(asize_c, a, b)
    ELSE
      differ = .FALSE.
    END IF
  END FUNCTION cmp_i8_arrays_a3d_a3d

  SUBROUTINE factorize(c, a, b)
    INTEGER, INTENT(in) :: c
    INTEGER, INTENT(out) :: a, b ! c = a*b

    INTEGER :: x0, i

    IF (c<1) CALL test_abort('factorize: invalid process space', &
         filename, __LINE__)
    IF (c <= 3 .OR. c == 5 .OR. c == 7) THEN
      a = c
      b = 1
      RETURN
    ENDIF

    ! simple approach, we try to be near c = (2*x) * x
    x0 = INT(SQRT(0.5 * REAL(c)) + 0.5)
    a = 2*x0
    f_loop: DO i = a, 1, -1
      IF (MOD(c,i) == 0) THEN
        a = i
        b = c/i
        EXIT f_loop
      ENDIF
    ENDDO f_loop

  END SUBROUTINE factorize

  !> computes n**(1/3)
  FUNCTION icbrt_i2(n) RESULT(icbrt)
    INTEGER(i2), INTENT(in) :: n
    INTEGER(i2) :: icbrt
    INTEGER, PARAMETER :: nbits = BIT_SIZE(n)-1
    INTEGER(i2) :: s
    INTEGER(i2) :: b, x

    x = ABS(n)
    icbrt = 0
    DO s = nbits, 0, -3
      icbrt = icbrt + icbrt
      b = 3_i2 * icbrt * (icbrt + 1_i2) + 1_i2
      IF (ISHFT(x, -s) >= b) THEN
        x = x - ISHFT(b, s)
        icbrt = icbrt + 1_i2
      END IF
    END DO
    icbrt = SIGN(icbrt, n)
  END FUNCTION icbrt_i2

  !> computes n**(1/3)
  FUNCTION icbrt_i4(n) RESULT(icbrt)
    INTEGER(i4), INTENT(in) :: n
    INTEGER(i4) :: icbrt
    INTEGER, PARAMETER :: nbits = BIT_SIZE(n)-1
    INTEGER(i4) :: s
    INTEGER(i4) :: b, x

    x = ABS(n)
    icbrt = 0
    DO s = nbits-1, 0, -3
      icbrt = icbrt + icbrt
      b = 3_i4 * icbrt * (icbrt + 1_i4) + 1_i4
      IF (ISHFT(x, -s) >= b) THEN
        x = x - ISHFT(b, s)
        icbrt = icbrt + 1_i4
      END IF
    END DO
    icbrt = SIGN(icbrt, n)
  END FUNCTION icbrt_i4

  !> computes n**(1/3)
  FUNCTION icbrt_i8(n) RESULT(icbrt)
    INTEGER(i8), INTENT(in) :: n
    INTEGER(i8) :: icbrt
    INTEGER, PARAMETER :: nbits = BIT_SIZE(n)-1
    INTEGER(i8) :: s
    INTEGER(i8) :: b, x

    x = ABS(n)
    icbrt = 0
    DO s = nbits, 0, -3
      icbrt = icbrt + icbrt
      b = 3_i8 * icbrt * (icbrt + 1_i8) + 1_i8
      IF (ISHFT(x, -s) >= b) THEN
        x = x - ISHFT(b, s)
        icbrt = icbrt + 1_i8
      END IF
    END DO
    icbrt = SIGN(icbrt, n)
  END FUNCTION icbrt_i8


  SUBROUTINE regular_deco(g_cn, c0, cn)
    INTEGER, INTENT(in) :: g_cn
    INTEGER, INTENT(out) :: c0(0:), cn(0:)

    ! convention: process space coords start at 0, grid point coords start at 1

    integer :: tn
    INTEGER :: d, m
    INTEGER :: it

    tn = SIZE(c0)
    IF (tn<0) CALL test_abort('(tn<0)', filename, __LINE__)
    IF (tn>g_cn) CALL test_abort('regular_deco: too many task for such a core&
         & region', &
         filename, __LINE__)

    d = g_cn/tn
    m = MOD(g_cn, tn)

    DO it = 0, m-1
      cn(it) = d + 1
    ENDDO
    DO it = m, tn-1
      cn(it) = d
    ENDDO

    c0(0)=0
    DO it = 1, tn-1
      c0(it) = c0(it-1) + cn(it-1)
    ENDDO
    IF (c0(tn-1)+cn(tn-1) /= g_cn) &
         CALL test_abort('regular_deco: internal error 1', filename, __LINE__)
  END SUBROUTINE regular_deco

  FUNCTION run_randomized_tests() RESULT(fully_random_tests)
    LOGICAL :: fully_random_tests
    CHARACTER(len=32) :: envval
    INTEGER :: envlen, envstat
    CALL get_environment_variable("YAXT_FULLY_RANDOM_TESTS", envval, envlen, &
         status=envstat)
    IF (envstat == 0 .AND. (envlen == 1 .OR. envlen == 3)) THEN
      IF (envlen == 1 .AND. (envval(1:1) == 'y' .OR. envval(1:1) == 'Y' &
           &                 .OR. envval(1:1) == '1')) THEN
        fully_random_tests = .TRUE.
      ELSE IF (str2lower(envval(1:3)) == 'yes') THEN
        fully_random_tests = .TRUE.
      ELSE
        fully_random_tests = .FALSE.
      END IF
    ELSE
      fully_random_tests = .FALSE.
    END IF
  END FUNCTION run_randomized_tests

  FUNCTION str2lower(s) RESULT(t)
    CHARACTER(len=*), INTENT(in) :: s
    CHARACTER(len=LEN(s)) :: t
    INTEGER, PARAMETER :: idel = ICHAR('a')-ICHAR('A')
    INTEGER :: i
    DO i = 1, LEN_TRIM(s)
      t(i:i) = CHAR( ICHAR(s(i:i)) &
           + MERGE(idel, 0,       ICHAR(s(i:i)) >= ICHAR('A') &
           &                .AND. ICHAR(s(i:i)) <= ICHAR('Z')))
    ENDDO
  END FUNCTION str2lower

  SUBROUTINE init_fortran_random(full_random)
    LOGICAL, INTENT(in) :: full_random
    INTEGER, ALLOCATABLE :: rseed(:)

    INTEGER :: rseed_size, i
    CHARACTER(len=32) :: fmt
    INTEGER :: tparts(8), timeseed
    INTEGER :: days_per_month(12), days_prefix
    INTEGER, PARAMETER :: tparts_mult(7) = (/ &
         365 * 24 * 60 * 60, & ! year
         0,                  & ! sum over days_per_month added to day
         24 * 60 * 60,       & ! day
         0,                  & ! ignore timezone offset
         60 * 60,            & ! hour of day
         60,                 & ! minute of hour
         1 /)                  ! seconnd
    CHARACTER(len=32) :: envval
    INTEGER :: envlen, envstat

    CALL random_seed(size=rseed_size)
    ALLOCATE(rseed(rseed_size))
    DO i = 1, rseed_size
      rseed(i) = 4711
    END DO
    IF (full_random) THEN

      CALL date_and_TIME(values=tparts)
      days_per_month( 1) = 31
      days_per_month( 2) = MERGE(28, 29, &
           MOD(tparts(1), 4) == 0 .AND. (     MOD(tparts(1), 100) /= 0 &
           &                             .OR. MOD(tparts(1), 400) == 0))
      days_per_month( 3) = 31
      days_per_month( 4) = 30
      days_per_month( 5) = 31
      days_per_month( 6) = 30
      days_per_month( 7) = 31
      days_per_month( 8) = 31
      days_per_month( 9) = 30
      days_per_month(10) = 31
      days_per_month(11) = 30
      days_per_month(12) = 31
      tparts(1) = tparts(1) - 1970
      days_prefix = SUM(days_per_month(1:tparts(2)-1))
      tparts(3) = tparts(3) + days_prefix - 1
      tparts(2) = 0
      timeseed = SUM(tparts(1:7) * tparts_mult)
      timeseed = IEOR(tparts(8), timeseed) ! mix in microseconds
      rseed(1) = timeseed
      CALL get_environment_VARIABLE("YAXT_RANDOM_SEED", envval, envlen, &
           status=envstat)
      IF (envstat == 0) THEN
        WRITE (fmt, '(a,i0,a)') '(i', DIGITS(rseed), ')'
        READ(envval(1:envlen), fmt) rseed(1)
      END IF
      WRITE(0, '(a,i0)') 'used extra seed=', rseed(1)
      FLUSH(0)
    END IF
    CALL random_seed(put=rseed)
  END SUBROUTINE init_fortran_random

END MODULE ftest_common
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
