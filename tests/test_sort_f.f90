PROGRAM test_sort
  USE iso_c_binding, ONLY: c_int, c_double, c_int32_t
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_finalize, &
       xt_sort_index, xt_sort_int, xt_sort_permutation
  USE xt_core, ONLY: i4, i8
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort, &
       random_fill, crc32, permute
  IMPLICIT NONE
  INTEGER, PARAMETER :: range_size=500000
  CHARACTER(len=*), PARAMETER :: filename = 'test_sort_f.f90'
  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  CALL test_int_sort()
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE test_int_sort()
    INTEGER(c_int), ALLOCATABLE :: a(:), permutation(:)
    REAL(c_double) :: rnd
    INTEGER :: n, i
    INTEGER(c_int32_t) :: crc_before, crc_after
    LOGICAL :: err_found
    CALL RANDOM_NUMBER(rnd)
    n = NINT(rnd * range_size)
    ALLOCATE(a(n), permutation(n))
    CALL random_fill(a)
    CALL xt_sort_int(a)
    err_found = .FALSE.
    DO i = 2, n
      err_found = err_found .OR. a(i-1) > a(i)
    END DO
    IF (err_found) CALL test_abort("discontinuity found!", filename, __LINE__)
    CALL random_fill(a)
    DO i = 1, n
      permutation(i) = i
    END DO
    crc_before = crc32(a)
    CALL xt_sort_permutation(a, permutation)
    DO i = 2, n
      err_found = err_found .OR. a(i-1) > a(i)
    END DO
    IF (err_found) CALL test_abort("discontinuity found!", filename, __LINE__)
    CALL permute(a, permutation)
    crc_after = crc32(a)
    IF (crc_after /= crc_before) &
         CALL test_abort("error in data reconstruction!", filename, __LINE__)

    DEALLOCATE(a, permutation)
  END SUBROUTINE test_int_sort

END PROGRAM test_sort
