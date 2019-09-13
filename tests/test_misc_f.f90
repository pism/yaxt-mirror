PROGRAM test_misc
  USE iso_c_binding, ONLY: c_long, c_int
  USE xt_core, ONLY: i2, i4, i8
  USE ftest_common, ONLY: icbrt, test_abort, run_randomized_tests, &
       init_fortran_random
  IMPLICIT NONE
  LOGICAL :: fully_random_tests
  fully_random_tests = run_randomized_tests()
  CALL init_fortran_random(fully_random_tests)
  CALL test_icbrt
CONTAINS
  SUBROUTINE test_icbrt
    CALL test_icbrt_i2
    CALL test_icbrt_i4
    CALL test_icbrt_i8
  END SUBROUTINE test_icbrt

  SUBROUTINE test_icbrt_i2
    INTEGER(i2), PARAMETER :: ulim = 31_i2
    INTEGER(i2) :: i, cubed, cbrt
    CHARACTER(len=132) :: msg
    DO i = -ulim, ulim, 1_i2
      cubed = i**3_i2
      cbrt = icbrt(cubed)
      IF (cbrt /= i) THEN
        WRITE (msg, '(4(a,i0))') &
             "integer cube root computation failed for ", i, &
             "**3 = ", cubed, ", but icbrt(", cubed, ") = ", cbrt
        CALL test_abort(TRIM(msg), &
             __FILE__, &
             __LINE__)
      END IF
    END DO
  END SUBROUTINE test_icbrt_i2

  SUBROUTINE test_icbrt_i4
    INTEGER(i4), PARAMETER :: ulim = 1290
    INTEGER(i4) :: i, cubed, cbrt, prev_cubed, other_cubed
    CHARACTER(len=132) :: msg
    DOUBLE PRECISION :: rnd
    i = icbrt(-8)
    prev_cubed = -HUGE(1_i4)
    DO i = -ulim, -1_i4, 1_i4
      cubed = i**3_i4
      cbrt = icbrt(cubed)
      IF (cbrt /= i) THEN
        WRITE (msg, '(4(a,i0))') &
             "integer cube root computation failed for ", i, &
             "**3 = ", cubed, ", but icbrt(", cubed, ") = ", cbrt
        CALL test_abort(TRIM(msg), &
             __FILE__, &
             __LINE__)
      END IF
      CALL random_number(rnd)
      other_cubed = prev_cubed + MAX(1_i4, INT(rnd * DBLE(cubed - prev_cubed), i4))
      cbrt = icbrt(other_cubed)
      IF (cbrt /= i) THEN
        WRITE (msg, '(4(a,i0))') &
             "integer cube root computation failed for ", other_cubed, &
             ", expected ", i, &
             ", but icbrt(", other_cubed, ") = ", cbrt
        CALL test_abort(TRIM(msg), &
             __FILE__, &
             __LINE__)
      END IF
      prev_cubed = cubed
    END DO
    prev_cubed = HUGE(1_i4)
    DO i = ulim, 0_i4, -1_i4
      cubed = i**3_i4
      cbrt = icbrt(cubed)
      IF (cbrt /= i) THEN
        WRITE (msg, '(4(a,i0))') &
             "integer cube root computation failed for ", i, &
             "**3 = ", cubed, ", but icbrt(", cubed, ") = ", cbrt
        CALL test_abort(TRIM(msg), &
             __FILE__, &
             __LINE__)
      END IF
      CALL random_number(rnd)
      other_cubed = cubed + MIN(INT(rnd * DBLE(prev_cubed - cubed), i4), -1_i4)
      cbrt = icbrt(cubed)
      IF (cbrt /= i) THEN
        WRITE (msg, '(4(a,i0))') &
             "integer cube root computation failed for ", other_cubed, &
             ", expected ", i, &
             ", but icbrt(", other_cubed, ") = ", cbrt
        CALL test_abort(TRIM(msg), &
             __FILE__, &
             __LINE__)
      END IF
      prev_cubed = cubed
    END DO
  END SUBROUTINE test_icbrt_i4

  SUBROUTINE test_icbrt_i8
    INTEGER(i8), PARAMETER :: ulim = 2097151_i8
    INTEGER(i8) :: i, cubed, cbrt
    CHARACTER(len=132) :: msg
    DO i = -ulim, ulim, 1_i8
      ! need to spell that out explcitely for crayftn 8.5.5 at least
#if defined _CRAYFTN && _REALEASE_MAJOR == 8 && _RELEASE_MINOR == 5
      cubed = i * i * i
#else
      cubed = i**3
#endif
      cbrt = icbrt(cubed)
      IF (cbrt /= i) THEN
        WRITE (msg, '(4(a,i0))') &
             "integer cube root computation failed for ", i, &
             "**3 = ", cubed, ", but icbrt(", cubed, ") = ", cbrt
        CALL test_abort(TRIM(msg), &
             __FILE__, &
             __LINE__)
      END IF
    END DO
  END SUBROUTINE test_icbrt_i8
END PROGRAM test_misc
