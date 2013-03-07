!>
!> @file xt_core_f.f90
!> @brief Fortran interface to yaxt core declarations
!>
!> @copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
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
MODULE xt_core
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_null_char, c_int, &
       c_long, c_short, c_long_long, c_ptr, c_null_ptr
  IMPLICIT NONE
  PRIVATE
  INTEGER, PUBLIC, PARAMETER :: xt_int_kind   = XT_INT_FC_KIND
  INTEGER, PUBLIC, PARAMETER :: pi2 = 4
  INTEGER, PUBLIC, PARAMETER :: pi4 = 9
  INTEGER, PUBLIC, PARAMETER :: pi8 = 14
  INTEGER, PUBLIC, PARAMETER :: i2 = SELECTED_INT_KIND(pi2)
  INTEGER, PUBLIC, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)
  INTEGER, PUBLIC, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)
  PUBLIC :: xt_initialize, xt_abort,  xt_get_default_comm, &
       xt_redist_c2f, xt_redist_f2c, &
       xt_xmap_c2f, xt_xmap_f2c, char

  INTEGER, PUBLIC, PARAMETER :: xt_mpi_fint_kind = XT_MPI_FINT_FC_KIND
  INTEGER(xt_int_kind), PARAMETER :: dummy = 0_xt_int_kind
  !> number of decimal places needed to print any variable of type
  !! INTEGER(xt_int_kind)
  INTEGER, PUBLIC, PARAMETER :: xt_int_dec_len &
       = CEILING(1.0 + REAL(DIGITS(dummy)) * LOG10(REAL(RADIX(dummy))))
  CHARACTER(9), PARAMETER :: xt_stripe_tag = 'xt_stripe'
  !> maximal length of string xt_stripe(a, b, c)
  INTEGER, PUBLIC, PARAMETER :: xt_stripe2s_len &
       = LEN(xt_stripe_tag) + 2 + 4 + 3 * xt_int_dec_len

  TYPE, BIND(C), PUBLIC :: xt_stripe
    INTEGER(xt_int_kind) :: start
    INTEGER(xt_int_kind) :: nstrides
    INTEGER(xt_int_kind) :: stride
  END TYPE xt_stripe

  ! note: this type must not be extended to contain any other
  ! components, its memory pattern has to match void * exactly, which
  ! it does because of C constraints
  TYPE, BIND(C), PUBLIC :: xt_xmap
    PRIVATE
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_xmap

  ! note: this type must not be extended to contain any other
  ! components, its memory pattern has to match void * exactly, which
  ! it does because of C constraints
  TYPE, BIND(C), PUBLIC :: xt_redist
    PRIVATE
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_redist

  TYPE, BIND(C), PUBLIC :: xt_bounds
    INTEGER(xt_int_kind) :: start, size
  END TYPE xt_bounds

  INTERFACE

    FUNCTION xt_get_default_comm() RESULT(comm) &
         BIND(c, name='xt_get_default_comm_f')
      IMPORT :: xt_mpi_fint_kind
      IMPLICIT NONE
      INTEGER(xt_mpi_fint_kind) :: comm
    END FUNCTION xt_get_default_comm

    SUBROUTINE xt_initialize(default_comm) BIND(C, name='xt_initialize_f')
      IMPORT:: c_int, xt_mpi_fint_kind
      IMPLICIT NONE
      INTEGER(xt_mpi_fint_kind), INTENT(in) :: default_comm
    END SUBROUTINE xt_initialize

    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_redist_f2c(redist) BIND(c, name='xt_redist_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_redist
      IMPLICIT NONE
      TYPE(xt_redist), INTENT(in) :: redist
      TYPE(c_ptr) :: p
    END FUNCTION xt_redist_f2c

    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_xmap_f2c(xmap) BIND(c, name='xt_xmap_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_xmap
      IMPLICIT NONE
      TYPE(xt_xmap), INTENT(in) :: xmap
      TYPE(c_ptr) :: p
    END FUNCTION xt_xmap_f2c

  END INTERFACE

  INTERFACE xt_abort
    MODULE PROCEDURE xt_abort4
    MODULE PROCEDURE xt_abort3
  END INTERFACE xt_abort

  INTERFACE
    SUBROUTINE xt_abort_c(comm, msg, source, line) BIND(c, name='xt_abort_f')
      IMPORT :: c_char, c_int, xt_mpi_fint_kind
      IMPLICIT NONE
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in):: comm
      CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: msg
      CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: source
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: line
    END SUBROUTINE xt_abort_c
  END INTERFACE

  INTERFACE char
    MODULE PROCEDURE xt_stripe2char
  END INTERFACE char

CONTAINS

  SUBROUTINE xt_abort4(comm, msg, source, line)
    INTEGER, INTENT(in):: comm
    CHARACTER(len=*), INTENT(in) :: msg
    CHARACTER(len=*), INTENT(in) :: source
    INTEGER, INTENT(in) :: line
    CALL xt_abort_c(comm, TRIM(msg)//c_null_char, &
         TRIM(source)//c_null_char, INT(line, c_int))
  END SUBROUTINE xt_abort4

  SUBROUTINE xt_abort3(msg, source, line)
    CHARACTER(len=*), INTENT(in) :: msg
    CHARACTER(len=*), INTENT(in) :: source
    INTEGER, INTENT(in) :: line
    CALL xt_abort_c(xt_get_default_comm(), TRIM(msg)//c_null_char, &
         TRIM(source)//c_null_char, INT(line, c_int))
  END SUBROUTINE xt_abort3

  FUNCTION xt_redist_c2f(redist) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: redist
    TYPE(xt_redist) :: p
    p%cptr = redist
  END FUNCTION xt_redist_c2f

  FUNCTION xt_xmap_c2f(xmap) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: xmap
    TYPE(xt_xmap) :: p
    p%cptr = xmap
  END FUNCTION xt_xmap_c2f

  ELEMENTAL FUNCTION xt_stripe2char(stripe) RESULT(str)
    CHARACTER(len=xt_stripe2s_len) :: str
    TYPE(xt_stripe), INTENT(in) :: stripe
    WRITE (str, '(2a,3(i0,a))') xt_stripe_tag, '(', stripe%start, ', ', &
         stripe%nstrides, ', ', stripe%stride, ')'
  END FUNCTION xt_stripe2char

END MODULE xt_core
