!>
!! @file xt_dmap_f.f90
!! @brief xt_dmap-related procedures of Fortran interface
!!
!! @copyright Copyright  (C)  2018 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!
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

!>
!! @example test_dmap_parallel_f.f90

MODULE xt_dmap_base
  USE xt_core, ONLY: xt_abort, xt_mpi_fint_kind, i2, i4, i8
  USE xt_xmap_abstract, ONLY: xt_xmap
  USE xt_redist_base, ONLY: xt_redist, xt_redist_c2f
  USE iso_c_binding, ONLY: c_int, c_null_ptr, c_ptr, c_loc, c_associated
  USE xt_mpi, ONLY: mpi_address_kind
  IMPLICIT NONE
  PRIVATE
  ! note: this type must not be extended to contain any other
  ! components, its memory pattern has to match void * exactly, which
  ! it does because of C constraints
  TYPE, BIND(C), PUBLIC :: xt_dmap
#ifndef __G95__
    PRIVATE
#endif
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_dmap

  ENUM, BIND( C )
    ENUMERATOR :: XT_REORDER_NONE, XT_REORDER_SEND_UP, XT_REORDER_RECV_UP
  END ENUM
  INTEGER, PARAMETER :: xt_reorder_type_kind = KIND(XT_REORDER_NONE)

  INTERFACE
    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_dmap_f2c(dmap) BIND(c, name='xt_dmap_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_dmap
      IMPLICIT NONE
      TYPE(xt_dmap), INTENT(in) :: dmap
      TYPE(c_ptr) :: p
    END FUNCTION xt_dmap_f2c
  END INTERFACE

  INTERFACE xt_dmap_delete
    MODULE PROCEDURE xt_dmap_delete_1
    MODULE PROCEDURE xt_dmap_delete_a1d
  END INTERFACE xt_dmap_delete

  INTERFACE
    SUBROUTINE xt_dmap_delete_c(dmap) &
         BIND(C, name='xt_dmap_delete')
      IMPORT :: c_ptr
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in) :: dmap
    END SUBROUTINE xt_dmap_delete_c
  END INTERFACE

  INTERFACE xt_is_null
    MODULE PROCEDURE xt_dmap_is_null
  END INTERFACE xt_is_null

  INTERFACE xt_dmap_repeat
    MODULE PROCEDURE xt_dmap_repeat_a1d
    MODULE PROCEDURE xt_dmap_repeat_i2_a1d
    MODULE PROCEDURE xt_dmap_repeat_i4_a1d
    MODULE PROCEDURE xt_dmap_repeat_i8_a1d
  END INTERFACE xt_dmap_repeat

  INTERFACE
    SUBROUTINE xt_dmap_repeat_f(dmap_f, num_repetitions, src_displacements, &
         dst_displacements) BIND(C, name='xt_dmap_repeat_f')
      IMPORT:: xt_dmap, c_int, xt_mpi_fint_kind
      TYPE(xt_dmap), INTENT(in) :: dmap_f
      INTEGER(c_int), VALUE, INTENT(in) :: num_repetitions
      INTEGER(xt_mpi_fint_kind), INTENT(in) :: src_displacements(num_repetitions)
      INTEGER(xt_mpi_fint_kind), INTENT(in) :: dst_displacements(num_repetitions)
    END SUBROUTINE xt_dmap_repeat_f
  END INTERFACE

  PUBLIC :: xt_dmap_c2f, xt_dmap_f2c, xt_is_null, xt_dmap_delete, &
       xt_dmap_offset_new, xt_dmap_generate_redist, xt_dmap_repeat, &
       xt_dmap_reorder, xt_reorder_type_kind, &
       XT_REORDER_NONE, XT_REORDER_SEND_UP, XT_REORDER_RECV_UP
CONTAINS

  FUNCTION xt_dmap_is_null(dmap) RESULT(p)
    TYPE(xt_dmap), INTENT(in) :: dmap
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(dmap%cptr)
  END FUNCTION xt_dmap_is_null

  FUNCTION xt_dmap_c2f(dmap) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: dmap
    TYPE(xt_dmap) :: p
    p%cptr = dmap
  END FUNCTION xt_dmap_c2f

  SUBROUTINE xt_dmap_delete_1(dmap)
    TYPE(xt_dmap), INTENT(inout) :: dmap
    CALL xt_dmap_delete_c(dmap%cptr)
    dmap%cptr = c_null_ptr
  END SUBROUTINE xt_dmap_delete_1

  SUBROUTINE xt_dmap_delete_a1d(dmaps)
    TYPE(xt_dmap), INTENT(inout) :: dmaps(:)
    INTEGER :: i, n
    n = SIZE(dmaps)
    DO i = 1, n
      CALL xt_dmap_delete_c(dmaps(i)%cptr)
      dmaps(i)%cptr = c_null_ptr
    END DO
  END SUBROUTINE xt_dmap_delete_a1d

  SUBROUTINE xt_dmap_repeat_a1d(dmap, src_offsets, dst_offsets)
    IMPLICIT NONE
    TYPE(xt_dmap), INTENT(in) :: dmap
    INTEGER, INTENT(in) :: src_offsets(:)
    INTEGER, INTENT(in) :: dst_offsets(:)
    INTEGER :: num_repetitions_src, num_repetitions_dst
    num_repetitions_src = SIZE(src_offsets)
    num_repetitions_dst = SIZE(dst_offsets)
    IF (num_repetitions_src /= num_repetitions_dst) &
         CALL xt_abort("invalid number of repetitions", &
         __FILE__, &
         __LINE__)
    IF (num_repetitions_src < 0_c_int .OR. &
        num_repetitions_src > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", &
         __FILE__, &
         __LINE__)
    CALL xt_dmap_repeat_f( &
      dmap, INT(num_repetitions_src, c_int), src_offsets, dst_offsets)
  END SUBROUTINE xt_dmap_repeat_a1d

  SUBROUTINE xt_dmap_repeat_i2_a1d(dmap, num_repetitions, src_offsets, &
                                   dst_offsets)
    IMPLICIT NONE
    TYPE(xt_dmap), INTENT(in) :: dmap
    INTEGER(i2), INTENT(in) :: num_repetitions
    INTEGER, INTENT(in) :: src_offsets(num_repetitions)
    INTEGER, INTENT(in) :: dst_offsets(num_repetitions)
    IF (num_repetitions < 0_c_int .OR. &
        num_repetitions > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", &
         __FILE__, &
         __LINE__)
    CALL xt_dmap_repeat_f( &
      dmap, INT(num_repetitions, c_int), src_offsets, dst_offsets)
  END SUBROUTINE xt_dmap_repeat_i2_a1d

  SUBROUTINE xt_dmap_repeat_i4_a1d(dmap, num_repetitions, src_offsets, &
                                   dst_offsets)
    IMPLICIT NONE
    TYPE(xt_dmap), INTENT(in) :: dmap
    INTEGER(i4), INTENT(in) :: num_repetitions
    INTEGER, INTENT(in) :: src_offsets(num_repetitions)
    INTEGER, INTENT(in) :: dst_offsets(num_repetitions)
    IF (num_repetitions < 0_c_int .OR. &
        num_repetitions > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", &
         __FILE__, &
         __LINE__)
    CALL xt_dmap_repeat_f( &
      dmap, INT(num_repetitions, c_int), src_offsets, dst_offsets)
  END SUBROUTINE xt_dmap_repeat_i4_a1d

  SUBROUTINE xt_dmap_repeat_i8_a1d(dmap, num_repetitions, src_offsets, &
                                   dst_offsets)
    IMPLICIT NONE
    TYPE(xt_dmap), INTENT(in) :: dmap
    INTEGER(i8), INTENT(in) :: num_repetitions
    INTEGER, INTENT(in) :: src_offsets(num_repetitions)
    INTEGER, INTENT(in) :: dst_offsets(num_repetitions)
    IF (num_repetitions < 0_c_int .OR. &
        num_repetitions > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", &
         __FILE__, &
         __LINE__)
    CALL xt_dmap_repeat_f( &
      dmap, INT(num_repetitions, c_int), src_offsets, dst_offsets)
  END SUBROUTINE xt_dmap_repeat_i8_a1d

  SUBROUTINE xt_dmap_reorder(dmap, reorder_type)
    IMPLICIT NONE
    TYPE(xt_dmap), INTENT(in) :: dmap
    INTEGER(xt_reorder_type_kind), INTENT(in) :: reorder_type

    INTERFACE
      SUBROUTINE xt_dmap_reorder_f(dmap, reorder_type) &
           BIND(C, name='xt_dmap_reorder_f')
        IMPORT:: xt_dmap, c_int
        TYPE(xt_dmap), INTENT(in) :: dmap
        INTEGER(c_int), VALUE, INTENT(in) :: reorder_type
      END SUBROUTINE xt_dmap_reorder_f
    END INTERFACE

    IF (reorder_type < 0_xt_reorder_type_kind .OR. &
         reorder_type > HUGE(1_c_int)) &
         CALL xt_abort("invalid reorder type", &
         __FILE__, &
         __LINE__)
    CALL xt_dmap_reorder_f(dmap, INT(reorder_type, c_int))
  END SUBROUTINE xt_dmap_reorder

  FUNCTION xt_dmap_generate_redist(dmap, datatype) RESULT(redist)
    IMPLICIT NONE
    TYPE(xt_dmap), INTENT(in) :: dmap
    INTEGER, VALUE, INTENT(in) :: datatype
    TYPE(xt_redist) :: redist

    INTERFACE
      FUNCTION xt_dmap_generate_redist_f(dmap_f, datatype) &
           BIND(c, name='xt_dmap_generate_redist_f') RESULT(redist)
        IMPORT :: xt_dmap, xt_mpi_fint_kind, c_ptr
        TYPE(xt_dmap), INTENT(in) :: dmap_f
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: redist
      END FUNCTION xt_dmap_generate_redist_f
    END INTERFACE

    redist = xt_redist_c2f(xt_dmap_generate_redist_f(dmap, datatype))
  END FUNCTION xt_dmap_generate_redist

  FUNCTION xt_dmap_offset_new(xmap, src_offsets, dst_offsets) RESULT(dmap)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: src_offsets(*)
    INTEGER, INTENT(in) :: dst_offsets(*)
    TYPE(xt_dmap) :: dmap

    INTERFACE
      FUNCTION xt_dmap_offset_new_f(xmap_f, src_offsets, dst_offsets) &
           BIND(C, name='xt_dmap_offset_new_f') RESULT(dmap)
        IMPORT :: xt_xmap, xt_mpi_fint_kind, c_ptr
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap_f
        INTEGER(xt_mpi_fint_kind), INTENT(in) :: src_offsets(*)
        INTEGER(xt_mpi_fint_kind), INTENT(in) :: dst_offsets(*)
        TYPE(c_ptr) :: dmap
      END FUNCTION xt_dmap_offset_new_f
    END INTERFACE

    dmap = xt_dmap_c2f(&
         xt_dmap_offset_new_f(xmap, src_offsets, dst_offsets))
  END FUNCTION xt_dmap_offset_new

END MODULE xt_dmap_base
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
