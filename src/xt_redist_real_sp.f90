!>
!! @file xt_redist_real_sp.f90
!! @brief convenience wrappers of xt_redist exchanges for Fortran data
!!
!! @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!

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
MODULE xt_redist_real_sp
  USE xt_redist_base, ONLY: xt_redist, xt_redist_s_exchange1
  USE iso_c_binding, ONLY: c_ptr, c_loc
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: ps = 6
  INTEGER, PARAMETER :: rs = 37
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(ps,rs)
  PUBLIC :: sp
  INTERFACE xt_redist_s_exchange
    MODULE PROCEDURE xt_redist_s_exchange_sp_1d
    MODULE PROCEDURE xt_redist_s_exchange_sp_2d
    MODULE PROCEDURE xt_redist_s_exchange_sp_3d
    MODULE PROCEDURE xt_redist_s_exchange_sp_4d
    MODULE PROCEDURE xt_redist_s_exchange_sp_5d
    MODULE PROCEDURE xt_redist_s_exchange_sp_6d
    MODULE PROCEDURE xt_redist_s_exchange_sp_7d
  END INTERFACE xt_redist_s_exchange
  PUBLIC :: xt_redist_s_exchange
CONTAINS

  ! see @ref xt_redist_s_exchange
  SUBROUTINE xt_redist_s_exchange_sp_1d_as(redist, src_size, src_data, &
       dst_size, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, INTENT(in) :: src_size, dst_size
    REAL(sp), TARGET, INTENT(in) :: src_data(src_size)
    REAL(sp), TARGET, INTENT(inout) :: dst_data(dst_size)
    TYPE(c_ptr) :: src_data_cptr, dst_data_cptr
    src_data_cptr = C_LOC(src_data)
    dst_data_cptr = C_LOC(dst_data)
    CALL xt_redist_s_exchange1(redist, src_data_cptr, dst_data_cptr)
  END SUBROUTINE xt_redist_s_exchange_sp_1d_as

  ! see @ref xt_redist_s_exchange
  SUBROUTINE xt_redist_s_exchange_sp_1d(redist, src_data, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    REAL(sp), TARGET, INTENT(in) :: src_data(:)
    REAL(sp), TARGET, INTENT(inout) :: dst_data(:)

    REAL(sp), POINTER :: src_p(:), dst_p(:)
    REAL(sp), TARGET :: dummy(1)
    INTEGER :: src_size, dst_size
    src_size = SIZE(src_data)
    dst_size = SIZE(dst_data)
    IF (src_size > 0) THEN
      src_p => src_data
    ELSE
      src_p => dummy
      src_size = 1
    END IF
    IF (dst_size > 0) THEN
      dst_p => dst_data
    ELSE
      dst_p => dummy
      dst_size = 1
    END IF
    CALL xt_redist_s_exchange_sp_1d_as(redist, src_size, src_p, dst_size, dst_p)
  END SUBROUTINE xt_redist_s_exchange_sp_1d

  ! see @ref xt_redist_s_exchange
  SUBROUTINE xt_redist_s_exchange_sp_2d(redist, src_data, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    REAL(sp), TARGET, INTENT(in) :: src_data(:,:)
    REAL(sp), TARGET, INTENT(inout) :: dst_data(:,:)

    REAL(sp), POINTER :: src_p(:,:), dst_p(:,:)
    REAL(sp), TARGET :: dummy(1,1)
    INTEGER :: src_size, dst_size
    src_size = SIZE(src_data)
    dst_size = SIZE(dst_data)
    IF (src_size > 0) THEN
      src_p => src_data
    ELSE
      src_p => dummy
      src_size = 1
    END IF
    IF (dst_size > 0) THEN
      dst_p => dst_data
    ELSE
      dst_p => dummy
      dst_size = 1
    END IF
    CALL xt_redist_s_exchange_sp_1d_as(redist, src_size, src_p, dst_size, dst_p)
  END SUBROUTINE xt_redist_s_exchange_sp_2d

  ! see @ref xt_redist_s_exchange
  SUBROUTINE xt_redist_s_exchange_sp_3d(redist, src_data, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    REAL(sp), TARGET, INTENT(in) :: src_data(:,:,:)
    REAL(sp), TARGET, INTENT(inout) :: dst_data(:,:,:)

    REAL(sp), POINTER :: src_p(:,:,:), dst_p(:,:,:)
    REAL(sp), TARGET :: dummy(1,1,1)
    INTEGER :: src_size, dst_size
    src_size = SIZE(src_data)
    dst_size = SIZE(dst_data)
    IF (src_size > 0) THEN
      src_p => src_data
    ELSE
      src_p => dummy
      src_size = 1
    END IF
    IF (dst_size > 0) THEN
      dst_p => dst_data
    ELSE
      dst_p => dummy
      dst_size = 1
    END IF
    CALL xt_redist_s_exchange_sp_1d_as(redist, src_size, src_p, dst_size, dst_p)
  END SUBROUTINE xt_redist_s_exchange_sp_3d

  ! see @ref xt_redist_s_exchange
  SUBROUTINE xt_redist_s_exchange_sp_4d(redist, src_data, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    REAL(sp), TARGET, INTENT(in) :: src_data(:,:,:,:)
    REAL(sp), TARGET, INTENT(inout) :: dst_data(:,:,:,:)

    REAL(sp), POINTER :: src_p(:,:,:,:), dst_p(:,:,:,:)
    REAL(sp), TARGET :: dummy(1,1,1,1)
    INTEGER :: src_size, dst_size
    src_size = SIZE(src_data)
    dst_size = SIZE(dst_data)
    IF (src_size > 0) THEN
      src_p => src_data
    ELSE
      src_p => dummy
      src_size = 1
    END IF
    IF (dst_size > 0) THEN
      dst_p => dst_data
    ELSE
      dst_p => dummy
      dst_size = 1
    END IF
    CALL xt_redist_s_exchange_sp_1d_as(redist, src_size, src_p, dst_size, dst_p)
  END SUBROUTINE xt_redist_s_exchange_sp_4d

  ! see @ref xt_redist_s_exchange
  SUBROUTINE xt_redist_s_exchange_sp_5d(redist, src_data, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    REAL(sp), TARGET, INTENT(in) :: src_data(:,:,:,:,:)
    REAL(sp), TARGET, INTENT(inout) :: dst_data(:,:,:,:,:)

    REAL(sp), POINTER :: src_p(:,:,:,:,:), dst_p(:,:,:,:,:)
    REAL(sp), TARGET :: dummy(1,1,1,1,1)
    INTEGER :: src_size, dst_size
    src_size = SIZE(src_data)
    dst_size = SIZE(dst_data)
    IF (src_size > 0) THEN
      src_p => src_data
    ELSE
      src_p => dummy
      src_size = 1
    END IF
    IF (dst_size > 0) THEN
      dst_p => dst_data
    ELSE
      dst_p => dummy
      dst_size = 1
    END IF
    CALL xt_redist_s_exchange_sp_1d_as(redist, src_size, src_p, dst_size, dst_p)
  END SUBROUTINE xt_redist_s_exchange_sp_5d

  ! see @ref xt_redist_s_exchange
  SUBROUTINE xt_redist_s_exchange_sp_6d(redist, src_data, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    REAL(sp), TARGET, INTENT(in) :: src_data(:,:,:,:,:,:)
    REAL(sp), TARGET, INTENT(inout) :: dst_data(:,:,:,:,:,:)

    REAL(sp), POINTER :: src_p(:,:,:,:,:,:), dst_p(:,:,:,:,:,:)
    REAL(sp), TARGET :: dummy(1,1,1,1,1,1)
    INTEGER :: src_size, dst_size
    src_size = SIZE(src_data)
    dst_size = SIZE(dst_data)
    IF (src_size > 0) THEN
      src_p => src_data
    ELSE
      src_p => dummy
      src_size = 1
    END IF
    IF (dst_size > 0) THEN
      dst_p => dst_data
    ELSE
      dst_p => dummy
      dst_size = 1
    END IF
    CALL xt_redist_s_exchange_sp_1d_as(redist, src_size, src_p, dst_size, dst_p)
  END SUBROUTINE xt_redist_s_exchange_sp_6d

  ! see @ref xt_redist_s_exchange
  SUBROUTINE xt_redist_s_exchange_sp_7d(redist, src_data, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    REAL(sp), TARGET, INTENT(in) :: src_data(:,:,:,:,:,:,:)
    REAL(sp), TARGET, INTENT(inout) :: dst_data(:,:,:,:,:,:,:)

    REAL(sp), POINTER :: src_p(:,:,:,:,:,:,:), dst_p(:,:,:,:,:,:,:)
    REAL(sp), TARGET :: dummy(1,1,1,1,1,1,1)
    INTEGER :: src_size, dst_size
    src_size = SIZE(src_data)
    dst_size = SIZE(dst_data)
    IF (src_size > 0) THEN
      src_p => src_data
    ELSE
      src_p => dummy
      src_size = 1
    END IF
    IF (dst_size > 0) THEN
      dst_p => dst_data
    ELSE
      dst_p => dummy
      dst_size = 1
    END IF
    CALL xt_redist_s_exchange_sp_1d_as(redist, src_size, src_p, dst_size, dst_p)
  END SUBROUTINE xt_redist_s_exchange_sp_7d
END MODULE xt_redist_real_sp
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! mode: f90
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
