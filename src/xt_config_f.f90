!>
!! @file xt_config_f.f90
!! @brief Fortran interface to yaxt configuration object
!!
!! @copyright Copyright  (C)  2020 Jörg Behrens <behrens@dkrz.de>
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
MODULE xt_config_f
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_null_ptr
  IMPLICIT NONE
  PRIVATE
  ! note: this type must not be extended to contain any other
  ! components, its memory pattern has to match void * exactly, which
  ! it does because of C constraints
  TYPE, BIND(C), PUBLIC :: xt_config
#ifndef __G95__
    PRIVATE
#endif
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_config
  PUBLIC :: xt_config_new, xt_config_delete
CONTAINS

  FUNCTION xt_config_new() RESULT(config)
    TYPE(xt_config) :: config
    INTERFACE
      FUNCTION xt_config_new_c() RESULT(config) &
           BIND(c, name='xt_config_new')
        IMPORT :: c_ptr
        IMPLICIT NONE
        TYPE(c_ptr) :: config
      END FUNCTION xt_config_new_c
    END INTERFACE
    config%cptr = xt_config_new_c()
  END FUNCTION xt_config_new

  SUBROUTINE xt_config_delete(config)
    TYPE(xt_config), INTENT(in) :: config
    INTERFACE
      SUBROUTINE xt_config_delete_c(config) BIND(c, name='xt_config_delete')
        IMPORT :: c_ptr
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: config
      END SUBROUTINE xt_config_delete_c
    END INTERFACE
    CALL xt_config_delete_c(config%cptr)
  END SUBROUTINE xt_config_delete

END MODULE xt_config_f
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
