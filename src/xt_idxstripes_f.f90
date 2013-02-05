!>
!> @file yaxt.f90
!> @brief Fortran interface to yaxt implementation
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
MODULE xt_idxstripes
  USE xt_core, ONLY: xt_int_kind, xt_stripe
  USE xt_idxlist_abstract, ONLY: xt_idxlist, xt_idxlist_c2f
  USE iso_c_binding, ONLY: c_ptr
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xt_idxstripes_new
CONTAINS
  FUNCTION xt_idxstripes_new(stripes, num_stripes) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_stripe), INTENT(in) :: stripes(*)
    INTEGER(xt_int_kind), VALUE, INTENT(in) :: num_stripes
    TYPE(xt_idxlist) :: res

    INTERFACE
      FUNCTION xt_idxstripes_new_f(stripes, num_stripes) &
           BIND(C, name='xt_idxstripes_new') RESULT(res_ptr)
        IMPORT:: xt_idxlist, xt_stripe, xt_int_kind, c_ptr
        IMPLICIT NONE
        TYPE(xt_stripe), INTENT(in) :: stripes(*)
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: num_stripes
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_idxstripes_new_f
    END INTERFACE

    res = xt_idxlist_c2f(xt_idxstripes_new_f(stripes, num_stripes))

  END FUNCTION xt_idxstripes_new
END MODULE xt_idxstripes
