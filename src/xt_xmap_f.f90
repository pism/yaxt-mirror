!>
!! @file xt_xmap_f.f90
!! @brief Fortran interface to yaxt xmap declarations
!!
!! @copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
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
!
MODULE xt_xmap_abstract
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_ptr, c_null_ptr, &
       c_associated
  USE xt_core, ONLY: xt_mpi_fint_kind
  USE xt_idxlist_abstract, ONLY: xt_idxlist
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xt_xmap_c2f, xt_xmap_f2c, xt_is_null
  PUBLIC :: xt_xmap_delete, xt_xmap_get_num_destinations, &
       xt_xmap_get_num_sources, xt_xmap_get_destination_ranks, &
       xt_xmap_get_source_ranks
  PUBLIC :: xt_xmap_all2all_new, xt_xmap_dist_dir_new

  ! note: this type must not be extended to contain any other
  ! components, its memory pattern has to match void * exactly, which
  ! it does because of C constraints
  TYPE, BIND(C), PUBLIC :: xt_xmap
    PRIVATE
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_xmap

  INTERFACE
    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_xmap_f2c(xmap) BIND(c, name='xt_xmap_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_xmap
      IMPLICIT NONE
      TYPE(xt_xmap), INTENT(in) :: xmap
      TYPE(c_ptr) :: p
    END FUNCTION xt_xmap_f2c

    SUBROUTINE xt_xmap_delete_c(xmap) BIND(C, name='xt_xmap_delete')
      IMPORT :: c_ptr
      IMPLICIT NONE
      TYPE(c_ptr), VALUE :: xmap
    END SUBROUTINE xt_xmap_delete_c

  END INTERFACE

  INTERFACE xt_xmap_delete
    MODULE PROCEDURE xt_xmap_delete_1
    MODULE PROCEDURE xt_xmap_delete_a1d
  END INTERFACE xt_xmap_delete

  INTERFACE xt_is_null
    MODULE PROCEDURE xt_xmap_is_null
  END INTERFACE xt_is_null

CONTAINS

  FUNCTION xt_xmap_is_null(xmap) RESULT(p)
    TYPE(xt_xmap), INTENT(in) :: xmap
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(xmap%cptr)
  END FUNCTION xt_xmap_is_null


  FUNCTION xt_xmap_c2f(xmap) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: xmap
    TYPE(xt_xmap) :: p
    p%cptr = xmap
  END FUNCTION xt_xmap_c2f

  SUBROUTINE xt_xmap_delete_1(xmap)
    TYPE(xt_xmap), INTENT(inout) :: xmap
    CALL xt_xmap_delete_c(xt_xmap_f2c(xmap))
    xmap%cptr = c_null_ptr
  END SUBROUTINE xt_xmap_delete_1

  SUBROUTINE xt_xmap_delete_a1d(xmaps)
    TYPE(xt_xmap), INTENT(inout) :: xmaps(:)
    INTEGER :: i, n
    n = SIZE(xmaps)
    DO i = 1, n
      CALL xt_xmap_delete_c(xt_xmap_f2c(xmaps(i)))
      xmaps(i)%cptr = c_null_ptr
    END DO
  END SUBROUTINE xt_xmap_delete_a1d

  FUNCTION xt_xmap_all2all_new(src_idxlist, dst_idxlist, comm) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist
    TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
    INTEGER, VALUE, INTENT(in) :: comm
    TYPE(xt_xmap) :: res

    INTERFACE
      FUNCTION xt_xmap_all2all_new_f(src_idxlist, dst_idxlist, comm) &
           BIND(C, name='xt_xmap_all2all_new_f') RESULT(res_ptr)
        IMPORT :: xt_idxlist, xt_xmap, xt_mpi_fint_kind, c_ptr
        IMPLICIT NONE
        TYPE(Xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(Xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_xmap_all2all_new_f
    END INTERFACE

    res = xt_xmap_c2f(xt_xmap_all2all_new_f(src_idxlist, dst_idxlist, comm))
  END FUNCTION xt_xmap_all2all_new

  FUNCTION xt_xmap_dist_dir_new(src_idxlist, dst_idxlist, comm) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist
    TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
    INTEGER, VALUE, INTENT(in) :: comm
    TYPE(xt_xmap) :: res

    INTERFACE
      FUNCTION xt_xmap_dist_dir_new_f(src_idxlist, dst_idxlist, comm) &
           BIND(C, name='xt_xmap_dist_dir_new_f') RESULT(res_ptr)
        IMPORT :: xt_idxlist, xt_xmap, xt_mpi_fint_kind, c_ptr
        IMPLICIT NONE
        TYPE(Xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(Xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_xmap_dist_dir_new_f
    END INTERFACE

    res = xt_xmap_c2f(xt_xmap_dist_dir_new_f(src_idxlist, dst_idxlist, comm))
  END FUNCTION xt_xmap_dist_dir_new

  FUNCTION xt_xmap_get_num_destinations(xmap) RESULT(num)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER :: num
    INTERFACE
      FUNCTION xt_xmap_get_num_destinations_c(xmap) RESULT(num) &
           BIND(c, name='xt_xmap_get_num_destinations')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int) :: num
      END FUNCTION xt_xmap_get_num_destinations_c
    END INTERFACE
    num = INT(xt_xmap_get_num_destinations_c(xmap%cptr))
  END FUNCTION xt_xmap_get_num_destinations

  FUNCTION xt_xmap_get_num_sources(xmap) RESULT(num)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER :: num
    INTERFACE
      FUNCTION xt_xmap_get_num_sources_c(xmap) RESULT(num) &
           BIND(c, name='xt_xmap_get_num_sources')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int) :: num
      END FUNCTION xt_xmap_get_num_sources_c
    END INTERFACE
    num = INT(xt_xmap_get_num_sources_c(xmap%cptr))
  END FUNCTION xt_xmap_get_num_sources

  SUBROUTINE xt_xmap_get_destination_ranks(xmap, ranks)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(out) :: ranks(*)
    INTERFACE
      SUBROUTINE xt_xmap_get_destination_ranks_c(xmap, ranks) &
           BIND(c, name='xt_xmap_get_destination_ranks')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int), INTENT(out) :: ranks(*)
      END SUBROUTINE xt_xmap_get_destination_ranks_c
    END INTERFACE
    CALL xt_xmap_get_destination_ranks_c(xmap%cptr, ranks)
  END SUBROUTINE xt_xmap_get_destination_ranks

  SUBROUTINE xt_xmap_get_source_ranks(xmap, ranks)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(out) :: ranks(*)
    INTERFACE
      SUBROUTINE xt_xmap_get_source_ranks_c(xmap, ranks) &
           BIND(c, name='xt_xmap_get_source_ranks')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int), INTENT(out) :: ranks(*)
      END SUBROUTINE xt_xmap_get_source_ranks_c
    END INTERFACE
    CALL xt_xmap_get_source_ranks_c(xmap%cptr, ranks)
  END SUBROUTINE xt_xmap_get_source_ranks

END MODULE xt_xmap_abstract
