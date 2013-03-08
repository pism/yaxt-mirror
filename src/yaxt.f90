!>
!> @file yaxt.f90
!> @brief Fortran interface to yaxt implementation
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

MODULE yaxt
  !
  ! Fortran interface to yaxt implementation
  !
  USE xt_core, ONLY: xt_int_kind, &
       xt_abort, xt_initialize, xt_get_default_comm, xt_stripe, xt_bounds, &
       xt_redist, xt_redist_c2f, xt_redist_f2c, &
       xt_mpi_fint_kind, char
  USE xt_idxlist_abstract, ONLY: &
       xt_idxlist_f2c, xt_idxlist_c2f, xt_idxlist, &
       xt_idxlist_delete, xt_idxlist_get_pack_size, &
       xt_idxlist_pack, xt_idxlist_unpack, xt_idxlist_copy, &
       xt_idxlist_get_index_at_position, xt_idxlist_get_indices_at_positions, &
       xt_idxlist_get_position_of_index, xt_idxlist_get_position_of_index_off, &
       xt_idxlist_get_positions_of_indices, xt_idxlist_get_index_stripes, &
       xt_idxlist_get_bounding_box, xt_idxlist_get_intersection, &
       xt_idxlist_get_num_indices, &
       xt_idxlist_get_indices, xt_idxlist_get_indices_const
  USE xt_idxvec, ONLY: xt_idxvec_from_stripes_new, xt_idxvec_new
  USE xt_idxstripes, ONLY: xt_idxstripes_new
  USE xt_idxsection, ONLY: xt_idxsection_new, xt_idxfsection_new
  USE xt_xmap_abstract, ONLY: xt_xmap, xt_xmap_c2f, xt_xmap_f2c, &
       xt_xmap_delete, xt_xmap_get_num_destinations, xt_xmap_get_num_sources, &
       xt_xmap_get_destination_ranks, xt_xmap_get_source_ranks, &
       xt_xmap_all2all_new
  USE iso_c_binding, ONLY: c_int, c_ptr, c_null_ptr, c_loc

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: xt_initialize, xt_finalize, xt_abort, xt_get_default_comm, &
       xt_int_kind, xt_stripe, xt_bounds, &
       xt_idxlist, xt_idxlist_delete,  &
       xt_idxlist_c2f, xt_idxlist_f2c, xt_idxlist_copy, &
       xt_idxlist_get_pack_size, xt_idxlist_pack, xt_idxlist_unpack, &
       xt_idxlist_get_intersection, xt_idxlist_get_positions_of_indices, &
       xt_idxlist_get_index_at_position, xt_idxlist_get_indices_at_positions, &
       xt_idxlist_get_position_of_index, xt_idxlist_get_position_of_index_off, &
       xt_idxlist_get_num_indices, &
       xt_idxlist_get_indices, xt_idxlist_get_indices_const, &
       xt_idxlist_get_index_stripes, xt_idxlist_get_bounding_box, &
       xt_idxvec_new, xt_idxvec_from_stripes_new, xt_idxstripes_new, &
       xt_idxlist_collection_new, xt_idxsection_new, xt_idxfsection_new, &
       xt_idxmod_new, xt_xmap, xt_xmap_c2f, xt_xmap_f2c, &
       xt_xmap_all2all_new, xt_xmap_delete, xt_xmap_get_num_destinations, &
       xt_xmap_get_num_sources, xt_xmap_get_destination_ranks, &
       xt_xmap_get_source_ranks, &
       xt_redist, xt_redist_f2c, xt_redist_c2f, &
       xt_redist_p2p_off_new, xt_redist_p2p_new, &
       &    xt_redist_p2p_blocks_off_new, xt_redist_p2p_blocks_new, &
       &    xt_redist_delete, &
       &    xt_redist_s_exchange1, &
       xt_idxempty_new, &
       OPERATOR(==), OPERATOR(/=), char

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE xt_bounds_eq
  END INTERFACE OPERATOR(==)

  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE xt_bounds_ne
  END INTERFACE OPERATOR(/=)

  TYPE, BIND(C), PUBLIC :: xt_modifier
    TYPE(xt_idxlist) :: extract
    TYPE(xt_idxlist) :: subst
    INTEGER(c_int) :: mask
  END TYPE xt_modifier

  INTERFACE
    SUBROUTINE xt_redist_delete(redist) &
         BIND(C, name='xt_redist_delete_f')
      IMPORT :: xt_redist
      IMPLICIT NONE
      TYPE(xt_redist), INTENT(inout) :: redist
    END SUBROUTINE xt_redist_delete
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_redist_s_exchange1(redist, src_data_cptr, dst_data_cptr) &
         BIND(C, name='xt_redist_s_exchange1_f')
      IMPORT:: xt_redist, c_ptr
      TYPE(xt_redist), INTENT(in) :: redist
      TYPE(c_ptr), VALUE :: src_data_cptr
      TYPE(c_ptr), VALUE :: dst_data_cptr
    END SUBROUTINE xt_redist_s_exchange1
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_finalize() BIND(C, name='xt_finalize')
    END SUBROUTINE xt_finalize
  END INTERFACE

CONTAINS

  FUNCTION xt_idxempty_new() RESULT(res)
    IMPLICIT NONE
    TYPE(Xt_idxlist) :: res

    INTERFACE
      FUNCTION xt_idxempty_new_c() &
           BIND(C, name='xt_idxempty_new') RESULT(res_ptr)
        IMPORT :: c_ptr
        IMPLICIT NONE
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_idxempty_new_c
    END INTERFACE

    res = xt_idxlist_c2f(xt_idxempty_new_c())

  END FUNCTION xt_idxempty_new

  FUNCTION xt_idxlist_collection_new(idxlists) RESULT(res)
    TYPE(xt_idxlist), INTENT(in) :: idxlists(:)
    TYPE(xt_idxlist) :: res

    INTERFACE
      FUNCTION xt_idxlist_collection_new_c(idxlists, num_lists) &
           BIND(C, name='xt_idxlist_collection_new') RESULT(res_ptr)
        IMPORT:: c_ptr, c_int, xt_idxlist
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: idxlists(*)
        INTEGER(c_int), VALUE, INTENT(in) :: num_lists
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_idxlist_collection_new_c
    END INTERFACE

    IF (SIZE(idxlists) > INT(HUGE(1_c_int))) &
         CALL xt_abort(xt_get_default_comm(), "idxlists array too large", &
         __FILE__, &
         __LINE__)
    res = xt_idxlist_c2f(xt_idxlist_collection_new_c(idxlists, &
         INT(SIZE(idxlists), c_int)))
  END FUNCTION xt_idxlist_collection_new

  FUNCTION xt_redist_p2p_blocks_off_new(xmap, src_block_offsets, &
       src_block_sizes, src_block_num, &
       dst_block_offsets, dst_block_sizes, dst_block_num, &
       datatype) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(in) :: src_block_offsets(*)
    INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
    INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
    INTEGER(c_int), INTENT(in) :: dst_block_offsets(*)
    INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
    INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
    INTEGER, VALUE, INTENT(in) :: datatype
    TYPE(xt_redist) :: res
    INTERFACE
      FUNCTION xt_redist_p2p_blocks_off_new_f(xmap, src_block_offsets, &
           src_block_sizes, src_block_num, &
           dst_block_offsets, dst_block_sizes, dst_block_num, &
           datatype) BIND(C, name='xt_redist_p2p_blocks_off_new_f') &
           RESULT(res_ptr)
        IMPORT :: xt_xmap, xt_redist, xt_mpi_fint_kind, c_int, c_ptr
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(c_int), INTENT(in) :: src_block_offsets(*)
        INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
        INTEGER(c_int), INTENT(in) :: dst_block_offsets(*)
        INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_blocks_off_new_f
    END INTERFACE

    res = xt_redist_c2f(&
         xt_redist_p2p_blocks_off_new_f(xmap, src_block_offsets, &
           src_block_sizes, src_block_num, &
           dst_block_offsets, dst_block_sizes, dst_block_num, &
           datatype))

  END FUNCTION xt_redist_p2p_blocks_off_new

  FUNCTION xt_redist_p2p_blocks_new(xmap, src_block_sizes, src_block_num, &
       &                                  dst_block_sizes, dst_block_num, &
       &                                  datatype) &
       RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
    INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
    INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
    INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
    INTEGER, VALUE, INTENT(in) :: datatype
    TYPE(xt_redist) :: res
    INTERFACE
      FUNCTION xt_redist_p2p_blocks_new_f(xmap, &
           &                              src_block_sizes, src_block_num, &
           &                              dst_block_sizes, dst_block_num, &
           &                              datatype) &
           BIND(C, name='xt_redist_p2p_blocks_new_f') RESULT(res_ptr)
        IMPORT :: xt_xmap, xt_mpi_fint_kind, xt_redist, c_int, c_ptr
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
        INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_blocks_new_f
    END INTERFACE

    res = xt_redist_c2f(&
         xt_redist_p2p_blocks_new_f(xmap, &
         &                          src_block_sizes, src_block_num, &
         &                          dst_block_sizes, dst_block_num, &
         &                          datatype))
  END FUNCTION xt_redist_p2p_blocks_new

  FUNCTION xt_redist_p2p_off_new(xmap, src_offsets, dst_offsets, datatype) &
       RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: src_offsets(*)
    INTEGER, INTENT(in) :: dst_offsets(*)
    INTEGER, VALUE, INTENT(in) :: datatype
    TYPE(xt_redist) :: res

    INTERFACE
      FUNCTION xt_redist_p2p_off_new_f(xmap, src_offsets, dst_offsets, &
           datatype) BIND(C, name='xt_redist_p2p_off_new_f') RESULT(res_ptr)
        IMPORT :: xt_xmap, xt_redist, c_int, c_ptr, xt_mpi_fint_kind
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(xt_mpi_fint_kind), INTENT(in) :: src_offsets(*)
        INTEGER(xt_mpi_fint_kind), INTENT(in) :: dst_offsets(*)
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_off_new_f
    END INTERFACE

    res = xt_redist_c2f(&
         xt_redist_p2p_off_new_f(xmap, src_offsets, dst_offsets, datatype))

  END FUNCTION xt_redist_p2p_off_new

  FUNCTION xt_redist_p2p_new(xmap, datatype) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, VALUE, INTENT(in) :: datatype
    TYPE(xt_redist) :: res

    INTERFACE
      FUNCTION xt_redist_p2p_new_f(xmap, datatype) &
           BIND(C, name='xt_redist_p2p_new_f') RESULT(res_ptr)
        IMPORT:: xt_xmap, xt_redist, c_int, c_ptr, xt_mpi_fint_kind
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_new_f
    END INTERFACE

    res = xt_redist_c2f(xt_redist_p2p_new_f(xmap, datatype))

  END FUNCTION xt_redist_p2p_new

  ELEMENTAL FUNCTION xt_bounds_eq(a, b) RESULT(a_equals_b)
    TYPE(xt_bounds), INTENT(in) :: a, b
    LOGICAL :: a_equals_b
    a_equals_b = a%size == b%size .AND. a%start == b%start
  END FUNCTION xt_bounds_eq

  ELEMENTAL FUNCTION xt_bounds_ne(a, b) RESULT(a_equals_b)
    TYPE(xt_bounds), INTENT(in) :: a, b
    LOGICAL :: a_equals_b
    a_equals_b = a%size /= b%size .OR. a%start /= b%start
  END FUNCTION xt_bounds_ne

  FUNCTION xt_idxmod_new(patch, modifier, modifier_num, mstate) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: patch
    TYPE(xt_modifier), INTENT(in) :: modifier(*)
    INTEGER, INTENT(in) :: modifier_num
    INTEGER(c_int), OPTIONAL, TARGET, INTENT(inout) :: mstate(*)
    TYPE(Xt_idxlist) :: res

    INTERFACE
      FUNCTION xt_idxmod_new_f(patch, modifier, modifier_num, mstate_ptr) &
           BIND(C, name='xt_idxmod_new_f') RESULT(res)
        IMPORT :: xt_idxlist, xt_int_kind, xt_modifier, c_ptr, c_int
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: patch
        TYPE(xt_modifier), INTENT(in) :: modifier(*)
        INTEGER(c_int), VALUE, INTENT(in) :: modifier_num
        TYPE(c_ptr), VALUE :: mstate_ptr
        TYPE(c_ptr) :: res
      END FUNCTION xt_idxmod_new_f
    END INTERFACE

    IF (PRESENT(mstate)) THEN
      res = xt_idxlist_c2f(xt_idxmod_new_f(patch, modifier, &
           INT(modifier_num, c_int), C_LOC(mstate)))
    ELSE
      res = xt_idxlist_c2f(xt_idxmod_new_f(patch, modifier, &
           INT(modifier_num, c_int), c_null_ptr))
    ENDIF

  END FUNCTION xt_idxmod_new

END MODULE yaxt
