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

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_null_char, c_int, &
       c_long, c_short, c_long_long, c_ptr, c_null_ptr, c_f_pointer, c_intptr_t

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: xt_initialize, xt_finalize, xt_abort, &
       xt_idxlist, xt_idxlist_delete,  &
       xt_idxlist_c2f, xt_idxlist_f2c, xt_idxlist_copy, &
       xt_idxlist_get_pack_size, xt_idxlist_pack, xt_idxlist_unpack, &
       xt_idxlist_get_intersection, xt_idxlist_get_num_indices, &
       xt_idxlist_get_positions_of_indices, &
       xt_idxlist_get_index_at_position, xt_idxlist_get_indices_at_positions, &
       xt_idxlist_get_position_of_index, xt_idxlist_get_position_of_index_off, &
       xt_idxlist_get_index_stripes, xt_idxlist_get_bounding_box, &
       xt_idxvec_new, xt_idxvec_from_stripes_new, xt_idxstripes_new, &
       xt_idxlist_collection_new, xt_idxsection_new, &
       &    xt_idxlist_get_indices, xt_modifier, xt_idxmod_new, &
       &    xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       &    xt_redist, xt_redist_p2p_off_new, xt_redist_p2p_new, &
       &    xt_redist_p2p_blocks_off_new, xt_redist_p2p_blocks_new, &
       &    xt_redist_delete, &
       &    xt_redist_s_exchange1, &
       xt_get_default_comm, xt_idxempty_new, &
       OPERATOR(==), OPERATOR(/=)

  INTEGER, PARAMETER, PUBLIC :: xt_int_kind   = XT_INT_FC_KIND

  TYPE, BIND(C), PUBLIC :: xt_stripe
    INTEGER(xt_int_kind) :: start
    INTEGER(xt_int_kind) :: nstrides
    INTEGER(xt_int_kind) :: stride
  END TYPE xt_stripe

  TYPE, BIND(C) :: Xt_idxlist
    PRIVATE
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE Xt_idxlist

  TYPE, BIND(C) :: xt_xmap
    PRIVATE
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_xmap

  TYPE, BIND(C) :: xt_redist
    PRIVATE
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_redist

  TYPE, BIND(C), PUBLIC :: xt_bounds
    INTEGER(xt_int_kind) :: start, size
  END TYPE xt_bounds

  TYPE, BIND(C) :: xt_modifier
    TYPE(Xt_idxlist) :: extract
    TYPE(Xt_idxlist) :: subst
    INTEGER(c_int) :: mask
  END TYPE xt_modifier

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE xt_bounds_eq
  END INTERFACE OPERATOR(==)

  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE xt_bounds_ne
  END INTERFACE OPERATOR(/=)

  INTERFACE

    SUBROUTINE xt_initialize(default_comm) BIND(C, name='xt_initialize_f')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER, INTENT(in) :: default_comm
    END SUBROUTINE xt_initialize

    FUNCTION xt_get_default_comm() RESULT(comm) &
         BIND(c, name='xt_get_default_comm_f')
      INTEGER :: comm
    END FUNCTION xt_get_default_comm

  END INTERFACE

  INTERFACE
    FUNCTION xt_idxlist_f2c(idxlist) BIND(c, name='xt_idxlist_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_idxlist
      IMPLICIT NONE
      TYPE(xt_idxlist), INTENT(in) :: idxlist
      TYPE(c_ptr) :: p
    END FUNCTION xt_idxlist_f2c

    SUBROUTINE xt_idxlist_delete(idxlist) BIND(C, name='xt_idxlist_delete_f')
      IMPORT :: Xt_idxlist
      IMPLICIT NONE
      TYPE(Xt_idxlist), INTENT(inout) :: idxlist
    END SUBROUTINE xt_idxlist_delete

    FUNCTION xt_idxlist_get_pack_size(idxlist, comm) &
         BIND(c, name='xt_idxlist_get_pack_size_f2c') RESULT(pack_size)
      IMPORT :: xt_idxlist
      IMPLICIT NONE
      TYPE(xt_idxlist), INTENT(in) :: idxlist
      INTEGER, VALUE, INTENT(in) :: comm
      INTEGER :: pack_size
    END FUNCTION xt_idxlist_get_pack_size

  END INTERFACE

  ! xt_idxlist_pack(type(xt_idxlist), intent(out) :: idxlist,
  !                 type<*>, intent(inout) :: buffer, integer,
  !                 intent(in) :: buffer_size,
  !                 integer, intent(inout) :: position, integer,
  !                 intent(in) :: comm)
  EXTERNAL :: xt_idxlist_pack

  ! xt_idxlist_unpack(type(xt_idxlist), intent(out) :: idxlist,
  !                   type<*>, intent(in) :: buffer,
  !                   integer, intent(in) :: buffer_size,
  !                   integer, intent(inout) :: position, integer,
  !                   intent(in) :: comm)
  EXTERNAL :: xt_idxlist_unpack

  INTERFACE
    SUBROUTINE xt_xmap_delete(xmap) BIND(C, name='xt_xmap_delete_f')
      IMPORT :: xt_xmap
      IMPLICIT NONE
      TYPE(xt_xmap), INTENT(inout) :: xmap
    END SUBROUTINE xt_xmap_delete
  END INTERFACE

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
      TYPE(c_ptr), INTENT(in) :: src_data_cptr
      TYPE(c_ptr), INTENT(in) :: dst_data_cptr
    END SUBROUTINE xt_redist_s_exchange1
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_finalize() BIND(C, name='xt_finalize')
    END SUBROUTINE xt_finalize
  END INTERFACE

  INTERFACE
    INTEGER(xt_int_kind) FUNCTION xt_idxlist_get_num_indices(idxlist) &
         &          BIND(C, name='xt_idxlist_get_num_indices_f') RESULT(res)
      IMPORT :: xt_int_kind, xt_idxlist
      TYPE(xt_idxlist), INTENT(in) :: idxlist
    END FUNCTION xt_idxlist_get_num_indices
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_idxlist_get_indices(idxlist, indices) &
         &          BIND(C, name='xt_idxlist_get_indices_f')
      IMPORT :: xt_int_kind, xt_idxlist
      TYPE(xt_idxlist), INTENT(in) :: idxlist
      INTEGER(xt_int_kind), INTENT(out) :: indices(*)
    END SUBROUTINE xt_idxlist_get_indices
  END INTERFACE

CONTAINS

  SUBROUTINE xt_abort(comm, msg, source, line)
    INTEGER, INTENT(in):: comm
    CHARACTER(len=*), INTENT(in) :: msg
    CHARACTER(len=*), INTENT(in) :: source
    INTEGER, INTENT(in) :: line
    INTERFACE
      SUBROUTINE xt_abort_f(comm, msg, source, line) BIND(c, name='xt_abort_f')
        IMPORT :: c_char
        INTEGER, VALUE, INTENT(in):: comm
        CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: msg
        CHARACTER(kind=c_char), DIMENSION(*), INTENT(in) :: source
        INTEGER, VALUE, INTENT(in) :: line
      END SUBROUTINE xt_abort_f
    END INTERFACE
    CALL xt_abort_f(comm, TRIM(msg)//c_null_char, &
         TRIM(source)//c_null_char, line)
  END SUBROUTINE xt_abort

  FUNCTION xt_idxvec_new(idxvec, num_indices) RESULT(res)
    IMPLICIT NONE
    INTEGER(xt_int_kind), INTENT(in) :: idxvec(*)
    INTEGER(xt_int_kind), VALUE, INTENT(in) :: num_indices
    TYPE(Xt_idxlist) :: res

    INTERFACE
      FUNCTION xt_idxvec_new_f(idxvec, num_indices) &
           BIND(C, name='xt_idxvec_new') RESULT(res_ptr)
        IMPORT :: xt_int_kind, c_ptr
        IMPLICIT NONE
        INTEGER(xt_int_kind), INTENT(in) :: idxvec(*)
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: num_indices
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_idxvec_new_f
    END INTERFACE

    res%cptr = xt_idxvec_new_f(idxvec, num_indices)

  END FUNCTION xt_idxvec_new

  FUNCTION xt_idxvec_from_stripes_new(stripes) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    TYPE(Xt_idxlist) :: res

    INTERFACE
      FUNCTION xt_idxvec_from_stripes_new_c(stripes, num_stripes) &
           BIND(C, name='xt_idxvec_from_stripes_new') RESULT(res_ptr)
        IMPORT :: xt_stripe, xt_int_kind, c_ptr
        IMPLICIT NONE
        TYPE(xt_stripe), INTENT(in) :: stripes(*)
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: num_stripes
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_idxvec_from_stripes_new_c
    END INTERFACE

    res%cptr = xt_idxvec_from_stripes_new_c(stripes, &
         INT(SIZE(stripes), xt_int_kind))
  END FUNCTION xt_idxvec_from_stripes_new

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

    res%cptr = xt_idxempty_new_c()

  END FUNCTION xt_idxempty_new

  FUNCTION xt_idxsection_new(start, global_size, local_size, local_start) &
       RESULT(idxsection)
    INTEGER(xt_int_kind), INTENT(in) :: start, global_size(:), local_size(:), &
         local_start(:)
    TYPE(xt_idxlist) :: idxsection
    INTEGER :: num_dimensions
    INTERFACE
      FUNCTION xt_idxsection_new_c(start, num_dimensions, global_size, &
           local_size, local_start) BIND(c, name='xt_idxsection_new') &
           RESULT(idxsection)
        IMPORT :: c_int, c_ptr, xt_idxlist, xt_int_kind
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: start
        INTEGER(c_int), VALUE, INTENT(in) :: num_dimensions
        INTEGER(xt_int_kind), INTENT(in) :: global_size(*), local_size(*), &
             local_start(*)
        TYPE(c_ptr) :: idxsection
      END FUNCTION xt_idxsection_new_c
    END INTERFACE

    num_dimensions = SIZE(global_size)
    IF (SIZE(local_size) /= num_dimensions &
         .OR. SIZE(local_start) /= num_dimensions) &
         CALL xt_abort(xt_get_default_comm(), "non-matching array sizes", &
         __FILE__, &
         __LINE__)
    idxsection%cptr = xt_idxsection_new_c(start, INT(num_dimensions, c_int), &
         global_size, local_size, local_start)
  END FUNCTION xt_idxsection_new

  FUNCTION xt_idxstripes_new(stripes, num_stripes) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_stripe), INTENT(in) :: stripes(*)
    INTEGER(xt_int_kind), VALUE, INTENT(in) :: num_stripes
    TYPE(Xt_idxlist) :: res

    INTERFACE
      FUNCTION xt_idxstripes_new_f(stripes, num_stripes) &
           BIND(C, name='xt_idxstripes_new') RESULT(res_ptr)
        IMPORT:: Xt_idxlist, xt_stripe, xt_int_kind, c_ptr
        IMPLICIT NONE
        TYPE(xt_stripe), INTENT(in) :: stripes(*)
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: num_stripes
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_idxstripes_new_f
    END INTERFACE

    res%cptr = xt_idxstripes_new_f(stripes, num_stripes)

  END FUNCTION xt_idxstripes_new

  FUNCTION xt_idxlist_collection_new(idxlists) RESULT(res)
    TYPE(xt_idxlist), INTENT(in) :: idxlists(:)
    TYPE(xt_idxlist) :: res

    TYPE(c_ptr) :: p(SIZE(idxlists))
    INTEGER(c_intptr_t) :: i
    INTERFACE
      FUNCTION xt_idxlist_collection_new_c(idxlists, num_lists) &
           BIND(C, name='xt_idxlist_collection_new') RESULT(res_ptr)
        IMPORT:: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), INTENT(in) :: idxlists(*)
        INTEGER(c_int), VALUE, INTENT(in) :: num_lists
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_idxlist_collection_new_c
    END INTERFACE

    IF (SIZE(idxlists) > INT(HUGE(1_c_int))) &
         CALL xt_abort(xt_get_default_comm(), "idxlists array too large", &
         __FILE__, &
         __LINE__)
    p = idxlists(:)%cptr
    IF (SIZE(idxlists) > 5000) WRITE (0, *) "p=", TRANSFER(p, i, SIZE(p))
    res%cptr = xt_idxlist_collection_new_c(p, &
         INT(SIZE(idxlists), c_int))

  END FUNCTION xt_idxlist_collection_new

  FUNCTION xt_xmap_all2all_new(src_idxlist, dst_idxlist, comm) RESULT(res)
    IMPLICIT NONE
    TYPE(Xt_idxlist), INTENT(in) :: src_idxlist
    TYPE(Xt_idxlist), INTENT(in) :: dst_idxlist
    INTEGER, VALUE, INTENT(in) :: comm
    TYPE(xt_xmap) :: res

    INTERFACE
      FUNCTION xt_xmap_all2all_new_f(src_idxlist, dst_idxlist, comm) &
           BIND(C, name='xt_xmap_all2all_new_f') RESULT(res_ptr)
        IMPORT :: Xt_idxlist, xt_xmap, c_ptr
        IMPLICIT NONE
        TYPE(Xt_idxlist), INTENT(in) :: src_idxlist
        TYPE(Xt_idxlist), INTENT(in) :: dst_idxlist
        INTEGER, VALUE, INTENT(in) :: comm
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_xmap_all2all_new_f
    END INTERFACE

    res%cptr = xt_xmap_all2all_new_f(src_idxlist, dst_idxlist, comm)
  END FUNCTION xt_xmap_all2all_new

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
        IMPORT :: xt_xmap, xt_redist, c_int, c_ptr
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(c_int), INTENT(in) :: src_block_offsets(*)
        INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
        INTEGER(c_int), INTENT(in) :: dst_block_offsets(*)
        INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
        INTEGER, VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_blocks_off_new_f
    END INTERFACE

    res%cptr = xt_redist_p2p_blocks_off_new_f(xmap, src_block_offsets, &
           src_block_sizes, src_block_num, &
           dst_block_offsets, dst_block_sizes, dst_block_num, &
           datatype)

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
        IMPORT :: xt_xmap, xt_redist, c_int, c_ptr
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
        INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
        INTEGER, VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_blocks_new_f
    END INTERFACE

    res%cptr = xt_redist_p2p_blocks_new_f(xmap, &
         &                                src_block_sizes, src_block_num, &
         &                                dst_block_sizes, dst_block_num, &
         &                                datatype)
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
        IMPORT :: xt_xmap, xt_redist, c_int, c_ptr
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER, INTENT(in) :: src_offsets(*)
        INTEGER, INTENT(in) :: dst_offsets(*)
        INTEGER, VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_off_new_f
    END INTERFACE

    res%cptr = xt_redist_p2p_off_new_f(xmap, src_offsets, dst_offsets, datatype)

  END FUNCTION xt_redist_p2p_off_new

  FUNCTION xt_redist_p2p_new(xmap, datatype) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, VALUE, INTENT(in) :: datatype
    TYPE(xt_redist) :: res

    INTERFACE
      FUNCTION xt_redist_p2p_new_f(xmap, datatype) &
           BIND(C, name='xt_redist_p2p_new_f') RESULT(res_ptr)
        IMPORT:: xt_xmap, xt_redist, c_int, c_ptr
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER, VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_new_f
    END INTERFACE

    res%cptr = xt_redist_p2p_new_f(xmap, datatype)

  END FUNCTION xt_redist_p2p_new

  FUNCTION xt_idxlist_c2f(idxlist) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: idxlist
    TYPE(xt_idxlist) :: p
    p%cptr = idxlist
  END FUNCTION xt_idxlist_c2f

  FUNCTION xt_idxlist_get_index_at_position(idxlist, position, idx) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER, VALUE, INTENT(in) :: position
    INTEGER(xt_int_kind), INTENT(out) :: idx
    LOGICAL :: res

    INTERFACE
      FUNCTION xt_idxlist_get_index_at_position_c(idxlist, position, idx) &
           BIND(c, name='xt_idxlist_get_index_at_position') RESULT(res)
        IMPORT :: c_ptr, c_int, xt_int_kind
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(c_int), VALUE, INTENT(in) :: position
        INTEGER(xt_int_kind), INTENT(out) :: idx
        INTEGER(c_int) :: res
      END FUNCTION xt_idxlist_get_index_at_position_c
    END INTERFACE

    res = xt_idxlist_get_index_at_position_c(idxlist%cptr, &
         INT(position, c_int), idx) /= 0
  END FUNCTION xt_idxlist_get_index_at_position

  FUNCTION xt_idxlist_get_indices_at_positions(idxlist, positions, indices, &
       undef_idx) RESULT(num_subst)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER, INTENT(in) :: positions(:)
    INTEGER(xt_int_kind), INTENT(out) :: indices(:)
    INTEGER(xt_int_kind), INTENT(in) :: undef_idx
    INTEGER :: num_subst

    INTERFACE
      FUNCTION xt_idxlist_get_indices_at_positions_c(idxlist, positions, &
           num_pos, indices, undef_idx) &
           BIND(c, name='xt_idxlist_get_indices_at_positions') RESULT(num_subst)
        IMPORT :: c_ptr, c_int, xt_int_kind
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(c_int), INTENT(in) :: positions(*)
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: num_pos, undef_idx
        INTEGER(xt_int_kind), INTENT(out) :: indices(*)
        INTEGER(c_int) :: num_subst
      END FUNCTION xt_idxlist_get_indices_at_positions_c
    END INTERFACE

    num_subst = xt_idxlist_get_indices_at_positions_c(idxlist%cptr, &
         INT(positions, c_int), INT(SIZE(positions), xt_int_kind), &
         indices, undef_idx)
  END FUNCTION xt_idxlist_get_indices_at_positions

  FUNCTION xt_idxlist_get_position_of_index(idxlist, idx, position) &
       RESULT(notfound)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), VALUE, INTENT(in) :: idx
    INTEGER, INTENT(out) :: position
    LOGICAL :: notfound
    INTEGER(c_int) :: position_c

    INTERFACE
      FUNCTION xt_idxlist_get_position_of_index_c(idxlist, idx, position) &
           BIND(c, name='xt_idxlist_get_position_of_index') RESULT(res)
        IMPORT :: Xt_idxlist, xt_int_kind, c_int, c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: idx
        INTEGER(c_int), INTENT(out) :: position
        INTEGER(c_int) :: res
      END FUNCTION xt_idxlist_get_position_of_index_c
    END INTERFACE

    notfound = xt_idxlist_get_position_of_index_c(idxlist%cptr, idx, &
         position_c) /= 0
    position = INT(position_c)
  END FUNCTION xt_idxlist_get_position_of_index

  FUNCTION xt_idxlist_get_position_of_index_off(idxlist, idx, position, &
       offset) RESULT(notfound)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), VALUE, INTENT(in) :: idx
    INTEGER, INTENT(out) :: position
    INTEGER, INTENT(in) :: offset
    LOGICAL :: notfound
    INTEGER(c_int) :: position_c

    INTERFACE
      FUNCTION xt_idxlist_get_position_of_index_off_c(idxlist, idx, position, &
           offset) BIND(c, name='xt_idxlist_get_position_of_index_off') &
           RESULT(res)
        IMPORT :: Xt_idxlist, xt_int_kind, c_int, c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: idx
        INTEGER(c_int), INTENT(out) :: position
        INTEGER(c_int), VALUE, INTENT(in) :: offset
        INTEGER(c_int) :: res
      END FUNCTION xt_idxlist_get_position_of_index_off_c
    END INTERFACE

    notfound = xt_idxlist_get_position_of_index_off_c(idxlist%cptr, idx, &
         position_c, INT(offset, c_int)) /= 0
    position = INT(position_c)
  END FUNCTION xt_idxlist_get_position_of_index_off

  FUNCTION xt_idxlist_get_positions_of_indices(idxlist, indices, positions, &
       single_match_only) RESULT(num_missing)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(in) :: indices(:)
    INTEGER, INTENT(out) :: positions(:)
    LOGICAL, INTENT(in) :: single_match_only
    INTEGER :: num_missing

    INTERFACE
      FUNCTION xt_idxlist_get_positions_of_indices_c(idxlist, indices, &
           num_indices, positions, single_match_only) &
           BIND(c, name='xt_idxlist_get_positions_of_indices') &
           RESULT(num_missing)
        IMPORT :: Xt_idxlist, xt_int_kind, c_int, c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(xt_int_kind), INTENT(in) :: indices(*)
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: num_indices
        INTEGER(c_int), INTENT(out) :: positions(*)
        INTEGER(c_int), VALUE, INTENT(in) :: single_match_only
        INTEGER(c_int) :: num_missing
      END FUNCTION xt_idxlist_get_positions_of_indices_c
    END INTERFACE
    IF (SIZE(positions) < SIZE(indices)) THEN
      CALL xt_abort(xt_get_default_comm(), "positions array too small", &
           __FILE__, &
           __LINE__)
    END IF
    num_missing = xt_idxlist_get_positions_of_indices_c(idxlist%cptr, indices, &
         INT(SIZE(indices), xt_int_kind), positions, &
         INT(MERGE(1, 0, single_match_only), c_int))
  END FUNCTION xt_idxlist_get_positions_of_indices

  SUBROUTINE xt_idxlist_get_index_stripes(idxlist, stripes)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    TYPE(xt_stripe), ALLOCATABLE, INTENT(out) :: stripes(:)

    INTERFACE
      SUBROUTINE xt_idxlist_get_index_stripes_c(idxlist, stripes,&
           num_stripes) BIND(c, name='xt_idxlist_get_index_stripes')
        IMPORT :: c_ptr, xt_int_kind
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        TYPE(c_ptr), INTENT(out) :: stripes
        INTEGER(xt_int_kind) :: num_stripes
      END SUBROUTINE xt_idxlist_get_index_stripes_c

      SUBROUTINE free_c(p) BIND(c, name='free')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE :: p
      END SUBROUTINE free_c
    END INTERFACE
    TYPE(c_ptr) :: stripes_c_ptr
    INTEGER(xt_int_kind) :: num_stripes
    TYPE(xt_stripe), POINTER :: stripes_f_ptr(:)
    INTEGER :: stripes_shape(1)
    CALL xt_idxlist_get_index_stripes_c(idxlist%cptr, stripes_c_ptr,&
         num_stripes)
    stripes_shape(1) = INT(num_stripes)
    IF (num_stripes > 0) THEN
      ALLOCATE(stripes(INT(num_stripes)))
      CALL C_F_POINTER(stripes_c_ptr, stripes_f_ptr, stripes_shape)
      stripes = stripes_f_ptr
    END IF
    CALL free_c(stripes_c_ptr)
  END SUBROUTINE xt_idxlist_get_index_stripes

  FUNCTION xt_idxlist_get_bounding_box(idxlist, global_size, &
       global_start_index) RESULT(bounds)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(in) :: global_size(:)
    INTEGER(xt_int_kind), INTENT(in) :: global_start_index
    TYPE(xt_bounds) :: bounds(SIZE(global_size))

    INTERFACE
      SUBROUTINE xt_idxlist_get_bounding_box_c(idxlist, ndim, global_size, &
           global_start_index, bounds) &
           BIND(c, name='xt_idxlist_get_bounding_box')
        IMPORT :: c_int, c_ptr, xt_int_kind, xt_bounds
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(c_int), VALUE, INTENT(in) :: ndim
        INTEGER(xt_int_kind), INTENT(in) :: global_size(ndim)
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: global_start_index
        TYPE(xt_bounds), INTENT(out) :: bounds(ndim)
      END SUBROUTINE xt_idxlist_get_bounding_box_c
    END INTERFACE

    CALL xt_idxlist_get_bounding_box_c(xt_idxlist_f2c(idxlist), &
         INT(SIZE(global_size), c_int), global_size, global_start_index, bounds)
  END FUNCTION xt_idxlist_get_bounding_box

  FUNCTION xt_idxlist_get_intersection(idxlist_src, idxlist_dst) &
       RESULT(intersection)
    TYPE(xt_idxlist), INTENT(in) :: idxlist_src, idxlist_dst
    TYPE(xt_idxlist) :: intersection

    INTERFACE
      FUNCTION xt_idxlist_get_intersection_c(idxlist_src, idxlist_dst) &
           BIND(c, name='xt_idxlist_get_intersection') RESULT(intersection)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist_src, idxlist_dst
        TYPE(c_ptr) :: intersection
      END FUNCTION xt_idxlist_get_intersection_c
    END INTERFACE

    intersection%cptr = xt_idxlist_get_intersection_c(&
         xt_idxlist_f2c(idxlist_src), xt_idxlist_f2c(idxlist_dst))
  END FUNCTION xt_idxlist_get_intersection

  FUNCTION xt_idxlist_copy(idxlist) RESULT(copy)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    TYPE(xt_idxlist) :: copy

    INTERFACE
      FUNCTION xt_idxlist_copy_c(idxlist) BIND(c, name='xt_idxlist_copy') &
           RESULT(copy)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(IN) :: idxlist
        TYPE(c_ptr) :: copy
      END FUNCTION xt_idxlist_copy_c
    END INTERFACE

    copy%cptr = xt_idxlist_copy_c(idxlist%cptr)

  END FUNCTION xt_idxlist_copy

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
    INTEGER(xt_int_kind), INTENT(in) :: modifier_num
    INTEGER(xt_int_kind), INTENT(inout) :: mstate(*)
    TYPE(Xt_idxlist) :: res

    INTERFACE
      FUNCTION xt_idxmod_new_f(patch, modifier, modifier_num, mstate) &
           BIND(C, name='xt_idxmod_new_f') RESULT(res)
        IMPORT :: xt_idxlist, xt_int_kind, xt_modifier, c_ptr, c_int
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: patch
        TYPE(xt_modifier), INTENT(in) :: modifier(*)
        INTEGER(c_int), VALUE, INTENT(in) :: modifier_num
        INTEGER(c_int), INTENT(in) :: mstate(*)
        TYPE(c_ptr) :: res
      END FUNCTION xt_idxmod_new_f
    END INTERFACE

    res%cptr = xt_idxmod_new_f(patch, modifier, modifier_num, mstate)

  END FUNCTION xt_idxmod_new

END MODULE yaxt
