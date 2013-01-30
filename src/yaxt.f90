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
       c_long, c_short, c_long_long, c_ptr, c_null_ptr, c_f_pointer

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: xt_abort, xt_idxlist, xt_idxlist_delete, xt_idxvec_new,  &
       xt_idxlist_f2c, &
       xt_idxlist_get_pack_size, xt_idxlist_pack, xt_idxlist_unpack, &
       xt_idxlist_get_intersection, &
       xt_idxlist_get_index_stripes, xt_idxlist_get_bounding_box, &
       &    xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       &    xt_redist, xt_redist_p2p_off_new, xt_redist_p2p_new, &
       &    xt_redist_p2p_blocks_off_new, xt_redist_p2p_blocks_new, &
       &    xt_redist_delete, &
       &    xt_redist_s_exchange1, xt_idxstripes_new, xt_initialize, &
       xt_get_default_comm, xt_idxempty_new

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
    SUBROUTINE xt_redist_delete(redist) BIND(C, name='xt_redist_delete_f')
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

CONTAINS

  SUBROUTINE  xt_abort(comm, msg, source, line)
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
        IMPORT :: Xt_idxlist, xt_int_kind, c_ptr
        IMPLICIT NONE
        INTEGER(xt_int_kind), INTENT(in) :: idxvec(*)
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: num_indices
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_idxvec_new_f
    END INTERFACE

    res%cptr = xt_idxvec_new_f(idxvec, num_indices)

  END FUNCTION xt_idxvec_new

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
        TYPE(c_ptr), INTENT(in) :: idxlist
        INTEGER(c_int), INTENT(in) :: position
        INTEGER(xt_int_kind), INTENT(out) :: idx
        INTEGER(c_int) :: res
      END FUNCTION xt_idxlist_get_index_at_position_c
    END INTERFACE

    res = xt_idxlist_get_index_at_position_c(idxlist%cptr, &
         INT(position, c_int), idx) /= 0
  END FUNCTION xt_idxlist_get_index_at_position

  FUNCTION xt_idxlist_get_position_of_index(idxlist, idx, position) &
       RESULT(notfound)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), VALUE, INTENT(in) :: idx
    INTEGER, INTENT(out) :: position
    LOGICAL :: notfound

    INTERFACE
      FUNCTION xt_idxlist_get_position_of_index_c(idxlist, idx, &
           position) BIND(c, name='xt_idxlist_get_position_of_index') &
           RESULT(res)
        IMPORT :: Xt_idxlist, xt_int_kind, c_int, c_ptr
        TYPE(c_ptr) :: idxlist
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: idx
        INTEGER(c_int), INTENT(out) :: position
        INTEGER(c_int) :: res
      END FUNCTION xt_idxlist_get_position_of_index_c
    END INTERFACE

    notfound = xt_idxlist_get_position_of_index_c(idxlist%cptr, idx, position) &
         /= 0
  END FUNCTION xt_idxlist_get_position_of_index

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

END MODULE yaxt
