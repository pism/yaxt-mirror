MODULE yaxt
  !
  ! Fortran interface to yaxt implementation
  !

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char, c_null_char, c_int,  c_ptr, c_null_ptr, c_loc

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: xt_abort, Xt_idxlist, xt_idxlist_delete, xt_idxvec_new,  &
       &    xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       &    xt_redist, xt_redist_p2p_off_new, xt_redist_p2p_new, &
       &    xt_redist_p2p_blocks_off_new, xt_redist_p2p_blocks_new, xt_redist_delete, &
       &    xt_redist_s_exchange1, xt_idxstripes_new, c_loc, xt_initialize

  INTEGER, PARAMETER, PUBLIC :: xt_idx_kind   = C_INT ! \todo xt_idx_kind and xt_count_kind should come out of configure
  INTEGER, PARAMETER, PUBLIC :: xt_count_kind = C_INT


  TYPE, BIND(C), PUBLIC :: xt_stripe
    INTEGER(xt_idx_kind) :: start
    INTEGER(xt_count_kind) :: nstrides
    INTEGER(xt_idx_kind) :: stride
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

  INTERFACE
    SUBROUTINE xt_abort_cmsl_f(comm_f, msg, source, line)  BIND(C, name='xt_abort_cmsl_f')
      IMPORT:: c_char, c_int
      INTEGER, INTENT(in):: comm_f
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: msg
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: source
      INTEGER(C_INT), VALUE, INTENT(in) :: line
    END SUBROUTINE xt_abort_cmsl_f
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_abort_msl_f(msg, source, line)  BIND(C, name='xt_abort_msl_f')
      IMPORT:: c_char, c_int
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: msg
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: source
      INTEGER(C_INT), VALUE, INTENT(in) :: line
    END SUBROUTINE xt_abort_msl_f
  END INTERFACE

  INTERFACE xt_abort
    MODULE PROCEDURE xt_abort_cmsl
    MODULE PROCEDURE xt_abort_msl
  END INTERFACE xt_abort

  INTERFACE
    SUBROUTINE xt_initialize(default_comm) BIND(C, name='xt_initialize_f')
      IMPORT:: C_INT
      IMPLICIT NONE
      INTEGER, INTENT(in) :: default_comm
    END SUBROUTINE xt_initialize
  END INTERFACE

  INTERFACE
    FUNCTION xt_idxvec_new(idxvec, num_indices) BIND(C, name='xt_idxvec_new_f') RESULT(res)
      IMPORT :: Xt_idxlist, xt_idx_kind, xt_count_kind
      IMPLICIT NONE
      INTEGER(xt_idx_kind), INTENT(in) :: idxvec(*)
      INTEGER(xt_count_kind), VALUE, INTENT(in) :: num_indices
      TYPE(Xt_idxlist) :: res
    END FUNCTION xt_idxvec_new
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_idxlist_delete(idxlist) BIND(C, name='xt_idxlist_delete_f')
      IMPORT :: Xt_idxlist
      IMPLICIT NONE
      TYPE(Xt_idxlist), INTENT(inout) :: idxlist
    END SUBROUTINE xt_idxlist_delete
  END INTERFACE

  INTERFACE
    FUNCTION xt_xmap_all2all_new(src_idxlist, dst_idxlist, comm) BIND(C, name='xt_xmap_all2all_new_f') RESULT(res)
      IMPORT :: Xt_idxlist, xt_xmap
      IMPLICIT NONE
      TYPE(Xt_idxlist), INTENT(in) :: src_idxlist
      TYPE(Xt_idxlist), INTENT(in) :: dst_idxlist
      INTEGER, INTENT(in) :: comm
      TYPE(xt_xmap) :: res
    END FUNCTION xt_xmap_all2all_new
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_xmap_delete(xmap) BIND(C, name='xt_xmap_delete_f')
      IMPORT :: xt_xmap
      IMPLICIT NONE
      TYPE(xt_xmap), INTENT(inout) :: xmap
    END SUBROUTINE xt_xmap_delete
  END INTERFACE

  INTERFACE
    FUNCTION xt_redist_p2p_blocks_off_new(xmap, src_block_offsets, src_block_sizes, src_block_num, &
         &                                  dst_block_offsets, dst_block_sizes, dst_block_num, &
         &                                  datatype) BIND(C, name='xt_redist_p2p_blocks_off_new_f') RESULT(res)
      IMPORT :: xt_xmap, xt_redist, c_int
      IMPLICIT NONE
      TYPE(xt_xmap), INTENT(in) :: xmap
      INTEGER(c_int), INTENT(in) :: src_block_offsets(*)
      INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
      INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
      INTEGER(c_int), INTENT(in) :: dst_block_offsets(*)
      INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
      INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
      INTEGER, INTENT(in) :: datatype
      TYPE(xt_redist) :: res
    END FUNCTION xt_redist_p2p_blocks_off_new
  END INTERFACE

  INTERFACE
    FUNCTION xt_redist_p2p_blocks_new(xmap, src_block_sizes, src_block_num, &
         &                                  dst_block_sizes, dst_block_num, &
         &                                  datatype) BIND(C, name='xt_redist_p2p_blocks_new_f') RESULT(res)
      IMPORT :: xt_xmap, xt_redist, c_int
      IMPLICIT NONE
      TYPE(xt_xmap), INTENT(in) :: xmap
      INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
      INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
      INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
      INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
      INTEGER, INTENT(in) :: datatype
      TYPE(xt_redist) :: res
    END FUNCTION xt_redist_p2p_blocks_new
  END INTERFACE

  INTERFACE
    FUNCTION xt_redist_p2p_off_new(xmap, src_offsets, dst_offsets, datatype) BIND(C, name='xt_redist_p2p_off_new_f') RESULT(res)
      IMPORT :: xt_xmap, xt_redist, c_int
      IMPLICIT NONE
      TYPE(xt_xmap), INTENT(in) :: xmap
      INTEGER(c_int), INTENT(in) :: src_offsets(*)
      INTEGER(c_int), INTENT(in) :: dst_offsets(*)
      INTEGER, INTENT(in) :: datatype
      TYPE(xt_redist) :: res
    END FUNCTION xt_redist_p2p_off_new
  END INTERFACE

  INTERFACE
    FUNCTION xt_redist_p2p_new(xmap, comm) BIND(C, name='xt_redist_p2p_new_f') RESULT(res)
      IMPORT:: xt_xmap, xt_redist, c_int
      IMPLICIT NONE
      TYPE(xt_xmap), INTENT(in) :: xmap
      INTEGER, INTENT(in) :: comm
      TYPE(xt_redist) :: res
    END FUNCTION xt_redist_p2p_new
  END INTERFACE

  INTERFACE
    SUBROUTINE xt_redist_delete(redist) BIND(C, name='xt_redist_delete_f')
      IMPORT :: xt_redist
      IMPLICIT NONE
      TYPE(xt_redist), INTENT(inout) :: redist
    END SUBROUTINE xt_redist_delete
  END INTERFACE


  INTERFACE
    SUBROUTINE xt_redist_s_exchange1(redist, src_data_cptr, dst_data_cptr) BIND(C, name='xt_redist_s_exchange1_f')
      IMPORT:: xt_redist, c_ptr
      TYPE(xt_redist), INTENT(in) :: redist
      TYPE(c_ptr), INTENT(in) :: src_data_cptr
      TYPE(c_ptr), INTENT(in) :: dst_data_cptr
    END SUBROUTINE xt_redist_s_exchange1
  END INTERFACE

  INTERFACE
    FUNCTION xt_idxstripes_new(stripes, num_stripes) BIND(C, name='xt_idxstripes_new_f') RESULT(res)
      IMPORT:: Xt_idxlist, xt_stripe, xt_count_kind
      IMPLICIT NONE
      TYPE(xt_stripe), INTENT(in) :: stripes(*)
      INTEGER(xt_count_kind), VALUE, INTENT(in) :: num_stripes
      TYPE(Xt_idxlist) :: res
    END FUNCTION xt_idxstripes_new
  END INTERFACE

CONTAINS

  SUBROUTINE  xt_abort_cmsl(comm_f, msg, source, line)
    INTEGER, INTENT(in):: comm_f
    CHARACTER(len=*), INTENT(in) :: msg
    CHARACTER(len=*), INTENT(in) :: source
    INTEGER, INTENT(in) :: line
    CALL xt_abort_cmsl_f(comm_f, TRIM(msg)//c_null_char, source, line)
  END SUBROUTINE Xt_abort_cmsl

  SUBROUTINE  xt_abort_msl(msg, source, line)
    CHARACTER(len=*), INTENT(in) :: msg
    CHARACTER(len=*), INTENT(in) :: source
    INTEGER, INTENT(in) :: line
    CALL xt_abort_msl_f(TRIM(msg)//c_null_char, source, line)
  END SUBROUTINE Xt_abort_msl

END MODULE yaxt
