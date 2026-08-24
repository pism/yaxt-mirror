!>
!! @file xt_exchanger_f.f90
!! @brief Fortran interface to yaxt xmap declarations
!!
!! @copyright Copyright  (C)  2024 Jörg Behrens <behrens@dkrz.de>
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
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#include "fc_feature_defs.inc"
MODULE xt_exchanger_abstract
  USE iso_c_binding, ONLY: c_int, c_ptr, c_null_ptr, &
       c_associated
  USE xt_core, ONLY: xt_mpi_fint_kind, xt_pos_ext, i2, i4, i8
  USE xt_mpi, ONLY: mpi_address_kind
  USE xt_config_f, ONLY: xt_config, xt_config_f2c
  USE xt_requests, ONLY: xt_request
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xt_exchanger_c2f, xt_exchanger_f2c, xt_is_null
  PUBLIC :: xt_exchanger_iter_f2c, &
       xt_exchanger_iter_c2f
  PUBLIC :: xt_exchanger_copy, xt_exchanger_delete, &
       xt_exchanger_get_num_msg
  PUBLIC :: xt_exchanger_new
  PUBLIC :: xt_exchanger_s_exchange, xt_exchanger_a_exchange
  PUBLIC :: xt_exchanger_get_iterator, &
       xt_exchanger_iter_next, &
       xt_exchanger_iter_get_rank, &
       xt_exchanger_iter_get_direction, &
       xt_exchanger_iter_get_mpi_datatype, &
       xt_exchanger_iter_get_displacement, &
       xt_exchanger_iter_get_count, &
       xt_exchanger_iter_delete


  ! note: this type must not be extended to contain any other
  ! components, its memory pattern has to match void * exactly, which
  ! it does because of C constraints
  TYPE, BIND(C), PUBLIC :: xt_exchanger
#ifndef __G95__
    PRIVATE
#endif
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_exchanger

  TYPE, BIND(c), PUBLIC :: xt_exchanger_iter
#ifndef __G95__
    PRIVATE
#endif
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_exchanger_iter

  INTERFACE
    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_exchanger_f2c(exchanger) BIND(c, name='xt_exchanger_f2c') &
         RESULT(p)
      IMPORT :: c_ptr, xt_exchanger
      IMPLICIT NONE
      TYPE(xt_exchanger), INTENT(in) :: exchanger
      TYPE(c_ptr) :: p
    END FUNCTION xt_exchanger_f2c

    SUBROUTINE xt_exchanger_delete_c(exchanger) BIND(C, name='xt_exchanger_delete')
      IMPORT :: c_ptr
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in) :: exchanger
    END SUBROUTINE xt_exchanger_delete_c

    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_exchanger_iter_f2c(iter) BIND(c, name='xt_exchanger_iter_f2c') &
         RESULT(p)
      IMPORT :: c_ptr, xt_exchanger_iter
      IMPLICIT NONE
      TYPE(xt_exchanger_iter), INTENT(in) :: iter
      TYPE(c_ptr) :: p
    END FUNCTION xt_exchanger_iter_f2c

    SUBROUTINE xt_exchanger_iter_delete_c(iter) &
         BIND(C, name='xt_exchanger_iter_delete')
      IMPORT :: c_ptr
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in) :: iter
    END SUBROUTINE xt_exchanger_iter_delete_c

  END INTERFACE

  INTERFACE xt_exchanger_delete
    MODULE PROCEDURE xt_exchanger_delete_1
    MODULE PROCEDURE xt_exchanger_delete_a1d
  END INTERFACE xt_exchanger_delete

  INTERFACE xt_is_null
    MODULE PROCEDURE xt_exchanger_is_null
    MODULE PROCEDURE xt_exchanger_iter_is_null
  END INTERFACE xt_is_null

  INTERFACE
    FUNCTION xt_exchanger_iter_get_mpi_datatype(iter) RESULT(dt) &
         BIND(c, name='xt_exchanger_iter_get_MPI_datatype_f2c')
      IMPORT :: xt_mpi_fint_kind, xt_exchanger_iter
      TYPE(xt_exchanger_iter), INTENT(in) :: iter
      INTEGER(xt_mpi_fint_kind) :: dt
    END FUNCTION xt_exchanger_iter_get_mpi_datatype
  END INTERFACE

  TYPE, BIND(c), PUBLIC :: xt_msg_param
    INTEGER(mpi_address_kind) :: displacement
    INTEGER(xt_mpi_fint_kind) :: datatype, rank, count
  END TYPE xt_msg_param

  CHARACTER(len=*), PARAMETER :: filename = 'xt_exchanger_f.f90'
CONTAINS

  FUNCTION xt_exchanger_is_null(exchanger) RESULT(p)
    TYPE(xt_exchanger), INTENT(in) :: exchanger
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(exchanger%cptr)
  END FUNCTION xt_exchanger_is_null


  FUNCTION xt_exchanger_c2f(exchanger) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: exchanger
    TYPE(xt_exchanger) :: p
    p%cptr = exchanger
  END FUNCTION xt_exchanger_c2f

  FUNCTION xt_exchanger_copy(exchanger) RESULT(exchanger_copy)
    TYPE(xt_exchanger), INTENT(in) :: exchanger
    TYPE(xt_exchanger) :: exchanger_copy
    INTERFACE
      FUNCTION xt_exchanger_copy_c(exchanger) &
           BIND(C, name='xt_exchanger_copy') RESULT(res_ptr)
        IMPORT :: xt_exchanger, c_ptr
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: exchanger
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_exchanger_copy_c
    END INTERFACE
    exchanger_copy%cptr = xt_exchanger_copy_c(exchanger%cptr)
  END FUNCTION xt_exchanger_copy

  SUBROUTINE xt_exchanger_delete_1(exchanger)
    TYPE(xt_exchanger), INTENT(inout) :: exchanger
    CALL xt_exchanger_delete_c(exchanger%cptr)
    exchanger%cptr = c_null_ptr
  END SUBROUTINE xt_exchanger_delete_1

  SUBROUTINE xt_exchanger_delete_a1d(exchangers)
    TYPE(xt_exchanger), INTENT(inout) :: exchangers(:)
    INTEGER :: i, n
    n = SIZE(exchangers)
    DO i = 1, n
      CALL xt_exchanger_delete_c(exchangers(i)%cptr)
      exchangers(i)%cptr = c_null_ptr
    END DO
  END SUBROUTINE xt_exchanger_delete_a1d

  FUNCTION xt_exchanger_new(nsend, nrecv, send_msgs, recv_msgs, &
       comm, tag_offset, config) RESULT(exchanger)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nsend, nrecv, comm, tag_offset
    TYPE(xt_msg_param), INTENT(in) :: send_msgs(nsend), recv_msgs(nrecv)
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_exchanger) :: exchanger

    INTEGER(c_int) :: nsend_c, nrecv_c
    TYPE(c_ptr) :: config_cptr

    INTERFACE
      FUNCTION xt_exchanger_new_f2c(nsend, nrecv, send_msgs, recv_msgs, &
           comm, tag_offset, config) &
           BIND(C, name='xt_exchanger_new_f2c') RESULT(exchanger_ptr)
        IMPORT :: c_int, xt_exchanger, xt_mpi_fint_kind, c_ptr, &
             xt_msg_param
        IMPLICIT NONE
        INTEGER(c_int), VALUE, INTENT(in) :: nsend, nrecv
        TYPE(xt_msg_param), INTENT(in) :: send_msgs(nsend), recv_msgs(nrecv)
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm, tag_offset
        TYPE(c_ptr), VALUE, INTENT(in) :: config
        TYPE(c_ptr) :: exchanger_ptr
      END FUNCTION xt_exchanger_new_f2c
    END INTERFACE
    nsend_c = INT(nsend, c_int)
    nrecv_c = INT(nrecv, c_int)
    config_cptr = xt_config_f2c(config)
    exchanger%cptr = xt_exchanger_new_f2c(nsend_c, nrecv_c, &
         send_msgs, recv_msgs, comm, tag_offset, config_cptr)
  END FUNCTION xt_exchanger_new

  FUNCTION xt_exchanger_get_num_msg(exchanger, direction) RESULT(num)
    TYPE(xt_exchanger), INTENT(in) :: exchanger
    INTEGER, INTENT(in) :: direction
    INTEGER :: num
    INTERFACE
      FUNCTION xt_exchanger_get_num_msg_c(exchanger, direction) RESULT(num) &
           BIND(c, name='xt_exchanger_get_num_msg')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: exchanger
        INTEGER(c_int), VALUE, INTENT(in) :: direction
        INTEGER(c_int) :: num
      END FUNCTION xt_exchanger_get_num_msg_c
    END INTERFACE
    num = INT(xt_exchanger_get_num_msg_c(exchanger%cptr, direction))
  END FUNCTION xt_exchanger_get_num_msg

  SUBROUTINE xt_exchanger_s_exchange(exchanger, src_data, dst_data)
    TYPE(xt_exchanger), INTENT(inout) :: exchanger
    TYPE(c_ptr) :: src_data, dst_data
    INTERFACE
      SUBROUTINE xt_exchanger_s_exchange_c(exchanger, src_data, dst_data) &
           BIND(c, name='xt_exchanger_s_exchange')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: exchanger
        TYPE(c_ptr), VALUE :: src_data, dst_data
      END SUBROUTINE xt_exchanger_s_exchange_c
    END INTERFACE
    CALL xt_exchanger_s_exchange_c(exchanger%cptr, src_data, dst_data)
  END SUBROUTINE xt_exchanger_s_exchange

  SUBROUTINE xt_exchanger_a_exchange(exchanger, src_data, dst_data, request)
    TYPE(xt_exchanger), INTENT(inout) :: exchanger
    TYPE(xt_request) :: request
    TYPE(c_ptr) :: src_data, dst_data
    INTERFACE
      SUBROUTINE xt_exchanger_a_exchange_c(exchanger, src_data, dst_data, &
           request) &
           BIND(c, name='xt_exchanger_a_exchange')
        IMPORT :: c_ptr, xt_request
        TYPE(c_ptr), VALUE, INTENT(in) :: exchanger
        TYPE(c_ptr), VALUE :: src_data, dst_data
        TYPE(xt_request), INTENT(out) :: request
      END SUBROUTINE xt_exchanger_a_exchange_c
    END INTERFACE
    CALL xt_exchanger_a_exchange_c(exchanger%cptr, src_data, dst_data, request)
  END SUBROUTINE xt_exchanger_a_exchange

  FUNCTION xt_exchanger_get_iterator(exchanger, config) RESULT(iter)
    TYPE(xt_exchanger), INTENT(in) :: exchanger
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_exchanger_iter) :: iter
    INTERFACE
      FUNCTION xt_exchanger_get_iterator_c(exchanger, config) RESULT(cptr) &
           BIND(c, name='xt_exchanger_get_iterator')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: exchanger, config
        TYPE(c_ptr) :: cptr
      END FUNCTION xt_exchanger_get_iterator_c
    END INTERFACE
    iter%cptr = xt_exchanger_get_iterator_c(exchanger%cptr, &
         xt_config_f2c(config))
  END FUNCTION xt_exchanger_get_iterator

  FUNCTION xt_exchanger_iter_c2f(iter) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: iter
    TYPE(xt_exchanger_iter) :: p
    p%cptr = iter
  END FUNCTION xt_exchanger_iter_c2f

  FUNCTION xt_exchanger_iter_is_null(iter) RESULT(p)
    TYPE(xt_exchanger_iter), INTENT(in) :: iter
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(iter%cptr)
  END FUNCTION xt_exchanger_iter_is_null

  FUNCTION xt_exchanger_iter_next(iter) RESULT(avail)
    TYPE(xt_exchanger_iter), INTENT(inout) :: iter
    LOGICAL :: avail
    INTERFACE
      FUNCTION xt_exchanger_iter_next_c(iter) RESULT(avail) &
           BIND(c, name='xt_exchanger_iter_next')
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        INTEGER(c_int) :: avail
      END FUNCTION xt_exchanger_iter_next_c
    END INTERFACE
    avail = xt_exchanger_iter_next_c(iter%cptr) /= 0
  END FUNCTION xt_exchanger_iter_next

  FUNCTION xt_exchanger_iter_get_rank(iter) RESULT(rank)
    TYPE(xt_exchanger_iter), INTENT(in) :: iter
    INTEGER :: rank
    INTERFACE
      FUNCTION xt_exchanger_iter_get_rank_c(iter) RESULT(rank) &
           BIND(c, name='xt_exchanger_iter_get_rank')
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        INTEGER(c_int) :: rank
      END FUNCTION xt_exchanger_iter_get_rank_c
    END INTERFACE
    rank = INT(xt_exchanger_iter_get_rank_c(iter%cptr))
  END FUNCTION xt_exchanger_iter_get_rank

  FUNCTION xt_exchanger_iter_get_count(iter) RESULT(msg_count)
    TYPE(xt_exchanger_iter), INTENT(in) :: iter
    INTEGER :: msg_count
    INTERFACE
      FUNCTION xt_exchanger_iter_get_count_c(iter) RESULT(msg_count) &
           BIND(c, name='xt_exchanger_iter_get_count')
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        INTEGER(c_int) :: msg_count
      END FUNCTION xt_exchanger_iter_get_count_c
    END INTERFACE
    msg_count = INT(xt_exchanger_iter_get_count_c(iter%cptr))
  END FUNCTION xt_exchanger_iter_get_count

  FUNCTION xt_exchanger_iter_get_direction(iter) RESULT(direction)
    TYPE(xt_exchanger_iter), INTENT(in) :: iter
    INTEGER :: direction
    INTERFACE
      FUNCTION xt_exchanger_iter_get_direction_c(iter) RESULT(direction) &
           BIND(c, name='xt_exchanger_iter_get_direction')
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        INTEGER(c_int) :: direction
      END FUNCTION xt_exchanger_iter_get_direction_c
    END INTERFACE
    direction = INT(xt_exchanger_iter_get_direction_c(iter%cptr))
  END FUNCTION xt_exchanger_iter_get_direction

  FUNCTION xt_exchanger_iter_get_displacement(iter) RESULT(displacement)
    TYPE(xt_exchanger_iter), INTENT(in) :: iter
    INTEGER(mpi_address_kind) :: displacement
    INTERFACE
      FUNCTION xt_exchanger_iter_get_displacement_c(iter) RESULT(displacement) &
           BIND(c, name='xt_exchanger_iter_get_displacement')
        IMPORT :: c_ptr, mpi_address_kind
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        INTEGER(mpi_address_kind) :: displacement
      END FUNCTION xt_exchanger_iter_get_displacement_c
    END INTERFACE
    displacement = xt_exchanger_iter_get_displacement_c(iter%cptr)
  END FUNCTION xt_exchanger_iter_get_displacement

  SUBROUTINE xt_exchanger_iter_delete(iter)
    TYPE(xt_exchanger_iter), INTENT(inout) :: iter
    INTERFACE
      SUBROUTINE xt_exchanger_iter_delete_c(iter) &
           BIND(c, name='xt_exchanger_iter_delete')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
      END SUBROUTINE xt_exchanger_iter_delete_c
    END INTERFACE
    CALL xt_exchanger_iter_delete_c(iter%cptr)
    iter%cptr = c_null_ptr
  END SUBROUTINE xt_exchanger_iter_delete
END MODULE xt_exchanger_abstract

!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
