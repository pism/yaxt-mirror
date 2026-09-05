!>
!! @file test_exchanger_parallel_f.f90
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
#define FP_EQ(a, b) (ABS((a)-(b)) < EPSILON((a)))
#ifdef VERBOSE
#define PUT_ERR(msg) WRITE (0, '(a)') msg ; err_count__ = err_count__ + 1
#else
#define PUT_ERR(msg) err_count__ = err_count__ + 1
#endif
PROGRAM test_exchanger_parallel
  USE iso_c_binding, ONLY: c_ptr, c_loc, c_int, c_null_ptr
  USE mpi
  USE yaxt
  USE xt_exchanger_abstract
  USE ftest_common
  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER :: filename = 'test_exchanger_parallel_f.f90'
  TYPE test_message
    INTEGER :: rank       !< rank of communication partner
    INTEGER, ALLOCATABLE :: pos(:)   !< positions to be sent/received
  END TYPE test_message

  TYPE(xt_config) :: config
  INTEGER :: test_freq = 3
  INTEGER, PARAMETER :: num_mt_modes = 1 &
#ifdef _OPENMP
       + 1 &
#endif
       + 0

  INTEGER, PARAMETER :: mt_configs_2_test(num_mt_modes) = (/ XT_MT_NONE &
#ifdef _OPENMP
    , XT_MT_OPENMP &
#endif
    /)
#if XT_ENABLE_ALTERNATIVE_EXCHANGER_SEND_MODES
  INTEGER, PARAMETER :: num_send_modes = 3
#else
  INTEGER, PARAMETER :: num_send_modes = 1
#endif
  INTEGER, PARAMETER :: send_comm_modes_2_test(num_send_modes) = (/ &
       xt_exch_send_standard &
#if XT_ENABLE_ALTERNATIVE_EXCHANGER_SEND_MODES
       , xt_exch_send_buffered, &
       xt_exch_send_synchronous &
#endif
       /)
  INTEGER :: comm
  INTEGER :: m, i, sm, split_rank, world_size_2
#ifdef _OPENMP
  INTEGER :: thread_support_provided
#endif
  INTEGER, ALLOCATABLE :: exchanger_code(:)
  INTEGER, ALLOCATABLE :: split_comm(:)
  INTEGER :: world_size, world_rank, buf_size, buf_size_, ierror
  INTEGER :: no_exch_id
  INTEGER(c_int), ALLOCATABLE :: buf(:)
  TYPE(c_ptr) :: bufp

  comm = mpi_comm_world

  CALL init_mpi
  CALL xt_initialize(comm)

  ! check that non-existent exchanger results in expected return
  no_exch_id = xt_exchanger_id_by_name("non-existant fantasy exchanger")
  IF (no_exch_id /= -1) &
       CALL test_abort("invalid exchanger", filename, __LINE__)
  CALL parse_options(exchanger_code, test_freq)
  config = xt_config_new()
  ! todo: replace computation with query of exchanger
  CALL mpi_comm_size(comm, world_size, ierror)
  CALL mpi_comm_rank(comm, world_rank, ierror)
  buf_size = (4 + mpi_bsend_overhead) * world_size
  ALLOCATE(buf((buf_size+3)/4))
  CALL mpi_buffer_attach(buf, buf_size, ierror)
  world_size_2 = world_size/2
  ALLOCATE(split_comm(0:world_size_2))
  split_comm(0) = mpi_comm_world
  DO split_rank = 1, world_size_2
    CALL mpi_comm_split(mpi_comm_world, &
         MERGE(1, 0, world_rank < split_rank), world_rank, &
         split_comm(split_rank), ierror)
  END DO
  DO m = 1, num_mt_modes
#ifdef _OPENMP
    CALL mpi_query_thread(thread_support_provided, ierror)
    IF (mt_configs_2_test(m) == XT_MT_OPENMP &
         .AND. thread_support_provided /= MPI_THREAD_MULTIPLE) CYCLE
#endif
    DO split_rank = 0, world_size_2
      comm = split_comm(split_rank)
      CALL xt_config_set_redist_mthread_mode(config, mt_configs_2_test(m))
      DO sm = 1, num_send_modes
        CALL xt_config_set_exchange_send_comm_mode(config, &
             send_comm_modes_2_test(sm))
        DO i = 1, SIZE(exchanger_code)
          CALL xt_config_set_exchange_method(config, exchanger_code(i))

          CALL test_bcast(comm, config)

          CALL test_gather(comm, config)

          CALL test_all2all(comm, config)

          CALL test_rr(comm, config)

          CALL test_intercomm_all2all(comm, config)
        END DO
      END DO
    END DO
  END DO
  DO split_rank = 1, world_size_2
    CALL mpi_comm_free(split_comm(split_rank), ierror)
  END DO

  CALL mpi_buffer_detach(bufp, buf_size_, ierror)
  CALL xt_config_delete(config)
  DEALLOCATE(exchanger_code, buf, split_comm)

  CALL xt_finalize
  CALL finish_mpi

CONTAINS

  SUBROUTINE parse_options(exchanger_codes, test_freq)
    INTEGER, ALLOCATABLE, INTENT(out) :: exchanger_codes(:)
    INTEGER, INTENT(inout) :: test_freq
    INTEGER, PARAMETER :: num_exchanger = 3 &
#ifdef XT_ENABLE_DDT_EXCHANGER
         + 1 &
#endif
#ifdef XT_CAN_USE_MPI_NEIGHBOR_ALLTOALL
         + 1 &
#endif
         + 1
    INTEGER, PARAMETER :: exchanger_table(num_exchanger) = (/ &
         xt_exchanger_irecv_send, &
         xt_exchanger_irecv_isend, &
         xt_exchanger_irecv_isend_packed, &
#ifdef XT_ENABLE_DDT_EXCHANGER
         xt_exchanger_irecv_isend_ddt_packed, &
#endif
#ifdef XT_CAN_USE_MPI_NEIGHBOR_ALLTOALL
         xt_exchanger_neigh_alltoall, &
#endif
         xt_exchanger_mix_isend_irecv /)
    LOGICAL :: selected_exchanger(num_exchanger)
    INTEGER :: j, i, num_cmd_args, arg_len, opt_len, &
         exchanger_new_id, freq, ierror, num_sel_exch
#ifndef XT_ENABLE_DDT_EXCHANGER
    INTEGER :: my_world_rank
#endif
    INTEGER, PARAMETER :: max_opt_arg_len = 80
    CHARACTER(max_opt_arg_len) :: optarg, opt

    selected_exchanger = .FALSE.
    selected_exchanger(xt_exchanger_mix_isend_irecv) = .TRUE.
    num_cmd_args = COMMAND_ARGUMENT_COUNT()
    i = 1
    DO WHILE (i < num_cmd_args)
      CALL GET_COMMAND_ARGUMENT(i, opt, opt_len)
      IF (opt_len == 2 .AND. opt(1:1) == '-') THEN
        IF (i < num_cmd_args) THEN
          CALL GET_COMMAND_ARGUMENT(i + 1, optarg, arg_len)
        ELSE
          WRITE (0, '(4a)') &
               'unexpected command-line argument parsing error: ', &
               opt(1:opt_len), ' ', optarg(1:arg_len)
          FLUSH(0)
          CALL test_abort('unexpected command-line argument', &
               filename, __LINE__)
        END IF
        SELECT CASE (opt(2:2))
        CASE('m')
          IF (arg_len > max_opt_arg_len) &
               CALL test_abort('incorrect argument to command-line option -m', &
               filename, __LINE__)
          exchanger_new_id &
               = xt_exchanger_id_by_name(optarg(1:LEN_TRIM(optarg)))
          IF (exchanger_new_id == -1) THEN
            WRITE (0, '(2a)') &
                 'unexpected command-line argument parsing error: ', &
                 optarg(1:arg_len)
            FLUSH(0)
            CALL test_abort('unexpected command-line argument', &
                 filename, __LINE__)
#ifndef XT_CAN_USE_MPI_NEIGHBOR_ALLTOALL
          ELSE IF (exchanger_new_id == xt_exchanger_neigh_alltoall) THEN
            WRITE (0, '(a)') &
                 "xt_exchanger_neigh_alltoall_new requires MPI version 3.0 or &
                 &higher"
            i = i + 2
            CYCLE
#endif
#ifndef XT_ENABLE_DDT_EXCHANGER
          ELSE IF (exchanger_new_id == xt_exchanger_irecv_isend_ddt_packed) THEN
            CALL mpi_comm_rank(mpi_comm_world, my_world_rank, ierror)
            IF (my_world_rank == 0) THEN
              WRITE (0, '(a)') &
                   "notice: xt_exchanger_irecv_isend_ddt_packed disabled."
            END IF
            i = i + 2
            CYCLE
#endif
          ELSE
            DO j = 1, num_exchanger
              IF (exchanger_table(j) == exchanger_new_id) THEN
                selected_exchanger(j) = .TRUE.
                EXIT
              END IF
            END DO
          END IF
        CASE('-s')
          READ (optarg(1:arg_len), fmt=*, iostat=ierror) freq
          IF (ierror /= 0) THEN
            WRITE (0, '(3a,i0)') &
                 'unexpected command-line argument parsing error: ', &
                 optarg(1:arg_len), ' ', ierror
            FLUSH(0)
            CALL test_abort('unexpected command-line argument', &
                 filename, __LINE__)
          END IF
          test_freq = freq
        END SELECT
      ELSE
        WRITE (0, '(2a)') 'unknown argument: ', opt
        CALL test_abort('unexpected command-line argument', &
             filename, __LINE__)
      END IF
      i = i + 2
    END DO
    num_sel_exch = COUNT(selected_exchanger)
    ALLOCATE(exchanger_codes(num_sel_exch))
    exchanger_codes = PACK(exchanger_table, selected_exchanger)
  END SUBROUTINE parse_options

  SUBROUTINE data_exchange(exchanger, src_data_p, dst_data_p, async)
    TYPE(xt_exchanger), INTENT(inout) :: exchanger
    TYPE(c_ptr), INTENT(in) :: src_data_p, dst_data_p
    LOGICAL, INTENT(in) :: async

    TYPE(xt_request) :: request
    LOGICAL :: request_finished

    IF (async) THEN
      CALL xt_exchanger_a_exchange(exchanger, src_data_p, &
           dst_data_p, request)
      CALL xt_request_test(request, request_finished)
      CALL xt_request_wait(request)
      CALL xt_request_test(request, request_finished)
      IF (.NOT. request_finished) &
           CALL test_abort("invalid flag result", filename, __LINE__)
    ELSE
      CALL xt_exchanger_s_exchange(exchanger, src_data_p, dst_data_p)
    END IF
  END SUBROUTINE data_exchange

  SUBROUTINE test_bcast(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config

    TYPE(xt_msg_param), ALLOCATABLE :: send_msgs(:)
    TYPE(xt_msg_param) :: recv_msgs(1)
    TYPE(xt_exchanger) :: exchanger, exchanger_copy
    TYPE(c_ptr) :: src_data_p, dst_data_p
    INTEGER, TARGET :: src_data(1), dst_data(1)
    INTEGER :: my_rank, comm_size, ierror
    INTEGER :: incr, i, j, nsend, nrecv, exchange_method, &
         test_async, test_async_ubound

    CALL mpi_comm_rank(comm, my_rank, ierror)
    CALL mpi_comm_size(comm, comm_size, ierror)
    exchange_method = xt_config_get_exchange_method(config)
    src_data_p = C_LOC(src_data)
#if defined __GNUC__ && ( __GNUC__ >= 5 && ( __GNUC__ <= 7 \
    || __GNUC__ == 8 && __GNUC_MINOR__ < 4 ) \
    || __GNUC__ == 4 && __GNUC_MINOR__ == 6 )
    CALL xt_slice_c_loc(dst_data, dst_data_p)
#else
    dst_data_p = C_LOC(dst_data)
#endif

    ! bcast pattern
    incr = comm_size/test_freq + MERGE(1, 0, comm_size < test_freq)
    DO i = 0, comm_size-1, incr
      ! setup
      nsend = 0
      nrecv = MERGE(1, 0, my_rank /= i)

      IF (my_rank == i) THEN
        nsend = comm_size - 1
        ALLOCATE(send_msgs(nsend))
        DO j = 0, i-1
          send_msgs(j+1) = xt_msg_param(rank=j, &
               datatype=mpi_integer, displacement=0_mpi_address_kind, count=1)
        END DO
        DO j = i, nsend-1
          send_msgs(j+1) = xt_msg_param(rank=MOD(1+j, comm_size), &
               datatype=mpi_integer, displacement=0_mpi_address_kind, count=1)
        END DO
        recv_msgs(1) = xt_msg_param(rank=-1, datatype=mpi_datatype_null, &
             displacement=0_mpi_address_kind, count=1)
      ELSE
        ALLOCATE(send_msgs(0))
        recv_msgs(1) = xt_msg_param(rank=i, datatype=mpi_integer, &
           displacement=0_mpi_address_kind, count=1)
      END IF

      exchanger = xt_exchanger_new(nsend, nrecv, send_msgs, &
           recv_msgs, comm, 0, config)

      IF (xt_is_null(exchanger)) THEN
        PUT_ERR("unexpected NULL result")
      END IF
      ! run the test
      test_async_ubound &
           = MERGE(1, 2, exchange_method == xt_exchanger_irecv_send)
      DO test_async = 1, test_async_ubound
        src_data(1) = MERGE(4711, -1, my_rank == i)
        dst_data(1) = MERGE(4711, -1, my_rank == i)
        CALL data_exchange(exchanger, src_data_p, &
               dst_data_p, test_async > 1)
        IF (dst_data(1) /= 4711) &
             CALL test_abort("invalid data", filename, __LINE__)
        exchanger_copy = xt_exchanger_copy(exchanger, comm, 0)
        CALL xt_exchanger_delete(exchanger)
        exchanger = exchanger_copy
      END DO

      ! cleanup
      DEALLOCATE(send_msgs)
      CALL xt_exchanger_delete(exchanger)
    END DO
  END SUBROUTINE test_bcast

  SUBROUTINE test_gather(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config

    INTEGER, PARAMETER :: num_exchanges=2
    TYPE(xt_exchanger) :: exchanger(2)
    TYPE(xt_msg_param) :: send_msgs(1)
    TYPE(xt_msg_param), ALLOCATABLE :: recv_msgs(:,:)
    TYPE(c_ptr) :: src_data_p, dst_data_p
    INTEGER, ALLOCATABLE :: dt_by_ofs(:)
    REAL, ALLOCATABLE, TARGET :: dst_data(:, :)

    INTEGER :: my_rank, comm_size, ierror
    INTEGER :: incr, i, j, nsend, nrecv, exchange_method, &
         test_async, test_async_ubound
    LOGICAL :: no_mismatch
    REAL, TARGET :: src_data(1)

    CALL mpi_comm_rank(comm, my_rank, ierror)
    CALL mpi_comm_size(comm, comm_size, ierror)
    exchange_method = xt_config_get_exchange_method(config)
    src_data_p = C_LOC(src_data)

    ! gather pattern
    ! prepare datatypes outside of loop for load-balance
    ALLOCATE(dt_by_ofs(comm_size))
    dt_by_ofs(1) = MPI_REAL
    DO j = 2, comm_size
      CALL mpi_type_indexed(1, (/ 1 /), (/ j-1 /), MPI_REAL, dt_by_ofs(j), &
           ierror)
      CALL mpi_type_commit(dt_by_ofs(j), ierror)
    END DO
    ALLOCATE(dst_data(comm_size-1, 2))
    incr = comm_size/test_freq + MERGE(1, 0, comm_size < test_freq)
    DO i = 0, comm_size-1, incr
      ! setup
      nsend = MERGE(1, 0, i /= my_rank)
      nrecv = 0


      send_msgs(1) = xt_msg_param(rank=i, datatype=MPI_REAL, &
           displacement=0_mpi_address_kind, count=1)
      ALLOCATE(recv_msgs(MERGE(comm_size-1, 1, my_rank == i), 2))
      IF (my_rank == i) THEN
        nrecv = comm_size - 1
        DO j = 1, nrecv
          recv_msgs(j, 1)%rank = MOD(i + j, comm_size)
          recv_msgs(j, 1)%datatype = dt_by_ofs(j)
          recv_msgs(j, 1)%displacement = 0_mpi_address_kind
          recv_msgs(j, 1)%count=1
          recv_msgs(j, 2)%rank &
               = MOD(comm_size - j + i, comm_size)
          recv_msgs(j, 2)%datatype = dt_by_ofs(j)
          recv_msgs(j, 2)%displacement = 0
          recv_msgs(j, 2)%count=1
        END DO
      END IF

      DO j = 1, num_exchanges
        exchanger(j) = xt_exchanger_new(nsend, nrecv, send_msgs, &
             recv_msgs(1:nrecv, j), comm, 0, config)
      END DO

      ! run test
      test_async_ubound &
           = MERGE(1, 2, exchange_method == xt_exchanger_irecv_send)
      dst_data_p = c_null_ptr
      DO test_async = 1, test_async_ubound
        src_data(1) = REAL(MOD(my_rank+comm_size-i, comm_size))
        dst_data = -1.0
        DO j = 1, num_exchanges
          IF (nrecv > 0) dst_data_p = C_LOC(dst_data(1, j))
          CALL data_exchange(exchanger(j), src_data_p, &
               dst_data_p, test_async > 1)
        END DO
      END DO

      no_mismatch = .TRUE.
      DO j = 1, nrecv
        no_mismatch = no_mismatch .AND. FP_EQ(dst_data(j, 1), REAL(j)) &
             .AND. FP_EQ(dst_data(j, 2), REAL(comm_size - j))
      END DO
      IF (no_mismatch) THEN
        CONTINUE
      ELSE
        CALL test_abort("invalid data", filename, __LINE__)
      END IF
      ! cleanup
      DEALLOCATE(recv_msgs)
      CALL xt_exchanger_delete(exchanger)
    END DO

    DEALLOCATE(dst_data)
    DO j = 2, comm_size
      CALL mpi_type_free(dt_by_ofs(j), ierror)
    END DO
    DEALLOCATE(dt_by_ofs)
  END SUBROUTINE test_gather

  SUBROUTINE test_all2all(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config

    TYPE(xt_msg_param), ALLOCATABLE :: send_msgs(:), recv_msgs(:)
    TYPE(xt_exchanger) :: exchanger
    DOUBLE PRECISION, ALLOCATABLE, TARGET :: dst_data(:)
    DOUBLE PRECISION, TARGET :: src_data
    TYPE(c_ptr) :: src_data_p, dst_data_p
    INTEGER(mpi_address_kind) :: extent_of_dbl, lb
    INTEGER :: my_rank, comm_size, ierror
    INTEGER :: incr, i, j, k, nsend, nrecv, exchange_method, &
         test_async, test_async_ubound, ofs
    LOGICAL :: no_mismatch

    CALL mpi_comm_rank(comm, my_rank, ierror)
    CALL mpi_comm_size(comm, comm_size, ierror)
    exchange_method = xt_config_get_exchange_method(config)
    src_data_p = C_LOC(src_data)
    ALLOCATE(dst_data(comm_size))
    dst_data_p = C_LOC(dst_data)
    src_data = DBLE(my_rank)

    ! all-to-all pattern
    ! setup
    nsend = comm_size - 1
    nrecv = comm_size - 1

    ALLOCATE(send_msgs(nsend), recv_msgs(nrecv))
    DO j = 1, nsend
      send_msgs(j)%datatype = mpi_double_precision
      send_msgs(j)%displacement = 0_mpi_address_kind
      send_msgs(j)%count = 1
    END DO
    CALL mpi_type_get_extent(mpi_double_precision, lb, extent_of_dbl, ierror)
    incr = comm_size/test_freq + MERGE(1, 0, comm_size < test_freq)
    DO i = 0, comm_size - 2, incr
      DO j = 1, nsend
        ofs = my_rank + i + j
        send_msgs(j)%rank = MOD(ofs + MERGE(1, 0, ofs >= comm_size + my_rank),&
             comm_size)
      END DO
      DO j = 0, comm_size - 2, incr
        DO k = 1, nrecv
          ofs = MOD(i + j + k-1, comm_size - 1)
          ofs = ofs + MERGE(1, 0, ofs >= my_rank)
          recv_msgs(k)%rank = ofs
          recv_msgs(k)%datatype = mpi_double_precision
          recv_msgs(k)%displacement &
               = extent_of_dbl * INT(ofs, mpi_address_kind)
          recv_msgs(k)%count = 1
        END DO
        exchanger = xt_exchanger_new(nsend, nrecv, send_msgs, &
             recv_msgs, comm, 0, config)
        ! test
        test_async_ubound &
           = MERGE(1, 2, exchange_method == xt_exchanger_irecv_send)
        DO test_async = 1, test_async_ubound
          dst_data(:) = DBLE(my_rank)

          CALL data_exchange(exchanger, src_data_p, &
               dst_data_p, test_async > 1)
          no_mismatch = .TRUE.
          DO k = 1, comm_size
            no_mismatch = no_mismatch .AND. FP_EQ(dst_data(k), DBLE(k-1))
          END DO
          IF (no_mismatch) THEN
            CONTINUE
          ELSE
            CALL test_abort("invalid data", filename, __LINE__)
          END IF
        END DO
        ! cleanup
        CALL xt_exchanger_delete(exchanger)
      END DO
    END DO
    DEALLOCATE(dst_data, recv_msgs, send_msgs)
  END SUBROUTINE test_all2all

  SUBROUTINE test_rr(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config

    INTEGER, PARAMETER :: nrecv = 1, nsend = 1
    TYPE(xt_msg_param) :: send_msgs(nsend), recv_msgs(nrecv)
    TYPE(xt_exchanger) :: exchanger
    TYPE(xt_exchanger_iter) :: msg_iter
    DOUBLE PRECISION, TARGET :: src_data(1), dst_data(1)
    TYPE(c_ptr) :: src_data_p, dst_data_p
    INTEGER :: my_rank, comm_size, ierror
    INTEGER :: incr, i, j, exchange_method, &
         test_async, test_async_ubound
    INTEGER :: iter_dt, msg_rank, msg_count, direction, nmsg_ref(0:2), nmsg(0:2)
#if defined __INTEL_COMPILER && (__INTEL_COMPILER_BUILD_DATE == 20110811 \
        || __INTEL_COMPILER < 1400)
    INTEGER(mpi_address_kind) :: addr_dummy
#endif

    CALL mpi_comm_rank(comm, my_rank, ierror)
    CALL mpi_comm_size(comm, comm_size, ierror)
    exchange_method = xt_config_get_exchange_method(config)

    src_data(1) = DBLE(my_rank)
    src_data_p = C_LOC(src_data)
#if defined __GNUC__ && ( __GNUC__ >= 5 && ( __GNUC__ <= 7 \
    || __GNUC__ == 8 && __GNUC_MINOR__ < 4 ) \
    || __GNUC__ == 4 && __GNUC_MINOR__ == 6 )
    CALL xt_slice_c_loc(dst_data, dst_data_p)
#else
    dst_data_p = C_LOC(dst_data)
#endif

    ! round robin pattern
    incr = comm_size/test_freq + MERGE(1, 0, comm_size < test_freq)
    DO i = 1, comm_size, incr
      ! setup
      send_msgs(1) = xt_msg_param(rank=MOD(my_rank+i, comm_size), &
           datatype=mpi_double_precision, displacement = 0_mpi_address_kind, &
           count=1)
      recv_msgs(1) = xt_msg_param(rank=MOD(my_rank+comm_size-i, comm_size), &
           datatype=mpi_double_precision, displacement=0_mpi_address_kind, &
           count=1)

      exchanger = xt_exchanger_new(nsend, nrecv, send_msgs, recv_msgs, &
           comm, 0, config)

      ! test
      test_async_ubound &
           = MERGE(1, 2, exchange_method == xt_exchanger_irecv_send)
      DO test_async = 1, test_async_ubound
        dst_data(1) = DBLE(-1)

        CALL data_exchange(exchanger, src_data_p, &
             dst_data_p, test_async > 1)
#if defined __INTEL_COMPILER && (__INTEL_COMPILER_BUILD_DATE == 20110811 \
        || __INTEL_COMPILER < 1400)
        ! the following call is only needed to make the optimizer
        ! in old Intel compilers less optimistic in their (mistaken)
        ! ability to predict changes to dst_data
        CALL mpi_get_address(dst_data, addr_dummy, ierror)
#endif
        IF (FP_EQ(dst_data(1), DBLE(MOD(my_rank + comm_size - i, comm_size)))) THEN
          CONTINUE
        ELSE
          CALL test_abort("invalid data", filename, __LINE__)
        END IF
      END DO

      ! test queries
      IF (i == 1) THEN
        msg_iter = xt_exchanger_get_iterator(exchanger, config)
        IF (.NOT. xt_is_null(msg_iter)) THEN
          nmsg_ref(xt_direction_send) = nsend
          nmsg_ref(xt_direction_recv) = nrecv
          nmsg_ref(2) = nsend + nrecv
          DO j = 0,2
            nmsg(j) = xt_exchanger_get_num_msg(exchanger, j)
          END DO
          IF (ANY(nmsg /= nmsg_ref)) THEN
            PUT_ERR("unexpected message counts")
          END IF
          DO
            iter_dt = xt_exchanger_iter_get_MPI_Datatype(msg_iter)
            IF (iter_dt /= mpi_double_precision) THEN
              PUT_ERR("unexpected datatype")
            END IF
            msg_rank = xt_exchanger_iter_get_rank(msg_iter)
            msg_count = xt_exchanger_iter_get_count(msg_iter)
            IF (msg_count /= 1) THEN
              PUT_ERR("unexpected count")
            END IF
            direction = xt_exchanger_iter_get_direction(msg_iter)
            IF (direction /= 0 .AND. direction /= 1) THEN
              PUT_ERR("unexpected direction")
            END IF
            IF ((direction == 0 .AND. msg_rank /= MOD(my_rank + i, comm_size)) &
                 .OR. (direction == 1 &
                 .AND. msg_rank /= MOD(my_rank + comm_size - i, comm_size))) &
                 THEN
              PUT_ERR("unexpected communication")
            END IF
            IF (.NOT. xt_exchanger_iter_next(msg_iter)) EXIT
          END DO
          IF (xt_exchanger_iter_get_direction(msg_iter) /= -1) THEN
            PUT_ERR("unexpected direction")
          END IF
          CALL xt_exchanger_iter_delete(msg_iter)
        END IF
      END IF

      ! cleanup
      CALL xt_exchanger_delete(exchanger)
    END DO
  END SUBROUTINE test_rr

  SUBROUTINE test_intercomm_all2all(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config

    TYPE(xt_msg_param), ALLOCATABLE :: send_msgs(:), recv_msgs(:)
    TYPE(xt_exchanger) :: exchanger
    DOUBLE PRECISION, TARGET, ALLOCATABLE :: dst_data(:)
    DOUBLE PRECISION, TARGET :: src_data(1)
    DOUBLE PRECISION :: dummy(2)
    INTEGER(mpi_address_kind) :: base, next, ofs
    INTEGER :: dst_data_offset
    TYPE(c_ptr) :: src_data_p, dst_data_p
    INTEGER :: my_rank, comm_size, ierror, &
         intra_group_comm, intra_rank, inter_comm, &
         local_size, remote_size
    INTEGER :: i, nsend, nrecv, exchange_method, &
         test_async, test_async_ubound, split_rank, group
    LOGICAL :: no_mismatch

    exchange_method = xt_config_get_exchange_method(config)
#ifdef XT_CAN_USE_MPI_NEIGHBOR_ALLTOALL
    ! inter-communicators are not defined for virtual topologies, which are
    ! used by xt_exchanger_neigh_alltoall_new
    IF (exchange_method == xt_exchanger_neigh_alltoall) RETURN
#endif
    CALL mpi_comm_rank(comm, my_rank, ierror)
    CALL mpi_comm_size(comm, comm_size, ierror)

    IF (comm_size == 1) RETURN

    ! generate intercomm
    split_rank = (comm_size * 2) / 3
    group = MERGE(1, 0, my_rank >= split_rank)
    ! split communicator with 2 to 1 ratio
    CALL mpi_comm_split(comm, group, 0, intra_group_comm, ierror)
    CALL mpi_intercomm_create(intra_group_comm, 0, comm, &
         MERGE(0, split_rank, group > 0), 0, inter_comm, ierror)

    CALL mpi_comm_rank(inter_comm, intra_rank, ierror)
    CALL mpi_comm_size(inter_comm, local_size, ierror)
    CALL mpi_comm_remote_size(inter_comm, remote_size, ierror)

    ! all-to-all pattern
    ! setup
    nsend = remote_size
    nrecv = remote_size
    src_data(1) = DBLE(my_rank)
    CALL mpi_get_address(dummy(1), base, ierror)
    CALL mpi_get_address(dummy(2), next, ierror)
    ofs = next - base

    ALLOCATE(send_msgs(nsend), recv_msgs(nrecv))
    DO i = 1, nsend
      send_msgs(i)%rank = i-1
      send_msgs(i)%datatype = mpi_double_precision
      send_msgs(i)%displacement = 0_mpi_address_kind
      send_msgs(i)%count = 1
    END DO

    ALLOCATE(dst_data(nrecv))
    DO i = 1, nrecv
      recv_msgs(i)%rank = i-1
      recv_msgs(i)%datatype = mpi_double_precision
      recv_msgs(i)%displacement = ofs * INT(i-1, mpi_address_kind)
      recv_msgs(i)%count = 1
    END DO

    dst_data_p = C_LOC(dst_data)
    src_data_p = C_LOC(src_data)

    exchanger = xt_exchanger_new(nsend, nrecv, send_msgs, recv_msgs, &
         inter_comm, 0, config)
    ! test
    test_async_ubound &
         = MERGE(1, 2, exchange_method == xt_exchanger_irecv_send)
    DO test_async = 1, test_async_ubound
      dst_data(:) = -1.0d0
      CALL data_exchange(exchanger, src_data_p, &
           dst_data_p, test_async > 1)
      dst_data_offset = MERGE(-1, split_rank-1, my_rank >= split_rank)
      no_mismatch = .TRUE.
      DO i = 1, nrecv
        no_mismatch = no_mismatch &
             .AND. FP_EQ(dst_data(i), DBLE(i+dst_data_offset))
      END DO
      IF (no_mismatch) THEN
        CONTINUE
      ELSE
        CALL test_abort("invalid data", filename, __LINE__)
      END IF
    END DO

    ! cleanup
    CALL xt_exchanger_delete(exchanger)

    DEALLOCATE(dst_data, recv_msgs, send_msgs)
    CALL mpi_comm_free(inter_comm, ierror)
    CALL mpi_comm_free(intra_group_comm, ierror)
  END SUBROUTINE test_intercomm_all2all

END PROGRAM test_exchanger_parallel

!
! Local Variables:
! c-basic-offset: 2
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
