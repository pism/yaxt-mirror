!>
!! @file test_redist_collection_parallel_f.f90
!! @brief parallel Fortran test of redist_collection class
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
#ifdef VERBOSE
#define PUT_ERR(msg) WRITE (0, '(a)') msg ; err_count__ = err_count__ + 1
#else
#define PUT_ERR(msg) err_count__ = err_count__ + 1
#endif
#define XT_UNUSED(x) IF (SIZE( (/(x)/) ) < 0) CONTINUE
#if ! defined __PGI || ( __PGIC__ >= 11 )
#  define WITH_COLL_OF_COLL_TEST
#endif
#ifdef WITH_COLL_OF_COLL_TEST
MODULE test_coll_of_coll_par
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_size_t, c_ptr, &
       c_funptr, c_funloc, c_char, c_null_char, c_null_ptr, c_f_pointer
  USE mpi
  USE yaxt, ONLY: xt_mpi_fint_kind, xt_config, xt_redist, xt_xmap, &
       xt_idxlist, xt_idxvec_new, xt_idxlist_delete, xt_int_kind, &
       xt_xmap_dist_dir_new, &
       xt_redist_p2p_new, xt_xmap_delete, xt_redist_collection_custom_new, &
       xt_redist_copy, xt_redist_delete, xt_redist_get_MPI_Comm
  USE ftest_common, ONLY: err_count__
  USE test_redist_common, ONLY: &
       sync_mode_test_all, communicators_are_congruent
#include "xt_slice_c_loc.inc"
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: nvalues = 5, nselect = (nvalues+1)/2
  PUBLIC :: test_het_rank_combination

  INTERFACE
    SUBROUTINE check_redist_(redist, sync_mode, num_redists, &
         src, dst_num_elems, dst, dst_prep, dst_prep_info, &
         ref_dst_data, dst_data_dt, ref_dst_data_dt, &
         file, line) BIND(c, name='check_redist_f2c')
      IMPORT :: c_ptr, c_int, c_size_t, c_funptr, &
           xt_mpi_fint_kind, c_char, xt_redist
      TYPE(xt_redist), INTENT(inout) :: redist
      INTEGER(c_int), VALUE :: sync_mode, num_redists
      TYPE(c_ptr), INTENT(in) :: src(num_redists)
      INTEGER(c_size_t), INTENT(in) :: dst_num_elems(num_redists)
      TYPE(c_ptr), INTENT(in) :: dst(num_redists)
      TYPE(c_funptr), INTENT(in) :: dst_prep(num_redists)
      TYPE(c_ptr), VALUE, INTENT(in) :: dst_prep_info
      TYPE(c_ptr), INTENT(in) ::  ref_dst_data(num_redists)
      INTEGER(xt_mpi_fint_kind), INTENT(in) :: dst_data_dt(num_redists), &
           ref_dst_data_dt(num_redists)
      CHARACTER(kind=c_char) :: file(*)
      INTEGER(c_int), VALUE :: line
    END SUBROUTINE check_redist_
  END INTERFACE


  CHARACTER(len=*), PARAMETER :: &
       filename = 'test_redist_collection_parallel_f.f90'
CONTAINS
  SUBROUTINE fill_array_dcmplx(dst, dst_prep_info, dst_num_elems) &
    BIND(c, name='trc_fill_array_real')
    INTEGER(c_size_t), VALUE :: dst_num_elems
    TYPE(c_ptr), VALUE :: dst, dst_prep_info
    DOUBLE COMPLEX, POINTER :: dst_(:)

    INTEGER :: i, n

    XT_UNUSED(dst_prep_info)
    n = INT(dst_num_elems)
    CALL C_F_POINTER(dst, dst_, (/ n /))
    DO i = 1, n
      dst_(i) = (-1.0d0, -1.0d0)
    END DO
  END SUBROUTINE fill_array_dcmplx

  !
  !  test aggregation of redist collections that are themselves
  !  composed of redists which address different ranks.
  !
  SUBROUTINE test_het_rank_combination(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    INTEGER, PARAMETER :: elems_per_rank = 2, num_p2p_redists = 2, &
         num_coll_redists = 1
    INTEGER, PARAMETER :: &
         num_redist_components = num_p2p_redists * num_coll_redists
    TYPE(xt_idxlist) :: src_indices, dst_indices
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redists(num_p2p_redists), redist, redist_coll_of_coll, &
         redist_copy
    DOUBLE COMPLEX, TARGET :: src_data(elems_per_rank), &
         results_(1, num_redist_components), &
         ref_dst_data_(num_redist_components)
    TYPE(c_ptr) :: results(num_redist_components), &
         input(num_redist_components), &
         ref_dst_data(num_redist_components)
    TYPE(c_funptr) :: prep(num_p2p_redists, num_coll_redists)
    INTEGER(xt_int_kind) :: src_indices_(elems_per_rank), dst_indices_(1,2), &
         temp
    INTEGER :: dt(num_redist_components)
    INTEGER(c_size_t), PARAMETER :: &
         num_dst_elems(num_redist_components) = 1_c_size_t
    INTEGER :: comm_rank, comm_size, ierror, i
    CHARACTER(kind=c_char) :: filename_c(LEN(filename)+1)

    DO i = 1, LEN(filename)
      filename_c(i) = filename(i:i)
    END DO
    filename_c(LEN(filename)+1) = c_null_char

    CALL mpi_comm_rank(comm, comm_rank, ierror)
    CALL mpi_comm_size(comm, comm_size, ierror)

    ! redist test with two redists that do a round robin exchange in
    ! different directions


    DO i = 1, elems_per_rank
      src_indices_(i) = INT(comm_rank * elems_per_rank + i - 1, xt_int_kind)
    END DO
    dst_indices_(1,1) = src_indices_(2)
    temp = INT(comm_rank * elems_per_rank - elems_per_rank, xt_int_kind)
    dst_indices_(1,2) &
         = MERGE(INT(temp + comm_size * elems_per_rank, xt_int_kind), &
         &       temp, temp < 0)

    src_indices = xt_idxvec_new(src_indices_)
    DO i = 1, 2
      dst_indices = xt_idxvec_new(dst_indices_(1:1, i:i))
      xmap = xt_xmap_dist_dir_new(src_indices, dst_indices, comm)
      CALL xt_idxlist_delete(dst_indices)
      redists(i) = xt_redist_p2p_new(xmap, mpi_double_complex)
      CALL xt_xmap_delete(xmap)
    END DO
    CALL xt_idxlist_delete(src_indices)
    redist = xt_redist_collection_custom_new(redists, -1, comm, config)
    redist_coll_of_coll &
         = xt_redist_collection_custom_new((/redist/), -1, comm, config)

    ! test communicator of redist
    IF (.NOT. communicators_are_congruent(xt_redist_get_MPI_Comm(redist), &
         &                                comm)) THEN
      PUT_ERR("error in xt_redist_get_MPI_Comm\n")
    END IF

    CALL xt_redist_delete(redists)
    CALL xt_redist_delete(redist)

    DO i = 1, elems_per_rank
      src_data(i) = CMPLX(comm_rank * elems_per_rank, i-1, KIND(1.0d0))
    END DO
    ref_dst_data_(1) = CMPLX(dst_indices_(1,1), 0.0d0, KIND(1.0d0))
    ref_dst_data_(2) = CMPLX(dst_indices_(1,2), 0.0d0, KIND(1.0d0))
    DO i = 1, num_redist_components
      XT_SLICE_C_LOC(results_(:, i), results(i))
      XT_SLICE_C_LOC(src_data, input(i))
      prep(i,1) = C_FUNLOC(fill_array_dcmplx)
      XT_SLICE_C_LOC(ref_dst_data_(i), ref_dst_data(i))
      dt(i) = mpi_double_complex
    END DO

    CALL check_redist_(redist_coll_of_coll, sync_mode_test_all, &
         num_redist_components, input, INT(num_dst_elems, c_size_t), results, &
         prep, c_null_ptr, ref_dst_data, dt, dt, filename_c, __LINE__)

    redist_copy = xt_redist_copy(redist_coll_of_coll)
    CALL xt_redist_delete(redist_coll_of_coll)
    CALL check_redist_(redist_copy, sync_mode_test_all, &
         num_redist_components, input, INT(num_dst_elems, c_size_t), &
         results, prep, c_null_ptr, ref_dst_data, dt, dt, filename_c, __LINE__)
    ! clean up
    CALL xt_redist_delete(redist_copy)
  END SUBROUTINE test_het_rank_combination

END MODULE test_coll_of_coll_par
#endif

PROGRAM test_redist_collection_parallel
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort, icbrt
  USE test_idxlist_utils, ONLY: test_err_count
  USE xt_core, ONLY: i8
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, xi => xt_int_kind, &
       xt_idxlist, xt_idxlist_delete, xt_stripe, xt_idxvec_new, &
       xt_idxsection_new, xt_idxlist_collection_new, xt_idxstripes_new, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_collection_new, &
       xt_redist_copy, xt_redist_delete, xt_redist_s_exchange, &
       xt_idxlist_get_indices, xt_int_mpidt, &
       xt_request, xt_redist_a_exchange, xt_config, xt_config_delete
  ! older PGI compilers do not handle generic interface correctly
#if defined __PGI && (__PGIC__ < 12 || (__PGIC__ ==  12 && __PGIC_MINOR__ <= 10))
  USE xt_redist_base, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
#endif
  USE test_redist_common, ONLY: check_wait_request, redist_exchanger_option
#ifdef WITH_COLL_OF_COLL_TEST
  USE test_coll_of_coll_par, ONLY: test_het_rank_combination
#endif
  USE iso_c_binding, ONLY: c_loc, c_ptr
#include "xt_slice_c_loc.inc"
  IMPLICIT NONE
  INTEGER :: rank, world_size, ierror
  CHARACTER(len=*), PARAMETER :: &
       filename = 'test_redist_collection_parallel_f.f90'
  CHARACTER(len=*), PARAMETER :: err_msg(2) = &
       (/ "error in xt_redist_s_exchange", "error in xt_redist_a_exchange" /)
  TYPE(xt_config) :: config

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  config = redist_exchanger_option()

  CALL mpi_comm_rank(mpi_comm_world, rank, ierror)
  IF (ierror /= MPI_SUCCESS) &
       CALL test_abort('mpi_comm_rank failed', filename, __LINE__)
  CALL mpi_comm_size(mpi_comm_world, world_size, ierror)
  IF (ierror /= MPI_SUCCESS) &
       CALL test_abort('mpi_comm_size failed', filename, __LINE__)

#ifdef WITH_COLL_OF_COLL_TEST
  CALL test_het_rank_combination(mpi_comm_world, config)
#endif
  IF (world_size > 1) THEN
    CALL test_4redist(mpi_comm_world, config)
    CALL test_rr_exchange(mpi_comm_world, config)
  END IF

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)
  CALL xt_config_delete(config)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE build_idxlists(indices_a, indices_b, indices_all)
    TYPE(xt_idxlist), INTENT(out) :: indices_a, indices_b, indices_all

    TYPE(xt_idxlist) :: indices_a_(2)
    INTEGER :: i
    INTEGER(xt_int_kind), PARAMETER :: start = 0
    INTEGER(xt_int_kind) :: global_size(2), local_start(2, 2)
    INTEGER :: local_size(2)

    TYPE(xt_stripe) :: stripe

    global_size(1) = INT(2 * world_size, xi)
    global_size(2) = INT(world_size**2, xi)
    local_size = world_size
    local_start = RESHAPE((/ 0_xi, INT(rank*world_size, xi), &
         INT(world_size, xi), &
         INT((world_size-(rank+1))*world_size, xi) /), (/ 2, 2 /))

    DO i = 1, 2
      indices_a_(i) = xt_idxsection_new(start, global_size, local_size, &
           local_start(:, i))
    END DO
    indices_a = xt_idxlist_collection_new(indices_a_)

    CALL xt_idxlist_delete(indices_a_(1))
    CALL xt_idxlist_delete(indices_a_(2))

    stripe = xt_stripe(INT(rank * 2 * world_size**2, xi), 1_xi, 2*world_size**2)
    indices_b = xt_idxstripes_new(stripe)

    stripe = xt_stripe(0_xi, 1_xi, 2*world_size**3)
    indices_all = xt_idxstripes_new(stripe)
  END SUBROUTINE build_idxlists

  SUBROUTINE test_4redist(comm, config)
    ! redist test with four different redists
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    INTEGER, PARAMETER :: num_tx = 4
    TYPE(xt_idxlist) :: indices_a, indices_b, indices_all
    INTEGER(xt_int_kind), ALLOCATABLE :: index_vector_a(:), &
         index_vector_b(:)
    TYPE(xt_xmap) :: xmaps(num_tx)
    TYPE(xt_redist) :: redists(num_tx), redist, redist_copy
    INTEGER :: i, vec_size

    IF (world_size &
         > icbrt((HUGE(1_xi)-MOD(HUGE(1_xi),2_xi))/2_xi)) &
         CALL test_abort('communicator too large for test', filename, __LINE__)

    vec_size = 2*world_size**2
    ALLOCATE(index_vector_a(vec_size), index_vector_b(vec_size))
    CALL build_idxlists(indices_a, indices_b, indices_all)

    xmaps(1) = xt_xmap_all2all_new(indices_a, indices_b, comm)
    xmaps(2) = xt_xmap_all2all_new(indices_b, indices_a, comm)
    xmaps(3) = xt_xmap_all2all_new(indices_a, indices_all, comm)
    xmaps(4) = xt_xmap_all2all_new(indices_b, indices_all, comm)

    CALL xt_idxlist_get_indices(indices_a, index_vector_a)
    CALL xt_idxlist_get_indices(indices_b, index_vector_b)

    CALL xt_idxlist_delete(indices_a)
    CALL xt_idxlist_delete(indices_b)
    CALL xt_idxlist_delete(indices_all)

    DO i = 1, num_tx
      redists(i) = xt_redist_p2p_new(xmaps(i), xt_int_mpidt)
      CALL xt_xmap_delete(xmaps(i))
    END DO

    redist = xt_redist_collection_new(redists, num_tx, -1, comm, config)

    ! test communicator of redist
    ! if (!test_communicator(xt_redist_get_MPI_Comm(redist), COMM))
    !   PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    CALL xt_redist_delete(redists)

    CALL exchange_4redist(redist, index_vector_a, index_vector_b)
    redist_copy = xt_redist_copy(redist)
    CALL xt_redist_delete(redist)
    CALL exchange_4redist(redist_copy, index_vector_a, index_vector_b)

    ! clean up
    CALL xt_redist_delete(redist_copy)
  END SUBROUTINE test_4redist

  SUBROUTINE exchange_4redist(redist, index_vector_a, index_vector_b)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(xt_int_kind), INTENT(in) :: index_vector_a(2*world_size**2), &
         index_vector_b(2*world_size**2)
    INTEGER(xt_int_kind), TARGET, ALLOCATABLE :: buf(:)
    INTEGER(xt_int_kind), POINTER :: results_1(:), &
         results_2(:), results_3(:), results_4(:)
    INTEGER :: result_sizes(4), buf_size, ofs
    INTEGER, PARAMETER :: result_spacing(4) = (/ 2, 14, 5, 8 /)
    INTEGER :: iexch

    result_sizes(1) = 2*world_size**2
    result_sizes(2) = 2*world_size**2
    result_sizes(3) = 2*world_size**3
    result_sizes(4) = 2*world_size**3

    buf_size = SUM(result_spacing) + SUM(result_sizes)
    ALLOCATE(buf(buf_size))
    DO iexch = 1, 2
      buf(:) = -1_xt_int_kind
      ofs = result_spacing(1)
      results_1 => buf(ofs+1:ofs+result_sizes(1))
      ofs = ofs + result_sizes(1) + result_spacing(2)
      results_2 => buf(ofs+1:ofs+result_sizes(2))
      ofs = ofs + result_sizes(2) + result_spacing(3)
      results_3 => buf(ofs+1:ofs+result_sizes(3))
      ofs = ofs + result_sizes(3) + result_spacing(4)
      results_4 => buf(ofs+1:ofs+result_sizes(4))

      CALL do_4redist(redist, index_vector_a, index_vector_b, &
           results_1, results_2, results_3, results_4, iexch)

      CALL check_4redist_results(results_1, results_2, results_3, results_4, &
           index_vector_a, index_vector_b, iexch)
      buf(:) = -1_xt_int_kind
      ! shift addresses around
      IF (rank == 0) THEN
        ofs = SUM(result_spacing(1:2)) + SUM(result_sizes(1:2))
        results_3 => buf(ofs+1:ofs+result_sizes(3))
      END IF

      CALL do_4redist(redist, index_vector_a, index_vector_b, &
           results_1, results_2, results_3, results_4, iexch)

      CALL check_4redist_results(results_1, results_2, results_3, results_4, &
           index_vector_a, index_vector_b, iexch)
    ENDDO
    ! clean up
    DEALLOCATE(buf)
  END SUBROUTINE exchange_4redist

  SUBROUTINE do_4redist(redist, index_vector_a, index_vector_b, &
       results_1, results_2, results_3, results_4, iexch)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(xt_int_kind), INTENT(in), TARGET :: &
         index_vector_a(*), index_vector_b(*)
    INTEGER(xt_int_kind), INTENT(inout), TARGET :: &
         results_1(*), results_2(*), results_3(*), results_4(*)
    INTEGER, INTENT(in) :: iexch

    TYPE(c_ptr) :: results(4), input(4)
    TYPE(xt_request) :: request

    results(1) = C_LOC(results_1)
    results(2) = C_LOC(results_2)
    results(3) = C_LOC(results_3)
    results(4) = C_LOC(results_4)

    input(1) = C_LOC(index_vector_a)
    input(2) = C_LOC(index_vector_b)
    input(3) = C_LOC(index_vector_a)
    input(4) = C_LOC(index_vector_b)
    IF (iexch == 1) THEN
      CALL xt_redist_s_exchange(redist, 4_i8, input, results)
    ELSE
      CALL xt_redist_a_exchange(redist, 4_i8, input, results, request)
      CALL check_wait_request(request, filename, __LINE__)
    ENDIF
  END SUBROUTINE do_4redist

  SUBROUTINE check_4redist_results(results_1, results_2, results_3, results_4, &
       index_vector_a, index_vector_b, iexch)
    INTEGER(xt_int_kind), INTENT(in) :: index_vector_a(:), index_vector_b(:), &
         results_1(:), results_2(:), results_3(0:), results_4(0:)
    INTEGER, INTENT(in) :: iexch
    INTEGER(xt_int_kind) :: i, n
    LOGICAL :: p_3, p_4

    IF (ANY(results_1 /= index_vector_b)) &
         CALL test_abort(err_msg(iexch), filename, __LINE__)

    IF (ANY(results_2 /= index_vector_a)) &
         CALL test_abort(err_msg(iexch), filename, __LINE__)

    n = INT(SIZE(results_3), xt_int_kind)
    p_3 = .FALSE.
    p_4 = .FALSE.
    DO i = 0_xi, n - 1_xi
      p_3 = p_3 .OR. results_3(i) /= i
      p_4 = p_4 .OR. results_4(i) /= i
    END DO
    IF (p_3 .OR. p_4) CALL test_abort(err_msg(iexch), filename, __LINE__)
  END SUBROUTINE check_4redist_results


  ! redist test with two redists that do a round robin exchange in
  ! different directions
  SUBROUTINE test_rr_exchange(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config

    TYPE(xt_idxlist) :: src_indices, dst_indices(2)
    INTEGER(xt_int_kind) :: src_indices_(5)
    INTEGER(xt_int_kind) :: i, temp, dst_indices_(5, 2)
    TYPE(xt_xmap) :: xmaps(2)
    TYPE(xt_redist) :: redists(2), redist, redist_copy

    IF (world_size > (HUGE(1_xi)-MOD(HUGE(1_xi),5_xi))/5_xi) &
      CALL test_abort('communicator too large for test', filename, __LINE__)

    DO i = 1_xi, 5_xi
      src_indices_(i) = INT(rank, xi) * 5_xi + (i - 1_xi)
      dst_indices_(i, 1) = MOD(src_indices_(i) + 1_xi, &
           &                   INT(world_size, xi) * 5_xi)
      temp = src_indices_(i) - 1_xi
      dst_indices_(i, 2) = MERGE(INT(world_size, xi) * 5_xi - 1_xi, &
           &                     temp, temp < 0_xi)
    END DO

    src_indices = xt_idxvec_new(src_indices_, 5)
    dst_indices(1) = xt_idxvec_new(dst_indices_(:, 1))
    dst_indices(2) = xt_idxvec_new(dst_indices_(:, 2))

    xmaps(1) = xt_xmap_all2all_new(src_indices, dst_indices(1), comm)
    xmaps(2) = xt_xmap_all2all_new(src_indices, dst_indices(2), comm)

    CALL xt_idxlist_delete(src_indices)
    CALL xt_idxlist_delete(dst_indices)

    redists(1) = xt_redist_p2p_new(xmaps(1), xt_int_mpidt)
    redists(2) = xt_redist_p2p_new(xmaps(2), xt_int_mpidt)

    CALL xt_xmap_delete(xmaps)

    redist = xt_redist_collection_new(redists, 2, -1, comm, config)

    ! test communicator of redist
    ! IF (!test_communicator(xt_redist_get_MPI_Comm(redist), comm))
    !     PUT_ERR("error in xt_redist_get_MPI_Comm\n");

    CALL xt_redist_delete(redists)

    CALL rr_exchange(redist, src_indices_, dst_indices_)
    redist_copy = xt_redist_copy(redist)
    CALL xt_redist_delete(redist)
    CALL rr_exchange(redist_copy, src_indices_, dst_indices_)

    ! clean up
    CALL xt_redist_delete(redist_copy)
  END SUBROUTINE test_rr_exchange

  SUBROUTINE rr_exchange(redist, src_indices_, ref_dst_indices_)
#if defined __GNUC__ && __GNUC__ >= 5 && ( __GNUC__ <= 7 \
    || __GNUC__ == 8 && __GNUC_MINOR__ < 4 )
    ! gcc versions 5.x to 8.x have a bug that lets them evaluate the
    ! ANY test too early if results never gets passed to some external
    ! routine directly, 9.x is not only fixed again, but requires to
    ! have the explicit escaping of a pointer to results via C_LOC
    USE yaxt, ONLY: xt_slice_c_loc
#undef XT_SLICE_C_LOC
#define XT_SLICE_C_LOC(slice, cptr) CALL xt_slice_c_loc(slice, cptr)
#endif
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, PARAMETER :: nredist = 2
    INTEGER(xt_int_kind), TARGET, INTENT(in) :: src_indices_(5)
    INTEGER(xt_int_kind), INTENT(in) :: ref_dst_indices_(5, nredist)

    INTEGER(xt_int_kind), TARGET :: results(5,nredist)
    TYPE(c_ptr) :: results_p(nredist), input(nredist)
    INTEGER :: iexch, i
    TYPE(xt_request) :: request

    DO i = 1, nredist
      XT_SLICE_C_LOC(results(:,i), results_p(i))
      input(i) = C_LOC(src_indices_)
    END DO

    DO iexch = 1, 2
      results = -1

      IF (iexch == 1) THEN
        CALL xt_redist_s_exchange(redist, input, results_p)
      ELSE
        CALL xt_redist_a_exchange(redist, input, results_p, request)
        CALL check_wait_request(request, filename, __LINE__)
      ENDIF

      ! check results
      IF (ANY(results /= ref_dst_indices_)) &
           CALL test_abort(err_msg(iexch), filename, __LINE__)
    ENDDO
  END SUBROUTINE rr_exchange

END PROGRAM test_redist_collection_parallel
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
