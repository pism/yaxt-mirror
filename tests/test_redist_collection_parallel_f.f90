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
PROGRAM test_redist_collection_parallel
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, xi => xt_int_kind, &
       xt_idxlist, xt_idxlist_delete, xt_stripe, xt_idxvec_new, &
       xt_idxsection_new, xt_idxlist_collection_new, xt_idxstripes_new, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_new, xt_redist_collection_new, &
       xt_redist_copy, xt_redist_delete, xt_redist_s_exchange, &
       xt_idxlist_get_indices, xt_int_mpidt
  ! older PGI compilers do not handle generic interface correctly
#if defined __PGI && (__PGIC__ < 12 || (__PGIC__ ==  12 && __PGIC_MINOR__ <= 10))
  USE xt_redist_base, ONLY: xt_redist_s_exchange
#endif
  USE iso_c_binding, ONLY: c_loc, c_ptr
  IMPLICIT NONE
  INTEGER :: rank, world_size, ierror
  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  CALL mpi_comm_rank(mpi_comm_world, rank, ierror)
  IF (ierror /= MPI_SUCCESS) &
       CALL test_abort('mpi_comm_rank failed', &
       __FILE__, &
       __LINE__)
  CALL mpi_comm_size(mpi_comm_world, world_size, ierror)
  IF (ierror /= MPI_SUCCESS) &
       CALL test_abort('mpi_comm_size failed', &
       __FILE__, &
       __LINE__)

  IF (world_size > 1) THEN
    CALL test_4redist
    CALL test_rr_exchange
  END IF

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)
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

  SUBROUTINE test_4redist
    ! redist test with four different redists
    INTEGER, PARAMETER :: num_tx = 4
    TYPE(xt_idxlist) :: indices_a, indices_b, indices_all
    INTEGER(xt_int_kind), ALLOCATABLE :: index_vector_a(:), &
         index_vector_b(:)
    TYPE(xt_xmap) :: xmaps(num_tx)
    TYPE(xt_redist) :: redists(num_tx), redist, redist_copy
    INTEGER :: i

    ALLOCATE(index_vector_a(2*world_size**2), index_vector_b(2*world_size**2))
    CALL build_idxlists(indices_a, indices_b, indices_all)

    xmaps(1) = xt_xmap_all2all_new(indices_a, indices_b, mpi_comm_world)
    xmaps(2) = xt_xmap_all2all_new(indices_b, indices_a, mpi_comm_world)
    xmaps(3) = xt_xmap_all2all_new(indices_a, indices_all, mpi_comm_world)
    xmaps(4) = xt_xmap_all2all_new(indices_b, indices_all, mpi_comm_world)

    CALL xt_idxlist_get_indices(indices_a, index_vector_a)
    CALL xt_idxlist_get_indices(indices_b, index_vector_b)

    CALL xt_idxlist_delete(indices_a)
    CALL xt_idxlist_delete(indices_b)
    CALL xt_idxlist_delete(indices_all)

    DO i = 1, num_tx
      redists(i) = xt_redist_p2p_new(xmaps(i), xt_int_mpidt)
      CALL xt_xmap_delete(xmaps(i))
    END DO

    redist = xt_redist_collection_new(redists, num_tx, -1, mpi_comm_world)

    ! test communicator of redist
    ! if (!test_communicator(xt_redist_get_MPI_Comm(redist), MPI_COMM_WORLD))
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

    result_sizes(1) = 2*world_size**2
    result_sizes(2) = 2*world_size**2
    result_sizes(3) = 2*world_size**3
    result_sizes(4) = 2*world_size**3

    buf_size = SUM(result_spacing) + SUM(result_sizes)
    ALLOCATE(buf(buf_size))

    ofs = result_spacing(1)
    results_1 => buf(ofs+1:ofs+result_sizes(1))
    ofs = ofs + result_sizes(1) + result_spacing(2)
    results_2 => buf(ofs+1:ofs+result_sizes(2))
    ofs = ofs + result_sizes(2) + result_spacing(3)
    results_3 => buf(ofs+1:ofs+result_sizes(3))
    ofs = ofs + result_sizes(3) + result_spacing(4)
    results_4 => buf(ofs+1:ofs+result_sizes(4))

    CALL do_4redist(redist, index_vector_a, index_vector_b, &
         results_1, results_2, results_3, results_4)

    CALL check_4redist_results(results_1, results_2, results_3, results_4, &
       index_vector_a, index_vector_b)

    ! shift addresses around
    IF (rank == 0) THEN
      ofs = SUM(result_spacing(1:2)) + SUM(result_sizes(1:2))
      results_3 => buf(ofs+1:ofs+result_sizes(3))
    END IF

    CALL do_4redist(redist, index_vector_a, index_vector_b, &
         results_1, results_2, results_3, results_4)

    CALL check_4redist_results(results_1, results_2, results_3, results_4, &
       index_vector_a, index_vector_b)
    ! clean up
    DEALLOCATE(buf)
  END SUBROUTINE exchange_4redist

  SUBROUTINE do_4redist(redist, index_vector_a, index_vector_b, &
       results_1, results_2, results_3, results_4)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(xt_int_kind), INTENT(in), TARGET :: &
         index_vector_a(*), index_vector_b(*)
    INTEGER(xt_int_kind), INTENT(inout), TARGET :: &
         results_1(*), results_2(*), results_3(*), results_4(*)

    TYPE(c_ptr) :: results(4), input(4)
    results(1) = C_LOC(results_1)
    results(2) = C_LOC(results_2)
    results(3) = C_LOC(results_3)
    results(4) = C_LOC(results_4)

    input(1) = C_LOC(index_vector_a)
    input(2) = C_LOC(index_vector_b)
    input(3) = C_LOC(index_vector_a)
    input(4) = C_LOC(index_vector_b)

    CALL xt_redist_s_exchange(redist, 4, input, results)

  END SUBROUTINE do_4redist

  SUBROUTINE check_4redist_results(results_1, results_2, results_3, results_4, &
       index_vector_a, index_vector_b)
    INTEGER(xt_int_kind), INTENT(in) :: index_vector_a(:), index_vector_b(:), &
         results_1(:), results_2(:), results_3(0:), results_4(0:)
    INTEGER(xt_int_kind) :: i, n
    LOGICAL :: p

    IF (ANY(results_1 /= index_vector_b)) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    IF (ANY(results_2 /= index_vector_a)) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    n = INT(SIZE(results_3), xt_int_kind)
    p = .FALSE.
    DO i = 0_xi, n - 1_xi
      p = p .OR. results_3(i) /= i
    END DO
    IF (p) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)

    DO i = 0_xi, n - 1_xi
      p = p .OR. results_4(i) /= i
    END DO
    IF (p) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)
  END SUBROUTINE check_4redist_results


  ! redist test with two redists that do a round robin exchange in
  ! different directions
  SUBROUTINE test_rr_exchange
    TYPE(xt_idxlist) :: src_indices, dst_indices(2)
    INTEGER(xt_int_kind) :: src_indices_(5)
    INTEGER(xt_int_kind) :: i, temp, dst_indices_(5, 2)
    TYPE(xt_xmap) :: xmaps(2)
    TYPE(xt_redist) :: redists(2), redist, redist_copy

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

    xmaps(1) = xt_xmap_all2all_new(src_indices, dst_indices(1), mpi_comm_world)
    xmaps(2) = xt_xmap_all2all_new(src_indices, dst_indices(2), mpi_comm_world)

    CALL xt_idxlist_delete(src_indices)
    CALL xt_idxlist_delete(dst_indices)

    redists(1) = xt_redist_p2p_new(xmaps(1), xt_int_mpidt)
    redists(2) = xt_redist_p2p_new(xmaps(2), xt_int_mpidt)

    CALL xt_xmap_delete(xmaps)

    redist = xt_redist_collection_new(redists, 2, -1, mpi_comm_world)

    ! test communicator of redist
    ! IF (!test_communicator(xt_redist_get_MPI_Comm(redist), MPI_COMM_WORLD))
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
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(xt_int_kind), TARGET, INTENT(in) :: src_indices_(5)
    INTEGER(xt_int_kind), INTENT(in) :: ref_dst_indices_(5, 2)


    INTEGER(xt_int_kind), TARGET :: results(5,2)
    TYPE(c_ptr) :: results_p(2), input(2)

    results = -1

    CALL xt_slice_c_loc(results(:,1), results_p(1))
    CALL xt_slice_c_loc(results(:,2), results_p(2))

    input(1) = C_LOC(src_indices_)
    input(2) = C_LOC(src_indices_)

    CALL xt_redist_s_exchange(redist, input, results_p)

    ! check results
    IF (ANY(results /= ref_dst_indices_)) &
         CALL test_abort("error on xt_redist_s_exchange", &
         __FILE__, &
         __LINE__)
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
