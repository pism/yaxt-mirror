!>
!> @file test_idxvec_f.f90
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
PROGRAM test_idxvec
  USE mpi
  USE yaxt, ONLY: xt_initialize, xt_idxlist, xt_idxvec_new, &
       xt_idxvec_from_stripes_new, &
       xt_int_kind, xt_idxlist_get_pack_size, xt_stripe, &
       xt_idxlist_pack, xt_idxlist_unpack, xt_idxlist_get_intersection, &
       xt_idxlist_get_index_stripes, xt_idxlist_delete, &
       xt_idxlist_get_positions_of_indices, &
       xt_idxlist_get_index_at_position, xt_idxlist_get_indices_at_positions, &
       xt_idxlist_get_position_of_index, xt_idxlist_get_position_of_index_off, &
       xt_bounds, xt_idxlist_get_bounding_box, xt_idxlist_copy
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: check_idxlist, test_err_count, check_stripes, &
       check_offsets
  IMPLICIT NONE
  INTEGER(xt_int_kind), PARAMETER :: index_vector(7) &
       = (/ 1, 2, 3, 4, 5, 6, 7 /)


  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)

  CALL test_idxvec_pack_unpack
  CALL test_copying
  CALL test_repeated_equal_indices
  CALL test_positions
  CALL test_intersection_surjective
  CALL test_intersection_partial
  CALL test_intersection_inverse_partial
  CALL test_intersection_unsort_partial
  CALL test_intersection_unsort_inverse_partial
  CALL test_stripes1
  CALL test_stripes2
  CALL test_stripes3
  CALL test_stripes4
  CALL test_stripes5
  CALL test_stripes6
  CALL test_stripes7
  CALL test_get_indices_at_positions1
  CALL test_get_indices_at_positions2
  CALL test_get_indices_at_positions3
  CALL test_get_positions_of_indices
  CALL test_bounding_box1
  CALL test_bounding_box2
  CALL test_bounding_box3

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", &
       __FILE__, &
       __LINE__)
  CALL finish_mpi
CONTAINS
  SUBROUTINE test_idxvec_pack_unpack
    TYPE(xt_idxlist) :: idxvector, idxvector_copy, intersection
    INTEGER :: buffer_size, position, rank, ierror, request
    CHARACTER(1), ALLOCATABLE :: send_buffer(:), recv_buffer(:)
    TYPE(xt_stripe), PARAMETER :: ref_stripes(1) = (/ xt_stripe(1, 7, 1) /)
    CALL setup_idxvec(idxvector, index_vector)
    buffer_size = xt_idxlist_get_pack_size(idxvector, mpi_comm_world)
    ALLOCATE(send_buffer(buffer_size), recv_buffer(buffer_size))
    position = 0
    CALL xt_idxlist_pack(idxvector, send_buffer, buffer_size, position, &
         mpi_comm_world)
    CALL mpi_comm_rank(mpi_comm_world, rank, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("mpi_comm_rank failed", &
         __FILE__, &
         __LINE__)
    ! send the buffer
    CALL mpi_isend(send_buffer, buffer_size, mpi_packed, rank, 0, &
         mpi_comm_world, request, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("mpi_isend failed", &
         __FILE__, &
         __LINE__)

    CALL mpi_request_free(request, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("mpi_request_free failed", &
         __FILE__, &
         __LINE__)

    ! receive buffer
    CALL mpi_recv(recv_buffer, buffer_size, mpi_packed, rank, 0, &
         mpi_comm_world, mpi_status_ignore, ierror)
    IF (ierror /= mpi_success) &
         CALL test_abort("mpi_recv failed", &
         __FILE__, &
         __LINE__)

    position = 0
    CALL xt_idxlist_unpack(idxvector_copy, recv_buffer, buffer_size, &
         position, mpi_comm_world)

    CALL check_idxlist(idxvector_copy, index_vector)

    intersection = xt_idxlist_get_intersection(idxvector, idxvector_copy)

    CALL check_idxlist(intersection, index_vector)

    CALL compare_stripes(idxvector, ref_stripes)

    CALL xt_idxlist_delete(idxvector)
    CALL xt_idxlist_delete(idxvector_copy)
    CALL xt_idxlist_delete(intersection)

  END SUBROUTINE test_idxvec_pack_unpack

  SUBROUTINE test_copying
    TYPE(xt_idxlist) :: idxvector, idxvector_copy
    CALL setup_idxvec(idxvector, index_vector)
    idxvector_copy = xt_idxlist_copy(idxvector)
    CALL check_idxlist(idxvector_copy, index_vector)
    CALL xt_idxlist_delete(idxvector)
    CALL xt_idxlist_delete(idxvector_copy)
  END SUBROUTINE test_copying

  SUBROUTINE test_repeated_equal_indices
    INTEGER(xt_int_kind), PARAMETER :: index_vector(8) &
         = (/ 1, 2, 3, 7, 5, 6, 7, 7 /)
    TYPE(xt_idxlist) :: idxvector
    CALL setup_idxvec(idxvector, index_vector)
    CALL xt_idxlist_delete(idxvector)
  END SUBROUTINE test_repeated_equal_indices

  SUBROUTINE test_positions
    LOGICAL, PARAMETER :: single_match_only = .TRUE.
    INTEGER(xt_int_kind), PARAMETER :: index_vector(20) &
         = (/ 10, 15, 14, 13, 12, 15, 10, 11, 12, 13, &
         &    23, 18, 19, 20, 21, 31, 26, 27, 28, 29 /), &
         intersection_vector(13) &
         = (/ 12, 12, 13, 13, 14, 15, 15, 20, 21, 23, 28, 29, 31 /)
    INTEGER :: intersection_pos(SIZE(intersection_vector))
    INTEGER, PARAMETER :: ref_intersection_pos(SIZE(intersection_vector)) &
         = (/ 4, 8, 3, 9, 2, 1, 5, 13, 14, 10, 18, 19, 15 /)
    TYPE(xt_idxlist) :: idxvector
    INTEGER :: notfound
    CALL setup_idxvec(idxvector, index_vector)
    notfound = xt_idxlist_get_positions_of_indices(idxvector, &
         intersection_vector, intersection_pos, single_match_only)
    IF (notfound /= 0) &
         CALL test_abort('expected indices not found!', &
         __FILE__, &
         __LINE__)
    CALL check_offsets(intersection_pos, ref_intersection_pos)
    CALL xt_idxlist_delete(idxvector)
  END SUBROUTINE test_positions

  SUBROUTINE test_intersection(index_vector_a, index_vector_b, &
       ref_intersection_indices)
    INTEGER(xt_int_kind), INTENT(in) :: index_vector_a(:), index_vector_b(:), &
         ref_intersection_indices(:)
    ! note: instead of declaring this as the move intuitive idxvector(2),
    ! two distinct variables are used to remain compatible with NAG 5.2
    TYPE(xt_idxlist) :: idxvector1, idxvector2, intersection
    CALL setup_idxvec(idxvector1, index_vector_a)
    CALL setup_idxvec(idxvector2, index_vector_b)
    intersection = xt_idxlist_get_intersection(idxvector1, idxvector2)
    CALL check_idxlist(intersection, ref_intersection_indices)
    CALL xt_idxlist_delete(intersection)
    CALL xt_idxlist_delete(idxvector2)
    CALL xt_idxlist_delete(idxvector1)
  END SUBROUTINE test_intersection

  SUBROUTINE test_intersection_surjective
    INTEGER(xt_int_kind), PARAMETER :: index_vector(3, 2) &
         = RESHAPE((/ 1, 2, 3, 1, 2, 3 /), SHAPE(index_vector)), &
         ref_intersection_indices(3) = (/ 1, 2, 3 /)
    CALL test_intersection(index_vector(:, 1), index_vector(:, 2), &
         ref_intersection_indices)
  END SUBROUTINE test_intersection_surjective

  SUBROUTINE test_intersection_partial
    INTEGER(xt_int_kind), PARAMETER :: index_vector(3, 2) &
         = RESHAPE((/ 1, 2, 3, 2, 3, 4 /), SHAPE(index_vector)), &
         ref_intersection_indices(2) = (/ 2, 3 /)
    CALL test_intersection(index_vector(:, 1), index_vector(:, 2), &
         ref_intersection_indices)
  END SUBROUTINE test_intersection_partial

  SUBROUTINE test_intersection_inverse_partial
    INTEGER(xt_int_kind), PARAMETER :: index_vector(3, 2) &
         = RESHAPE((/ 2, 3, 4, 1, 2, 3 /), SHAPE(index_vector)), &
         ref_intersection_indices(2) = (/ 2, 3 /)
    CALL test_intersection(index_vector(:, 1), index_vector(:, 2), &
         ref_intersection_indices)
  END SUBROUTINE test_intersection_inverse_partial

  SUBROUTINE test_intersection_unsort_partial
    INTEGER(xt_int_kind), PARAMETER :: index_vector(3, 2) &
         = RESHAPE((/ 4, 2, 3, 3, 1, 2 /), SHAPE(index_vector)), &
         ref_intersection_indices(2) = (/ 2, 3 /)
    CALL test_intersection(index_vector(:, 1), index_vector(:, 2), &
         ref_intersection_indices)
  END SUBROUTINE test_intersection_unsort_partial

  SUBROUTINE test_intersection_unsort_inverse_partial
    INTEGER(xt_int_kind), PARAMETER :: index_vector(3, 2) &
         = RESHAPE((/ 3, 1, 2, 4, 2, 3 /), SHAPE(index_vector)), &
         ref_intersection_indices(2) = (/ 2, 3 /)
    CALL test_intersection(index_vector(:, 1), index_vector(:, 2), &
         ref_intersection_indices)
  END SUBROUTINE test_intersection_unsort_inverse_partial

  SUBROUTINE test_idxvec_from_stripes(stripes, ref_indices)
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    INTEGER(xt_int_kind), INTENT(in) :: ref_indices(:)

    TYPE(xt_idxlist) :: idxvec
    idxvec = xt_idxvec_from_stripes_new(stripes)
    CALL check_idxlist(idxvec, ref_indices)
    CALL xt_idxlist_delete(idxvec)
  END SUBROUTINE test_idxvec_from_stripes

  SUBROUTINE test_stripes1
    TYPE(xt_stripe), PARAMETER :: stripes(2) &
         = (/ xt_stripe(5, 5, 1), xt_stripe(4, 5, -1) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(10) &
         = (/ 5, 6, 7, 8, 9, 4, 3, 2, 1, 0 /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes1

  SUBROUTINE test_stripes2
    TYPE(xt_stripe), PARAMETER :: stripes(3) &
         = (/ xt_stripe(0, 5, 1), xt_stripe(2, 5, 1), xt_stripe(4, 5, 1) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 0, 1, 2, 3, 4, 2, 3, 4, 5, 6, 4, 5, 6, 7, 8 /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes2

  SUBROUTINE test_stripes3
    TYPE(xt_stripe), PARAMETER :: stripes(3) &
         = (/ xt_stripe(2, 5, 1), xt_stripe(0, 5, 1), xt_stripe(4, 5, 1) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 4, 5, 6, 7, 8 /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes3

  SUBROUTINE test_stripes4
    TYPE(xt_stripe), PARAMETER :: stripes(3) &
         = (/ xt_stripe(2, 5, 1), xt_stripe(4, 5, -1), xt_stripe(4, 5, 1) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 2, 3, 4, 5, 6, 4, 3, 2, 1, 0, 4, 5, 6, 7, 8 /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes4

  SUBROUTINE test_stripes5
    TYPE(xt_stripe), PARAMETER :: stripes(3) &
         = (/ xt_stripe(0, 5, 3), xt_stripe(1, 5, 3), xt_stripe(14, 5, -3) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 0, 3, 6, 9, 12, 1, 4, 7, 10, 13, 14, 11, 8, 5, 2 /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes5

  SUBROUTINE test_stripes6
    TYPE(xt_stripe), PARAMETER :: stripes(3) &
         = (/ xt_stripe(0, 5, 3), xt_stripe(2, 5, 3), xt_stripe(14, 5, -3) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(15) &
         = (/ 0, 3, 6, 9, 12, 2, 5, 8, 11, 14, 14, 11, 8, 5, 2 /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes6

  SUBROUTINE test_stripes7
    TYPE(xt_stripe), PARAMETER :: stripes(4) &
         = (/ xt_stripe(0, 5, -1), xt_stripe(1, 5, 1), xt_stripe(-5, 5, -1), &
         &    xt_stripe(6, 5, 1) /)
    INTEGER(xt_int_kind), PARAMETER :: ref_indices(20) &
         = (/ 0, -1, -2, -3, -4, 1, 2, 3, 4, 5, -5, -6, -7, -8, -9, &
         &    6, 7, 8, 9, 10 /)
    CALL test_idxvec_from_stripes(stripes, ref_indices)
  END SUBROUTINE test_stripes7

  SUBROUTINE test_get_indices_at_positions(indices, undef_idx, pos)
    INTEGER(xt_int_kind), INTENT(in) :: indices(:), undef_idx
    INTEGER, INTENT(in) :: pos(:)
    INTEGER(xt_int_kind) :: ref_sel_idx(SIZE(pos)), sel_idx(SIZE(pos))
    INTEGER :: i, num_pos, undef_count, ref_undef_count
    TYPE(xt_idxlist) :: idxvec
    CHARACTER(80) :: msg
    num_pos = SIZE(pos)
    idxvec = xt_idxvec_new(indices, INT(SIZE(indices), xt_int_kind))
    ref_undef_count = 0
    DO i = 1, num_pos
      IF (xt_idxlist_get_index_at_position(idxvec, pos(i), ref_sel_idx(i))) &
           THEN
        ref_sel_idx(i) = undef_idx
        ref_undef_count = ref_undef_count + 1
      END IF
    END DO
    undef_count = xt_idxlist_get_indices_at_positions(idxvec, pos, sel_idx, &
         undef_idx)
    IF (undef_count /= ref_undef_count) THEN
      CALL test_abort("test_idxvec_f.f90: (undef_count /= ref_undef_count)", &
           __FILE__, &
           __LINE__)
    END IF
    DO i = 1, num_pos
      IF (sel_idx(i) /= ref_sel_idx(i)) THEN
        WRITE (msg, '(2(a,i0),a)') "test_idxvec_f.f90: (sel_idx(", i, &
             ") /= ref_sel_idx(", i, "))"
        CALL test_abort(msg, &
             __FILE__, &
             __LINE__)
      END IF
    END DO
    CALL xt_idxlist_delete(idxvec)
  END SUBROUTINE test_get_indices_at_positions

  SUBROUTINE test_get_indices_at_positions1
    INTEGER(xt_int_kind), PARAMETER :: indices(16) &
         = (/ 0, 3, 6, 9, 12, 1, 4, 7, 10, 13, 14, 11, 8, 5, 2, 1 /)
    INTEGER(xt_int_kind), PARAMETER :: undef_idx = -HUGE(undef_idx)
    INTEGER, PARAMETER :: pos(13) &
         = (/ 0, 2, 7, 9, 11, 100, 11, 200, 9, 300, 7, 400, 5 /)
    CALL test_get_indices_at_positions(indices, undef_idx, pos)
  END SUBROUTINE test_get_indices_at_positions1

  SUBROUTINE test_get_indices_at_positions2
    INTEGER(xt_int_kind), PARAMETER :: indices(16) &
         = (/ 0, 3, 6, 9, 12, 1, 4, 7, 10, 13, 14, 11, 8, 5, 2, 1 /)
    INTEGER(xt_int_kind), PARAMETER :: undef_idx = -HUGE(undef_idx)
    INTEGER, PARAMETER :: pos(9) &
         = (/ 0, 2, 7, 9, 11, 11, 9, 7, 5 /)
    CALL test_get_indices_at_positions(indices, undef_idx, pos)
  END SUBROUTINE test_get_indices_at_positions2

  SUBROUTINE test_get_indices_at_positions3
    INTEGER(xt_int_kind), PARAMETER :: indices(16) &
         = (/ 0, 3, 6, 9, 12, 1, 4, 7, 10, 13, 14, 11, 8, 5, 2, 1 /)
    INTEGER(xt_int_kind), PARAMETER :: undef_idx = -HUGE(undef_idx)
    INTEGER, PARAMETER :: pos(9) &
         = (/ 100, 102, 107, 109, 1011, 1011, 109, 107, 105 /)
    CALL test_get_indices_at_positions(indices, undef_idx, pos)
  END SUBROUTINE test_get_indices_at_positions3

  SUBROUTINE test_get_positions_of_indices
    INTEGER(xt_int_kind), PARAMETER :: indices(2) = (/ 0, 2 /), &
         selection(3) = (/ 1, 2, 3 /)
    INTEGER :: i, position, positions(SIZE(selection))
    INTEGER, PARAMETER :: ref_positions(3) = (/ -1, 1, -1 /)
    TYPE(xt_idxlist) :: idxvec

    idxvec = xt_idxvec_new(indices, INT(SIZE(indices), xt_int_kind))
    IF (.NOT. xt_idxlist_get_position_of_index(idxvec, 1_xt_int_kind, &
         position)) &
         CALL test_abort('xt_idxlist_get_position_of_index did not return &
         &an error', &
         __FILE__, &
         __LINE__)
    IF (.NOT. xt_idxlist_get_position_of_index_off(idxvec, 1_xt_int_kind, &
         position, 0)) &
         CALL test_abort('xt_idxlist_get_position_of_index_off did not return &
         &an error', &
         __FILE__, &
         __LINE__)
    IF (.NOT. xt_idxlist_get_position_of_index_off(idxvec, 0_xt_int_kind, &
         position, 1)) &
         CALL test_abort('xt_idxlist_get_position_of_index_off did not return &
         &an error', &
         __FILE__, &
         __LINE__)
    IF (xt_idxlist_get_positions_of_indices(idxvec, selection, positions, &
         .FALSE.) /= 2) THEN
      CALL test_abort('xt_idxlist_get_positions_of_indices did not return &
           &correct number of matches', &
           __FILE__, &
           __LINE__)
    END IF
    DO i = 1, SIZE(selection)
      IF (positions(i) /= ref_positions(i)) &
           CALL test_abort('xt_idxlist_get_positions_of_indices returned &
           &incorrect position', &
           __FILE__, &
           __LINE__)
    END DO
    CALL xt_idxlist_delete(idxvec)
  END SUBROUTINE test_get_positions_of_indices

  SUBROUTINE test_bounding_box(indices, global_size, global_start_index, &
       ref_bounds)
    INTEGER(xt_int_kind), INTENT(in) :: indices(:), global_size(:), &
         global_start_index
    TYPE(xt_bounds), INTENT(in) :: ref_bounds(:)

    TYPE(xt_bounds) :: bounds(SIZE(global_size))
    INTEGER :: i, ndim
    TYPE(xt_idxlist) :: idxvec
    CHARACTER(80) :: msg

    ndim = SIZE(global_size)
    IF (SIZE(ref_bounds) /= ndim) &
         CALL test_abort('ERROR: inequal dimensions', &
         __FILE__, &
         __LINE__)

    idxvec = xt_idxvec_new(indices, INT(SIZE(indices), xt_int_kind))

    bounds = xt_idxlist_get_bounding_box(idxvec, global_size, &
         global_start_index)

    DO i = 1, ndim
      IF (bounds(i)%start /= ref_bounds(i)%start) THEN
        WRITE (0, '(2(a,i0))') 'bounds(', i, ')%start=', bounds(i)%start
        WRITE (0, '(2(a,i0))') 'ref_bounds(', i, ')%start=', ref_bounds(i)%start
        WRITE (msg, '(a,i0)') "ERROR: xt_idxlist_get_bounding_box inequal &
             &starts at i=", i
        CALL test_abort(msg, &
             __FILE__, &
             __LINE__)
      END IF
      IF (bounds(i)%size /= ref_bounds(i)%size) THEN
        WRITE (0, '(2(a,i0))') 'bounds(', i, ')%size=', bounds(i)%size
        WRITE (0, '(2(a,i0))') 'ref_bounds(', i, ')%size=', ref_bounds(i)%size
        WRITE (msg, '(a,i0)') "ERROR: xt_idxlist_get_bounding_box inequal &
             &size at i=", i
        CALL test_abort(msg, &
             __FILE__, &
             __LINE__)
      END IF
    END DO
  END SUBROUTINE test_bounding_box

  SUBROUTINE test_bounding_box1
    INTEGER(xt_int_kind), PARAMETER :: indices(2) = (/ 21, 42 /), &
         global_size(3) = 4, global_start_index = 0
    TYPE(xt_bounds), PARAMETER :: ref_bounds(3) = xt_bounds(1, 2)
    CALL test_bounding_box(indices, global_size, global_start_index, ref_bounds)
  END SUBROUTINE test_bounding_box1

  SUBROUTINE test_bounding_box2
    INTEGER(xt_int_kind), PARAMETER :: indices(5) = (/ 45, 35, 32, 48, 33 /), &
         global_size(3) = (/ 5, 4, 3 /), global_start_index = 1
    TYPE(xt_bounds), PARAMETER :: ref_bounds(3) &
         = (/ xt_bounds(2, 2), xt_bounds(2, 2), xt_bounds(1, 2) /)
    CALL test_bounding_box(indices, global_size, global_start_index, ref_bounds)
  END SUBROUTINE test_bounding_box2

  SUBROUTINE test_bounding_box3
    INTEGER(xt_int_kind), PARAMETER :: indices(1) = (/ -1 /), &
         global_size(3) = 4, global_start_index = 0
    TYPE(xt_bounds), PARAMETER :: ref_bounds(3) = xt_bounds(0, 0)
    CALL test_bounding_box(indices(1:0), global_size, global_start_index, &
         ref_bounds)
  END SUBROUTINE test_bounding_box3

  SUBROUTINE setup_idxvec(idxlist, index_vector)
    TYPE(xt_idxlist), INTENT(out) :: idxlist
    INTEGER(xt_int_kind), INTENT(in) :: index_vector(:)

    idxlist = xt_idxvec_new(index_vector, INT(SIZE(index_vector), xt_int_kind))
    CALL check_idxlist(idxlist, index_vector)
  END SUBROUTINE setup_idxvec

  SUBROUTINE compare_stripes(idxlist, ref_stripes)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    TYPE(xt_stripe), INTENT(in) :: ref_stripes(:)

    TYPE(xt_stripe), ALLOCATABLE :: stripes(:)

    CALL xt_idxlist_get_index_stripes(idxlist, stripes)

    CALL check_stripes(stripes, ref_stripes)

  END SUBROUTINE compare_stripes

END PROGRAM test_idxvec
