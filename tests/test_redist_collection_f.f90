!>
!! @file test_redist_collection_f.f90
!! @brief Fortran test of redist_collection class
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
MODULE test_coll_of_coll
  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_size_t, c_ptr, &
       c_funptr, c_funloc, c_char, c_null_char, c_null_ptr, c_f_pointer
  USE mpi
  USE yaxt, ONLY: xt_mpi_fint_kind, xt_config, xt_redist, xt_xmap, &
       xt_redist_p2p_new, xt_xmap_delete, xt_redist_collection_new, &
       xt_redist_copy, xt_redist_delete, xt_redist_get_MPI_Comm
  USE ftest_common, ONLY: err_count__
  USE test_redist_common, ONLY: build_odd_selection_xmap, &
       sync_mode_test_all, communicators_are_congruent
#include "xt_slice_c_loc.inc"
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: nvalues = 5, nselect = (nvalues+1)/2
  PUBLIC :: test_collection_of_collection
  TYPE :: src_databuf
    REAL :: data_real(nvalues)
    INTEGER(xt_mpi_fint_kind) :: data_int(nvalues)
  END TYPE src_databuf
  TYPE :: src_databuf_displaced
    REAL :: data_real(nvalues+1)
    INTEGER(xt_mpi_fint_kind) :: data_int(nvalues)
  END TYPE src_databuf_displaced
  TYPE :: dst_databuf
    REAL :: data_real(nselect)
    INTEGER(xt_mpi_fint_kind) :: data_int(nselect)
  END TYPE dst_databuf
  INTEGER, PARAMETER :: num_redist_coll = 2, num_redist_p2p = 2

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


  CHARACTER(len=*), PARAMETER :: filename = 'test_redist_collection_f.f90'
CONTAINS
  SUBROUTINE fill_array_real(dst, dst_prep_info, dst_num_elems) &
    BIND(c, name='trc_fill_array_real')
    INTEGER(c_size_t), VALUE :: dst_num_elems
    TYPE(c_ptr), VALUE :: dst, dst_prep_info
    REAL, POINTER :: dst_(:)

    INTEGER :: i, n

    XT_UNUSED(dst_prep_info)
    n = INT(dst_num_elems)
    CALL c_f_POINTER(dst, dst_, (/ n /))
    DO i = 1, n
      dst_(i) = -1.0
    END DO
  END SUBROUTINE fill_array_real

  SUBROUTINE fill_array_int(dst, dst_prep_info, dst_num_elems) &
       BIND(c, name='trc_fill_array_int')
    INTEGER(c_size_t), VALUE :: dst_num_elems
    TYPE(c_ptr), VALUE :: dst, dst_prep_info
    INTEGER(xt_mpi_fint_kind), POINTER :: dst_(:)

    INTEGER :: i, n

    XT_UNUSED(dst_prep_info)
    n = INT(dst_num_elems)
    CALL C_F_POINTER(dst, dst_, (/ n /))
    DO i = 1, n
      dst_(i) = -1
    END DO
  END SUBROUTINE fill_array_int

  SUBROUTINE check_redist_coll2(redist_coll2, src_data, dst_data, &
       ref_dst_data)
    TYPE(xt_redist), INTENT(inout) :: redist_coll2
    TYPE(src_databuf), TARGET, INTENT(in) :: src_data
    TYPE(dst_databuf), TARGET, INTENT(out) :: dst_data(num_redist_coll)
    TYPE(dst_databuf), TARGET, INTENT(in) :: ref_dst_data
    TYPE(c_ptr) :: src_data_p(num_redist_p2p, num_redist_coll), &
         dst_data_p(num_redist_p2p, num_redist_coll), &
         ref_dst_data_p(num_redist_p2p, num_redist_coll)
    INTEGER(c_size_t) :: dst_num_elem(num_redist_p2p, num_redist_coll)
    TYPE(c_funptr) :: prep(num_redist_p2p, num_redist_coll)
    INTEGER(xt_mpi_fint_kind) :: dt(num_redist_p2p, num_redist_coll)
    CHARACTER(kind=c_char) :: filename_c(LEN(filename)+1)
    TYPE(src_databuf_displaced), TARGET :: src_data_disp
    type(xt_redist) :: redist_coll2_copy
    INTEGER :: i

    DO i = 1, num_redist_coll
      XT_SLICE_C_LOC(src_data%data_real, src_data_p(1, i))
      XT_SLICE_C_LOC(src_data%data_int, src_data_p(2, i))
      XT_SLICE_C_LOC(dst_data(i)%data_real, dst_data_p(1, i))
      XT_SLICE_C_LOC(dst_data(i)%data_int, dst_data_p(2, i))
      XT_SLICE_C_LOC(ref_dst_data%data_real, ref_dst_data_p(1, i))
      XT_SLICE_C_LOC(ref_dst_data%data_int, ref_dst_data_p(2, i))
      dt(1, i) = mpi_real
      dt(2, i) = mpi_integer
      prep(1, i) = C_FUNLOC(fill_array_real)
      prep(2, i) = C_FUNLOC(fill_array_int)
    END DO
    dst_num_elem = INT(nselect, c_size_t)
    DO i = 1, LEN(filename)
      filename_c(i) = filename(i:i)
    END DO
    filename_c(LEN(filename)+1) = c_null_char
    CALL check_redist_(redist_coll2, sync_mode_test_all, &
         num_redist_p2p*num_redist_coll, &
         src_data_p, dst_num_elem, dst_data_p, prep, c_null_ptr, &
         ref_dst_data_p, dt, dt, filename_c, __LINE__)
    DO i = 1, SIZE(src_data%data_real)
      src_data_disp%data_real(i) = src_data%data_real(i)
      src_data_disp%data_int(i) = src_data%data_int(i)
    END DO
    src_data_disp%data_real(SIZE(src_data_disp%data_real)) = -2.0
    ! check the same with the same data at different relative offset
    DO i = 1, num_redist_coll
      XT_SLICE_C_LOC(src_data_disp%data_real, src_data_p(1, i))
      XT_SLICE_C_LOC(src_data_disp%data_int, src_data_p(2, i))
    END DO
    CALL check_redist_(redist_coll2, sync_mode_test_all, &
         num_redist_p2p*num_redist_coll, &
         src_data_p, dst_num_elem, dst_data_p, prep, c_null_ptr, &
         ref_dst_data_p, dt, dt, filename_c, __LINE__)
    redist_coll2_copy = xt_redist_copy(redist_coll2)
    CALL check_redist_(redist_coll2, sync_mode_test_all, &
         num_redist_p2p*num_redist_coll, &
         src_data_p, dst_num_elem, dst_data_p, prep, c_null_ptr, &
         ref_dst_data_p, dt, dt, filename_c, __LINE__)
    CALL xt_redist_delete(redist_coll2_copy)
  END SUBROUTINE check_redist_coll2

  SUBROUTINE test_collection_of_collection(comm, config)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    ! general test WITH one redist
    ! set up data
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redists(num_redist_p2p), redist_coll(num_redist_coll), &
         redist_coll_coll
    TYPE(src_databuf), PARAMETER :: src_data &
         = src_databuf( (/ 1.0,2.0,3.0,4.0,5.0 /), (/ 1,2,3,4,5 /) )
    TYPE(dst_databuf) :: dst_data(num_redist_coll)
    TYPE(dst_databuf), PARAMETER :: ref_dst_data &
         = dst_databuf( (/ 1.0,3.0,5.0 /), (/ 1,3,5 /) )
    INTEGER :: i

    xmap = build_odd_selection_xmap(nvalues, comm)

    redists(1) = xt_redist_p2p_new(xmap, MPI_REAL)
    redists(2) = xt_redist_p2p_new(xmap, MPI_INTEGER)

    CALL xt_xmap_delete(xmap)

    ! generate redist_collection
    redist_coll(1) &
         = xt_redist_collection_new(redists, -1, comm)
    DO i = 2, num_redist_coll
      redist_coll(i) = redist_coll(1)
    END DO
    CALL xt_redist_delete(redists)

    redist_coll_coll &
         = xt_redist_collection_new(redist_coll, 2, comm, config)

    ! test communicator of redist
    IF (.NOT. communicators_are_congruent(&
         xt_redist_get_MPI_Comm(redist_coll_coll), comm)) THEN
      PUT_ERR("error in xt_redist_get_MPI_Comm")
    END IF

    CALL xt_redist_delete(redist_coll(1));

    ! test exchange
    CALL check_redist_coll2(redist_coll_coll, src_data, dst_data, &
         ref_dst_data)

    ! clean up
    CALL xt_redist_delete(redist_coll_coll)

  END SUBROUTINE test_collection_of_collection
END MODULE test_coll_of_coll
#endif
PROGRAM test_redist_collection
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort, cmp_arrays
  USE test_idxlist_utils, ONLY: test_err_count
  USE xt_core, ONLY: i2
  USE yaxt, ONLY: xt_initialize, xt_finalize, &
       xt_xmap, xt_xmap_delete, xt_xmap_all2all_new, &
       xt_redist, xt_redist_p2p_new, xt_redist_collection_new, &
       xt_redist_delete, xt_redist_copy, &
       xt_redist_s_exchange, xt_redist_a_exchange, &
       xt_idxempty_new, xt_idxlist_delete, &
       xt_idxlist, xt_request, xt_config, xt_config_delete
#if !defined HAVE_FC_LOGICAL_INTEROP || !defined(__GNUC__) || __GNUC__ > 4 \
  || (__GNUC__ == 4 && __GNUC_MINOR__ > 8)
#else
  USE yaxt, ONLY: xt_slice_c_loc
#endif
  ! older PGI compilers do not handle generic interface correctly
#if defined __PGI && (__PGIC__ < 12 || (__PGIC__ ==  12 && __PGIC_MINOR__ <= 10))
  USE xt_redist_base, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
#endif
  USE test_redist_common, ONLY: build_odd_selection_xmap, check_redist, &
       check_wait_request, redist_exchanger_option
  USE iso_c_binding, ONLY: c_loc, c_ptr
  USE redist_collection_displace, ONLY: test_displacement_variations
#ifdef WITH_COLL_OF_COLL_TEST
  USE test_coll_of_coll, ONLY: test_collection_of_collection
#endif
#include "xt_slice_c_loc.inc"
  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER :: filename = 'test_redist_collection_f.f90'
  CHARACTER(len=*), PARAMETER :: err_msg(2) = &
       (/ "error on xt_redist_s_exchange", "error on xt_redist_a_exchange" /)
  TYPE(xt_config) :: config

  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  config = redist_exchanger_option()

  CALL simple_test(mpi_comm_world, config)
  CALL simple_test2(mpi_comm_world, config)
  CALL test_empty_redist(mpi_comm_world, config)
  CALL test_repeated_redist(mpi_comm_world, config, -1)
  CALL test_repeated_redist(mpi_comm_world, config, 0)
  CALL test_displacement_variations(mpi_comm_world, config, 3)
  CALL test_displacement_variations(mpi_comm_world, config, 40)
#ifdef WITH_COLL_OF_COLL_TEST
  CALL test_collection_of_collection(mpi_comm_world, config)
#endif
  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)
  CALL xt_config_delete(config)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE simple_test(comm, config)
    ! general test with one redist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    ! set up data
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_coll, redist_copy
    INTEGER, PARAMETER :: src_slice_len = 5, dst_slice_len = 3
    DOUBLE PRECISION, PARAMETER :: &
         ref_dst_data(dst_slice_len) = (/ 1.0d0, 3.0d0, 5.0d0 /), &
         src_data(src_slice_len) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /)
    DOUBLE PRECISION :: dst_data(dst_slice_len)


    xmap = build_odd_selection_xmap(src_slice_len, comm)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)
    CALL xt_xmap_delete(xmap)
    redist_copy = xt_redist_copy(redist)
    CALL xt_redist_delete(redist)
    redist = redist_copy

    ! generate redist_collection
    redist_coll = xt_redist_collection_new((/ redist /), 1, -1, comm, config)

    CALL xt_redist_delete(redist)

    ! test exchange
    CALL check_redist(redist_coll, src_data, dst_data, ref_dst_data)

    ! clean up
    CALL xt_redist_delete(redist_coll)
  END SUBROUTINE simple_test

  SUBROUTINE simple_test2(comm, config)
    ! general test with one redist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    ! set up data
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist_coll, redist_copy, &
         redist_components(2)
    INTEGER, PARAMETER :: src_slice_len = 5, dst_slice_len = 3
    TYPE src_data_collection
      DOUBLE PRECISION :: dp(src_slice_len)
      LOGICAL :: l(src_slice_len)
    END TYPE src_data_collection
    TYPE dst_data_collection
      DOUBLE PRECISION :: dp(dst_slice_len)
      LOGICAL :: l(dst_slice_len)
    END TYPE dst_data_collection
    TYPE(src_data_collection), SAVE, TARGET :: src_data = src_data_collection(&
         (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /), &
         (/ .TRUE., .FALSE., .TRUE., .FALSE., .TRUE. /))
    TYPE(dst_data_collection), PARAMETER :: &
         ref_dst_data = dst_data_collection((/ 1.0d0, 3.0d0, 5.0d0 /), &
         (/ .TRUE., .TRUE., .TRUE. /))
    TYPE(dst_data_collection), TARGET :: dst_data
    TYPE(c_ptr) :: src_data_p(2), dst_data_p(2)

    xmap = build_odd_selection_xmap(src_slice_len, comm)

    redist_components(1) = xt_redist_p2p_new(xmap, mpi_double_precision)
    redist_components(2) = xt_redist_p2p_new(xmap, mpi_logical)
    CALL xt_xmap_delete(xmap)

    ! generate redist_collection
    redist_coll = xt_redist_collection_new(redist_components, comm, config)
    CALL xt_redist_delete(redist_components)
    redist_copy = xt_redist_copy(redist_coll)
    CALL xt_redist_delete(redist_coll)
    redist_coll = redist_copy

    ! test exchange
    ! GNU Fortran versions up to 4.8 cannot call c_loc for type components,
    ! instant ICE, and
#if !defined(__GNUC__) || __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 8)
#  define COMP_C_LOC(v, p) p = C_LOC(v)
#else
#  define COMP_C_LOC(v, p) CALL xt_slice_c_loc(v, p)
#endif
    ! some compilers won't create c_ptr's to LOGICALs and gfortran
    ! version 8 and below  has a bug that prevents it from noticing the escaped
    ! pointer (and computes ANY(dst_data%l .NEQV. ref_dst_data%l early
    ! for that reason).
#if defined HAVE_FC_LOGICAL_INTEROP && ( !defined __GNUC__ || __GNUC__ > 8 )
#  define L_COMP_C_LOC(v, p) p = C_LOC(v)
#else
#  define L_COMP_C_LOC(v, p) CALL xt_slice_c_loc(v, p)
#endif

    COMP_C_LOC(src_data%dp, src_data_p(1))
    L_COMP_C_LOC(src_data%l, src_data_p(2))
    dst_data%dp = -1.0d0
    dst_data%l = .FALSE.
    COMP_C_LOC(dst_data%dp, dst_data_p(1))
    L_COMP_C_LOC(dst_data%l, dst_data_p(2))
    CALL xt_redist_s_exchange(redist_coll, src_data_p, dst_data_p)
    IF (ANY(dst_data%l .NEQV. ref_dst_data%l)) &
         CALL test_abort("error in xt_redist_s_exchange", filename, __LINE__)
    IF (cmp_arrays(dst_data%dp, ref_dst_data%dp)) &
         CALL test_abort("error in xt_redist_s_exchange", filename, __LINE__)

    ! clean up
    CALL xt_redist_delete(redist_coll)
  END SUBROUTINE simple_test2

  SUBROUTINE test_empty_redist(comm, config)
    ! general test with empty redist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    ! set up data
    TYPE(xt_idxlist) :: src_idxlist, dst_idxlist
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redist, redist_coll, redist_copy


    src_idxlist = xt_idxempty_new()
    dst_idxlist = xt_idxempty_new()

    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, comm)

    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

    redist = xt_redist_p2p_new(xmap, mpi_double_precision)
    CALL xt_xmap_delete(xmap)
    redist_copy = xt_redist_copy(redist)
    CALL xt_redist_delete(redist)
    redist = redist_copy

    ! generate redist_collection
    redist_coll = xt_redist_collection_new((/ redist /), 1, -1, comm, config)

    CALL xt_redist_delete(redist)

    ! clean up
    CALL xt_redist_delete(redist_coll)
  END SUBROUTINE test_empty_redist

  SUBROUTINE test_repeated_redist_ds(redist_coll, src_data, permutation)
    TYPE(xt_redist), INTENT(in) :: redist_coll
    DOUBLE PRECISION, INTENT(in), TARGET :: src_data(5, 3)
    INTEGER, INTENT(in) :: permutation(3)
    INTEGER :: i, j
    DOUBLE PRECISION, PARAMETER :: ref_dst_data(3, 3) &
         = RESHAPE((/ ((DBLE(i + j), i = 1,5,2), j = 0,10,5) /), (/ 3, 3 /))
    DOUBLE PRECISION, TARGET :: dst_data(3, 3)
    TYPE(c_ptr) :: src_data_p(3), dst_data_p(3)

    INTEGER :: iexch
    TYPE(xt_request) :: request

    DO i = 1, 3
      XT_SLICE_C_LOC(src_data(:, permutation(i)), src_data_p(i))
      XT_SLICE_C_LOC(dst_data(:, permutation(i)), dst_data_p(i))
    END DO
    DO iexch = 1, 2
      dst_data = -1.0d0
      IF (iexch == 1) THEN
        CALL xt_redist_s_exchange(redist_coll, 3_i2, src_data_p, dst_data_p)
      ELSE
        CALL xt_redist_a_exchange(redist_coll, 3_i2, src_data_p, dst_data_p, &
             request)
        CALL check_wait_request(request, filename, __LINE__)
      ENDIF
      IF (cmp_arrays(ref_dst_data, dst_data)) &
           CALL test_abort(err_msg(iexch), filename, __LINE__)
    ENDDO
  END SUBROUTINE test_repeated_redist_ds

  SUBROUTINE test_repeated_redist(comm, config, cache_size)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    INTEGER, INTENT(in) :: cache_size
    ! test with one redist used three times (with two different input data
    ! displacements -> test of cache) (with default cache size)
    ! set up data
    INTEGER, PARAMETER :: num_slice = 3
    INTEGER, PARAMETER :: src_slice_len = 5
    INTEGER :: i
    DOUBLE PRECISION, PARAMETER :: src_data(src_slice_len, num_slice) &
         = RESHAPE((/ (DBLE(i), i = 1, 15)/), &
         &         (/ src_slice_len, num_slice /))
    TYPE(xt_xmap) :: xmap
    TYPE(xt_redist) :: redists(num_slice), redist_coll, redist_coll_copy
    INTEGER, PARAMETER :: permutation(3, 3) &
         = RESHAPE((/ 1, 2, 3, 2, 1, 3, 1, 2, 3 /), (/ 3, 3 /))

    xmap = build_odd_selection_xmap(src_slice_len, comm)

    redists = xt_redist_p2p_new(xmap, mpi_double_precision)

    CALL xt_xmap_delete(xmap)

    ! generate redist_collection

    redist_coll = xt_redist_collection_new(redists, 3, cache_size, &
         comm, config)

    CALL xt_redist_delete(redists(1))

    DO i = 1, 3
      ! test exchange, first with simple sequence, then permuted
      ! offsets and then original offsets again
      CALL test_repeated_redist_ds(redist_coll, src_data, permutation(:, i))
    END DO

    ! and the copy
    redist_coll_copy = xt_redist_copy(redist_coll)
    CALL xt_redist_delete(redist_coll)
    DO i = 1, 3
      ! test exchange, first with simple sequence, then permuted
      ! offsets and then original offsets again
      CALL test_repeated_redist_ds(redist_coll_copy, src_data, permutation(:, i))
    END DO

    ! clean up
    CALL xt_redist_delete(redist_coll_copy)
  END SUBROUTINE test_repeated_redist

END PROGRAM test_redist_collection
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
