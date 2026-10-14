!>
!! @file test_xmap_bucket_tiler_f.f90
!! @brief Fortran test of distributed directory
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
#define XT_UNUSED(x) IF (SIZE( (/(x)/) ) < 0) CONTINUE
MODULE test_xmap_bucket_tiler
  USE yaxt, ONLY: xt_is_null, xt_idxlist, xt_idxlist_c2f, xt_idxlist_delete, &
       xt_idxlist_is_stripe_conversion_profitable, &
       xt_config, xt_config_new, xt_config_get_xmdd_bucket_gen, &
       xt_config_set_xmdd_bucket_gen, &
       xt_idxfsection_new, &
       xt_xmdd_bucket_gen, xt_xmdd_bucket_gen_define_interface, &
       xt_xmdd_bucket_gen_new, xt_xmdd_bucket_gen_delete, &
       xt_xmdd_bucket_gen_comms, &
       xt_xmdd_bucket_gen_get_init_params, &
       xt_xmdd_bucket_gen_set_extended_init_params, &
       xt_com_list, xi => xt_int_kind
#if defined HAVE_FC_C_FUNLOC && defined HAVE_FC_ABSTRACT_INTERFACE && !defined NDEBUG
  USE yaxt, ONLY: xt_xmdd_bucket_gen_init_state, &
       xt_xmdd_bucket_gen_get_init, xt_abort
#endif
  USE iso_c_binding, ONLY: c_ptr, c_long_long, c_f_pointer, c_int, &
       c_null_ptr, c_size_t, c_char, c_null_char
#ifdef HAVE_FC_C_SIZEOF
  USE iso_c_binding, ONLY: c_sizeof
#endif
#if defined __PGI && ( __PGIC__ <= 11 )
  USE xt_idxlist_abstract, ONLY: xt_is_null
#endif
#if defined __PGI && __PGIC__ >= 22 && __PGIC__ <= 24
  USE iso_c_binding, ONLY: c_funloc, c_funptr
  USE yaxt, ONLY: xt_xmdd_bucket_gen_f2c, xt_config_c2f
#endif

#if ! defined HAVE_GFORTRAN_SIZEOF
  USE mpi
#endif
  IMPLICIT NONE
  INTEGER, PARAMETER :: max_arank = 1
  TYPE, BIND(c) :: bucket_gen_1d_tile_param
    INTEGER(xi) :: global_sizes(max_arank)
  END TYPE bucket_gen_1d_tile_param

  TYPE bucket_gen_1d_tile_state
    INTEGER(xi) :: global_sizes(max_arank)
    INTEGER :: deco(max_arank)
    INTEGER :: intra_comm_size, last_rank_generated, prev_rank_generated
    TYPE(xt_idxlist) :: prev_tile
  END TYPE bucket_gen_1d_tile_state
  CHARACTER(*), PARAMETER :: filename = "test_xmap_bucket_tiler_f.f90"

  PUBLIC
CONTAINS
    ! determine which MPI rank owns sub-section a particular coordinate falls into
  ! note: assumes row-major distribution of ranks
  FUNCTION rank_of_coord(gsz, deco, coord)
    INTEGER(xi), INTENT(in) :: gsz(max_arank), coord(max_arank)
    INTEGER, INTENT(in) :: deco(max_arank)
    INTEGER :: rank_of_coord, i, prod

    prod = 1
    rank_of_coord = 0
    DO i = 1, max_arank
      rank_of_coord = rank_of_coord + &
           INT((INT(coord(i)-1, c_long_long) &
           &    * INT(deco(i), c_long_long) &
           &    + INT(deco(i)-1, c_long_long)) &
           &              / INT(gsz(i), c_long_long)) * prod
      prod = prod * deco(i)
    END DO
    rank_of_coord = rank_of_coord
  END FUNCTION rank_of_coord

  ! nvfortran 22 to 24 exhibit incorrect handling of directly passed
  ! TYPE(x), VALUE arguments where x contains default initialized
  ! members, i.e. the value passed in from the caller is overwritten
  ! on function entry
  FUNCTION bucket_gen_1d_tile_init(gen_state_, &
       src_idxlist, dst_idxlist, config, comms, init_params_) &
       BIND(c) RESULT(convert_to_stripes)
    TYPE(c_ptr), VALUE :: gen_state_
#if defined __PGI && __PGIC__ >= 22 && __PGIC__ <= 24
    TYPE(c_ptr), VALUE :: src_idxlist, dst_idxlist, config
#else
    TYPE(xt_idxlist), VALUE :: src_idxlist, dst_idxlist
    TYPE(xt_config), VALUE :: config
#endif
    TYPE(xt_xmdd_bucket_gen_comms), INTENT(in) :: comms
    TYPE(c_ptr), VALUE :: init_params_
    TYPE(xt_idxlist) :: init_prev_tile

    TYPE(bucket_gen_1d_tile_state), POINTER :: gen_state
    TYPE(bucket_gen_1d_tile_param), POINTER :: init_params
    INTEGER :: ierror
    INTEGER(c_int) :: convert_to_stripes

#if defined __PGI && __PGIC__ >= 22 && __PGIC__ <= 24
    TYPE(xt_idxlist) :: src_idxlist_, dst_idxlist_
    TYPE(xt_config) :: config_
    src_idxlist_ = xt_idxlist_c2f(src_idxlist)
    dst_idxlist_ = xt_idxlist_c2f(dst_idxlist)
    config_ = xt_config(config)
#define src_idxlist src_idxlist_
#define dst_idxlist dst_idxlist_
#define config config_
#endif
    CALL C_F_POINTER(gen_state_, gen_state)
    CALL MPI_Comm_size(comms%intra_comm, gen_state%intra_comm_size, ierror)
    gen_state%prev_tile = init_prev_tile
    IF (xt_idxlist_is_stripe_conversion_profitable(&
         src_idxlist, config) > 0) THEN
      convert_to_stripes = 1_c_int
    ELSE IF (xt_idxlist_is_stripe_conversion_profitable(&
         dst_idxlist, config) > 0) THEN
      convert_to_stripes = 1_c_int
    ELSE
      convert_to_stripes = 0_c_int
    END IF
    ! todo: query idxlists for minimum and maximum to assess
    !  range of intersecting bucket ranks */
    ! establish local ranges
    CALL C_F_POINTER(init_params_, init_params)
    gen_state%global_sizes = init_params%global_sizes
    gen_state%last_rank_generated = gen_state%intra_comm_size-1
    gen_state%prev_rank_generated = -1
    ! todo: distribute ranks to more than one dimension
    gen_state%deco(1) = gen_state%intra_comm_size
#if defined __PGI && __PGIC__ >= 22 && __PGIC__ <= 24
#undef src_idxlist
#undef dst_idxlist
#undef config
#endif
  END FUNCTION bucket_gen_1d_tile_init

  SUBROUTINE bucket_gen_1d_tile_destroy(gen_state_) BIND(c)
    TYPE(c_ptr), VALUE :: gen_state_
    TYPE(bucket_gen_1d_tile_state), POINTER :: gen_state

    CALL C_F_POINTER(gen_state_, gen_state)

    IF (.NOT. xt_is_null(gen_state%prev_tile)) &
         CALL xt_idxlist_delete(gen_state%prev_tile)
  END SUBROUTINE bucket_gen_1d_tile_destroy

  FUNCTION bucket_gen_1d_tile_get_intersect_max_num(&
       gen_state_, bucket_type) BIND(c) RESULT(max_num)
    TYPE(c_ptr), VALUE :: gen_state_
    INTEGER(c_int), VALUE :: bucket_type
    INTEGER(c_int) :: max_num
    TYPE(bucket_gen_1d_tile_state), POINTER :: gen_state

    XT_UNUSED(bucket_type)
    CALL C_F_POINTER(gen_state_, gen_state)
    max_num = INT(gen_state%last_rank_generated &
         &        - gen_state%prev_rank_generated, c_int)
  END FUNCTION bucket_gen_1d_tile_get_intersect_max_num

  FUNCTION bucket_gen_1d_tile_next(gen_state_, bucket_type) &
       BIND(C) RESULT(bucket)
    TYPE(c_ptr), VALUE :: gen_state_
    INTEGER(c_int), VALUE :: bucket_type

    TYPE(xt_com_list) :: bucket

    TYPE(bucket_gen_1d_tile_state), POINTER :: gen_state
    INTEGER :: next_rank, lsz(max_arank), rank_mod, rank_idx, i
    INTEGER(xi) :: lst(max_arank)
    LOGICAL :: all_lsz_non_zero

    ! todo: properly handle bucket_type
    XT_UNUSED(bucket_type)

    CALL C_F_POINTER(gen_state_, gen_state)
    IF (.NOT. xt_is_null(gen_state%prev_tile)) &
         CALL xt_idxlist_delete(gen_state%prev_tile)
    next_rank = gen_state%prev_rank_generated
    DO next_rank = next_rank+1, gen_state%last_rank_generated
      gen_state%prev_rank_generated = next_rank
      rank_idx = next_rank
      all_lsz_non_zero = .TRUE.
      DO i = 1, max_arank
        rank_mod = MOD(rank_idx, gen_state%deco(i))
        rank_idx = rank_idx / gen_state%deco(i)
        lst(i) &
             = INT((INT(gen_state%global_sizes(i), c_long_long) &
             &      * INT(rank_mod, c_long_long)) &
             &     / INT(gen_state%deco(i), c_long_long) + 1_c_long_long, xi)
        lsz(i) &
             = INT(INT((INT(gen_state%global_sizes(i), c_long_long) &
             &         * INT(rank_mod+1, c_long_long)) &
             &         / INT(gen_state%deco(i), c_long_long), xi) &
             &     - (lst(i)-1_xi))
        all_lsz_non_zero = all_lsz_non_zero .AND. lsz(i) > 0
      END DO
      IF (all_lsz_non_zero) THEN
        bucket%list = xt_idxfsection_new(0_xi, &
             gen_state%global_sizes, &
             lsz, lst)
        bucket%rank = next_rank
        gen_state%prev_tile = bucket%list
        RETURN
      END IF
    END DO
    bucket%list = xt_idxlist_c2f(c_null_ptr)
    bucket%rank = -1
  END FUNCTION bucket_gen_1d_tile_next

  FUNCTION tile_bucket_config_new(ndims, global_dims, num_part) &
       RESULT(custom_bucket_config)
    INTEGER, INTENT(in) :: ndims
    INTEGER(xi), INTENT(in) :: global_dims(ndims)
    INTEGER, INTENT(in) :: num_part(ndims)
    TYPE(Xt_config) :: custom_bucket_config

    TYPE(Xt_xmdd_bucket_gen) :: bucket_gen
    TYPE(bucket_gen_1d_tile_param), POINTER :: bucket_params
    TYPE(bucket_gen_1d_tile_state) :: state_dummy(2)
#ifdef FC_C_SIZEOF
    TYPE(bucket_gen_1d_tile_param) :: bucket_params_dummy
#else
    TYPE(bucket_gen_1d_tile_param) :: bucket_params_dummy(2)
#endif
    TYPE(c_ptr) :: bucket_params_c
    INTEGER(c_size_t) :: bucket_params_c_size, state_c_size
#if ! defined HAVE_GFORTRAN_SIZEOF
    INTEGER(mpi_address_kind) :: dummy_addr(2)
    INTEGER :: ierror
#endif
#if defined __PGI && __PGIC__ >= 22 && __PGIC__ <= 24
    TYPE(c_funptr) :: bucket_gen_1d_tile_init_, bucket_gen_1d_tile_destroy_, &
         bucket_gen_1d_tile_get_intersect_max_num_, bucket_gen_1d_tile_next_
#endif
    INTEGER :: i
    CHARACTER(len=1, kind=c_char) :: filename_c(LEN(filename)+1)
    INTERFACE
      FUNCTION xt_xmalloc(char_size, source, line) BIND(c, name='Xt_xmalloc') &
           RESULT(ptr)
        IMPORT :: c_size_t, c_ptr, c_char, c_int
        INTEGER(c_size_t), VALUE, INTENT(in) :: char_size
        CHARACTER(kind=c_char, len=1), INTENT(in) :: source(*)
        INTEGER(c_int), VALUE, INTENT(in) :: line
        TYPE(c_ptr) :: ptr
      END FUNCTION xt_xmalloc
    END INTERFACE
    XT_UNUSED(num_part)
    custom_bucket_config = xt_config_new()
    bucket_gen = xt_xmdd_bucket_gen_new()
    ! todo: generalize for arbitrary number of dimensions
    IF (ndims /= 1) STOP

#ifdef FC_C_SIZEOF
    bucket_params_c_size = FC_C_SIZEOF(bucket_params_dummy)
#else
    DO i = 1, 2
      CALL mpi_get_address(bucket_params_dummy(i), dummy_addr(i), ierror)
    END DO
    bucket_params_c_size = INT(dummy_addr(2) - dummy_addr(1), c_size_t)
#endif
    DO i = 1, LEN(filename)
      filename_c(i) = filename(i:i)
    END DO
    filename_c(LEN(filename)+1) = c_null_char
    bucket_params_c = xt_xmalloc(bucket_params_c_size, filename_c, __LINE__)
    CALL C_F_POINTER(bucket_params_c, bucket_params)
    bucket_params%global_sizes(1:ndims) = global_dims
#if defined HAVE_GFORTRAN_SIZEOF
    state_c_size = SIZEOF(state_dummy(1))
#else
    DO i = 1, 2
      CALL mpi_get_address(state_dummy(i), dummy_addr(i), ierror)
    END DO
    state_c_size = INT(dummy_addr(2) - dummy_addr(1), c_size_t)
#endif
    ! unfortunately, for the above
    ! work-around, additional ISO C binding sauce is needed here
#if defined __PGI && __PGIC__ >= 22 && __PGIC__ <= 24
    bucket_gen_1d_tile_init_ = C_FUNLOC(bucket_gen_1d_tile_init)
    bucket_gen_1d_tile_destroy_ = C_FUNLOC(bucket_gen_1d_tile_destroy)
    bucket_gen_1d_tile_get_intersect_max_num_ &
         = C_FUNLOC(bucket_gen_1d_tile_get_intersect_max_num)
    bucket_gen_1d_tile_next_ = C_FUNLOC(bucket_gen_1d_tile_next)
    CALL xt_xmdd_bucket_gen_define_interface(&
         xt_xmdd_bucket_gen_f2c(bucket_gen), &
         bucket_gen_1d_tile_init_, bucket_gen_1d_tile_destroy_, &
         bucket_gen_1d_tile_get_intersect_max_num_, &
         bucket_gen_1d_tile_next_, state_c_size, bucket_params_c)
#else
    CALL xt_xmdd_bucket_gen_define_interface(bucket_gen, &
         bucket_gen_1d_tile_init, bucket_gen_1d_tile_destroy, &
         bucket_gen_1d_tile_get_intersect_max_num, &
         bucket_gen_1d_tile_next, state_c_size, bucket_params_c)
#endif
    CALL xt_xmdd_bucket_gen_set_extended_init_params(&
         bucket_gen, bucket_params_c_size)

    CALL xt_config_set_xmdd_bucket_gen(custom_bucket_config, bucket_gen)
    CALL xt_xmdd_bucket_gen_delete(bucket_gen)
  END FUNCTION tile_bucket_config_new

  SUBROUTINE tile_bucket_config_set_global_dims(custom_bucket_config, &
       global_dims)
    TYPE(xt_config) :: custom_bucket_config
    INTEGER(xi), INTENT(in) :: global_dims(*)
#if defined HAVE_FC_C_FUNLOC && defined HAVE_FC_ABSTRACT_INTERFACE \
    && !defined NDEBUG
    PROCEDURE(xt_xmdd_bucket_gen_init_state), POINTER :: ini
#endif

    TYPE(xt_xmdd_bucket_gen) ::  bucket_gen
    TYPE(c_ptr) :: init_params_c
    TYPE(bucket_gen_1d_tile_param), POINTER :: init_params

    bucket_gen = xt_config_get_xmdd_bucket_gen(custom_bucket_config)
#if defined HAVE_FC_C_FUNLOC && defined HAVE_FC_ABSTRACT_INTERFACE \
    && !defined NDEBUG
    ini => xt_xmdd_bucket_gen_get_init(bucket_gen)
    IF (.NOT. ASSOCIATED(ini, bucket_gen_1d_tile_init)) &
         CALL xt_abort("unexpected function pointer query result", &
         filename, __LINE__)
#endif
    init_params_c = xt_xmdd_bucket_gen_get_init_params(bucket_gen)
    CALL C_F_POINTER(init_params_c, init_params)
    init_params%global_sizes = global_dims(1:1)
  END SUBROUTINE tile_bucket_config_set_global_dims

END MODULE test_xmap_bucket_tiler
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
