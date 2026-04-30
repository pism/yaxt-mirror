!>
!! @file xt_xmap_dist_dir_bucket_gen2f.f90
!! @brief Fortran interface to yaxt bucket generator definition
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
MODULE xt_xmap_dist_dir_bucket_gen2
  USE iso_c_binding, ONLY: c_int, c_ptr, c_size_t
#ifdef HAVE_FC_C_FUNLOC
  USE iso_c_binding, ONLY: c_funptr
#endif
#if defined HAVE_FC_C_FUNLOC && defined HAVE_FC_ABSTRACT_INTERFACE
  USE iso_c_binding, ONLY: c_f_procpointer
#endif
  USE xt_core, ONLY: xt_mpi_fint_kind
  USE xt_config_f, ONLY: xt_config
  USE xt_xmap_intersection, ONLY: xt_com_list
  USE xt_idxlist_abstract, ONLY: xt_idxlist
  USE xt_xmap_dist_dir_bucket_gen, ONLY: xt_xmdd_bucket_gen, &
       xt_xmdd_bucket_gen_f2c
  IMPLICIT NONE
  PRIVATE

  INTEGER(c_int), PUBLIC, PARAMETER :: &
       xt_dist_dir_bucket_gen_type_send = 1, &
       xt_dist_dir_bucket_gen_type_recv = 2, &
       xt_dist_dir_bucket_gen_type_sendrecv = 3

  TYPE, BIND(c), PUBLIC :: xt_xmdd_bucket_gen_comms
    INTEGER(xt_mpi_fint_kind) :: intra_comm, inter_comm, &
         tag_offset_intra, tag_offset_inter
  END TYPE xt_xmdd_bucket_gen_comms

#ifdef HAVE_FC_ABSTRACT_INTERFACE
  ABSTRACT INTERFACE
    !  The init function sets up the generator state.
    FUNCTION xt_xmdd_bucket_gen_init_state(gen_state, src_idxlist, &
         dst_idxlist, config, comms, init_params) BIND(c) RESULT(stripify)
      IMPORT :: c_int, c_ptr, xt_idxlist, xt_config, &
           xt_xmdd_bucket_gen_comms
      INTEGER(c_int) :: stripify
      TYPE(c_ptr), VALUE :: gen_state
      TYPE(xt_idxlist), VALUE, INTENT(in) :: src_idxlist,dst_idxlist
      TYPE(xt_config), VALUE, INTENT(in) :: config
      TYPE(xt_xmdd_bucket_gen_comms), INTENT(in) :: comms
      TYPE(c_ptr), VALUE, INTENT(in) :: init_params
    END FUNCTION xt_xmdd_bucket_gen_init_state
  END INTERFACE
  PUBLIC :: xt_xmdd_bucket_gen_init_state
#endif

  ! Interface of bucket generator
  !
  ! Essentially, the generator needs to be able to enumerate all
  ! buckets used to form intersections. This means the generator only
  ! needs to produce buckets that actually can intersect and,
  ! consequently, buckets that won't intersect the requested type of
  ! list can be skipped.
  PUBLIC :: xt_xmdd_bucket_gen_define_interface
  INTERFACE xt_xmdd_bucket_gen_define_interface
    SUBROUTINE xt_xmdd_bucket_gen_define_interface_f(gen, init, destroy, &
         get_intersect_max_num, next, gen_state_size, init_params)
      IMPORT :: xt_xmdd_bucket_gen, xt_config, xt_idxlist, c_int, &
           c_ptr, c_size_t, xt_mpi_fint_kind, xt_com_list, &
           xt_xmdd_bucket_gen_comms
#ifdef HAVE_FC_ABSTRACT_INTERFACE
      IMPORT :: xt_xmdd_bucket_gen_init_state
#endif
      TYPE(xt_xmdd_bucket_gen), INTENT(inout) :: gen
#ifdef HAVE_FC_ABSTRACT_INTERFACE
      PROCEDURE(xt_xmdd_bucket_gen_init_state) :: init
      INTERFACE
#else
      INTERFACE
        !  The init function sets up the generator state.
        FUNCTION init(gen_state, src_idxlist, &
             dst_idxlist, config, comms, init_params) BIND(c) RESULT(stripify)
          IMPORT :: c_int, c_ptr, xt_idxlist, xt_config, &
               xt_xmdd_bucket_gen_comms
          INTEGER(c_int) :: stripify
          TYPE(c_ptr), VALUE :: gen_state
          TYPE(xt_idxlist), VALUE, INTENT(in) :: src_idxlist,dst_idxlist
          TYPE(xt_config), VALUE, INTENT(in) :: config
          TYPE(xt_xmdd_bucket_gen_comms), intent(in) :: comms
          TYPE(c_ptr), VALUE, INTENT(in) :: init_params
        END FUNCTION init
#endif
        !  The destroy function clean up the generator state. Can be zero
        !  if no cleaning is needed.
        SUBROUTINE destroy(gen_state) BIND(c)
          IMPORT :: c_ptr
          TYPE(c_ptr), VALUE, INTENT(in) :: gen_state
        END SUBROUTINE destroy

        FUNCTION get_intersect_max_num(gen_state, &
             bucket_type) BIND(C) RESULT(max_num)
          IMPORT :: c_ptr, c_int
          TYPE(c_ptr), VALUE :: gen_state
          INTEGER(c_int), VALUE :: bucket_type
          INTEGER(c_int) :: max_num
        END FUNCTION get_intersect_max_num

        ! The next function returns the next bucket and corresponding rank
        ! (ranks can be skipped when the intersection will be empty
        ! anyway).
        ! Any previously returned buckets become invalid.
        FUNCTION next(gen_state, bucket_type) &
             BIND(C) RESULT(bucket)
          IMPORT :: c_ptr, c_int, xt_com_list
          TYPE(c_ptr), VALUE :: gen_state
          INTEGER(c_int), VALUE :: bucket_type
          TYPE(xt_com_list) :: bucket
        END FUNCTION next
      END INTERFACE
          ! gen_state_size is the size of the generator state
      !
      ! The distributed directory will provide memory aligned to
      ! TYPE(c_ptr) variables of this size.
      INTEGER(c_size_t), INTENT(in) :: gen_state_size
      TYPE(c_ptr), INTENT(in) :: init_params
    END SUBROUTINE xt_xmdd_bucket_gen_define_interface_f
#ifdef HAVE_FC_C_FUNLOC
    SUBROUTINE xt_xmdd_bucket_gen_def_if_f2c(gen, init, destroy, &
         get_intersect_max_num, next, gen_state_size, init_params) &
         BIND(c, name='xt_xmdd_bucket_gen_def_if_f2c')
      IMPORT :: c_funptr, c_size_t, c_ptr
      IMPLICIT NONE
      TYPE(c_ptr), VALUE :: gen
      TYPE(c_funptr), VALUE :: init, destroy, get_intersect_max_num, next
      INTEGER(c_size_t), VALUE :: gen_state_size
      TYPE(c_ptr), INTENT(in) :: init_params
    END SUBROUTINE xt_xmdd_bucket_gen_def_if_f2c
#endif
  END INTERFACE xt_xmdd_bucket_gen_define_interface
  PUBLIC :: xt_xmdd_bucket_gen_get_init_params, &
       xt_xmdd_bucket_gen_set_init_params

  INTERFACE xt_xmdd_bucket_gen_set_extended_init_params
    SUBROUTINE xt_xmddbg_set_ext_init_params_f(gen, &
         iparm_size, destroy, cpy)
      IMPORT :: xt_xmdd_bucket_gen, c_ptr, c_size_t
      TYPE(xt_xmdd_bucket_gen), INTENT(inout) :: gen
      INTEGER(c_size_t), INTENT(in) :: iparm_size
      INTERFACE
        SUBROUTINE destroy(init_params, init_params_size)
          IMPORT :: c_ptr, c_size_t
          TYPE(c_ptr), VALUE, INTENT(in) :: init_params
          INTEGER(c_size_t), VALUE, INTENT(in) :: init_params_size
        END SUBROUTINE destroy
        SUBROUTINE cpy(dest_init_params, src_init_params, init_params_size)
          IMPORT :: c_ptr, c_size_t
          TYPE(c_ptr), VALUE, INTENT(in) :: dest_init_params, src_init_params
          INTEGER(c_size_t), VALUE, INTENT(in) :: init_params_size
        END SUBROUTINE cpy
      END INTERFACE
      OPTIONAL :: destroy, cpy
    END SUBROUTINE xt_xmddbg_set_ext_init_params_f
  END INTERFACE xt_xmdd_bucket_gen_set_extended_init_params
  PUBLIC :: xt_xmdd_bucket_gen_set_extended_init_params

  INTERFACE xt_xmdd_bucket_gen_set_init
    SUBROUTINE xt_xmdd_bucket_gen_set_init_f(gen, init)
      IMPORT :: xt_xmdd_bucket_gen, xt_xmdd_bucket_gen_comms
#ifdef HAVE_FC_ABSTRACT_INTERFACE
      IMPORT :: xt_xmdd_bucket_gen_init_state
#else
      IMPORT :: c_int, c_ptr, xt_idxlist, xt_config
#endif
      TYPE(xt_xmdd_bucket_gen), INTENT(in) :: gen
#ifdef HAVE_FC_ABSTRACT_INTERFACE
      PROCEDURE(xt_xmdd_bucket_gen_init_state) :: init
#else
      INTERFACE
        FUNCTION init(gen_state, src_idxlist, &
             dst_idxlist, config, comms, init_params) BIND(c) RESULT(stripify)
          IMPORT :: c_int, c_ptr, xt_idxlist, xt_config, &
               xt_xmdd_bucket_gen_comms
          INTEGER(c_int) :: stripify
          TYPE(c_ptr), VALUE :: gen_state
          TYPE(xt_idxlist), VALUE, INTENT(in) :: src_idxlist,dst_idxlist
          TYPE(xt_config), VALUE, INTENT(in) :: config
          TYPE(xt_xmdd_bucket_gen_comms), INTENT(in) :: comms
          TYPE(c_ptr), VALUE, INTENT(in) :: init_params
        END FUNCTION init
      END INTERFACE
#endif
    END SUBROUTINE xt_xmdd_bucket_gen_set_init_f
  END INTERFACE xt_xmdd_bucket_gen_set_init
  PUBLIC :: xt_xmdd_bucket_gen_set_init

#if defined HAVE_FC_C_FUNLOC && defined HAVE_FC_ABSTRACT_INTERFACE
  PUBLIC :: xt_xmdd_bucket_gen_get_init
#endif

CONTAINS


  FUNCTION xt_xmdd_bucket_gen_get_init_params(gen) &
       RESULT(init_params)
    TYPE(xt_xmdd_bucket_gen), INTENT(in) :: gen
    TYPE(c_ptr) :: init_params
    INTERFACE
      FUNCTION xt_xmdd_bg_get_init_params_c(gen) &
           BIND(c, name='xt_xmdd_bucket_gen_get_init_params') &
           RESULT(init_params)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: gen
        TYPE(c_ptr) :: init_params
      END FUNCTION xt_xmdd_bg_get_init_params_c
    END INTERFACE
    init_params = xt_xmdd_bg_get_init_params_c(&
         xt_xmdd_bucket_gen_f2c(gen))
  END FUNCTION xt_xmdd_bucket_gen_get_init_params

  SUBROUTINE xt_xmdd_bucket_gen_set_init_params(gen, init_params)
    TYPE(xt_xmdd_bucket_gen), INTENT(in) :: gen
    TYPE(c_ptr), INTENT(in) :: init_params
    INTERFACE
      SUBROUTINE xt_xmdd_bg_set_init_params_c(gen, init_params) &
           BIND(c, name='xt_xmdd_bucket_gen_set_init_params')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: gen
        TYPE(c_ptr), VALUE, INTENT(in) :: init_params
      END SUBROUTINE xt_xmdd_bg_set_init_params_c
    END INTERFACE
    CALL xt_xmdd_bg_set_init_params_c(xt_xmdd_bucket_gen_f2c(gen), &
         init_params)
  END SUBROUTINE xt_xmdd_bucket_gen_set_init_params

#if defined HAVE_FC_C_FUNLOC && defined HAVE_FC_ABSTRACT_INTERFACE
  FUNCTION xt_xmdd_bucket_gen_get_init(gen) RESULT(ini)
    TYPE(Xt_xmdd_bucket_gen), INTENT(in) :: gen
    PROCEDURE(xt_xmdd_bucket_gen_init_state), POINTER :: ini
    TYPE(c_funptr) :: ini_c

    INTERFACE
      FUNCTION xt_xmdd_bucket_gen_get_init_c(gen) &
           BIND(c, name='xt_xmdd_bucket_gen_get_init') RESULT(ini_ptr)
        IMPORT :: c_funptr, c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: gen
        TYPE(c_funptr) :: ini_ptr
      END FUNCTION xt_xmdd_bucket_gen_get_init_c
    END INTERFACE

    ini_c = xt_xmdd_bucket_gen_get_init_c(xt_xmdd_bucket_gen_f2c(gen))
    CALL C_F_PROCPOINTER(ini_c, ini)
  END FUNCTION xt_xmdd_bucket_gen_get_init
#endif


END MODULE xt_xmap_dist_dir_bucket_gen2
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
