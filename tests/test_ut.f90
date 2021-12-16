!>
!! @file test_ut.f90
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
PROGRAM test_ut
  USE iso_c_binding, ONLY: C_LOC
  USE mpi
  USE yaxt, ONLY: xt_finalize, xt_idxlist, xt_idxvec_new, xt_xmap, &
       xt_xmap_all2all_new, xt_redist, xt_redist_p2p_new, &
       xt_redist_p2p_off_new, xt_redist_s_exchange1, xt_idxlist_delete, &
       xt_xmap_delete, xt_redist_delete, xt_int_kind, xt_initialize
  USE ftest_common, ONLY: finish_mpi, icmp, id_map, factorize, regular_deco, &
       test_abort
  IMPLICIT NONE

  ! global extents including halos
  INTEGER, PARAMETER :: g_ie = 8, g_je = 4
  LOGICAL, PARAMETER :: verbose = .FALSE.
  INTEGER, PARAMETER :: nlev = 3
  INTEGER, PARAMETER :: undef_int = g_ie * g_je * nlev + 1
  INTEGER(xt_int_kind), PARAMETER :: undef_index = -1_xt_int_kind
  INTEGER, PARAMETER :: nhalo = 1 ! 1dim. halo border size
  CHARACTER(len=*), PARAMETER :: filename = 'test_ut.f90'

  INTEGER :: ie, je ! local extents, including halos
  INTEGER :: p_ioff, p_joff ! offsets within global domain
  INTEGER :: nprocx, nprocy ! process space extents
  INTEGER :: nprocs ! == nprocx*nprocy
  ! process rank, process coords within (0:, 0:) process space
  INTEGER :: mype, mypx, mypy
  LOGICAL :: lroot ! true only for proc 0

  INTEGER(xt_int_kind) :: g_id(g_ie, g_je) ! global id

  ! global "tripolar-like" toy bounds exchange
  INTEGER(xt_int_kind) :: g_tpex(g_ie, g_je)

  INTEGER(xt_int_kind), ALLOCATABLE :: loc_id(:,:), loc_tpex(:,:)
  INTEGER, ALLOCATABLE :: fval(:,:), gval(:,:)
  INTEGER, ALLOCATABLE :: gval3d(:,:,:)
  INTEGER, ALLOCATABLE :: id_pos(:,:), pos3d_surf(:,:)

  ! yaxt data:
  TYPE(xt_idxlist) :: start_idxlist, end_idxlist
  TYPE(xt_xmap) :: tpex_xmap
  TYPE(xt_redist) :: tpex_redist, surf_tpex_redist

  ! mpi & decomposition & allocate mem:
  CALL init_all

  ! full global index space:
  CALL id_map(g_id)

  ! local window of global index space:
  CALL get_window(g_id, loc_id)

  ! define bounds exchange for full global index space
  CALL def_exchange(g_id, g_tpex)

  ! local window of global bounds exchange:
  CALL get_window(g_tpex, loc_tpex)

  ! gen xmap: loc_id -> loc_tpex:
  start_idxlist = xt_idxvec_new(loc_id)
  end_idxlist = xt_idxvec_new(loc_tpex)
  tpex_xmap = xt_xmap_all2all_new(start_idxlist, end_idxlist, MPI_COMM_WORLD)

  ! gen transposition with default offsets:
  tpex_redist = xt_redist_p2p_new(tpex_xmap,  MPI_INTEGER)

  ! test 2d-to-2d bounds exchange:
  fval = INT(loc_id)
  gval = undef_int
  CALL exchange2d_2d(tpex_redist, fval, gval)
  CALL icmp('2d to 2d check', gval, INT(loc_tpex), mype)

  ! define positions of surface elements within (i,k,j) array
  CALL gen_id_pos(id_pos)
  CALL gen_id_pos(pos3d_surf)
  CALL gen_pos3d_surf(pos3d_surf)

  ! generate surface redist for (i,j)-source data and (i,k,j)-target data:
  surf_tpex_redist = xt_redist_p2p_off_new(tpex_xmap, id_pos(:,:)-1, pos3d_surf(:,:)-1, MPI_INTEGER)

  ! test surface bounds exchange:
  gval3d = -1
  CALL exchange2d_3d(surf_tpex_redist, fval, gval3d)
  CALL icmp('surface check', gval3d(:,1,:), INT(loc_tpex), mype)

  ! check sub surface:
  CALL icmp('sub surface check', gval3d(:,2,:), INT(loc_tpex)*0-1, mype)

  ! cleanup:
  CALL xt_idxlist_delete(start_idxlist)
  CALL xt_idxlist_delete(end_idxlist)
  CALL xt_xmap_delete(tpex_xmap)
  CALL xt_redist_delete(tpex_redist)

  DEALLOCATE(fval, gval, loc_id, loc_tpex, id_pos, gval3d, pos3d_surf)

  ! finalize:
  CALL xt_finalize
  CALL finish_mpi

CONTAINS

  SUBROUTINE gen_pos3d_surf(pos)
    INTEGER, INTENT(inout) :: pos(:,:)
    ! positions for zero based arrays (ECHAM grid point dim order)
    ! old pos = i + j*ie
    ! new pos = i + k*ie + j*ie*nlev
    INTEGER :: ii,jj, i,j,k, p,q

    k = 0 ! surface
    DO jj=1,je
      DO ii=1,ie
        p = pos(ii,jj) - 1 ! shift to 0-based index
        j = p/ie
        i = MOD(p,ie)
        q = i + k*ie + j*ie*nlev
        pos(ii,jj) = q + 1 ! shift to 1-based index
      ENDDO
    ENDDO

  END SUBROUTINE gen_pos3d_surf

  SUBROUTINE init_all
    CHARACTER(len=*), PARAMETER :: context = 'init_all: '
    INTEGER :: ierror

    CALL MPI_INIT(ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_INIT failed', &
         filename, __LINE__)

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_COMM_SIZE failed', &
         filename, __LINE__)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_COMM_RANK failed', &
         filename, __LINE__)
    IF (mype==0) THEN
      lroot = .true.
    ELSE
      lroot = .FALSE.
    ENDIF

    CALL factorize(nprocs, nprocx, nprocy)
    IF (verbose .AND. lroot) WRITE(0,*) 'nprocx, nprocy=',nprocx, nprocy
    mypy = mype / nprocx
    mypx = MOD(mype, nprocx)

    CALL xt_initialize(MPI_COMM_WORLD)

    CALL deco

    ALLOCATE(fval(ie,je), gval(ie,je))
    ALLOCATE(loc_id(ie,je), loc_tpex(ie,je))
    ALLOCATE(id_pos(ie,je), gval3d(ie,nlev,je), pos3d_surf(ie,je))

    fval = undef_int
    gval = undef_int
    loc_id = INT(undef_int, xt_int_kind)
    loc_tpex = INT(undef_int, xt_int_kind)
    id_pos = undef_int
    gval3d = undef_int
    pos3d_surf = undef_int

  END SUBROUTINE init_all

  SUBROUTINE gen_id_pos(pos)
    INTEGER, INTENT(out) :: pos(:,:)

    INTEGER :: i,j
    INTEGER :: p

    p = 0
    DO j = 1, SIZE(pos,2)
      DO i = 1, SIZE(pos,1)
        p = p + 1
        pos(i,j) = p
      ENDDO
    ENDDO

  END SUBROUTINE gen_id_pos

  SUBROUTINE exchange2d_2d(redist, start_data, end_data)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, TARGET, INTENT(in) :: start_data(:,:)
    INTEGER, TARGET, INTENT(out) :: end_data(:,:)

    CALL xt_redist_s_exchange1(redist, C_LOC(start_data), C_LOC(end_data))

  END SUBROUTINE exchange2d_2d

  SUBROUTINE exchange2d_3d(redist, src_data, target_data)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, TARGET, INTENT(in) :: src_data(:,:)
    INTEGER, TARGET, INTENT(out) :: target_data(:,:,:)

    CALL xt_redist_s_exchange1(redist, C_LOC(src_data), C_LOC(target_data))

  END SUBROUTINE exchange2d_3d

  SUBROUTINE get_window(gval, win)
    INTEGER(xt_int_kind), INTENT(in) :: gval(:,:)
    INTEGER(xt_int_kind), INTENT(out) :: win(:,:)

    INTEGER :: i, j, ig, jg

    DO j = 1, je
      jg = p_joff + j
      DO i = 1, ie
        ig = p_ioff + i
        win(i,j) =  gval(ig,jg)
      ENDDO
    ENDDO

  END SUBROUTINE get_window

  SUBROUTINE def_exchange(id_in, id_out)
    INTEGER(xt_int_kind), INTENT(in) :: id_in(:,:)
    INTEGER(xt_int_kind), INTENT(out) :: id_out(:,:)

    LOGICAL, PARAMETER :: increased_north_halo = .FALSE.
    LOGICAL, PARAMETER :: with_north_halo = .true.
    INTEGER :: i, j
    INTEGER :: g_core_is, g_core_ie, g_core_js, g_core_je
    INTEGER :: north_halo

    ! global core domain:
    g_core_is = nhalo + 1
    g_core_ie = g_ie-nhalo
    g_core_js = nhalo + 1
    g_core_je = g_je-nhalo

    ! global tripolar boundsexchange:
    id_out = undef_index
    id_out(g_core_is:g_core_ie, g_core_js:g_core_je) &
         = id_in(g_core_is:g_core_ie, g_core_js:g_core_je)

    IF (with_north_halo) THEN

      ! north inversion, (maybe with increased north halo)
      IF (increased_north_halo) THEN
        north_halo = nhalo+1
      ELSE
        north_halo = nhalo
      ENDIF

      IF (2*north_halo > g_core_je) &
           CALL test_abort('def_exchange: grid too small (or halo too large)&
           & for tripolar north exchange', &
           filename, __LINE__)
      DO j = 1, north_halo
        DO i = g_core_is, g_core_ie
          id_out(i,j) = id_out(g_core_ie + (g_core_is-i), 2*north_halo + (1-j))
        ENDDO
      ENDDO

    ELSE

      DO j = 1, nhalo
        DO i = nhalo+1, g_ie-nhalo
          id_out(i,j) = id_in(i,j)
        ENDDO
      ENDDO

    ENDIF

    ! south: no change
    DO j = g_core_je+1, g_je
      DO i = nhalo+1, g_ie-nhalo
        id_out(i,j) = id_in(i,j)
      ENDDO
    ENDDO

    ! PBC
    DO j = 1, g_je
      DO i = 1, nhalo
        id_out(g_core_is-i,j) = id_out(g_core_ie+(1-i),j)
      ENDDO
      DO i = 1, nhalo
        id_out(g_core_ie+i,j) = id_out(nhalo+i,j)
      ENDDO
    ENDDO

    CALL check_g_idx(id_out)

  END SUBROUTINE def_exchange

  SUBROUTINE check_g_idx(gidx)
    INTEGER(xt_int_kind), INTENT(in) :: gidx(:,:)

    IF (ANY(gidx == undef_index)) THEN
      CALL test_abort('check_g_idx: check failed', filename, __LINE__)
    ENDIF
  END SUBROUTINE check_g_idx

  SUBROUTINE deco
    INTEGER :: cx0(0:nprocx-1), cxn(0:nprocx-1)
    INTEGER :: cy0(0:nprocy-1), cyn(0:nprocy-1)

    CALL regular_deco(g_ie-2*nhalo, cx0, cxn)
    CALL regular_deco(g_je-2*nhalo, cy0, cyn)

    ! process local deco variables:
    ie = cxn(mypx) + 2*nhalo
    je = cyn(mypy) + 2*nhalo
    p_ioff = cx0(mypx)
    p_joff = cy0(mypy)

  END SUBROUTINE deco

END PROGRAM test_ut
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
