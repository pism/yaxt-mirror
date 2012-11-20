!>
!> @file test_ut.f90
!>
!> @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
!>                                 Moritz Hanke <hanke@dkrz.de>
!>
!> @author Jörg Behrens <behrens@dkrz.de>
!>         Moritz Hanke <hanke@dkrz.de>
!>

!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
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

PROGRAM test_ut
  USE mpi
  USE xt_ut, ONLY: ut_abort, ut_init, ut_init_decomposition, ut_init_oneway_transposition_template, &
       & ut_destroy_transposition_template, ut_mode_dt_p2p, ut_init_transposition, &
       & comm_forward, ut_transpose, ut_destroy_decomposition, ut_destroy_transposition, &
       & xt_idx_kind
  IMPLICIT NONE

  ! global extents including halos
  INTEGER, PARAMETER :: g_ie = 8, g_je = 4
  LOGICAL, PARAMETER :: verbose = .FALSE.
  INTEGER, PARAMETER :: nlev = 3
  INTEGER, PARAMETER :: undef_int = HUGE(undef_int)/2 - 1
  INTEGER(xt_idx_kind), PARAMETER :: undef_index = -1
  INTEGER, PARAMETER :: nhalo = 1 ! 1dim. halo border size

  INTEGER :: ie, je ! local extents, including halos
  INTEGER :: p_ioff, p_joff ! offsets within global domain
  INTEGER :: nprocx, nprocy ! process space extents
  INTEGER :: nprocs ! == nprocx*nprocy
  INTEGER :: mype, mypx, mypy ! process rank, process coords within (0:, 0:) process space
  LOGICAL :: lroot ! true only for proc 0

  INTEGER(xt_idx_kind) :: g_id(g_ie, g_je) ! global id

  INTEGER(xt_idx_kind) :: g_tpex(g_ie, g_je) ! global "tripolar-like" toy bounds exchange
  INTEGER :: template_tpex, trans_tpex

  INTEGER(xt_idx_kind), ALLOCATABLE :: loc_id(:,:), loc_tpex(:,:)
  INTEGER, ALLOCATABLE :: fval(:,:), gval(:,:)
  INTEGER, ALLOCATABLE :: gval3d(:,:,:)
  INTEGER, ALLOCATABLE :: id_pos(:,:), pos3d_surf(:,:)
  INTEGER :: surf_trans_tpex

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

  ! template: loc_id -> loc_tpex
  CALL gen_template(loc_id, loc_tpex, template_tpex)

  ! transposition: loc_id:data -> loc_tpex:data
  CALL gen_trans(template_tpex, MPI_INTEGER, MPI_INTEGER, trans_tpex)

  ! test 2d-to-2d transposition:
  fval = loc_id
  CALL exchange2d_2d(trans_tpex, fval, gval)
  CALL icmp('2d to 2d check', gval, INT(loc_tpex))

  ! define positions of surface elements within (i,k,j) array
  CALL gen_id_pos(id_pos)
  CALL gen_id_pos(pos3d_surf)
  CALL gen_pos3d_surf(pos3d_surf)

  ! generate surface transposition:
  CALL gen_off_trans(template_tpex, MPI_INTEGER, id_pos(:,:)-1, MPI_INTEGER, pos3d_surf(:,:)-1, surf_trans_tpex)

  ! 2d to surface boundsexchange:
  gval3d = -1
  CALL exchange2d_3d(surf_trans_tpex, fval, gval3d)
  CALL icmp('surface check', gval3d(:,1,:), INT(loc_tpex))

  ! check sub surface:
  CALL icmp('sub surface check', gval3d(:,2,:), INT(loc_tpex)*0-1)

  ! cleanup:
  CALL ut_destroy_transposition_template(template_tpex)
  CALL ut_destroy_transposition(trans_tpex)
  CALL ut_destroy_transposition(surf_trans_tpex)

  CALL exit_all

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

  SUBROUTINE icmp(label, f,g)
    CHARACTER(len=*), PARAMETER :: context = 'test_ut::icmp: '
    CHARACTER(len=*), INTENT(in)  :: label
    INTEGER, INTENT(in)  :: f(:,:)
    INTEGER, INTENT(in)  :: g(:,:)

    INTEGER :: i, j, n1, n2

    n1 = SIZE(f,1)
    n2 = SIZE(f,2)
    IF (SIZE(g,1) /= n1 .OR. SIZE(g,2) /= n2) CALL ut_abort(context//'internal error', __FILE__, __LINE__)

    DO j = 1, n2
      DO i = 1, n1
        IF (f(i,j) /= g(i,j)) THEN
          WRITE(0,*) context,'test failed: i, j, f(i,j), g(i,j) =',i, j, f(i,j), g(i,j)
          CALL ut_abort(context//'test failed', __FILE__, __LINE__)
        ENDIF
      ENDDO
    ENDDO
    IF (verbose) WRITE(0,*) mype,':',context//label//' passed'
  END SUBROUTINE icmp

  SUBROUTINE init_all
    CHARACTER(len=*), PARAMETER :: context = 'init_all: '
    INTEGER :: ierror

    CALL MPI_INIT(ierror)
    IF (ierror /= MPI_SUCCESS) CALL ut_abort(context//'MPI_INIT failed', __FILE__, __LINE__)

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
    IF (ierror /= MPI_SUCCESS) CALL ut_abort(context//'MPI_COMM_SIZE failed', __FILE__, __LINE__)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierror)
    IF (ierror /= MPI_SUCCESS) CALL ut_abort(context//'MPI_COMM_RANK failed', __FILE__, __LINE__)
    IF (mype==0) THEN
      lroot = .true.
    ELSE
      lroot = .FALSE.
    ENDIF

    CALL factorize(nprocs, nprocx, nprocy)
    IF (verbose .AND. lroot) WRITE(0,*) 'nprocx, nprocy=',nprocx, nprocy
    mypy = mype / nprocx
    mypx = MOD(mype, nprocx)

    CALL ut_init(decomp_size=30, comm_tmpl_size=30, comm_size=30, &
         &       debug_lvl=0, mode=ut_mode_dt_p2p, debug_unit=0)

    CALL deco

    ALLOCATE(fval(ie,je), gval(ie,je))
    ALLOCATE(loc_id(ie,je), loc_tpex(ie,je))
    ALLOCATE(id_pos(ie,je), gval3d(ie,nlev,je), pos3d_surf(ie,je))

    fval = undef_int
    gval = undef_int
    loc_id = undef_int
    loc_tpex = undef_int
    id_pos = undef_int
    gval3d = undef_int
    pos3d_surf = undef_int

  END SUBROUTINE init_all

  SUBROUTINE exit_all
    CHARACTER(len=*), PARAMETER :: context = 'exit_all: '
    INTEGER :: ierror

    CALL MPI_FINALIZE(ierror)
    IF (ierror /= MPI_SUCCESS) CALL ut_abort(context//'MPI_FINALIZE failed', __FILE__, __LINE__)

  END SUBROUTINE exit_all

  SUBROUTINE id_map(map)
    INTEGER(xt_idx_kind), INTENT(out) :: map(:,:)

    INTEGER :: i,j
    INTEGER(xt_idx_kind) :: p

    p = 0
    DO j = 1, SIZE(map,2)
      DO i = 1, SIZE(map,1)
        p = p + 1
        map(i,j) = p
      ENDDO
    ENDDO

  END SUBROUTINE id_map

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


  SUBROUTINE exchange2d_2d(itrans, f, g)
    INTEGER, INTENT(in) :: itrans
    INTEGER, TARGET, INTENT(in) :: f(:,:)
    INTEGER, TARGET, INTENT(out) :: g(:,:)

    INTEGER, POINTER :: p_in, p_out

    p_in => f(1,1)
    p_out=> g(1,1)

    CALL ut_transpose(p_in, itrans, comm_forward, p_out)

  END SUBROUTINE exchange2d_2d

  SUBROUTINE exchange2d_3d(itrans, f, g)
    INTEGER, INTENT(in) :: itrans
    INTEGER, TARGET, INTENT(in) :: f(:,:)
    INTEGER, TARGET, INTENT(out) :: g(:,:,:)

    INTEGER, POINTER :: p_in, p_out

    p_in => f(1,1)
    p_out=> g(1,1,1)

    CALL ut_transpose(p_in, itrans, comm_forward, p_out)

  END SUBROUTINE exchange2d_3d

  SUBROUTINE gen_trans(itemp, send_dt, recv_dt, itrans)
    INTEGER,INTENT(in) :: itemp, send_dt, recv_dt
    INTEGER,INTENT(out) :: itrans

    INTEGER :: dt

    IF (send_dt /= recv_dt) CALL ut_abort('gen_trans: (send_dt /= recv_dt) unsupported', __FILE__, __LINE__)
    dt = send_dt
    CALL ut_init_transposition(itemp, dt, itrans)

  END SUBROUTINE gen_trans

  SUBROUTINE gen_off_trans(itemp, send_dt, send_off, recv_dt, recv_off, itrans)
    INTEGER,INTENT(in) :: itemp, send_dt, recv_dt
    INTEGER,INTENT(in) :: send_off(:,:), recv_off(:,:)
    INTEGER,INTENT(out) :: itrans

    INTEGER :: send_offsets(SIZE(send_off)), recv_offsets(SIZE(recv_off))

    send_offsets = RESHAPE(send_off, (/SIZE(send_off)/) )
    recv_offsets = RESHAPE(recv_off, (/SIZE(recv_off)/) )

    CALL ut_init_transposition(itemp, send_offsets, recv_offsets, send_dt, recv_dt, itrans)

  END SUBROUTINE gen_off_trans

  SUBROUTINE get_window(gval, win)
    INTEGER(xt_idx_kind), INTENT(in) :: gval(:,:)
    INTEGER(xt_idx_kind), INTENT(out) :: win(:,:)

    INTEGER :: i, j, ig, jg

    DO j = 1, je
      jg = p_joff + j
      DO i = 1, ie
        ig = p_ioff + i
        win(i,j) =  gval(ig,jg)
      ENDDO
    ENDDO

  END SUBROUTINE get_window

  SUBROUTINE gen_template(local_src_idx, local_dst_idx, ihandle)
    INTEGER(xt_idx_kind), INTENT(in) :: local_src_idx(:,:)
    INTEGER(xt_idx_kind), INTENT(in) :: local_dst_idx(:,:)
    INTEGER, INTENT(out) :: ihandle

    INTEGER(xt_idx_kind) :: src(SIZE(local_src_idx)), dst(size(local_dst_idx))

    INTEGER :: src_handle, dst_handle

    src = RESHAPE( local_src_idx, (/SIZE(local_src_idx)/) )
    dst = RESHAPE( local_dst_idx, (/SIZE(local_dst_idx)/) )

    CALL ut_init_decomposition(src, INT(g_ie, xt_idx_kind) &
         * INT(g_je, xt_idx_kind), src_handle)
    CALL ut_init_decomposition(dst, INT(g_ie, xt_idx_kind) &
         * INT(g_je, xt_idx_kind), dst_handle)

    CALL ut_init_oneway_transposition_template(src_handle, dst_handle, MPI_COMM_WORLD, ihandle)

    CALL ut_destroy_decomposition(dst_handle)
    CALL ut_destroy_decomposition(src_handle)

  END SUBROUTINE gen_template

  SUBROUTINE def_exchange(id_in, id_out)
    INTEGER(xt_idx_kind), INTENT(in) :: id_in(:,:)
    INTEGER(xt_idx_kind), INTENT(out) :: id_out(:,:)

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
    id_out(g_core_is:g_core_ie, g_core_js:g_core_je) = id_in(g_core_is:g_core_ie, g_core_js:g_core_je)

    IF (with_north_halo) THEN

      ! north inversion, (maybe with increased north halo)
      IF (increased_north_halo) THEN
        north_halo = nhalo+1
      ELSE
        north_halo = nhalo
      ENDIF

      IF (2*north_halo > g_core_je) CALL ut_abort('def_exchange: grid too small (or halo too large) for tripolar north exchange', &
           &  __FILE__, __LINE__)
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
    INTEGER(xt_idx_kind), INTENT(in) :: gidx(:,:)
    INTEGER :: i,j

    DO j=1,g_je
      DO i=1,g_ie
        IF (gidx(i,j) == undef_index) THEN
          CALL ut_abort('check_g_idx: check failed', __FILE__, __LINE__)
        ENDIF
      ENDDO
    ENDDO

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

  SUBROUTINE regular_deco(g_cn, c0, cn)
    INTEGER, INTENT(in) :: g_cn
    INTEGER, INTENT(out) :: c0(0:), cn(0:)

    ! convention: process space coords start at 0, grid point coords start at 1

    integer :: tn
    INTEGER :: d, m
    INTEGER :: it

    tn = SIZE(c0)
    IF (tn<0) CALL ut_abort('(tn<0)', __FILE__, __LINE__)
    IF (tn>g_cn) CALL ut_abort('regular_deco: too many task for such a core region', __FILE__, __LINE__)

    d = g_cn/tn
    m = MOD(g_cn, tn)

    DO it = 0, m-1
      cn(it) = d + 1
    ENDDO
    DO it = m, tn-1
      cn(it) = d
    ENDDO

    c0(0)=0
    DO it = 1, tn-1
      c0(it) = c0(it-1) + cn(it-1)
    ENDDO
    IF (c0(tn-1)+cn(tn-1) /= g_cn) CALL ut_abort('regular_deco: internal error 1', __FILE__, __LINE__)

  END SUBROUTINE regular_deco

  SUBROUTINE factorize(c, a, b)
    INTEGER, INTENT(in) :: c
    INTEGER, INTENT(out) :: a, b ! c = a*b

    INTEGER :: x0, i

    IF (c<1) CALL ut_abort('factorize: invalid process space', __FILE__, __LINE__)
    IF (c <= 3 .OR. c == 5 .OR. c == 7) THEN
      a = c
      b = 1
      RETURN
    ENDIF

    ! simple approach, we try to be near c = (2*x) * x
    x0 = INT(SQRT(0.5*c) + 0.5)
    a = 2*x0
    f_loop: DO i = a, 1, -1
      IF (MOD(c,a) == 0) THEN
        b = c/a
        EXIT f_loop
      ENDIF
      a = a - 1
    ENDDO f_loop

  END SUBROUTINE factorize

END PROGRAM test_ut
