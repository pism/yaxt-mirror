
PROGRAM test_perf_stripes
  USE mpi
  USE yaxt, ONLY: xt_abort, Xt_idxlist, xt_idxlist_delete, xt_idxvec_new, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, xt_redist, &
       xt_redist_p2p_new, xt_redist_delete, xt_redist_p2p_off_new, &
       xt_redist_p2p_blocks_off_new, xt_redist_s_exchange1, &
       xt_stripe, xt_idxstripes_new, xt_idx_kind, xt_count_kind

  USE ISO_C_BINDING, ONLY: c_int, c_loc

  IMPLICIT NONE
  ! global extents including halos:

  INTEGER, PARAMETER :: nlev = 20
  INTEGER, PARAMETER :: undef_int = HUGE(undef_int)/2 - 1
  INTEGER(xt_idx_kind), PARAMETER :: undef_index = -1
  INTEGER, PARAMETER :: nhalo = 1 ! 1dim. halo border size
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 307)

  INTEGER, PARAMETER :: grid_kind_test = 1
  INTEGER, PARAMETER :: grid_kind_toy  = 2
  INTEGER, PARAMETER :: grid_kind_tp10 = 3
  INTEGER, PARAMETER :: grid_kind_tp04 = 4
  INTEGER, PARAMETER :: grid_kind_tp6M = 5
  INTEGER :: grid_kind = grid_kind_test

  CHARACTER(len=10) :: grid_label
  INTEGER :: g_ie, g_je ! global domain extents
  INTEGER :: ie, je, ke ! local extents, including halos
  INTEGER :: p_ioff, p_joff ! offsets within global domain
  INTEGER :: nprocx, nprocy ! process space extents
  INTEGER :: nprocs ! == nprocx*nprocy
  ! process rank, process coords within (0:, 0:) process space
  INTEGER :: mype, mypx, mypy
  LOGICAL :: lroot ! true only for proc 0

  INTEGER, ALLOCATABLE :: g_id(:,:) ! global id
  ! global "tripolar-like" toy bounds exchange
  INTEGER, ALLOCATABLE :: g_tpex(:, :)
  TYPE(xt_xmap) :: xmap_tpex_2d, xmap_tpex_3d, xmap_tpex_3d_ws
  TYPE(xt_redist) :: redist_tpex_2d, redist_surf_tpex_2d, redist_tpex_3d, &
       redist_tpex_3d_ws, redist_tpex_3d_wb
  TYPE(Xt_idxlist) :: loc_id_3d_ws, loc_tpex_3d_ws

  INTEGER(xt_idx_kind), ALLOCATABLE :: loc_id_2d(:,:), loc_tpex_2d(:,:)
  INTEGER(xt_idx_kind), ALLOCATABLE :: loc_id_3d(:,:,:), loc_tpex_3d(:,:,:)
  INTEGER, ALLOCATABLE :: fval_2d(:,:), gval_2d(:,:)
  INTEGER, ALLOCATABLE :: fval_3d(:,:,:), gval_3d(:,:,:)
  INTEGER, ALLOCATABLE :: id_pos(:,:), pos3d_surf(:,:)
  LOGICAL, PARAMETER :: full_test = .TRUE.
  LOGICAL, PARAMETER :: debug =.FALSE.
  LOGICAL :: verbose

  TYPE timer
    CHARACTER(len=30) :: label = 'undef'
    INTEGER  :: istate  = -1
    REAL(dp) :: t0      = 0.0_dp
    REAL(dp) :: dt_work = 0.0_dp
  END TYPE timer

  REAL(dp) :: sync_dt_sum = 0.0_dp

  TYPE(timer) :: t_all, t_surf_redist, t_exch_surf
  TYPE(timer) :: t_xmap_2d, t_redist_2d, t_exch_2d
  TYPE(timer) :: t_xmap_3d, t_redist_3d, t_exch_3d
  TYPE(timer) :: t_xmap_3d_ws, t_redist_3d_ws, t_exch_3d_ws, t_exch_3d_wb
  TYPE(timer) :: t_redist_3d_wb

  !WRITE(0,*) '(debug) test_perf_stripes: verbose=', verbose

  CALL treset(t_all, 'all')
  CALL treset(t_surf_redist, 'surf_redist')
  CALL treset(t_exch_surf, 'exch_surf')
  CALL treset(t_xmap_2d, 'xmap_2d')
  CALL treset(t_redist_2d, 'redist_2d')
  CALL treset(t_exch_2d, 'exch_2d')
  CALL treset(t_xmap_3d, 'xmap_3d')
  CALL treset(t_redist_3d, 'redist_3d')
  CALL treset(t_exch_3d, 'exch_3d')

  CALL treset(t_xmap_3d_ws, 'xmap_3d_ws')
  CALL treset(t_redist_3d_ws, 'redist_3d_ws')
  CALL treset(t_exch_3d_ws, 'exch_3d_ws')

  CALL treset(t_redist_3d_wb, 'redist_3d_wb')
  CALL treset(t_exch_3d_wb, 'exch_3d_wb')

  CALL init_mpi

  CALL tstart(t_all)

  ! mpi & decomposition & allocate mem:
  CALL init_all

  ALLOCATE(fval_3d(nlev,ie,je), gval_3d(nlev,ie,je))

  ! full global index space:
  CALL id_map(g_id)

  ! local window of global index space:
  CALL get_window(g_id, loc_id_2d)

  ! define bounds exchange for full global index space
  CALL def_exchange()
  !g_tpex = g_id

  ! local window of global bounds exchange:
  CALL get_window(g_tpex, loc_tpex_2d)

  IF (full_test) THEN
  ! xmap: loc_id_2d -> loc_tpex_2d
  CALL tstart(t_xmap_2d)
  CALL gen_xmap_2d(loc_id_2d, loc_tpex_2d, xmap_tpex_2d)
  CALL tstop(t_xmap_2d)

  ! transposition: loc_id_2d:data -> loc_tpex_2d:data
  CALL tstart(t_redist_2d)
  CALL gen_redist(xmap_tpex_2d, MPI_INTEGER, MPI_INTEGER, redist_tpex_2d)
  CALL tstop(t_redist_2d)

  ! test 2d-to-2d transposition:
  fval_2d = loc_id_2d
  CALL tstart(t_exch_2d)
  CALL exchange_int(redist_tpex_2d, fval_2d, gval_2d)
  CALL tstop(t_exch_2d)
  CALL icmp_2d('2d to 2d check', gval_2d, INT(loc_tpex_2d))

  ! define positions of surface elements within (k,i,j) array
  CALL id_map(id_pos)
  CALL id_map(pos3d_surf)
  CALL gen_pos3d_surf(pos3d_surf)

  ! generate surface transposition:
  CALL tstart(t_surf_redist)
  CALL gen_off_redist(xmap_tpex_2d, MPI_INTEGER, id_pos(:,:)-1, &
       MPI_INTEGER, pos3d_surf(:,:)-1, redist_surf_tpex_2d)
  CALL tstop(t_surf_redist)

  ! 2d to surface boundsexchange:
  gval_3d = -1
  CALL tstart(t_exch_surf)
  CALL exchange_int(redist_surf_tpex_2d, fval_2d, gval_3d)
  CALL tstop(t_exch_surf)

  ! check surface:
  CALL icmp_2d('surface check', gval_3d(1,:,:), INT(loc_tpex_2d))

  IF (nlev>1) THEN
    ! check sub surface:
    CALL icmp_2d('sub surface check', gval_3d(2,:,:), INT(loc_tpex_2d)*0-1)
  ENDIF
  endif

  ! inflate (i,j) -> (k,i,j)
  CALL inflate_idx(1, loc_id_2d, loc_id_3d)
  CALL inflate_idx(1, loc_tpex_2d, loc_tpex_3d)

  IF (full_test) THEN

  ! xmap: loc_id_3d -> loc_tpex_3d
  CALL tstart(t_xmap_3d)
  CALL gen_xmap_3d(loc_id_3d, loc_tpex_3d, xmap_tpex_3d)
  CALL tstop(t_xmap_3d)

  ! transposition: loc_id_3d:data -> loc_tpex_3d:data
  CALL tstart(t_redist_3d)
  CALL gen_redist(xmap_tpex_3d, MPI_INTEGER, MPI_INTEGER, redist_tpex_3d)
  CALL tstop(t_redist_3d)

  CALL xt_xmap_delete(xmap_tpex_3d)

  ! test 3d-to-3d transposition:
  !DEALLOCATE(gval_3d)
  fval_3d = loc_id_3d
  gval_3d = -1
  CALL tstart(t_exch_3d)
  CALL exchange_int(redist_tpex_3d, fval_3d, gval_3d)
  CALL tstop(t_exch_3d)

  CALL xt_redist_delete(redist_tpex_3d)

  ! check 3d exchange:
  CALL icmp_3d('3d to 3d check', gval_3d, INT(loc_tpex_3d))
  endif


  ! gen stripes, xmap, redist:
  CALL gen_stripes(loc_id_3d, loc_id_3d_ws)
  CALL gen_stripes(loc_tpex_3d, loc_tpex_3d_ws)

  xmap_tpex_3d_ws = xt_xmap_all2all_new(loc_id_3d_ws, loc_tpex_3d_ws, &
       MPI_COMM_WORLD)

  CALL tstart(t_redist_3d_ws)
  redist_tpex_3d_ws  = xt_redist_p2p_new(xmap_tpex_3d_ws, MPI_INTEGER)
  CALL tstop(t_redist_3d_ws)

  ! test redist_tpex_3d_ws:
  fval_3d = loc_id_3d
  gval_3d = -1

  CALL tstart(t_exch_3d_ws)
  CALL exchange_int(redist_tpex_3d_ws, fval_3d, gval_3d)
  CALL tstop(t_exch_3d_ws)
  if (full_test) then
  ! check 3d exchange:
    CALL icmp_3d('3d to 3d check (using stripes)', gval_3d, INT(loc_tpex_3d))
  endif

  CALL tstart(t_redist_3d_wb)
  CALL gen_redist_3d_wb(xmap_tpex_2d, MPI_INTEGER, redist_tpex_3d_wb)
  CALL tstop(t_redist_3d_wb)

  fval_3d = loc_id_3d
  gval_3d = -1
  CALL tstart(t_exch_3d_wb)
  CALL exchange_int(redist_tpex_3d_wb, fval_3d, gval_3d)
  CALL tstop(t_exch_3d_wb)

  CALL xt_redist_delete(redist_tpex_3d_wb)

  CALL icmp_3d('redist_tpex_3d_wb check', gval_3d, INT(loc_tpex_3d))

  ! cleanup:
  IF (full_test) THEN
    CALL xt_xmap_delete(xmap_tpex_2d)
    CALL xt_redist_delete(redist_tpex_2d)
    CALL xt_redist_delete(redist_surf_tpex_2d)
  ENDIF

  CALL tstop(t_all)

  IF (verbose) WRITE(0,*) 'timer report for nprocs=',nprocs

  CALL treport(t_all)
  CALL treport(t_surf_redist)
  CALL treport(t_exch_surf)
  CALL treport(t_xmap_2d)
  CALL treport(t_redist_2d)
  CALL treport(t_exch_2d)
  CALL treport(t_xmap_3d)
  CALL treport(t_redist_3d)
  CALL treport(t_exch_3d)

  CALL treport(t_xmap_3d_ws)
  CALL treport(t_redist_3d_ws)
  CALL treport(t_exch_3d_ws)

  CALL treport(t_redist_3d_wb)
  CALL treport(t_exch_3d_wb)

  CALL exit_mpi

CONTAINS

  SUBROUTINE gen_redist_3d_wb(xmap_2d, dt, redist_3d)
    TYPE(xt_xmap), INTENT(in) :: xmap_2d
    INTEGER, INTENT(in) :: dt
    TYPE(xt_redist), INTENT(out) :: redist_3d

    INTEGER :: block_disp(ie,je), block_size(ie,je)
    INTEGER :: i, j
    ! data(k,i,j)
    DO j = 1, je
      DO i = 1, ie
        block_disp(i,j) = ( (j-1) * ie + i - 1 ) * nlev
        block_size(i,j) =  nlev
      ENDDO
    ENDDO
    !WRITE(0,*) '(gen_redist_3d_wb) call redist with field sizes =',ie*je
    redist_3d = xt_redist_p2p_blocks_off_new(xmap_2d, block_disp, block_size, &
         SIZE(block_size), block_disp, block_size, SIZE(block_size),dt)

  END SUBROUTINE gen_redist_3d_wb

  SUBROUTINE msg(s)
    CHARACTER(len=*), INTENT(in) :: s
    IF (verbose) WRITE(0,*) s
  END SUBROUTINE msg

  SUBROUTINE treset(t, label)
    TYPE(timer), INTENT(inout) :: t
    CHARACTER(len=*), INTENT(in) :: label
    t%label   = label
    t%istate  = 0
    t%t0      = 0.0_dp
    t%dt_work = 0.0_dp
  END SUBROUTINE treset

  SUBROUTINE tstart(t)
    TYPE(timer), INTENT(inout) :: t
    IF (debug) WRITE(0,*) 'tstart: ',t%label
    CALL sync
    t%istate = 1
    t%t0 = work_time()
  END SUBROUTINE tstart

  SUBROUTINE tstop(t)
    TYPE(timer), INTENT(inout) :: t
    REAL(dp) :: t1
    IF (debug) WRITE(0,*) 'tstop: ',t%label
    t1 = work_time()
    t%dt_work = t%dt_work + (t1 - t%t0)
    t%istate = 0
    CALL sync

  END SUBROUTINE tstop

  SUBROUTINE treport(t)
    CHARACTER(len=*), PARAMETER :: context = 'treport: '
    TYPE(timer), INTENT(in) :: t

    REAL(dp) :: work_sum, work_max, work_avg, e
    REAL(dp) :: sbuf, rbuf(0:nprocs-1)
    INTEGER :: ierror

    sbuf = t%dt_work
    rbuf = -1.0_dp
    CALL MPI_GATHER(sbuf, 1, MPI_DOUBLE_PRECISION, &
         &  rbuf, 1, MPI_DOUBLE_PRECISION, &
         &  0, MPI_COMM_WORLD, ierror)


    IF (lroot) THEN
      IF (rbuf(0) /= sbuf) CALL xt_abort(context//'internal error (1)', &
           __FILE__, &
           __LINE__)
      IF (ANY(rbuf < 0.0_dp)) CALL xt_abort(context//'internal error (2)', &
           __FILE__, &
           __LINE__)
      work_sum = SUM(rbuf)
      work_max = MAXVAL(rbuf)
      work_avg = work_sum / nprocs
      e = work_avg / (work_max+1.e-20)

      IF (verbose) WRITE(0,'(A,I4,2X,A16,3F18.8)') &
           'nprocs, label, wmax, wavg, e =', &
           nprocs, TRIM(grid_label)//':'//t%label, &
           work_max, work_avg, e
    ENDIF

  END SUBROUTINE treport


  REAL(dp) FUNCTION work_time()
    work_time = MPI_WTIME() - sync_dt_sum
    RETURN
  END FUNCTION work_time

  SUBROUTINE sync
    CHARACTER(len=*), PARAMETER :: context = 'sync: '
    INTEGER :: ierror
    REAL(dp) :: t0, dt

    t0 = MPI_WTIME()

    CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    IF (ierror /= MPI_SUCCESS) CALL xt_abort(context//'MPI_BARRIER failed', &
         __FILE__, &
         __LINE__)

    dt = (MPI_WTIME() - t0)
    sync_dt_sum = sync_dt_sum + dt

  END SUBROUTINE sync

  SUBROUTINE inflate_idx(inflate_pos, idx_2d, idx_3d)
    CHARACTER(len=*), PARAMETER :: context = 'test_perf::inflate_idx: '
    INTEGER, INTENT(in) :: inflate_pos
    INTEGER(xt_idx_kind), INTENT(in) :: idx_2d(:,:)
    INTEGER(xt_idx_kind), ALLOCATABLE, INTENT(out) :: idx_3d(:,:,:)

    INTEGER :: i, j, k

    IF (ALLOCATED(idx_3d)) DEALLOCATE(idx_3d)

    SELECT CASE(inflate_pos)
    CASE(1)
      ALLOCATE(idx_3d(ke, ie, je))
        DO j=1,je
          DO i=1,ie
            DO k=1,ke
            idx_3d(k,i,j) = k + (idx_2d(i,j)-1) * ke
          ENDDO
        ENDDO
      ENDDO
    CASE(3)
      ALLOCATE(idx_3d(ie, je, ke))
      DO k=1,ke
        DO j=1,je
          DO i=1,ie
            idx_3d(i,j,k) = idx_2d(i,j) + (k-1) * g_ie * g_je
          ENDDO
        ENDDO
      ENDDO
    CASE DEFAULT
      CALL xt_abort(context//' unsupported inflate position', &
           __FILE__, &
           __LINE__)
    END SELECT

  END SUBROUTINE inflate_idx

  SUBROUTINE gen_pos3d_surf(pos)
    INTEGER, INTENT(inout) :: pos(:,:)

    ! positions for zero based arrays ([k,i,j] dim order):
    ! old pos = i + j*ie
    ! new pos = k + (i + j*ie)*nlev

    INTEGER :: ii,jj, i,j,k, p,q

    k = 0 ! surface
    DO jj=1,je
      DO ii=1,ie
        p = pos(ii,jj) - 1 ! shift to 0-based index
        j = p/ie
        i = MOD(p,ie)
        q = k +  (i + j*ie)*nlev
        pos(ii,jj) = q + 1 ! shift to 1-based index
      ENDDO
    ENDDO

  END SUBROUTINE gen_pos3d_surf

  SUBROUTINE icmp_2d(label, f,g)
    CHARACTER(len=*), PARAMETER :: context = 'test_perf::icmp_2d: '
    CHARACTER(len=*), INTENT(in)  :: label
    INTEGER, INTENT(in)  :: f(:,:)
    INTEGER, INTENT(in)  :: g(:,:)

    INTEGER :: i, j, n1, n2

    n1 = SIZE(f,1)
    n2 = SIZE(f,2)
    IF (SIZE(g,1) /= n1 .OR. SIZE(g,2) /= n2) &
         CALL xt_abort(context//'shape mismatch error', &
         __FILE__, &
         __LINE__)

    DO j = 1, n2
      DO i = 1, n1
        IF (f(i,j) /= g(i,j)) THEN
          WRITE(0,*) context,label,' test failed: i, j, f(i,j), g(i,j) =', &
               i, j, f(i,j), g(i,j)
          CALL xt_abort(context//label//' test failed', &
               __FILE__, &
               __LINE__)
        ENDIF
      ENDDO
    ENDDO
    IF (verbose) WRITE(0,*) mype,':',context//label//' passed'
  END SUBROUTINE icmp_2d

  SUBROUTINE icmp_3d(label, f,g)
    CHARACTER(len=*), PARAMETER :: context = 'test_perf::icmp_3d: '
    CHARACTER(len=*), INTENT(in)  :: label
    INTEGER, INTENT(in)  :: f(:,:,:)
    INTEGER, INTENT(in)  :: g(:,:,:)

    INTEGER :: i1, i2, i3, n1, n2, n3

    n1 = SIZE(f,1)
    n2 = SIZE(f,2)
    n3 = SIZE(f,3)
    IF (SIZE(g,1) /= n1 .OR. SIZE(g,2) /= n2 .OR. SIZE(g,3) /= n3) &
      CALL xt_abort(context//label//'shape mismatch', __FILE__, __LINE__)

    DO i3 = 1, n3
      DO i2 = 1, n2
        DO i1 = 1, n1
          IF (f(i1,i2,i3) /= g(i1,i2,i3)) THEN
            WRITE(0,*) context,label,&
                 ' test failed: i1, i2, i3, f(i1,i2,i3), g(i1,i2,i3) =', &
                 i1, i2, i3, f(i1,i2,i3), g(i1,i2,i3)
            CALL xt_abort(context//label//' test failed', __FILE__, __LINE__)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    IF (verbose) WRITE(0,*) mype,':',context//label//' passed'
  END SUBROUTINE icmp_3d

  SUBROUTINE init_mpi
    CHARACTER(len=*), PARAMETER :: context = 'init_mpi: '
    INTEGER :: ierror
    CALL MPI_INIT(ierror)
    IF (ierror /= MPI_SUCCESS) CALL xt_abort(context//'MPI_INIT failed', &
         __FILE__, __LINE__)
  END SUBROUTINE init_mpi

  SUBROUTINE init_all
    CHARACTER(len=*), PARAMETER :: context = 'init_all: '
    INTEGER :: ierror
    CHARACTER(len=20) :: grid_str

    CALL get_environment_variable('YAXT_TEST_PERF_GRID', grid_str)

    verbose = .TRUE.

    SELECT CASE (TRIM(ADJUSTL(grid_str)))
    CASE('TOY')
      grid_kind = grid_kind_toy
      grid_label = 'TOY'
      g_ie = 66
      g_je = 36
    CASE('TP10')
      grid_kind = grid_kind_tp10
      grid_label = 'TP10'
      g_ie = 362
      g_je = 192
    CASE('TP04')
      grid_kind = grid_kind_tp04
      grid_label = 'TP04'
      g_ie = 802
      g_je = 404
    CASE('TP6M')
      grid_kind = grid_kind_tp6m
      grid_label = 'TP6M'
      g_ie = 3602
      g_je = 2394
    CASE default
      grid_kind = grid_kind_test
      grid_label = 'TEST'
      g_ie = 32
      g_je = 12
      verbose = .FALSE.
    END SELECT

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
    IF (ierror /= MPI_SUCCESS) CALL xt_abort(context//'MPI_COMM_SIZE failed', &
         __FILE__, &
         __LINE__)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierror)
    IF (ierror /= MPI_SUCCESS) CALL xt_abort(context//'MPI_COMM_RANK failed', &
         __FILE__, &
         __LINE__)
    IF (mype==0) THEN
      lroot = .true.
    ELSE
      lroot = .FALSE.
      verbose = .FALSE.
    ENDIF

    CALL factorize(nprocs, nprocx, nprocy)
    IF (lroot .AND. verbose) WRITE(0,*) 'nprocx, nprocy=',nprocx, nprocy
    IF (lroot .AND. verbose) WRITE(0,*) 'g_ie, g_je=',g_ie, g_je
    mypy = mype / nprocx
    mypx = MOD(mype, nprocx)

    CALL deco
    ke = nlev

    ALLOCATE(g_id(g_ie, g_je), g_tpex(g_ie, g_je))

    ALLOCATE(fval_2d(ie,je), gval_2d(ie,je))
    ALLOCATE(loc_id_2d(ie,je), loc_tpex_2d(ie,je))
    ALLOCATE(id_pos(ie,je), pos3d_surf(ie,je))

    fval_2d = undef_int
    gval_2d = undef_int
    loc_id_2d = undef_int
    loc_tpex_2d = undef_int
    id_pos = undef_int
    pos3d_surf = undef_int

  END SUBROUTINE init_all

  SUBROUTINE exit_mpi
    CHARACTER(len=*), PARAMETER :: context = 'exit_mpi: '
    INTEGER :: ierror
    CALL MPI_FINALIZE(ierror)
    IF (ierror /= MPI_SUCCESS) CALL xt_abort(context//'MPI_FINALIZE failed', &
         __FILE__, &
         __LINE__)
  END SUBROUTINE exit_mpi

  SUBROUTINE id_map(map)
    INTEGER, INTENT(out) :: map(:,:)

    INTEGER :: i,j,p

    p = 0
    DO j = 1, SIZE(map,2)
      DO i = 1, SIZE(map,1)
        p = p + 1
        map(i,j) = p
      ENDDO
    ENDDO

  END SUBROUTINE id_map

  SUBROUTINE exchange_int(redist, f, g)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, TARGET, INTENT(in) :: f(*)
    INTEGER, TARGET, VOLATILE, INTENT(out) :: g(*)

    CALL xt_redist_s_exchange1(redist, c_loc(f), c_loc(g));

  END SUBROUTINE exchange_int

  SUBROUTINE gen_redist(xmap, send_dt, recv_dt, redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: send_dt, recv_dt
    TYPE(xt_redist),INTENT(out) :: redist

    INTEGER :: dt

    IF (send_dt /= recv_dt) &
         CALL xt_abort('gen_redist: (send_dt /= recv_dt) unsupported', &
         __FILE__, &
         __LINE__)
    dt = send_dt
    redist = xt_redist_p2p_new(xmap, dt)

  END SUBROUTINE gen_redist

  SUBROUTINE gen_off_redist(xmap, send_dt, send_off, recv_dt, recv_off, redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER,INTENT(in) :: send_dt, recv_dt
    INTEGER,INTENT(in) :: send_off(:,:), recv_off(:,:)
    TYPE(xt_redist),INTENT(out) :: redist

    INTEGER :: dt

    IF (send_dt /= recv_dt) &
         CALL xt_abort('gen_off_redist: (send_dt /= recv_dt) unsupported', &
         __FILE__, &
         __LINE__)
    dt = send_dt

    redist = xt_redist_p2p_off_new(xmap, send_off, recv_off, dt)
  END SUBROUTINE gen_off_redist

  SUBROUTINE get_window(gval, win)
    INTEGER, INTENT(in) :: gval(:,:)
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

  SUBROUTINE gen_stripes(local_idx, local_stripes)
    CHARACTER(len=*), PARAMETER :: context = 'gen_stripes: '

    INTEGER(xt_idx_kind), INTENT(in) :: local_idx(:,:,:)
    TYPE(Xt_idxlist), INTENT(out) :: local_stripes

    TYPE(xt_stripe), ALLOCATABLE :: stripes(:,:)
    INTEGER :: i, j, k, ni, nj, nk

    nk = SIZE(local_idx,1)
    ni = SIZE(local_idx,2)
    nj = SIZE(local_idx,3)

    ALLOCATE(stripes(ni,nj))

    DO j = 1, nj
      DO i = 1, ni
        ! start, nstrides, stride
        stripes(i,j) = xt_stripe(local_idx(1,i,j), nk, 1)
        DO k = 1, nk
          IF (local_idx(1,i,j)-1+k /= local_idx(k,i,j)) &
               CALL xt_abort(context//'stripe condition violated', &
               __FILE__, &
               __LINE__)
        ENDDO
      ENDDO
    ENDDO

    local_stripes = xt_idxstripes_new(stripes, SIZE(stripes))

  END SUBROUTINE gen_stripes

  SUBROUTINE gen_xmap_2d(local_src_idx, local_dst_idx, xmap)
    INTEGER(xt_idx_kind), INTENT(in) :: local_src_idx(:,:)
    INTEGER(xt_idx_kind), INTENT(in) :: local_dst_idx(:,:)
    TYPE(xt_xmap), INTENT(out) :: xmap

    TYPE(Xt_idxlist) :: src_idxlist, dst_idxlist

    src_idxlist = xt_idxvec_new(local_src_idx, SIZE(local_src_idx))
    dst_idxlist = xt_idxvec_new(local_dst_idx, SIZE(local_dst_idx))
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist,  MPI_COMM_WORLD)

    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

  END SUBROUTINE gen_xmap_2d

  SUBROUTINE gen_xmap_3d(local_src_idx, local_dst_idx, xmap)
    INTEGER(xt_idx_kind), INTENT(in) :: local_src_idx(:,:,:)
    INTEGER(xt_idx_kind), INTENT(in) :: local_dst_idx(:,:,:)
    TYPE(xt_xmap), INTENT(out) :: xmap

    TYPE(Xt_idxlist) :: src_idxlist, dst_idxlist

    src_idxlist = xt_idxvec_new(local_src_idx, SIZE(local_src_idx))
    dst_idxlist = xt_idxvec_new(local_dst_idx, SIZE(local_dst_idx))
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist,  MPI_COMM_WORLD)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)


  END SUBROUTINE gen_xmap_3d

  SUBROUTINE def_exchange()

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
    g_tpex = undef_index
    g_tpex(g_core_is:g_core_ie, g_core_js:g_core_je) &
         = g_id(g_core_is:g_core_ie, g_core_js:g_core_je)

    IF (with_north_halo) THEN

      ! north inversion, (maybe with increased north halo)
      IF (increased_north_halo) THEN
        north_halo = nhalo+1
      ELSE
        north_halo = nhalo
      ENDIF

      IF (2*north_halo > g_core_je) &
           CALL xt_abort('def_exchange: grid too small (or halo too large&
           &) for tripolar north exchange', &
           __FILE__, &
           __LINE__)
      DO j = 1, north_halo
        DO i = g_core_is, g_core_ie
          g_tpex(i,j) = g_tpex(g_core_ie + (g_core_is-i), 2*north_halo + (1-j))
        ENDDO
      ENDDO

    ELSE

      DO j = 1, nhalo
        DO i = nhalo+1, g_ie-nhalo
          g_tpex(i,j) = g_id(i,j)
        ENDDO
      ENDDO

    ENDIF

    ! south: no change
    DO j = g_core_je+1, g_je
      DO i = nhalo+1, g_ie-nhalo
        g_tpex(i,j) = g_id(i,j)
      ENDDO
    ENDDO

    ! PBC
    DO j = 1, g_je
      DO i = 1, nhalo
        g_tpex(g_core_is-i,j) = g_tpex(g_core_ie+(1-i),j)
      ENDDO
      DO i = 1, nhalo
        g_tpex(g_core_ie+i,j) = g_tpex(nhalo+i,j)
      ENDDO
    ENDDO

    CALL check_g_idx (g_tpex)

  END SUBROUTINE def_exchange

  SUBROUTINE check_g_idx(gidx)
    INTEGER,INTENT(in) :: gidx(:,:)
    INTEGER :: i,j

    DO j=1,g_je
      DO i=1,g_ie
        IF (gidx(i,j) == undef_index) THEN
          CALL xt_abort('check_g_idx: check failed', __FILE__, __LINE__)
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE check_g_idx

  SUBROUTINE deco
    INTEGER :: cx0(0:nprocx-1), cxn(0:nprocx-1)
    INTEGER :: cy0(0:nprocy-1), cyn(0:nprocy-1)

    cx0 = 0
    cxn = 0
    CALL regular_deco(g_ie-2*nhalo, cx0, cxn)

    cy0 = 0
    cyn = 0
    CALL regular_deco(g_je-2*nhalo, cy0, cyn)

    ! process local deco variables:
    ie = cxn(mypx) + 2*nhalo
    je = cyn(mypy) + 2*nhalo
    p_ioff = cx0(mypx)
    p_joff = cy0(mypy)

  END SUBROUTINE deco

  SUBROUTINE regular_deco(g_cn, c0, cn)
    INTEGER, INTENT(in) :: g_cn
    INTEGER, INTENT(inout) :: c0(0:), cn(0:)

    ! convention: process space coords start at 0, grid point coords start at 1

    INTEGER :: tn
    INTEGER :: d, m
    INTEGER :: it

    tn = SIZE(c0)
    IF (tn<0) CALL xt_abort('(tn<0)', __FILE__, __LINE__)
    IF (tn>g_cn) &
         CALL xt_abort('regular_deco: too many task for such a core region', &
         __FILE__, &
         __LINE__)

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
    IF (c0(tn-1)+cn(tn-1) /= g_cn) &
         CALL xt_abort('regular_deco: internal error 1', &
         __FILE__, &
         __LINE__)

  END SUBROUTINE regular_deco

  SUBROUTINE factorize(c, a, b)
    INTEGER, INTENT(in) :: c
    INTEGER, INTENT(out) :: a, b ! c = a*b

    INTEGER :: x0, i

    IF (c<1) CALL xt_abort('factorize: invalid process space', &
         __FILE__, &
         __LINE__)
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

END PROGRAM test_perf_stripes
