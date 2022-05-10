module const_mod

! Define some constants that can be used anywhere
use, intrinsic :: iso_fortran_env, only : real32, real64, int32, int64
integer(int32) :: total_procs, local_proc
integer(int32), parameter :: max_factor_iters = 10000

end module const_mod

!-----------------------------------------------------------------------

module data_mod
  use mpi
  use yaxt, ONLY: &
       xt_idxlist, xt_idxvec_new, xt_idxlist_delete, &
       xt_xmap, xt_xmap_dist_dir_new, xt_xmap_delete, &
       xt_redist, xt_redist_p2p_off_new, xt_xmap_all2all_new, &
       xt_redist_delete, xt_redist_s_exchange
  use const_mod
  implicit none

! Type to hold data and perform a halo exchange
  type, public :: data_type
    private
    type(xt_redist) :: redist
    real(real64), allocatable :: fort_array(:)
  contains
    procedure, public :: set_values
    procedure, public :: halo_exchange
    procedure, public :: print_values
    procedure, public :: clean
  end type

  interface data_type
    module procedure data_constructor
  end interface

contains

! Constructor for the data type. Setup up everything ready to halo exchange
function data_constructor(partition, halo, count_no) result(self)
  implicit none
  integer(int64), intent(in)  :: partition( : ), halo( : )
  integer(int64), intent(out) :: count_no(7)
  type(data_type) :: self

  integer(int64), allocatable :: src_indices(:), tgt_indices(:)
  integer(int32) :: i
  type(xt_idxlist) :: src_idxlist, tgt_idxlist
  type(xt_xmap) :: xmap
  integer(int32), allocatable :: src_offsets(:)
  integer(int32), allocatable :: tgt_offsets(:)

  call system_clock(count_no(1))
  allocate(self%fort_array(size(partition)+size(halo)))
  
  allocate(src_indices(size(partition)+size(halo)))
  src_indices = -1
  do i = 1, size(partition)
    src_indices(i) = partition(i)
  end do
  allocate(tgt_indices(size(halo)))
  do i = 1, size(halo)
    tgt_indices(i) = halo(i)
  end do

  call system_clock(count_no(2))
  ! create decomposition descriptors
  src_idxlist = xt_idxvec_new(src_indices, size(src_indices))
  tgt_idxlist = xt_idxvec_new(tgt_indices, size(tgt_indices))

  call system_clock(count_no(3))
  ! generate exchange map
  xmap = xt_xmap_dist_dir_new(src_idxlist, tgt_idxlist, mpi_comm_world)
  !xmap = xt_xmap_all2all_new(src_idxlist, tgt_idxlist, mpi_comm_world)

  call system_clock(count_no(4))
  allocate(src_offsets(size(partition)))
  allocate(tgt_offsets(size(halo)))
  do i = 0, size(partition) - 1
    src_offsets(i+1) = i
  end do
  do i = 0, size(halo),1
    tgt_offsets(i+1) = i + size(partition)
  end do

  call system_clock(count_no(5))
  self%redist = xt_redist_p2p_off_new(xmap, src_offsets, tgt_offsets, &
                                      mpi_double_precision)

  call system_clock(count_no(6))
  call xt_xmap_delete(xmap)
  call xt_idxlist_delete(tgt_idxlist)
  call xt_idxlist_delete(src_idxlist)
  deallocate(src_indices)
  deallocate(src_offsets)
  deallocate(tgt_offsets)
  call system_clock(count_no(7))
end function data_constructor

subroutine set_values(self, data_values)
  implicit none
  class( data_type ) :: self
  integer(int64), intent(in) :: data_values(:)
  integer i
  ! Fill the data points
  do i=1,size(data_values)
    self%fort_array(i) = data_values(i)*10.0
  end do
end subroutine set_values

subroutine halo_exchange(self)
  implicit none
  class(data_type) :: self
  call xt_redist_s_exchange(self%redist, self%fort_array, self%fort_array)
end subroutine halo_exchange

subroutine print_values(self,title)
  implicit none
  class(data_type) :: self
  character (len=*), intent(in) :: title

  write(*,"(i2,a,a)")local_proc," : ",title
  write(*,"(i2,a,3e15.3,'...',3e15.3)")local_proc," : ",self%fort_array(1:3),self%fort_array(size(self%fort_array)-2:)

end subroutine print_values

subroutine clean(self)
  implicit none
  class( data_type ) :: self
  ! clean up
  call xt_redist_delete(self%redist)
end subroutine clean

end module

!-----------------------------------------------------------------------

program YAXT_BiP_halo_exchange

! Example of how to perform a halo exchange on a BiPeriodic planar mesh
! (which is treated as if it is unstructured), so the data is held as a
! 1-d list of points owned by the local partition, followed by the halo points
! Can be run on virtually any number of MPI tasks and will partition accordingly

use mpi
use yaxt, only: xt_initialize, xt_finalize
use data_mod
use const_mod
implicit none

  integer(int32) :: x_size, y_size, nlayers, ndata, halo_depth
  namelist/YAXT_BiP/ x_size, y_size, nlayers, ndata, halo_depth

  integer(int32) :: column_sz
  integer(int32) :: start_factor
  integer(int32) :: end_factor
  integer(int32) :: fact_count
  logical(int32) :: found_factors
  integer(int32) :: x_proc, y_proc  !total no of proc in x- and y-
  integer(int32) :: i_proc, j_proc  !local process in x- and y-
  integer(int32) :: start_x, start_y, num_x, num_y
  integer(int32) :: left_halo_x, lower_halo_y, right_halo_x, upper_halo_y !lower-left and upper-right halo coordinates
  integer(int32) :: i, j, k, dep
  integer(int32) :: halo_index
  integer(int64) :: count_rate, count_no(5), data_count_no(7)
  real(real64)   :: clock_rate
  real(real64)   :: time_local

  type(data_type)             :: my_data
  integer(int64), allocatable :: partition( : ), halo( : )
  integer(int32)              :: num_local, num_halo
  integer(int32)              :: ierror

  ! Read namelist input for dimensions of domain
  open(unit=50,file="YAXT_BiP.nml",status="old")
  read(50,nml=YAXT_BiP)
  close(50)
  column_sz = nlayers * ndata

  ! Initialise MPI and YAXT
  call mpi_init(ierror)
  if (ierror /= mpi_success) stop "ERROR: MPI failed to initialise"
  call xt_initialize(mpi_comm_world)
  call mpi_comm_rank(mpi_comm_world, local_proc, ierror)
  if (ierror /= mpi_success) stop "ERROR: could not retrieve current rank"
  call mpi_comm_size (mpi_comm_world, total_procs, ierror)
  if (ierror /= mpi_success) stop "ERROR: could not retrieve num of ranks"

  ! Initialise timer
  call system_clock(count_rate=count_rate)
  clock_rate = real(count_rate, real64)
  call system_clock(count_no(1))

  ! Calculate the most even decomposition in x- and y-direction given total_procs
  start_factor  = nint(sqrt(real(total_procs, kind=real64)), kind=int32)
  end_factor    = max(1,(start_factor-max_factor_iters))
  found_factors = .false.
  do fact_count = start_factor, end_factor, -1
    if (mod(total_procs,fact_count) == 0) then
      found_factors = .true.
      exit
    end if
  end do
  if (found_factors) then
    x_proc = fact_count
    y_proc = total_procs/fact_count
  else
    stop "ERROR: could not automatically partition domain."
  end if

  ! Work out the start index and number of cells (in x- and y-dirn) for
  ! the local partition - this algorithm should spread out the number of
  ! cells each partition gets fairly evenly
  j_proc = ( local_proc / x_proc )
  i_proc = local_proc - ( j_proc * x_proc )
  start_x = ( (i_proc*x_size) / x_proc ) + 1
  num_x   = ( ((i_proc+1)*x_size) / x_proc ) - start_x + 1
  start_y = ( (j_proc*y_size) / y_proc ) + 1
  num_y   = ( ((j_proc+1)*y_size) / y_proc ) - start_y + 1

  ! Store the global indices of the owned cells in the partition
  allocate(partition(num_x*num_y*column_sz))
  do i = 1,num_x
    do j = 1, num_y
      do k = 1, column_sz
        partition( xy2lid(i,j,k) ) = xy2gid(i+start_x-1,j+start_y-1,k)
      end do
    end do
  end do

  ! Store the global indices of the halo cells in increasing halo depths
  allocate(halo( (4*(halo_depth**2) + (2*num_x+2*num_y)*halo_depth)*column_sz ))
  halo_index = 1
  left_halo_x  = start_x
  lower_halo_y = start_y
  right_halo_x = start_x + num_x - 1
  upper_halo_y = start_y + num_x - 1
  do dep = 1, halo_depth
    left_halo_x  = left_halo_x  - 1
    lower_halo_y = lower_halo_y - 1
    right_halo_x = right_halo_x + 1
    upper_halo_y = upper_halo_y + 1
    do k = 1, column_sz
      do i = 1, ((dep-1)*2) + 1 + num_x
        halo(halo_index)=xy2gid(left_halo_x+i-1,lower_halo_y,k)
        halo_index = halo_index + 1
      end do
      do j = 1, ((dep-1)*2) + 1 + num_y
        halo(halo_index)=xy2gid(right_halo_x,lower_halo_y+j-1,k)
        halo_index = halo_index + 1
      end do
      do i = 1, ((dep-1)*2) + 1 + num_x
        halo(halo_index)=xy2gid(right_halo_x-i+1,upper_halo_y,k)
        halo_index = halo_index + 1
      end do
      do j = 1, ((dep-1)*2) + 1 + num_y
        halo(halo_index)=xy2gid(left_halo_x,upper_halo_y-j+1,k)
        halo_index = halo_index + 1
      end do
    end do
  end do
  call system_clock(count_no(2))

  ! Create the data object (call its constructor)
  my_data=data_type(partition, halo, data_count_no)
  call system_clock(count_no(3))

  ! Fill the data values owned by this cell with data
  call my_data%set_values(partition)

  ! Print out what we have before the halo exchange
  !call my_data%print_values("Data before halo swap...")

  call system_clock(count_no(4))
  ! Do the halo exchange
  call my_data%halo_exchange()
  call system_clock(count_no(5))

  ! Print out what we have after the halo exchange
  !call my_data%print_values("Data after halo swap...")

  call my_data%clean()

  if ( local_proc == 0 ) then
    write(*,'("Overall:")')
    write(*,'("                  Mean      Max       Min")')
  end if
  time_local = (count_no(5)-count_no(1))/clock_rate
  call output_times( time_local, "          All: " )
  time_local = (count_no(2)-count_no(1))/clock_rate
  call output_times( time_local, " Partitioning: " )
  time_local = (count_no(3)-count_no(2))/clock_rate
  call output_times( time_local, "Setup Routing: " )
  time_local = (count_no(5)-count_no(4))/clock_rate
  call output_times( time_local, "Halo Exchange: " )

  if ( local_proc == 0 ) then
    write(*,'("---------------------")')
    write(*,'("Routing Construction:")')
    write(*,'("                  Mean      Max       Min")')
  end if
  time_local = (data_count_no(3)-data_count_no(1))/clock_rate
  call output_times( time_local, "  Set Indices: " )
  time_local = (data_count_no(4)-data_count_no(3))/clock_rate
  call output_times( time_local, "   Create Map: " )
  time_local = (data_count_no(5)-data_count_no(4))/clock_rate
  call output_times( time_local, "Create Offset: " )
  time_local = (data_count_no(6)-data_count_no(5))/clock_rate
  call output_times( time_local, "Create Redist: " )

  ! finalise
  call xt_finalize()
  call mpi_finalize(ierror)

contains

function xy2gid( gx, gy ,gz) result( gid )
  implicit none
  integer(int32), intent(in) :: gx   !global x
  integer(int32), intent(in) :: gy   !global y
  integer(int32), intent(in) :: gz   !global z
  integer(int64)             :: gid  !global index

  integer(int32) :: x_dom
  integer(int32) :: y_dom

  x_dom = gx
  if ( x_dom < 1 ) x_dom = x_size + gx
  if ( x_dom > x_size ) x_dom = gx - x_size
  y_dom = gy
  if ( y_dom < 1 ) y_dom = y_size + gy
  if ( y_dom > y_size ) y_dom = gy - y_size

  gid = (y_dom-1)*x_size*column_sz + (x_dom-1)*column_sz + gz
end function xy2gid

function xy2lid( lx, ly ,lz) result( lid )
  implicit none
  integer(int32), intent(in) :: lx   !local x
  integer(int32), intent(in) :: ly   !local y
  integer(int32), intent(in) :: lz   !local z
  integer(int64)             :: lid  !local index

  lid = (ly-1)*num_x*column_sz + (lx-1)*column_sz + lz
end function xy2lid

subroutine output_times( time_local, string )
  implicit none
  real(real64), intent(in)  :: time_local
  character(*), intent(in)  :: string
  real(real64)              :: time_sum, time_mean, time_max, time_min
  integer(int32)            :: ierror
  call mpi_allreduce( time_local, time_sum, 1, MPI_DOUBLE_PRECISION, &
                      mpi_sum, mpi_comm_world, ierror )
  time_mean = time_sum/total_procs
  call mpi_allreduce( time_local, time_max, 1, MPI_DOUBLE_PRECISION, &
                      mpi_max, mpi_comm_world, ierror )
  call mpi_allreduce( time_local, time_min, 1, MPI_DOUBLE_PRECISION, &
                      mpi_min, mpi_comm_world, ierror )
  if ( local_proc == 0 ) &
    write(*,'(a,3f10.5)') string, time_mean, time_max, time_min
  end subroutine output_times
end program YAXT_BiP_halo_exchange
