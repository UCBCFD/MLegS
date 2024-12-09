module mlegs_envir 
  !> module for global environment variables
  use MPI
  implicit none
  private
  
  !> Basic math
  !> 32-bit (single-precision) real/complex
  integer, public, parameter :: p4 = selected_real_kind(6, 37) 
  !> 64-bit (double-precision) real/complex
  integer, public, parameter :: p8 = selected_real_kind(15, 307)
  !> 32-bit (default) integer
  integer, public, parameter :: i4 = selected_int_kind(9)
  !> 64-bit (long) integer
  integer, public, parameter :: i8 = selected_int_kind(18)
  !> pi == 3.1415926535 ...
  real(p8), public, parameter :: pi = acos(-1.D0)
  !> pure imaginary unit (1i)
  complex(p8), public, parameter :: iu = (0.D0, 1.D0)
  !> string size (per number to be stored) for formatted file I/O
  integer(i4), public, parameter :: formatted_num_str_len = 24

  !> Timestep
  !> Current timestep
  integer(i4), public :: curr_n = 0
  !> Current simulation time
  real(p8), public :: curr_t = 0.D0

  !> MPI
  !> global communicator group
  integer, public :: comm_glb, rank_glb, nprocs_glb
  !> cartesian communicator groups (MPI_GROUP)
  integer, public :: comm_grps(2), rank_grps(2), nprocs_grps(2)
  !> scalar element data type for MPI communications
  integer, public:: scalar_element_type = MPI_Double_complex
  integer, public:: cp8_size = 0
  !> error flags
  integer, public :: MPI_err
  integer, public, parameter :: &
  error_flag_comm = 201, &    ! error related to MPI communicator
  error_flag_alloc = 202, &   ! error related to memory allocation
  error_flag_dealloc = 203, & ! error related to memory deallocation
  error_flag_file = 204, &    ! error related to file I/O
  error_flag_misc = 205, &    ! miscellaneous error
  error_flag_warning = 206    ! warning
  logical, public, parameter :: is_warning = .false.
  !> 1D MPI decompose
  interface
    module subroutine decompose(nsize,nprocs,proc_num,nsize_proc,index_proc)
      implicit none
      integer :: nsize, nprocs, proc_num, nsize_proc, index_proc
    end subroutine
  end interface
  public :: decompose
  !> Local data size after 1D decompose
  interface
    module function local_size(nsize,comm) result(l_size)
      implicit none
      integer :: nsize, comm
      integer :: l_size
    end function
  end interface
  public :: local_size
  !> Local index after 1D decompose
  interface
    module function local_index(nsize,comm) result(l_index)
      implicit none
      integer :: nsize, comm
      integer :: l_index
    end function
  end interface
  public :: local_index
  !> Local MPI rank
  interface
    module function local_proc(comm) result(l_proc)
      implicit none
      integer :: comm
      integer :: l_proc
    end function
  end interface
  public :: local_proc
  !> Number of MPI procs
  interface
    module function count_proc(comm) result(nprocs)
      implicit none
      integer :: comm
      integer :: nprocs
    end function
  end interface
  public :: count_proc

end module