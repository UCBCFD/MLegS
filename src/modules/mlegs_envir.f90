module mlegs_envir 
  !> module for global environment variables
  use MPI
  implicit none
  private
  !> basic math
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

  !> MPI
  !> global communicator group
  integer, public :: comm_glb, rank_glb, nprocs_glb, mpi_thread_mode
  !> cartesian communicator groups (MPI_GROUP)
  integer, public :: comm_grps(2), rank_grps(2), nprocs_grps(2)
  !> scalar element data type for MPI communications
  integer, public:: scalar_element_type = MPI_DOUBLE_COMPLEX
  integer, public:: cp8_size = 0
  !> error flags
  integer, public :: mpi_ierr
  integer, public, parameter :: &
  error_flag_comm = 201, & ! error related to MPI communicator
  error_flag_alloc = 202, & ! error related to memory allocation
  error_flag_dealloc = 203, & ! error related to memory deallocation
  error_flag_file = 204, & ! error related to file I/O
  error_flag_misc = 205, & ! miscellaneous error
  error_flag_warning = 206 ! warning

  !> I/O
  logical, public, parameter :: is_warning = .true.

end module