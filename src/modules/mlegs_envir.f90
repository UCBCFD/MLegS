module mlegs_envir
  implicit none
  private

  !> environment variables
  
  !> 1. basic math
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

  !> 2. MPI
  !> global communicator group
  integer, public :: comm_glb, rank_glb, nprocs_glb
  !> cartesian communicator groups (MPI_GROUP)
  integer, public :: comm_grps(2), rank_grps(2), nprocs_grps(2)
  !> scalar element data type for MPI communications
  integer, public:: scalar_element_type = MPI_DOUBLE_COMPLEX
  integer, public:: cp8_size = 0
  !> error flags
  integer, public :: mpi_ierr
  integer, public, parameter :: error_flag_comm = 201

  !> 3. I/O
  logical, public, parameter :: is_warning = .true.

end module