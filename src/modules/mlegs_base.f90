module mlegs_base
  implicit none
  private

  !> parameters
  !> ----------
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

  !> global variables
  !> ----------------
  !> discretization -- number of colloc. pts. / spec. elems. in the radial direction (2*p)
  integer(i4), public :: nr, nrchop
  !> discretization -- number of colloc. pts. / spec. elems. in the azimuthal direction (2^p 3^q 5^r)
  integer(i4), public :: np, npchop
  !> discretization -- number of colloc. pts. / spec. elems. in the axial direction (if 3D) (2^p 3^q 5^r)
  integer(i4), public :: nz, nzchop
  !> computational domain -- Legendre map parameter (0 < r < ell: high-resolution region)
  real(p8), public :: ell
  !> computational domain -- axial length (== longest axial wavelength size to be considered)
  real(p8), public :: zlen
  !> time stepping -- time step size
  real(p8), public :: dt
  !> time stepping -- start time
  real(p8), public :: ti
  !> time stepping -- end time (termination criterion)
  real(p8), public :: totaltime
  !> time stepping -- start step index
  integer(i4), public :: ni
  !> time stepping -- end step index (termination criterion)
  integer(i4), public :: totaln
  !> flow properties -- viscosity (== inverse of Reynolds number)
  real(p8), public :: visc
  !> flow properties -- hyperviscosity power (0 if no hyperviscosity, otherwise greater than 2 and even)
  integer(i4), public :: hyperpow
  !> flow properties -- hypervisocisty value
  real(p8), public :: hypervisc
  !> file I/O -- field output directory
  character(len=1024), public :: flddir
  !> file I/O -- log output directory
  character(len=1024), public :: logdir

end module