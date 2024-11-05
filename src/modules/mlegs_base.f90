module mlegs_base 
  !> module for base simulation setup parameters
  use mlegs_envir
  implicit none
  private
  !> discretization -- number of colloc. pts. / spec. elems. in the radial direction (must be even)
  integer(i4), public :: nr = 256, nrchop = 256
  !> discretization -- number of colloc. pts. / spec. elems. in the azimuthal direction (2^p 3^q 5^r) (rec. npchop = np/2+1)
  integer(i4), public :: np = 128, npchop = 65
  !> discretization -- number of colloc. pts. / spec. elems. in the axial direction (if 3D) (2^p 3^q 5^r) (rec. nzchop = nz/2+1)
  integer(i4), public :: nz = 128, nzchop = 65
  !> computational domain -- Legendre map parameter (0 < r < ell: high-resolution region)
  real(p8), public :: ell = 4.D0
  !> computational domain -- axial length (== longest axial wavelength size to be considered)
  real(p8), public :: zlen = 2.D0*pi
  !> time stepping -- time step size
  real(p8), public :: dt = 1.D-3
  !> time stepping -- start time
  real(p8), public :: ti = 0.D0
  !> time stepping -- end time (termination criterion)
  real(p8), public :: totaltime = 1.D2
  !> time stepping -- start step index
  integer(i4), public :: ni = 0
  !> time stepping -- end step index (termination criterion)
  integer(i4), public :: totaln = 100000
  !> flow properties -- viscosity (== inverse of Reynolds number)
  real(p8), public :: visc = 1.D-5
  !> flow properties -- hyperviscosity power (0 if no hyperviscosity, otherwise greater than 2 and even)
  integer(i4), public :: hyperpow = 0
  !> flow properties -- hypervisocisty value
  real(p8), public :: hypervisc = 0.D0
  !> file I/O -- field output directory
  character(len=256), public :: flddir = './output/fld/'
  !> file I/O -- enabling field output (T/F)
  logical, public :: isfldsav = .true.
  !> file I/O -- field output interval
  integer(i4), public :: fldsavintvl = 100
  !> file I/O -- log output directory
  character(len=256), public :: logdir = './output/log/'
  !> file I/O -- enabling log output (T/F)
  logical, public :: islogsav = .true.
  !> file I/O -- log output interval
  integer(i4), public :: logsavintvl = 1

end module