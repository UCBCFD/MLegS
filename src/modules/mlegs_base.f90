module mlegs_base 
  !> module for base simulation setup parameters
  use mlegs_envir
  implicit none
  private
  !> discretization -- number of colloc. pts. / spec. elems. in the radial direction (2*p)
  integer(i4), public :: nr = 128, nrchop = 128
  !> discretization -- number of colloc. pts. / spec. elems. in the azimuthal direction (2^p 3^q 5^r)
  integer(i4), public :: np = 64, npchop = 33
  !> discretization -- number of colloc. pts. / spec. elems. in the axial direction (if 3D) (2^p 3^q 5^r)
  integer(i4), public :: nz = 64, nzchop = 33
  !> computational domain -- Legendre map parameter (0 < r < ell: high-resolution region)
  real(p8), public :: ell = 4.D0
  !> computational domain -- axial length (== longest axial wavelength size to be considered)
  real(p8), public :: zlen = 2*pi
  !> time stepping -- time step size
  real(p8), public :: dt = 1.D-2
  !> time stepping -- start time
  real(p8), public :: ti = 0.D0
  !> time stepping -- end time (termination criterion)
  real(p8), public :: totaltime = 1.D2
  !> time stepping -- start step index
  integer(i4), public :: ni = 0
  !> time stepping -- end step index (termination criterion)
  integer(i4), public :: totaln = 1.D4
  !> flow properties -- viscosity (== inverse of Reynolds number)
  real(p8), public :: visc = 1.D-5
  !> flow properties -- hyperviscosity power (0 if no hyperviscosity, otherwise greater than 2 and even)
  integer(i4), public :: hyperpow = 0
  !> flow properties -- hypervisocisty value
  real(p8), public :: hypervisc = 0.D0
  !> file I/O -- field output directory
  character(len=1024), public :: flddir = './output/fld/'
  !> file I/O -- log output directory
  character(len=1024), public :: logdir = './output/log/'

end module