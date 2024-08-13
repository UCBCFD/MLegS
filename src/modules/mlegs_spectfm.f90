module mlegs_spectfm
  use mlegs_base; use mlegs_mpi; use mlegs_genmat; use mlegs_bndmat
  implicit none
  private

  !> transformation kit counter (to only allow one tfm_kit per program)
  integer(i4), public :: tfm_kit_counter = 0

  !> base transformation kit (r; m=0) for f(r)==sum_(n=0~N) a_n*P^0_L_n(r)
  type, public :: tfm_kit_base
    !> check if intialized (by default false)
    logical :: is_set = .false.
    !> gauss-legendre abscissa before mapping (-1 ~ 1)
    real(p8), dimension(:), allocatable :: x
    !> gauss-legendre weights
    real(p8), dimension(:), allocatable :: w
    !> mapped collocation points in r (0 ~ inf)
    real(p8), dimension(:), allocatable :: r
    !> -ln(1-X_I) (I = 1, 2,...)
    real(p8), dimension(:), allocatable :: ln
    !> number of MLeg basis elem. used |M| ~ nrchop (M = 0, 1, ...)
    integer(i4), dimension(:), allocatable :: chops
    !> log(norm.) factors of MLeg basis elem. (N, M)
    real(p8), dimension(:,:), allocatable :: lognorm
    !> MLeg basis vals for P^M_N(X_I) (I, N, M)
    real(p8), dimension(:,:,:), allocatable :: pf
    !> MLeg values at origin & infinity when M = 0 (N = 0, 1, ...)
    real(p8), dimension(:), allocatable :: at0, at1
    !> spectral workspace dimensions
    integer(i4) :: nrdim
    contains
      procedure :: set => tfm_kit_set
      procedure :: dealloc => tfm_kit_dealloc
  end type

  !> 2d spectral transformation kit (r, p; ak=0) for f(r,p)==[sum_(n=|m|~|m|+N)a_n*P^m_L_n(r)]e^(i*m*p)
  type, extends(tfm_kit_base), public :: tfm_kit_2d
    !> equispaced collocation points in p (0 ~ 2pi)
    real(p8), dimension(:), allocatable :: p
    !> accessible azimuthal wavenumbers 0 ~ npchop - 1
    integer(i4), dimension(:), allocatable :: m
    !> spectral workspace dimensions
    integer(i4) :: npdim
  end type

  !> 3d spectral transformation kit (r, p, z) for f(r,p,z)==[sum_(n=|m|~|m|+N)a_n*P^m_L_n(r)]e^(i*(m*p+k*z))
  type, extends(tfm_kit_2d), public :: tfm_kit_3d
    !> equispaced collocation points in z (0 ~ zlen)
    real(p8), dimension(:), allocatable :: z
    !> accessible axial wavenumbers 0 ~ nzchop - 1
    integer(i4), dimension(:), allocatable :: ak
    !> spectral workspace dimensions
    integer(i4) :: nzdim, zchopl, zchopu
  end type

  !> initialization
  interface
    module subroutine tfm_kit_set(this)
      implicit none
      class(tfm_kit_base), intent(inout) :: this
    end subroutine
    module subroutine tfm_kit_dealloc(this)
      implicit none
      class(tfm_kit_base), intent(inout) :: this
    end subroutine
  end interface

  !> spectral differential operators -- (1-x**2)d/dx (r*d/dr) 
  interface leg_xxdx
    module function leg_xxdx_with_tfm_kit(mval, n_input, n_output, tfm_kit) result(xxdx)
      implicit none
      integer(i4), intent(in) :: mval, n_input, n_output
      class(tfm_kit_base), intent(in) :: tfm_kit
      type(real_bndm) :: xxdx
    end function
  end interface
  public :: leg_xxdx

  !> spectral differential operators -- Del^2_H (Horizontal laplacian) (1/r*d/dr(r*d/dr) - m^2/r^2)
  interface leg_del2h
    module function leg_del2h_with_tfm_kit(mval, n_input, n_output, tfm_kit) result(del2h)
      implicit none
      integer(i4), intent(in) :: mval, n_input, n_output
      class(tfm_kit_base), intent(in) :: tfm_kit
      type(real_bndm) :: del2h
    end function
  end interface
  public :: leg_del2h

  !> spectral differential operators -- Del^2 (Laplacian)  (1/r*d/dr(r*d/dr) - m^2/r^2 - ak^2)
  interface leg_del2
    module function leg_del2_with_tfm_kit(mval, akval, n_input, n_output, tfm_kit) result(del2)
      implicit none
      integer(i4), intent(in) :: mval, n_input, n_output
      real(p8), intent(in) :: akval
      class(tfm_kit_base), intent(in) :: tfm_kit
      type(real_bndm) :: del2
    end function
  end interface
  public :: leg_del2

end module