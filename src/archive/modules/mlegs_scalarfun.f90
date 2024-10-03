module mlegs_scalarfun
  use mlegs_base; use mlegs_genmat; use mlegs_bndmat
  use mlegs_mpi; use mlegs_spectfm
  implicit none
  private

  !> bind all scalar types with the class name 'scalar'
  type, public :: scalar
    !> logarithmic term coeff. in order to describe O(ln(r)) as r -> infty. (P_ln == ln(r^2+L^2/2L^2))
    real(p8) :: ln
    contains
      procedure :: set => scalar_set
      procedure :: chopset => scalar_chopset
      procedure :: chopdo => scalar_chopdo
      procedure :: dealloc => scalar_dealloc
  end type

  !> 1d scalar (always serial; single-processor execution is recommended)
  type, extends(scalar), public :: scalar_1d 
    !> 1d scalar function value e(r) (Physical) <--> e(n) (Spectral)
    complex(p8), dimension(:), allocatable :: e
    !> space where the scalar is currently in - (radial) ('P' or 'F')
    character(len=1) :: space
    !> local radial chopping info. (# or radial spectral elements per m)
    integer(i4), dimension(:), allocatable :: local_chops
    contains
      procedure :: tfm => scalar_1d_tfm
  end type

  !> 2d scalar (slab decomposition for parallelism)
  type, extends(scalar), public :: scalar_2d
    !> 2d scalar function value e(p,r) (Physical) <--> e(m,n) (Spectral)
    complex(p8), dimension(:,:), allocatable :: e
    !> space where the scalar is currently in - (azimuthal, radial) ('P' or 'F')
    character(len=2) :: space
    !> local radial chopping info. (# or radial spectral elements per m)
    integer(i4), dimension(:), allocatable :: local_chops
    !> local azimuthal chopping info.
    integer(i4) :: local_chopp
    !> direction of a slab (i.e., the dimension completely residing in local memory; 1: p, 2: r)
    integer(i4) :: slab_dir
    contains
      procedure :: tfm => scalar_2d_tfm
  end type

  !> 3d scalar (pencil decomposition for parallelism)
  type, extends(scalar), public :: scalar_3d
    !> 2d scalar function value e(p,r,z) (Physical) <--> e(m,n,k) (Spectral)
    complex(p8), dimension(:,:,:), allocatable :: e
    !> space where the scalar is currently in - (azimuthal, radial, axial) ('P' or 'F')
    character(len=3) :: space
    !> local radial chopping info. (# or radial spectral elements per m)
    integer(i4), dimension(:), allocatable :: local_chops
    !> local azimuthal/axial chopping info.
    integer(i4) :: local_chopp, local_chopzl, local_chopzu
    !> direction of a pencil (i.e., the dimension completely residing in local memory; 1: p, 2: r, 3: z)
    integer(i4) :: pencil_dir
    contains
      procedure :: tfm => scalar_3d_tfm
  end type

  !> common initialization/deallocation procedures
  interface
    module subroutine scalar_set(this, kit)
      implicit none
      class(scalar), intent(inout) :: this
      class(tfm_kit), intent(in) :: kit
    end subroutine
    module subroutine scalar_chopset(this, iof, iofp, iofz)
      implicit none
      class(scalar), intent(inout) :: this
      integer(i4), intent(in) :: iof
      integer(i4), optional :: iofp, iofz
    end subroutine
    module subroutine scalar_chopdo(this)
      implicit none
      class(scalar), intent(inout) :: this
    end subroutine
    module subroutine scalar_dealloc(this)
      implicit none
      class(scalar), intent(inout) :: this
    end subroutine
  end interface  

  !> 1-d/2-d/3-d spectral transformations
  interface
    module subroutine scalar_1d_tfm(this, sp)
      class(scalar_1d), intent(inout) :: this
      character(len=1), intent(in) :: sp
    end subroutine
    module subroutine scalar_2d_tfm(this, sp)
      class(scalar_2d), intent(inout) :: this
      character(len=2), intent(in) :: sp
    end subroutine
    module subroutine scalar_3d_tfm(this, sp)
      class(scalar_3d), intent(inout) :: this
      character(len=3), intent(in) :: sp
    end subroutine
  end interface

end module