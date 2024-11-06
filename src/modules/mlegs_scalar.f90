module mlegs_scalar
  !> module for distributed scalar class
  use MPI
  use mlegs_envir
  use mlegs_base
  use mlegs_misc
  use mlegs_genmat
  use mlegs_bndmat
  use mlegs_spectfm

  implicit none
  private

  !> class: distributed scalar
  type, public :: scalar
    !> global size is distributed into local size on ea proc
    integer(i4) :: glb_sz(3), loc_sz(3), loc_st(3)
    
    !> ea axis' sub communicator index
    !> axis_comm(i) = 0: i-th axis is not distributed
    !> axis_comm(i) = j: i-th axis is distributed by communicator group comm_grps(j)
    integer(i4) :: axis_comm(3)
    
    !> local distributed data
    complex(p8), dimension(:,:,:), pointer :: e

    !> log term expression. used for toroidal-poloidal decomposition
    real(p8) :: ln = 0.D0

    !> spectral element chopping offset
    integer(i4) :: nrchop_offset = 0, npchop_offset = 0, nzchop_offset = 0

    !> field status indicator in the order of radial, azimuthal, and axial directions. 
    character(len=3) :: space ! 'PPP', 'PFP', 'FFP' or 'FFF'

    !> procedures
    contains
      !> initialization
      procedure :: init => scalar_init
      procedure :: dealloc => scalar_dealloc
      
      !> MPI distribution for parallel processing
      procedure :: exchange => scalar_exchange
      procedure :: assemble => scalar_assemble
      procedure :: disassemble => scalar_disassemble

      procedure :: chop_offset => scalar_chop_offset
      ! procedure :: calcat0 => scalar_calcat0
      ! procedure :: calcat1 => scalar_calcat1
      ! procedure :: zeroat1 => scalar_zeroat1
      ! procedure :: msave => scalar_msave
      ! procedure :: mload => scalar_mload
  end type

  interface !> type bound procedures
    !> initialize a scalar: glb_sz, loc_sz, loc_st, and axis_comm
    module subroutine scalar_init(this, glb_sz, axis_comm)
      implicit none
      class(scalar), intent(inout) :: this
      integer(i4), dimension(3), intent(in) :: glb_sz
      integer(i4), intent(in) :: axis_comm(:)
    end subroutine
    !> deallocate a scalar
    module subroutine scalar_dealloc(this)
      implicit none
      class(scalar), intent(inout) :: this
    end subroutine
    !> re-orient distributed data
    !> before: data is complete (not distributed) along axis_old for ea proc
    !> after : data is complete (not distributed) along axis_new for ea proc
    module subroutine scalar_exchange(this,axis_old,axis_new)
      implicit none
      class(scalar), intent(inout) :: this
      integer, intent(in) :: axis_old, axis_new
    end subroutine
    !> assemble/disassemble distributed data into/from a single proc
    module recursive function scalar_assemble(this,axis_input) result(array_glb)
      implicit none
      class(scalar), intent(in) :: this
      complex(P8), dimension(:,:,:), allocatable :: array_glb
      integer(i4), optional :: axis_input
    end function
    module recursive subroutine scalar_disassemble(this,array_glb,axis_input)
      implicit none
      class(scalar), intent(inout) :: this
      complex(P8), dimension(:,:,:), allocatable :: array_glb
      integer(i4), optional :: axis_input
    end subroutine
    !> set up new chopping offsets for a scalar field
    module subroutine scalar_chop_offset(this, iof1, iof2, iof3)
      implicit none
      class(scalar), intent(inout) :: this
      integer(i4), intent(in) :: iof1 ! offset of nrchop
      integer(i4), optional :: iof2 ! offset of npchop (default is 0)
      integer(i4), optional :: iof3 ! offset of nzchop (default is 0)
    end subroutine
    !> set up new chopping offsets for a scalar field
    module subroutine scalar_chop_do(this)
      implicit none
      class(scalar), intent(inout) :: this
    end subroutine
  end interface

  !> copy a scalar
  interface assignment(=)
    module subroutine scalar_copy(this,that)
      implicit none
      class(scalar), intent(inout) :: this
      class(scalar), intent(in) :: that
    end subroutine
  end interface
  public :: assignment(=)

  !> set up cartesian communicator groups
  interface set_comm_grps
    module subroutine subcomm_cart_2d(comm,dims)
      implicit none
      integer:: comm
      integer, optional:: dims(2)
    end subroutine
  end interface
  public :: set_comm_grps

  !> chop a scalar
  interface chop
    module subroutine chop(s, tfm)
      implicit none
      class(scalar) :: s
      class(tfm_kit) :: tfm
    end subroutine
  end interface
  public :: chop

  !> perform spectral transformation of a scalar
  interface trans
    module subroutine trans(s, space, tfm)
      implicit none
      class(scalar) :: s
      character(len=3) :: space
      class(tfm_kit) :: tfm
    end subroutine
  end interface
  public :: trans

end module