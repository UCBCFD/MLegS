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

      !> chopping index setup
      procedure :: chop_offset => scalar_chop_offset
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
      integer :: comm
      integer, optional :: dims(2)
    end subroutine
  end interface
  public :: set_comm_grps

  !> chop a scalar
  interface chop
    module subroutine chop(s, tfm)
      implicit none
      class(scalar), intent(inout) :: s
      class(tfm_kit), intent(in) :: tfm
    end subroutine
  end interface
  public :: chop

  !> perform spectral transformation of a scalar
  interface trans
    module subroutine trans(s, space, tfm)
      implicit none
      class(scalar), intent(inout) :: s
      character(len=3), intent(in) :: space
      class(tfm_kit), intent(in) :: tfm
    end subroutine
  end interface
  public :: trans

  !> calculate the values of scalar at origin
  interface calcat0
    module function calcat0(s, tfm) result(calc)
      implicit none
      class(scalar) :: s
      class(tfm_kit) :: tfm
      complex(p8), dimension(:), allocatable :: calc
    end function
  end interface
  public :: calcat0

  !> calculate the values of scalar at infinity
  interface calcat1
    module function calcat1(s, tfm) result(calc)
      implicit none
      class(scalar) :: s
      class(tfm_kit) :: tfm
      complex(p8), dimension(:), allocatable :: calc
    end function
  end interface
  public :: calcat1

  !> zero out the values of scalar at infinity
  interface zeroat1
    module subroutine zeroat1(s, tfm)
      implicit none
      class(scalar), intent(inout) :: s
      class(tfm_kit), intent(in) :: tfm
    end subroutine
  end interface
  public :: zeroat1

  !> del^2_perp
  interface delsqp
    module function delsqp(s, tfm) result(so)
      implicit none
      class(scalar) :: s
      class(tfm_kit) :: tfm
      type(scalar) :: so
    end function
  end interface
  public :: delsqp

  !> inverse of del^2_perp
  interface idelsqp
    module function idelsqp(s, tfm) result(so)
      implicit none
      class(scalar) :: s
      class(tfm_kit) :: tfm
      type(scalar) :: so
    end function
  end interface
  public :: idelsqp

  !> (1-x)^2*d()/dx (equivalent to r*d()/dr)
  interface xxdx
    module function xxdx(s, tfm) result(so)
      implicit none
      class(scalar) :: s
      class(tfm_kit) :: tfm
      type(scalar) :: so
    end function
  end interface
  public :: xxdx

  !> del^2
  interface del2
    module function del2(s, tfm) result(so)
      implicit none
      class(scalar) :: s
      class(tfm_kit) :: tfm
      type(scalar) :: so
    end function
  end interface
  public :: del2

  !> inverse of del^2 (predetermined ln, or ln from scalar entry (1,1,1))
  interface idel2
    module function idel2_preln(s, tfm, preln) result (so)
      implicit none
      class(scalar) :: s
      class(tfm_kit) :: tfm
      real(p8) :: preln
      type(scalar) :: so
    end function
    module function idel2_proln(s, tfm) result (so)
      implicit none
      class(scalar) :: s
      class(tfm_kit) :: tfm
      type(scalar) :: so
    end function
  end interface
  public :: idel2

  !> helmholtz (del^2 + alpha*identity)
  interface helm
    module function helm(s, alpha, tfm) result (so)
      implicit none
      class(scalar) :: s
      real(p8) :: alpha
      class(tfm_kit) :: tfm
      type(scalar) :: so
    end function
  end interface
  public :: helm

  !> inverse of helmholtz
  interface ihelm
    module function ihelm(s, alpha, tfm) result (so)
      implicit none
      class(scalar) :: s
      real(p8) :: alpha
      class(tfm_kit) :: tfm
      type(scalar) :: so
    end function
  end interface
  public :: ihelm

  !> powered helmholtz (del^p + beta*del^2 + alpha*identity)
  interface helmp
    module function helmp(s, power, alpha, beta, tfm) result (so)
      implicit none
      class(scalar) :: s
      integer(i4) :: power
      real(p8) :: alpha, beta
      class(tfm_kit) :: tfm
      type(scalar) :: so
    end function
  end interface
  public :: helmp

  !> inverse of powered helmholtz (del^p + beta*del^2 + alpha*identity)
  interface ihelmp
    module function ihelmp(s, power, alpha, beta, tfm) result (so)
      implicit none
      class(scalar) :: s
      integer(i4) :: power
      real(p8) :: alpha, beta
      class(tfm_kit) :: tfm
      type(scalar) :: so
    end function
  end interface
  public :: ihelmp

  !> forward euler - forward euler time advancement (1st explicit order)
  !> note: we solve the following:
  !> d(s)/dt = (s_rhs_nonlin) + visc*del^2(s) + hypervisc*del^p(s)
  !>           --advct.term--   ------- diffusion (stiff) --------
  !> FEFE: F(K+1) = F(K) + DT * (NLTERM(K) + LTERM(K))
  interface fefe
    module subroutine fefe(s, s_rhs_nonlin, dt, tfm)
      implicit none
      class(scalar), intent(inout) :: s
      class(scalar), intent(in) :: s_rhs_nonlin
      real(p8), intent(in) :: dt
      class(tfm_kit) :: tfm
    end subroutine
  end interface

  !> adams bashforth - adams bashforth time advancement (2nd explicit order)
  !> note: we solve the following:
  !> d(s)/dt = (s_rhs_nonlin) + visc*del^2(s) + hypervisc*del^p(s)
  !>           --advct.term--   ------- diffusion (stiff) --------
  !> ABAB: F(K+1) = F(K) + DT * [1.5D0*(NLTERM(K)+LTERM(K)) - 0.5D0*(NLTERM(K-1)+LTERM(K-1))]
  interface abab
    module subroutine abab(s, s_p, s_rhs_nonlin, s_rhs_nonlin_p, dt, tfm, is_2nd_svis_p)
      implicit none
      class(scalar), intent(inout) :: s
      class(scalar), intent(in) :: s_p
      class(scalar), intent(in) :: s_rhs_nonlin, s_rhs_nonlin_p
      real(p8), intent(in) :: dt
      class(tfm_kit) :: tfm
      logical, optional :: is_2nd_svis_p
    end subroutine
  end interface

  !> forward euler - backward euler time advancement (1st semi-implicit order)
  !> note: we solve the following:
  !> d(s)/dt = (s_rhs_nonlin) + visc*del^2(s) + hypervisc*del^p(s)
  !>           --advct.term--   ------- diffusion (stiff) --------
  !> FE : F(K+1/2) = F(K) + DT * NLTERM(K)
  !> BE : F(K+1) = F(K+1/2) + DT * LTERM(K+1)
  interface febe
    module subroutine febe(s, s_rhs_nonlin, dt, tfm)
      implicit none
      class(scalar), intent(inout) :: s
      class(scalar), intent(in) :: s_rhs_nonlin
      real(p8), intent(in) :: dt
      class(tfm_kit) :: tfm
    end subroutine
  end interface

  !> adams bashforth - crank nicolson time advancement (2nd semi-implicit order)
  !> note: we solve the following:
  !> d(s)/dt = (s_rhs_nonlin) + visc*del^2(s) + hypervisc*del^p(s)
  !>           --advct.term--   ------- diffusion (stiff) --------
  !> AB : F(K+1/2) = F(K) + DT * [1.5D0*NLTERM(K) - 0.5D0*NLTERM(K-1)]
  !> CN : F(K+1) = F(K+1/2) + DT/2 * [LTERM(K+1) + LTERM(K)]
  interface abcn
    module subroutine abcn(s, s_p, s_rhs_nonlin, s_rhs_nonlin_p, dt, tfm)
      implicit none
      class(scalar), intent(inout) :: s
      class(scalar), intent(in) :: s_p
      class(scalar), intent(in) :: s_rhs_nonlin, s_rhs_nonlin_p
      real(p8), intent(in) :: dt
      class(tfm_kit) :: tfm
    end subroutine
  end interface

  !> i/o save
  interface msave
    module subroutine msave_scalar(s, fn, is_binary, is_global)
      implicit none
      class(scalar),intent(in) :: s
      character(len=*), intent(in) :: fn
      logical, optional :: is_binary, is_global
    end subroutine
  end interface
  public :: msave

  !> i/o load
  interface mload
    module subroutine mload_scalar(fn, s, is_binary, is_global)
      implicit none
      character(len=*), intent(in) :: fn
      class(scalar),intent(inout) :: s
      logical, optional :: is_binary, is_global
    end subroutine
  end interface
  public :: mload

end module