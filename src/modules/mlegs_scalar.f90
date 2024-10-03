module mlegs_scalar
    use mlegs_envir
    use mlegs_misc, only: itoa
    use MPI

    implicit none
    private

    !> class: distributed scalar
    type, public :: scalar
        !> glboal size is distributed into local size on ea proc
        integer :: glb_sz(3), loc_sz(3), loc_st(3)
        !> ea axis' sub communicator index
        !> axis_comm(i) = 0: i-th axis is not distributed
        !> axis_comm(i) = j: i-th axis is distributed by communicator group comm_grps(j)
        integer :: axis_comm(3)
        !> local distributed data
        complex(p8), dimension(:,:,:), pointer :: e
        !> logarithmic term coeff. to describe O(ln(r)) as r -> infty. (P_ln == ln(r^2+L^2/2L^2))
        real(p8) :: ln

        !> methods
        contains
        ! procedure :: set => scalar_set
        procedure :: initialise => scalar_initialise
        procedure :: alloc => scalar_alloc
        procedure :: dealloc => scalar_dealloc

        ! procedure :: chopset => scalar_chopset
        ! procedure :: chopdo => scalar_chopdo

        ! MPI-related
        procedure :: exchange => scalar_exchange
        procedure :: assemble => scalar_assemble ! right now, only supports 3D: axis = 1 or 2
        procedure :: disassemble => scalar_disassemble ! right now, only supports 3D: axis = 1 or 2
    end type

    !> scalar: class methods
    interface
        ! module subroutine scalar_set(this, kit)
        ! implicit none
        ! class(scalar), intent(inout) :: this
        ! class(tfm_kit), intent(in) :: kit
        ! end subroutine

        ! !> adjusting spectral truncation settings
        ! module subroutine scalar_chopset(this, iof, iofp, iofz)
        ! implicit none
        ! class(scalar), intent(inout) :: this
        ! integer(i4), intent(in) :: iof
        ! integer(i4), optional :: iofp, iofz
        ! end subroutine

        ! !> perform spectral truncation
        ! module subroutine scalar_chopdo(this)
        ! implicit none
        ! class(scalar), intent(inout) :: this
        ! end subroutine
        
        !> initialise a scalar: glb_sz, loc_sz, loc_st, and axis_comm
        module subroutine scalar_initialise(this, glb_sz, axis_comm)
            implicit none
            class(scalar), intent(inout) :: this
            integer, intent(in) :: glb_sz(3)
            integer, allocatable, intent(in) :: axis_comm(:)
        end subroutine

        !> NOTE: add overloaded initialise & alloc that takes XXX_SPACE tag (public)
        !> OR create extended type for 3D FFT which directly create loc_sz and loc_st based on XXX_SPACE

        !> allocate/deallocate a scalr
        module subroutine scalar_alloc(this)
            implicit none
            class(scalar), intent(inout) :: this
        end subroutine
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
        module function scalar_assemble(this,axis) result(array_glb)
            implicit none
            class(scalar), intent(in) :: this
            complex(P8), dimension(:,:,:), allocatable :: array_glb
            integer :: axis
        end function
        module subroutine scalar_disassemble(this,axis,array_glb)
            implicit none
            class(scalar), intent(in) :: this
            complex(P8), dimension(:,:,:), allocatable :: array_glb
            integer :: axis
        end subroutine

    end interface  

    !> set up cartesian communicator groups
    interface set_comm_grps
        module subroutine subcomm_cart_2d(comm)
            implicit none
            integer:: comm
        end subroutine
    end interface
    public :: set_comm_grps

end module