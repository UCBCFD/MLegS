module mlegs_scalar
  !> module for scalar class
  use mlegs_envir
  use mlegs_misc, only: itoa
  use MPI

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
        
        !> logarithmic term coeff. to describe O(ln(r)) as r -> infty. (P_ln == ln(r^2+L^2/2L^2))
        real(p8) :: ln

        !> methods
        contains
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
            integer, intent(in) :: axis_comm(:)
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

end module