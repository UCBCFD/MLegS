module mlegs_bndmat
  use mlegs_envir
  implicit none
  private

  ! real banded matrix derived data type
  type, public :: real_bndm
    real(p8), dimension(:,:), allocatable :: e ! elements
    integer(i4) :: sublen = 0 ! # of subdiagonals
    integer(i4) :: suplen = 0 ! # of superdiagonals
    contains 
      procedure :: alloc => real_bndm_alloc
      procedure :: dealloc => real_bndm_dealloc
      procedure :: t => real_bndm_transpose
  end type

  interface
    module subroutine real_bndm_alloc(this, nj, sublen, suplen)
      class(real_bndm), intent(inout) :: this
      integer(i4) :: nj
      integer(i4), optional :: sublen, suplen
    end subroutine
    module subroutine real_bndm_dealloc(this)
      class(real_bndm), intent(inout) :: this
    end subroutine
    module function real_bndm_transpose(this) result(tp)
      class(real_bndm) :: this
      type(real_bndm) :: tp
    end function
  end interface

  ! complex banded matrix derived data type
  type, public :: cmpx_bndm
    complex(p8), dimension(:,:), allocatable :: e ! elements
    integer(i4) :: sublen = 0 ! # of subdiagonals
    integer(i4) :: suplen = 0 ! # of superdiagonals
    contains 
      procedure :: alloc => cmpx_bndm_alloc
      procedure :: dealloc => cmpx_bndm_dealloc
      procedure :: t => cmpx_bndm_transpose
      procedure :: real => cmpx_real
      procedure :: aimag => cmpx_aimag
  end type

  interface
    module subroutine cmpx_bndm_alloc(this, nj, sublen, suplen)
      implicit none
      class(cmpx_bndm), intent(inout) :: this
      integer(i4) :: nj
      integer(i4), optional :: sublen, suplen
    end subroutine
    module subroutine cmpx_bndm_dealloc(this)
      implicit none
      class(cmpx_bndm), intent(inout) :: this
    end subroutine
    module function cmpx_bndm_transpose(this) result(tp)
      class(cmpx_bndm) :: this
      type(cmpx_bndm) :: tp
    end function
    module function cmpx_real(this) result(real_part)
      implicit none
      class(cmpx_bndm) :: this
      type(real_bndm) :: real_part
    end function
    module function cmpx_aimag(this) result(imag_part)
      implicit none
      class(cmpx_bndm) :: this
      type(real_bndm) :: imag_part
    end function
  end interface

  !> bndmat assignment (conversion between real & cmpx like intrinsic floats)
  interface assignment(=)
    module subroutine realbnd_eq_cmpxbnd(a, b)
      type(real_bndm), intent(out) :: a
      type(cmpx_bndm), intent(in) :: b
    end subroutine
    module subroutine cmpxbnd_eq_realbnd(a, b)
      type(cmpx_bndm), intent(out) :: a
      type(real_bndm), intent(in) :: b
    end subroutine
  end interface
  public :: assignment(=)

  !> band-to-generic conversion
  interface fullmat
    module function fullmatr(bndm) result(genm)
      implicit none
      type(real_bndm) :: bndm
      real(p8), dimension(:,:), allocatable :: genm
    end function
    module function fullmatc(bndm) result(genm)
      implicit none
      type(cmpx_bndm) :: bndm
      complex(p8), dimension(:,:), allocatable :: genm
    end function
  end interface
  public :: fullmat

  !> generic-to-band conversion
  interface bandmat
    module function bandmatr(genm, sublen, suplen) result(bndm)
      implicit none
      real(p8), dimension(:,:), intent(in) :: genm
      integer(i4) :: sublen, suplen
      type(real_bndm) :: bndm
    end function
    module function bandmatc(genm, sublen, suplen) result(bndm)
      implicit none
      complex(p8), dimension(:,:), intent(in) :: genm
      integer(i4) :: sublen, suplen
      type(cmpx_bndm) :: bndm
    end function
  end interface
  public :: bandmat

  !> matrix multiplications
  interface operator(.mul.)
    module function mulccb(a, b) result(c)
      implicit none
      complex(p8), dimension(:,:), intent(in) :: a
      type(cmpx_bndm), intent(in) :: b
      complex(p8), dimension(:,:), allocatable :: c
    end function
    module function mulc1cb(a, b) result(c)
      implicit none
      complex(p8), dimension(:), intent(in) :: a
      type(cmpx_bndm), intent(in) :: b
      complex(p8), dimension(:,:), allocatable :: c
    end function
    module function mulcbc(a, b) result(c)
      implicit none
      type(cmpx_bndm), intent(in) :: a
      complex(p8), dimension(:,:), intent(in) :: b
      complex(p8), dimension(:,:), allocatable :: c
    end function
    module function mulcbc1(a, b) result(c)
      implicit none
      type(cmpx_bndm), intent(in) :: a
      complex(p8), dimension(:), intent(in) :: b
      complex(p8), dimension(:), allocatable :: c
    end function
    module function mulrcb(a, b) result(c)
      implicit none
      real(p8), dimension(:,:), intent(in) :: a
      type(cmpx_bndm), intent(in) :: b
      complex(p8), dimension(:,:), allocatable :: c
    end function
    module function mulr1cb(a, b) result(c)
      implicit none
      real(p8), dimension(:), intent(in) :: a
      type(cmpx_bndm), intent(in) :: b
      complex(p8), dimension(:,:), allocatable :: c
    end function
    module function mulcbr(a, b) result(c)
      implicit none
      type(cmpx_bndm), intent(in) :: a
      real(p8), dimension(:,:), intent(in) :: b
      complex(p8), dimension(:,:), allocatable :: c
    end function
    module function mulcbr1(a, b) result(c)
      implicit none
      type(cmpx_bndm), intent(in) :: a
      real(p8), dimension(:), intent(in) :: b
      complex(p8), dimension(:), allocatable :: c
    end function
    module function mulcrb(a, b) result(c)
      implicit none
      complex(p8), dimension(:,:), intent(in) :: a
      type(real_bndm), intent(in) :: b
      complex(p8), dimension(:,:), allocatable :: c
    end function
    module function mulc1rb(a, b) result(c)
      implicit none
      complex(p8), dimension(:), intent(in) :: a
      type(real_bndm), intent(in) :: b
      complex(p8), dimension(:,:), allocatable :: c
    end function
    module function mulrbc(a, b) result(c)
      implicit none
      type(real_bndm), intent(in) :: a
      complex(p8), dimension(:,:), intent(in) :: b
      complex(p8), dimension(:,:), allocatable :: c
    end function
    module function mulrbc1(a, b) result(c)
      implicit none
      type(real_bndm), intent(in) :: a
      complex(p8), dimension(:), intent(in) :: b
      complex(p8), dimension(:), allocatable :: c
    end function
    module function mulrrb(a, b) result(c)
      implicit none
      real(p8), dimension(:,:), intent(in) :: a
      type(real_bndm), intent(in) :: b
      real(p8), dimension(:,:), allocatable :: c
    end function
    module function mulr1rb(a, b) result(c)
      implicit none
      real(p8), dimension(:), intent(in) :: a
      type(real_bndm), intent(in) :: b
      real(p8), dimension(:,:), allocatable :: c
    end function
    module function mulrbr(a, b) result(c)
      implicit none
      type(real_bndm), intent(in) :: a
      real(p8), dimension(:,:), intent(in) :: b
      real(p8), dimension(:,:), allocatable :: c
    end function
    module function mulrbr1(a, b) result(c)
      implicit none
      type(real_bndm), intent(in) :: a
      real(p8), dimension(:), intent(in) :: b
      real(p8), dimension(:), allocatable :: c
    end function
    module function mulcbcb(a, b) result(c)
      implicit none
      type(cmpx_bndm), intent(in) :: a
      type(cmpx_bndm), intent(in) :: b
      type(cmpx_bndm) :: c
    end function
    module function mulrbrb(a, b) result(c)
      implicit none
      type(real_bndm), intent(in) :: a
      type(real_bndm), intent(in) :: b
      type(real_bndm) :: c
    end function
    module function mulrbcb(a, b) result(c)
      implicit none
      type(real_bndm), intent(in) :: a
      type(cmpx_bndm), intent(in) :: b
      type(cmpx_bndm) :: c
    end function
    module function mulcbrb(a, b) result(c)
      implicit none
      type(cmpx_bndm), intent(in) :: a
      type(real_bndm), intent(in) :: b
      type(cmpx_bndm) :: c
    end function
  end interface
  public :: operator(.mul.)

  !> LU decomposition
  interface lu
    module subroutine lurb(a, ipiv)
      type(real_bndm), intent(inout) :: a
      integer(i4), dimension(:), optional :: ipiv
    end subroutine
    module subroutine lucb(a, ipiv)
      type(cmpx_bndm), intent(inout) :: a
      integer(i4), dimension(:), optional :: ipiv
    end subroutine
  end interface
  public :: lu

  !> solve a linear system
  interface lsolve
    module subroutine solvecbc(a, b)
      type(cmpx_bndm), intent(inout) :: a
      complex(p8), dimension(:,:), intent(inout) :: b
    end subroutine
    module subroutine solvecbc1(a, b)
      type(cmpx_bndm), intent(inout) :: a
      complex(p8), dimension(:), intent(inout) :: b
    end subroutine
    module subroutine solverbc(a, b)
      type(real_bndm), intent(inout) :: a
      complex(p8), dimension(:,:), intent(inout) :: b
    end subroutine
    module subroutine solverbc1(a, b)
      type(real_bndm), intent(inout) :: a
      complex(p8), dimension(:), intent(inout) :: b
    end subroutine
    module subroutine solverbr(a, b)
      type(real_bndm), intent(inout) :: a
      real(p8), dimension(:,:), intent(inout) :: b
    end subroutine
    module subroutine solverbr1(a, b)
      type(real_bndm), intent(inout) :: a
      real(p8), dimension(:), intent(inout) :: b
    end subroutine
  end interface
  public :: lsolve

end module