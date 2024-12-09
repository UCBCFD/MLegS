module mlegs_genmat
  !> module for general matrix operations
  use mlegs_envir
  implicit none
  private

  !> matrix/vector multiplication
  interface operator(.mul.)
    module function mulrr(a, b) result(c)
      implicit none
      real(p8), dimension(:,:), intent(in) :: a, b
      real(p8), dimension(size(a,1),size(b,2)) :: c
    end function
    module function mulcc(a, b) result(c)
      implicit none
      complex(p8), dimension(:,:), intent(in) :: a, b
      complex(p8), dimension(size(a,1),size(b,2)) :: c
    end function
    module function mulrc(a, b) result(c)
      implicit none
      real(p8), dimension(:,:), intent(in) :: a
      complex(p8), dimension(:,:), intent(in) :: b
      complex(p8), dimension(size(a,1),size(b,2)) :: c
    end function
    module function mulcr(a, b) result(c)
      implicit none
      complex(p8), dimension(:,:), intent(in) :: a
      real(p8), dimension(:,:), intent(in) :: b
      complex(p8), dimension(size(a,1),size(b,2)) :: c
    end function
    module function mulrr1(a, b) result(c)
      implicit none
      real(p8), dimension(:,:), intent(in) :: a 
      real(p8), dimension(:), intent(in) :: b
      real(p8), dimension(size(a,1)) :: c
    end function
    module function mulcc1(a, b) result(c)
      implicit none
      complex(p8), dimension(:,:), intent(in) :: a
      complex(p8), dimension(:), intent(in) :: b
      complex(p8), dimension(size(a,1)) :: c
    end function
    module function mulrc1(a, b) result(c)
      implicit none
      real(p8), dimension(:,:), intent(in) :: a
      complex(p8), dimension(:), intent(in) :: b
      complex(p8), dimension(size(a,1)) :: c
    end function
    module function mulcr1(a, b) result(c)
      implicit none
      complex(p8), dimension(:,:), intent(in) :: a
      real(p8), dimension(:), intent(in) :: b
      complex(p8), dimension(size(a,1)) :: c
    end function
    module function mulr1r(a, b) result(c)
      implicit none
      real(p8), dimension(:), intent(in) :: a 
      real(p8), dimension(:,:), intent(in) :: b
      real(p8), dimension(size(a), size(b,2)) :: c
    end function
    module function mulc1c(a, b) result(c)
      implicit none
      complex(p8), dimension(:), intent(in) :: a
      complex(p8), dimension(:,:), intent(in) :: b
      complex(p8), dimension(size(a), size(b,2)) :: c
    end function
    module function mulr1c(a, b) result(c)
      implicit none
      real(p8), dimension(:), intent(in) :: a
      complex(p8), dimension(:,:), intent(in) :: b
      complex(p8), dimension(size(a), size(b,2)) :: c
    end function
    module function mulc1r(a, b) result(c)
      implicit none
      complex(p8), dimension(:), intent(in) :: a
      real(p8), dimension(:,:), intent(in) :: b
      complex(p8), dimension(size(a), size(b,2)) :: c
    end function
  end interface
  public :: operator(.mul.)

  !> LU decomposition
  interface lu
    module subroutine lur(a, ipiv)
      real(p8), dimension(:,:), intent(inout) :: a
      integer(i4), dimension(:), optional :: ipiv
    end subroutine
    module subroutine luc(a, ipiv)
      complex(p8), dimension(:,:), intent(inout) :: a
      integer(i4), dimension(:), optional :: ipiv
    end subroutine
  end interface
  public :: lu

  !> solve a linear system
  interface lsolve
    module subroutine solvecc(a, b)
      complex(p8), dimension(:,:), intent(inout) :: a
      complex(p8), dimension(:,:), intent(inout) :: b
    end subroutine
    module subroutine solvecc1(a, b)
      complex(p8), dimension(:,:), intent(inout) :: a
      complex(p8), dimension(:), intent(inout) :: b
    end subroutine
    module subroutine solverc(a, b)
      real(p8), dimension(:,:), intent(inout) :: a
      complex(p8), dimension(:,:), intent(inout) :: b
    end subroutine
    module subroutine solverc1(a, b)
      real(p8), dimension(:,:), intent(inout) :: a
      complex(p8), dimension(:), intent(inout) :: b
    end subroutine
    module subroutine solverr(a, b)
      real(p8), dimension(:,:), intent(inout) :: a
      real(p8), dimension(:,:), intent(inout) :: b
    end subroutine
    module subroutine solverr1(a, b)
      real(p8), dimension(:,:), intent(inout) :: a
      real(p8), dimension(:), intent(inout) :: b
    end subroutine
  end interface
  public :: lsolve

  !> matrix inversion
  interface inv
    module function invr(a) result(ainv)
      implicit none
      real(p8), dimension(:,:), intent(in) :: a
      real(p8), dimension(size(a,1), size(a,2)) :: ainv
    end function
    module function invc(a) result(ainv)
      implicit none
      complex(p8), dimension(:,:), intent(in) :: a
      complex(p8), dimension(size(a,1), size(a,2)) :: ainv
    end function
  end interface
  public :: inv

  !> obtain eigenvalues and eigenvectors
  interface eig
    module subroutine geigrr(a, b, ev, er, el) ! A(R) = B(R)D // (L**H)A = D(L**H)B
      implicit none
      real(p8), dimension(:,:), intent(in) :: a, b
      complex(p8), dimension(size(a,1)), intent(inout) :: ev
      complex(p8), dimension(size(a,1),size(a,1)), optional :: er, el
    end subroutine
    module subroutine steigr(a, ev, er, el) ! A(R) = (R)D // (L**H)A = D(L**H)
      implicit none
      real(p8), dimension(:,:), intent(in) :: a
      complex(p8), dimension(size(a,1)), intent(inout) :: ev
      complex(p8), dimension(size(a,1),size(a,1)), optional :: er, el
    end subroutine
    module subroutine geigcc(a, b, ev, er, el) ! same as geigrr, but complex
      implicit none
      complex(p8), dimension(:,:), intent(in) :: a, b
      complex(p8), dimension(size(a,1)), intent(inout) :: ev
      complex(p8), dimension(size(a,1),size(a,1)), optional :: er, el
    end subroutine
    module subroutine geigcr(a, b, ev, er, el)
      implicit none
      complex(p8), dimension(:,:), intent(in) :: a
      real(p8), dimension(:,:), intent(in) :: b
      complex(p8), dimension(size(a,1)), intent(inout) :: ev
      complex(p8), dimension(size(a,1),size(a,1)), optional :: er, el
    end subroutine
    module subroutine geigrc(a, b, ev, er, el)
      implicit none
      real(p8), dimension(:,:), intent(in) :: a
      complex(p8), dimension(:,:), intent(in) :: b
      complex(p8), dimension(size(a,1)), intent(inout) :: ev
      complex(p8), dimension(size(a,1),size(a,1)), optional :: er, el
    end subroutine
    module subroutine steigc(a, ev, er, el) ! same as steigr, but complex
      implicit none
      complex(p8), dimension(:,:), intent(in) :: a
      complex(p8), dimension(size(a,1)), intent(inout) :: ev
      complex(p8), dimension(size(a,1),size(a,1)), optional :: er, el
    end subroutine
  end interface
  public :: eig

end module