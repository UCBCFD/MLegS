module mlegs_misc
  !> module for miscellaneous procedures
  use mlegs_envir
  implicit none
  private

  !> start and end time (in ticks)
  integer(i8), public :: start_time, end_time
  !> clock rate (ticks per second)
  real(p8) :: clock_rate

  !> print real-time clock on display
  interface print_real_time
    module subroutine print_real_time()
      implicit none
    end subroutine
  end interface
  public :: print_real_time

  !> start the timer
  interface tic
    module subroutine tic()
      implicit none
    end subroutine
  end interface
  public :: tic

  !> stop the timer and get the elapsed time
  interface toc
    module function toc() result(elapsed_time)
      implicit none
      real(p8) :: elapsed_time
    end function
  end interface
  public :: toc

  !> save a matrix (or vector/array) into a file
  interface msave
    module subroutine msavedr(a, fn, is_binary)
      implicit none
      real(p8), dimension(:,:), intent(in) :: a
      character(len=*), intent(in) :: fn
      logical, optional :: is_binary
    end subroutine
    module subroutine msavedc(a, fn, is_binary)
      implicit none
      complex(p8), dimension(:,:), intent(in) :: a
      character(len=*), intent(in) :: fn
      logical, optional :: is_binary
    end subroutine
    module subroutine msave1r(a, fn, is_binary)
      implicit none
      real(p8), dimension(:), intent(in) :: a
      character(len=*), intent(in) :: fn
      logical, optional :: is_binary
    end subroutine
    module subroutine msave1c(a, fn, is_binary)
      implicit none
      complex(p8), dimension(:), intent(in) :: a
      character(len=*), intent(in) :: fn
      logical, optional :: is_binary
    end subroutine
    module subroutine msave3r(a, fn, is_binary)
      implicit none
      real(p8), dimension(:,:,:), intent(in) :: a
      character(len=*), intent(in) :: fn
      logical, optional :: is_binary
    end subroutine
    module subroutine msave3c(a, fn, is_binary)
      implicit none
      complex(p8), dimension(:,:,:), intent(in) :: a
      character(len=*), intent(in) :: fn
      logical, optional :: is_binary
    end subroutine
  end interface
  public :: msave

  !> load a matrix (or vector/array) into a file
  interface mload
    module subroutine mloaddr(fn, a, is_binary)
      implicit none
      character(len=*), intent(in) :: fn
      real(p8), dimension(:,:) :: a
      logical, optional :: is_binary
    end subroutine
    module subroutine mloaddc(fn, a, is_binary)
      implicit none
      character(len=*), intent(in) :: fn
      complex(p8), dimension(:,:) :: a
      logical, optional :: is_binary
    end subroutine
    module subroutine mload1r(fn, a, is_binary)
      implicit none
      character(len=*), intent(in) :: fn
      real(p8), dimension(:) :: a
      logical, optional :: is_binary
    end subroutine
    module subroutine mload1c(fn, a, is_binary)
      implicit none
      character(len=*), intent(in) :: fn
      complex(p8), dimension(:) :: a
      logical, optional :: is_binary
    end subroutine
    module subroutine mload3r(fn, a, is_binary)
      implicit none
      character(len=*), intent(in) :: fn
      real(p8), dimension(:,:,:) :: a
      logical, optional :: is_binary
    end subroutine
    module subroutine mload3c(fn, a, is_binary)
      implicit none
      character(len=*), intent(in) :: fn
      complex(p8), dimension(:,:,:) :: a
      logical, optional :: is_binary
    end subroutine
  end interface
  public :: mload

  !> view a matrix in an organized way
  interface mcat
    module subroutine mcatdr(a, width, precision)
      implicit none
      real(p8), dimension(:,:) :: a
      integer(i4), optional :: width, precision
    end subroutine
    module subroutine mcatdc(a, width, precision)
      implicit none
      complex(p8), dimension(:,:) :: a
      integer(i4), optional :: width, precision
    end subroutine
    module subroutine mcat1r(a, width, precision)
      implicit none
      real(p8), dimension(:) :: a
      integer(i4), optional :: width, precision
    end subroutine
    module subroutine mcat1c(a, width, precision)
      implicit none
      complex(p8), dimension(:) :: a
      integer(i4), optional :: width, precision
    end subroutine
    module subroutine mcat3r(a, width, precision)
      implicit none
      real(p8), dimension(:,:,:) :: a
      integer(i4), optional :: width, precision
    end subroutine
    module subroutine mcat3c(a, width, precision)
      implicit none
      complex(p8), dimension(:,:,:) :: a
      integer(i4), optional :: width, precision
    end subroutine
  end interface
  public :: mcat

  !> I/O utilities
  !> integer to string for print out
  interface itoa
    module function itoa(i) result(a)
      implicit none
      integer :: i
      character(len=36) :: a
    end function
  end interface
  public :: itoa

end module