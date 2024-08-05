module mlegs_misc
  use mlegs_base
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

end module