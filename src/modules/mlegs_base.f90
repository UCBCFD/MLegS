module mlegs_base
  implicit none
  private

  !> 32-bit (single-precision) real/complex
  integer, public, parameter :: p4 = selected_real_kind(6, 37) 
  !> 64-bit (double-precision) real/complex
  integer, public, parameter :: p8 = selected_real_kind(15, 307)
  !> 32-bit (default) integer
  integer, public, parameter :: i4 = selected_int_kind(9)
  !> 64-bit (long) integer
  integer, public, parameter :: i8 = selected_int_kind(18)
  !> pi == 3.1415926535 ...
  real(p8), public, parameter :: pi = acos(-1.D0)
  !> pure imaginary unit (1i)
  complex(p8), public, parameter :: iu = (0.D0, 1.D0)
  !> string size (per number to be stored) for formatted data save/load
  integer(i4), public, parameter :: formatted_num_str_len = 24

end module