submodule (mlegs_misc) mlegs_misc_timer
  implicit none

contains

  module procedure print_real_time
    character(len=10) ::  dum
    integer(p4), dimension(8) :: vals
    call date_and_time(dum, dum, dum, vals)
    write(*,101) vals(1), vals(2), vals(3), vals(5), vals(6), vals(7)
    write(*,*) ''
    101 format(' @ ',I0.4,'-',I0.2,'-',I0.2,' ',I0.2,':',I0.2,':',I0.2)
  end procedure

  module procedure tic
    call system_clock(count=start_time, count_rate=clock_rate)
  end procedure

  module procedure toc
    integer(i8) :: elapsed_ticks
    call system_clock(count=end_time)
    elapsed_ticks = end_time - start_time
    elapsed_time = elapsed_ticks / clock_rate
  end procedure

end submodule