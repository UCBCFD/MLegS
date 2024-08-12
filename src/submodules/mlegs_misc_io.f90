submodule (mlegs_misc) mlegs_misc_io
  implicit none

  !> file i/o parameters
  integer(i4) :: ie, issav, ii
  !> file i/o buffer
  character(len=formatted_num_str_len) :: buf = " "

contains

  module procedure msavedr
    implicit none
    logical :: is_binary_
    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif
    call msaver(fn, reshape(a, (/size(a,1), size(a,2), 1/)), size(a,1), size(a,2), 1, is_binary_)
    write(*, 102) size(a,1), size(a,2), trim(fn)
  102 format(I6, ' x ', I6, ' 2D real matrix data written to ', A)
  end procedure

  module procedure msavedc
    implicit none
    logical :: is_binary_
    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif
    call msavec(fn, reshape(a, (/size(a,1), size(a,2), 1/)), size(a,1), size(a,2), 1, is_binary_)
    write(*, 102) size(a,1), size(a,2), trim(fn)
  102 format(I6, ' x ', I6, ' 2D complex matrix data written to ', A)
  end procedure

  module procedure msave1r
    implicit none
    logical :: is_binary_
    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif
    call msaver(fn, reshape(a, (/size(a), 1, 1/)), size(a), 1, 1, is_binary_)
    write(*, 102) size(a), trim(fn)
  102 format(I6, ' x    1 1D real vector data written to ', A)
  end procedure

  module procedure msave1c
    implicit none
    logical :: is_binary_
    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif
    call msavec(fn, reshape(a, (/size(a), 1, 1/)), size(a), 1, 1, is_binary_)
    write(*, 102) size(a), trim(fn)
  102 format(I6, ' x    1 1D complex vector data written to ', A)
  end procedure

  module procedure msave3r
    implicit none
    logical :: is_binary_
    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif
    call msaver(fn, a, size(a,1), size(a,2), size(a,3), is_binary_)
    write(*, 102) size(a,1), size(a,2), size(a,3), trim(fn)
  102 format(I6, ' x ', I6, ' x ', I6, ' real 3D array data written to ', A)
  end procedure

  module procedure msave3c
    implicit none
    logical :: is_binary_
    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif
    call msavec(fn, a, size(a,1), size(a,2), size(a,3), is_binary_)
    write(*, 102) size(a,1), size(a,2), size(a,3), trim(fn)
  102 format(I6, ' x ', I6, ' x ', I6, ' complex 3D array data written to ', A)
  end procedure

  module procedure mloaddr
    implicit none
    real(p8), dimension(:,:,:), allocatable :: a_
    logical :: is_binary_
    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif
    allocate(a_(size(a,1),size(a,2),1))
    call mloadr(fn, a_, is_binary_)
    a = reshape(a_, (/size(a,1), size(a,2)/))
    write(*, 102) size(a,1), size(a,2), trim(fn)
    deallocate(a_)
  102 format(I6, ' x ', I6, ' 2D real matrix data loaded from ', A)
  end procedure

  module procedure mloaddc
    implicit none
    complex(p8), dimension(:,:,:), allocatable :: a_
    logical :: is_binary_
    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif
    allocate(a_(size(a,1),size(a,2),1))
    call mloadc(fn, a_, is_binary_)
    a = reshape(a_, (/size(a,1), size(a,2)/))
    write(*, 102) size(a,1), size(a,2), trim(fn)
    deallocate(a_)
  102 format(I6, ' x ', I6, ' 2D complex matrix data loaded from ', A)
  end procedure

  module procedure mload1r
    implicit none
    real(p8), dimension(:,:,:), allocatable :: a_
    logical :: is_binary_
    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif
    allocate(a_(size(a),1,1))
    call mloadr(fn, a_, is_binary_)
    a = reshape(a_, (/size(a)/))
    write(*, 102) size(a), trim(fn)
    deallocate(a_)
  102 format(I6, ' x      1 1D real vector data loaded from ', A)
  end procedure

  module procedure mload1c
    implicit none
    complex(p8), dimension(:,:,:), allocatable :: a_
    logical :: is_binary_
    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif
    allocate(a_(size(a),1,1))
    call mloadc(fn, a_, is_binary_)
    a = reshape(a_, (/size(a)/))
    write(*, 102) size(a), trim(fn)
    deallocate(a_)
  102 format(I6, ' x      1 1D complex vector data loaded from ', A)
  end procedure

  module procedure mload3r
    implicit none
    logical :: is_binary_
    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif
    call mloadr(fn, a, is_binary_)
    write(*, 102) size(a,1), size(a,2), size(a,3), trim(fn)
  102 format(I6, ' x ', I6, ' x ', I6, ' complex 3D array data loaded from ', A)
  end procedure

  module procedure mload3c
    implicit none
    logical :: is_binary_
    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif
    call mloadc(fn, a, is_binary_)
    write(*, 102) size(a,1), size(a,2), size(a,3), trim(fn)
  102 format(I6, ' x ', I6, ' x ', I6, ' complex 3D array data loaded from ', A)
  end procedure

  module procedure mcatdr
    implicit none
    integer(i4) :: width_, precision_
    if (present(precision) .eqv. .false.) then
      precision_ = 3
    else
      precision_ = max(precision, 0)
    endif
    if (present(width) .eqv. .false.) then
      width_ = 8
    else
      width_ = max(width, 1)
    endif
    call mcatr(reshape(a, (/ size(a,1), size(a,2), 1 /)), width_, precision_)
  end procedure

  module procedure mcatdc
    implicit none
    integer(i4) :: width_, precision_
    if (present(precision) .eqv. .false.) then
      precision_ = 3
    else
      precision_ = max(precision, 0)
    endif
    if (present(width) .eqv. .false.) then
      width_ = 4
    else
      width_ = max(width, 1)
    endif
    call mcatc(reshape(a, (/ size(a,1), size(a,2), 1 /)), width_, precision_)
  end procedure

  module procedure mcat1r
    implicit none
    integer(i4) :: width_, precision_
    if (present(precision) .eqv. .false.) then
      precision_ = 3
    else
      precision_ = max(precision, 0)
    endif
    if (present(width) .eqv. .false.) then
      width_ = 8
    else
      width_ = max(width, 1)
    endif
    call mcatr(reshape(a, (/ size(a), 1, 1 /)), width_, precision_)
  end procedure

  module procedure mcat1c
    implicit none
    integer(i4) :: width_, precision_
    if (present(precision) .eqv. .false.) then
      precision_ = 3
    else
      precision_ = max(precision, 0)
    endif
    if (present(width) .eqv. .false.) then
      width_ = 4
    else
      width_ = max(width, 1)
    endif
    call mcatc(reshape(a, (/ size(a), 1, 1 /)), width_, precision_)
  end procedure

  module procedure mcat3r
    implicit none
    integer(i4) :: width_, precision_
    if (present(precision) .eqv. .false.) then
      precision_ = 3
    else
      precision_ = max(precision, 0)
    endif
    if (present(width) .eqv. .false.) then
      width_ = 8
    else
      width_ = max(width, 1)
    endif
    call mcatr(a, width_, precision_)
  end procedure

  module procedure mcat3c
    implicit none
    integer(i4) :: width_, precision_
    if (present(precision) .eqv. .false.) then
      precision_ = 3
    else
      precision_ = max(precision, 0)
    endif
    if (present(width) .eqv. .false.) then
      width_ = 4
    else
      width_ = max(width, 1)
    endif
    call mcatc(a, width_, precision_)
  end procedure

! ======================================================================================================== !
! VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV INTERNAL (PRIVATE) SUBROUTINES/FUNCTIONS VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV !
! ======================================================================================================== !

  subroutine mcatr(a, width, precision)
    implicit none
    real(p8), dimension(:,:,:) :: a
    integer(i4) :: width, precision

    integer(i4) :: ni, nj ,nk, i, j, k
    integer(i4) :: i, j, jb, je, jj
    character(len=72) :: fmt1, fmt2, fmt3

    ni = size(a,1)
    nj = size(a,2)
    nk = size(a,3)

    fmt1 = '(I4, A)'
    fmt2 = '(5X, '//trim(itoa(width))//'I'//trim(itoa(precision+9))//')'
    fmt3 = '(I4, ":", 1P'//trim(itoa(width))//'E'//trim(itoa(precision+9))//'.'//trim(itoa(precision))&
                                                                                                     //'E3)'

    do k = 1, nk
    if (nk .gt. 1) write(*, fmt1) k, ' (Axis 3): '
      do j = 1, (nj+(width-1))/width
        jb = width*j - (width-1)
        je = min(width*j, nj)
        write(*,*) ''
        write(*, fmt2) (jj, jj = jb, je)
        do i = 1, ni
          write(*, fmt3) i, (a(i,jj,k), jj = jb, je)
        enddo
      enddo
    write(*,*) ''
    enddo
  end subroutine

! ======================================================================================================== !

  subroutine mcatc(a, width, precision)
    implicit none
    complex(p8), dimension(:,:,:) :: a
    integer(i4) :: width, precision

    integer(i4) :: ni, nj ,nk, i, j, k
    integer(i4) :: i, j, jb, je, jj
    character(len=72) :: fmt1, fmt2, fmt3

    ni = size(a,1)
    nj = size(a,2)
    nk = size(a,3)

    fmt1 = '(I4, A)'
    fmt2 = '(5X, '//trim(itoa(2*width))//'I'//trim(itoa(2*precision+9+8+1))//')'
    fmt3 = '(I4, ":", 1P'//trim(itoa(2*width))//'(S,E'//trim(itoa(precision+9))//'.'//trim(itoa(precision))&
                                //'E3,SP,E'//trim(itoa(precision+8))//'.'//trim(itoa(precision))//'E3,"i"))'

    do k = 1, nk
    if (nk .gt. 1) write(*, fmt1) k, ' (Axis 3): '
      do j = 1, (nj+(width-1))/width
        jb = width*j - (width-1)
        je = min(width*j, nj)
        write(*,*) ''
        write(*, fmt2) (jj, jj = jb, je)
        do i = 1, ni
          write(*, fmt3) i, (a(i,jj,k), jj = jb, je)
        enddo
      enddo
    write(*,*) ''
    enddo
  end subroutine

! ======================================================================================================== !

  subroutine msaver(fn, a, ni, nj, nk, is_binary)
    implicit none
    integer(i4) :: ni, nj, nk
    real(p8), dimension(ni, nj, nk) :: a
    character(len=*) :: fn
    logical :: is_binary

    integer(i4) :: i, j, k, fo = 11
    character(:), allocatable :: dum1, dum2, bsstr

    dum1 = itoa(formatted_num_str_len)
    dum2 = itoa(formatted_num_str_len - 8)
    bsstr = '(1PE' // trim(adjustl(dum1)) // '.' // trim(adjustl(dum2)) // 'E3)'

    if (is_binary) then
      open(newunit=fo, file=fn, form='unformatted', access='stream', status='replace')
      write(fo) ni, nj, nk
      write(fo) a
    else
      open(newunit=fo, file=fn, status='unknown')
      write(unit=fo, fmt='(3(1X,I10))', advance='no') ni, nj, nk
      write(fo, *) ''
      do k = 1, nk
        do i = 1, ni
          do j = 1, nj
            write(unit=fo, fmt=trim(bsstr), advance='no') a(i,j,k)
          enddo
          if (i .lt. ni) write(fo, *) ''
        enddo
        if (k .lt. nk) then 
          write(fo, *) ''
          write(fo, *) ''
        endif
      enddo
    endif
    close(fo)
  end subroutine

! ======================================================================================================== !

  subroutine mloadr(fn, a, is_binary)
    implicit none
    character(len=*), intent(in) :: fn
    real(p8), dimension(:,:,:), intent(inout) :: a
    logical :: is_binary

    integer(i4) :: ni, nj, nk
    integer(i4) :: i, j, k, is, fo = 11
    character(:), allocatable :: dum1, dum2, dum3, bsstr

    dum1 = itoa(formatted_num_str_len)
    dum2 = itoa(formatted_num_str_len - 8)
    dum3 = itoa(nj)
    bsstr = '('//trim(adjustl(dum1))//'(1PE'//trim(adjustl(dum1))//'.'//trim(adjustl(dum2))//'E3))'

    ni = size(a,1)
    nj = size(a,2)
    nk = size(a,3)

    if (is_binary) then
      open(unit=fo, file=fn, form='unformatted', access='stream', status='old', iostat=is, action='read')
      if (is .ne. 0) then
        write(*,*) 'mloadr: cannot open ', fn
        stop 
      endif
      read(fo) i, j, k
      if ((i .ne. ni) .or. (j .ne. nj) .or. (k .ne. nk)) then
        write(*,*) 'mloadr: size inconsistency between data and array'
        write(*,103) i, ni
        write(*,104) j, nj
        write(*,105) k, nk
        stop
      endif
      read(fo) a
    else
      open(newunit=fo, file=fn, status='old')
      read(unit=fo, fmt='(3(1X,I10))') i, j, k
      if ((i .ne. ni) .or. (j .ne. nj) .or. (k .ne. nk)) then
        write(*,*) 'mloadr: size inconsistency between data and array'
        write(*,103) i, ni
        write(*,104) j, nj
        write(*,105) k, nk
        stop
      endif
      do k = 1, nk
        do i = 1, ni
          read(unit=fo, fmt=bsstr) a(i,:,k)
        enddo
        if (k .lt. nk) read(fo, *)
      enddo
    endif
    close(fo)
  103 format('axis 1 (data): ', I6, '   ', 'axis 1 (array): ', I6)
  104 format('axis 2 (data): ', I6, '   ', 'axis 2 (array): ', I6)
  105 format('axis 3 (data): ', I6, '   ', 'axis 3 (array): ', I6)
  end subroutine

! ======================================================================================================== !

  subroutine msavec(fn, a, ni, nj, nk, is_binary)
    implicit none
    integer(i4) :: ni, nj, nk
    complex(p8), dimension(ni, nj, nk) :: a
    character(len=*) :: fn
    logical :: is_binary

    integer(i4) :: i, j, k, fo = 11
    character(:), allocatable :: dum1, dum2, bsstr

    dum1 = itoa(formatted_num_str_len)
    dum2 = itoa(formatted_num_str_len - 8)
    bsstr = '(1PE' // trim(adjustl(dum1)) // '.' // trim(adjustl(dum2)) // 'E3)'

    if (is_binary) then
      open(newunit=fo, file=fn, form='unformatted', access='stream', status='replace')
      write(fo) ni, nj, nk
      write(fo) a
    else
      open(newunit=fo, file=fn, status='replace')
      write(unit=fo, fmt='(3(1X,I10))', advance='no') ni, nj, nk
      write(fo, *) ' '
      do k = 1, nk
        do i = 1, ni
          do j = 1, nj
            write(unit=fo, fmt=trim(bsstr), advance='no') real(a(i,j,k))
            write(unit=fo, fmt=trim(bsstr), advance='no') aimag(a(i,j,k))
          enddo
          if (i .lt. ni) write(fo, *) ' '
        enddo
        if (k .lt. nk) then 
          write(fo, *) ' '
          write(fo, *) ' '
        endif
      enddo
    endif
    close(fo)
  end subroutine

! ======================================================================================================== !

  subroutine mloadc(fn, a, is_binary)
    implicit none
    character(len=*), intent(in) :: fn
    complex(p8), dimension(:,:,:), intent(inout) :: a
    logical :: is_binary

    integer(i4) :: ni, nj, nk
    integer(i4) :: i, j, k, is, fo = 11
    character(:), allocatable :: dum1, dum2, dum3, bsstr
    real(p8), dimension(:), allocatable :: real_imag_pairs

    dum1 = itoa(formatted_num_str_len)
    dum2 = itoa(formatted_num_str_len - 8)
    dum3 = itoa(2*nj)
    bsstr = '('//trim(adjustl(dum1))//'(1PE'//trim(adjustl(dum1))//'.'//trim(adjustl(dum2))//'E3))'

    ni = size(a,1)
    nj = size(a,2)
    nk = size(a,3)

    if (is_binary) then
      open(unit=fo, file=fn, form='unformatted', access='stream', status='old', iostat=is, action='read')
      if (is .ne. 0) then
        write(*,*) 'mloadc: cannot open ', fn
        stop 
      endif
      read(fo) i, j, k
      if ((i .ne. ni) .or. (j .ne. nj) .or. (k .ne. nk)) then
        write(*,*) 'mloadc: size inconsistency between data and array'
        write(*,103) i, ni
        write(*,104) j, nj
        write(*,105) k, nk
        stop
      endif
      read(fo) a
    else
      open(newunit=fo, file=fn, status='old')
      read(unit=fo, fmt='(3(1X,I10))') i, j, k
      if ((i .ne. ni) .or. (j .ne. nj) .or. (k .ne. nk)) then
        write(*,*) 'mloadc: size inconsistency between data and array'
        write(*,103) i, ni
        write(*,104) j, nj
        write(*,105) k, nk
        stop
      endif
      allocate(real_imag_pairs(2*nj))
      do k = 1, nk
        do i = 1, ni
          read(unit=fo, fmt=bsstr) real_imag_pairs
          do j = 1, nj
            a(i,j,k) = real_imag_pairs(2*j-1) + iu * real_imag_pairs(2*j)
          enddo
        enddo
        if (k .lt. nk) read(fo, *)
      enddo
      deallocate(real_imag_pairs)
    endif
    close(fo)
  103 format('axis 1 (data): ', I6, '   ', 'axis 1 (array): ', I6)
  104 format('axis 2 (data): ', I6, '   ', 'axis 2 (array): ', I6)
  105 format('axis 3 (data): ', I6, '   ', 'axis 3 (array): ', I6)
  end subroutine

! ======================================================================================================== !

  function itoa(i) result(a)
    implicit none
    integer(i4) :: i
    character(len=36) :: a

    write(a, '(i36)') i
    a = trim(adjustl(a))
  end function

end submodule