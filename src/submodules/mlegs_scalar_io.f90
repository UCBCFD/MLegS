submodule (mlegs_scalar) mlegs_scalar_io
  implicit none

contains

  module procedure msave_scalar
    implicit none
    complex(p8), dimension(:,:,:), allocatable :: a
    logical :: is_binary_, is_global_
    integer(i4) :: status, i, j, k, fo
    character(:), allocatable :: dum1, dum2, bsstr
    type(scalar) :: st

    dum1 = ntoa(formatted_num_str_len)
    dum2 = ntoa(formatted_num_str_len - 9)
    bsstr = '(1PE' // trim(adjustl(dum1)) // '.' // trim(adjustl(dum2)) // 'E3)'

    fo = 11 + rank_glb

    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif

    if (present(is_binary) .eqv. .false.) then 
      is_global_ = .true.
    else
      is_global_ = is_global
    endif

    if (is_global_) then
      if (s%axis_comm(3) .eq. 0) then
        call st%init(s%glb_sz, s%axis_comm); st = s
        call st%exchange(3, 1) ! because %assemble cannot work when dim3 resides locally
        a = st%assemble() ! collect all local arrays into a global array at proc 1
        call st%dealloc()
      else
        a = s%assemble() ! collect all local arrays into a global array at proc 1
      endif
      if (rank_glb .eq. 0) then
        if (is_binary_) then
          open(newunit=fo, file=fn, form='unformatted', access='stream', status='replace')
          write(fo) size(a,1), size(a,2), size(a,3)
          write(fo) a
          write(fo) s%ln, s%nrchop_offset, s%npchop_offset, s%nzchop_offset
          write(fo) s%space
        else
          open(newunit=fo, file=fn, status='unknown')
          write(unit=fo, fmt='(3(1X,I10))', advance='no') size(a,1), size(a,2), size(a,3) 
          write(fo, *) ''
          do k = 1, size(a,3)
            do i = 1, size(a,1)
              do j = 1, size(a,2)
                write(unit=fo, fmt=trim(bsstr), advance='no') real(a(i,j,k))
                write(unit=fo, fmt=trim(bsstr), advance='no') aimag(a(i,j,k))
              enddo
              if (i .lt. size(a,1)) write(fo, *) ''
            enddo
            if (k .lt. size(a,3)) then
              write(fo, *) ''
              write(fo, *) ''
            endif
          enddo
          write(fo ,*) ''
          write(fo, *) s%ln, s%nrchop_offset, s%npchop_offset, s%nzchop_offset
          write(fo, *) s%space
        endif
        close(fo)
      endif
      call MPI_barrier(comm_glb, MPI_err)
      deallocate( a )
    else
      if (is_binary_) then
        open(newunit=fo, file=trim(adjustl(fn))//'_'//trim(adjustl(ntoa(rank_glb))), &
             form='unformatted', access='stream', status='replace')
        write(fo) s%glb_sz, s%loc_sz, s%loc_st
        write(fo) s%e
        write(fo) s%ln, s%nrchop_offset, s%npchop_offset, s%nzchop_offset
        write(fo) s%space
      else
        open(newunit=fo, file=trim(adjustl(fn))//'_'//trim(adjustl(ntoa(rank_glb))), &
             status='unknown')
        write(unit=fo, fmt='(3(1X,I10))') s%glb_sz
        write(unit=fo, fmt='(3(1X,I10))') s%loc_sz
        write(unit=fo, fmt='(3(1X,I10))') s%loc_st 
        do k = 1, s%loc_sz(3)
          do i = 1, s%loc_sz(1)
            do j = 1, s%loc_sz(2)
              write(unit=fo, fmt=trim(bsstr), advance='no') real(s%e(i,j,k))
              write(unit=fo, fmt=trim(bsstr), advance='no') aimag(s%e(i,j,k))
            enddo
            if (i .lt. s%loc_sz(1)) write(fo, *) ''
          enddo
          if (k .lt. s%loc_sz(3)) then
            write(fo, *) ''
            write(fo, *) ''
          endif
        enddo
        write(fo, *) ''
        write(fo, *) s%ln, s%nrchop_offset, s%npchop_offset, s%nzchop_offset
        write(fo, *) s%space
      endif
      close(fo)
      call MPI_barrier(comm_glb, MPI_err)
    endif
  end procedure

  module procedure mload_scalar
    implicit none
    complex(p8), dimension(:,:,:), allocatable :: a
    logical :: is_binary_, is_global_
    integer(i4) :: status, i, j, k, is, fo
    character(:), allocatable :: dum1, dum2, bsstr
    real(p8), dimension(:), allocatable :: real_imag_pairs
    type(scalar) :: st

    dum1 = ntoa(formatted_num_str_len)
    dum2 = ntoa(formatted_num_str_len - 9)
    bsstr = '(1PE'//trim(adjustl(dum1))//'.'//trim(adjustl(dum2))//'E3))'

    fo = 11 + rank_glb

    if (present(is_binary) .eqv. .false.) then 
      is_binary_ = .false.
    else
      is_binary_ = is_binary
    endif

    if (present(is_binary) .eqv. .false.) then 
      is_global_ = .true.
    else
      is_global_ = is_global
    endif

    if (is_global_) then
      allocate( a(s%glb_sz(1), s%glb_sz(2), s%glb_sz(3)) )
      if (rank_glb .eq. 0) then
        if (is_binary_) then
          open(newunit=fo, file=fn, form='unformatted', access='stream', &
               status='old', iostat=is, action='read')
          if (is .ne. 0) then
            write(*,*) 'mload_scalar: cannot open ', fn
            stop
          endif
          read(fo) i, j, k
          if ((i .ne. s%glb_sz(1)) .or. (j .ne. s%glb_sz(2)) .or. (k .ne. s%glb_sz(3))) then
            write(*,*) 'mloadc: size inconsistency between data and array'
            write(*,103) i, s%glb_sz(1)
            write(*,104) j, s%glb_sz(2)
            write(*,105) k, s%glb_sz(3)
            stop
          endif
          read(fo) a
          read(fo) s%ln, s%nrchop_offset, s%npchop_offset, s%nzchop_offset
          read(fo) s%space
        else
          open(newunit=fo, file=fn, status='old')
          read(unit=fo, fmt='(3(1X,I10))') i, j, k
          if ((i .ne. s%glb_sz(1)) .or. (j .ne. s%glb_sz(2)) .or. (k .ne. s%glb_sz(3))) then
            write(*,*) 'mloadc: size inconsistency between data and array'
            write(*,103) i, s%glb_sz(1)
            write(*,104) j, s%glb_sz(2)
            write(*,105) k, s%glb_sz(3)
            stop
          endif
          allocate(real_imag_pairs(2*s%glb_sz(2)))
          do k = 1, s%glb_sz(3)
            do i = 1, s%glb_sz(1)
              read(unit=fo, fmt='('//trim(adjustl(ntoa(2*s%glb_sz(2))))//bsstr) real_imag_pairs
              do j = 1, s%glb_sz(2)
                a(i,j,k) = real_imag_pairs(2*j-1) + iu * real_imag_pairs(2*j)
              enddo
            enddo
            if (k .lt. size(a,3)) read(fo, *)
          enddo
          deallocate(real_imag_pairs)
          read(fo, *) s%ln, s%nrchop_offset, s%npchop_offset, s%nzchop_offset
          read(fo, *) s%space
        endif
        close(fo)
      endif
      call MPI_barrier(comm_glb, MPI_err)

      if (s%axis_comm(3) .eq. 0) then
        call st%init(s%glb_sz, s%axis_comm); st = s
        call st%exchange(3, 1) ! because %assemble cannot work when dim3 resides locally
        call st%disassemble(a) ! distribute the global array into a local array at each proc
        call st%exchange(1, 3); s = st
        call st%dealloc()
      else
        call s%disassemble(a) ! distribute the global array into a local array at each proc
      endif
      deallocate( a )
    else
      if (is_binary_) then
        open(newunit=fo, file=trim(adjustl(fn))//'_'//trim(adjustl(ntoa(rank_glb))), &
             form='unformatted', access='stream', status='old', iostat=is, action='read')
        if (is .ne. 0) then
          write(*,*) 'mload_scalar: cannot open ', fn
          stop
        endif
        i = s%loc_sz(1); j = s%loc_sz(2); k = s%loc_sz(3)
        read(fo) s%glb_sz, s%loc_sz, s%loc_st
        if ((i .ne. s%loc_sz(1)) .or. (j .ne. s%loc_sz(2)) .or. (k .ne. s%loc_sz(3))) then
          write(*,*) 'mloadc: size inconsistency between data and array'
          write(*,103) s%loc_sz(1), i
          write(*,104) s%loc_sz(2), j
          write(*,105) s%loc_sz(3), k
          stop
        endif
        read(fo) s%e
        read(fo) s%ln, s%nrchop_offset, s%npchop_offset, s%nzchop_offset
        read(fo) s%space
      else
        open(newunit=fo, file=trim(adjustl(fn))//'_'//trim(adjustl(ntoa(rank_glb))), &
             status='old')
        i = s%loc_sz(1); j = s%loc_sz(2); k = s%loc_sz(3)
        read(unit=fo, fmt='(3(1X,I10))') s%glb_sz
        read(unit=fo, fmt='(3(1X,I10))') s%loc_sz
        read(unit=fo, fmt='(3(1X,I10))') s%loc_st 
        if ((i .ne. s%loc_sz(1)) .or. (j .ne. s%loc_sz(2)) .or. (k .ne. s%loc_sz(3))) then
          write(*,*) 'mloadc: size inconsistency between data and array'
          write(*,103) s%loc_sz(1), i
          write(*,104) s%loc_sz(2), j
          write(*,105) s%loc_sz(3), k
          stop
        endif
        allocate(real_imag_pairs(2*s%loc_sz(2)))
        do k = 1, s%loc_sz(3)
          do i = 1, s%loc_sz(1)
            read(unit=fo, fmt='('//trim(adjustl(ntoa(2*s%loc_sz(2))))//bsstr) real_imag_pairs
            do j = 1, s%loc_sz(2)
              s%e(i,j,k) = real_imag_pairs(2*j-1) + iu * real_imag_pairs(2*j)
            enddo
          enddo
          if (k .lt. size(a,3)) read(fo, *)
        enddo
        deallocate(real_imag_pairs)
        write(fo, *) s%ln, s%nrchop_offset, s%npchop_offset, s%nzchop_offset
        write(fo, *) s%space
      endif
      close(fo)
      call MPI_barrier(comm_glb, MPI_err)
    endif
  103 format('axis 1 (data): ', I6, '   ', 'axis 1 (array): ', I6)
  104 format('axis 2 (data): ', I6, '   ', 'axis 2 (array): ', I6)
  105 format('axis 3 (data): ', I6, '   ', 'axis 3 (array): ', I6)
  end procedure

end submodule