submodule (mlegs_scalar) mlegs_scalar_ops
  implicit none

contains

  module procedure chop
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    integer(i4) :: mm, kk

    select type(tfm)
      type is (tfm_kit_1d)
        nrdim = tfm%nrdim; nrc = tfm%chops(1); nrcs = tfm%chops
        npdim = 1        ; npc = 1;
        nzdim = 1        ; nzc = 1; nzcu = 1
        m = (/ 0 /)
      type is (tfm_kit_2d)
        nrdim = tfm%nrdim; nrc = tfm%chops(1); nrcs = tfm%chops
        npdim = tfm%npdim; npc = tfm%chopp
        nzdim = 1        ; nzc = 1; nzcu = 1
        m = tfm%m
      type is (tfm_kit_3d)
        nrdim = tfm%nrdim; nrc = tfm%chops(1); nrcs = tfm%chops
        npdim = tfm%npdim; npc = tfm%chopp
        nzdim = tfm%nzdim; nzc = tfm%chopzl; nzcu = tfm%chopzu
        m = tfm%m 
    end select

    if ( (is_warning) .and. &
         ((s%glb_sz(1) .ne. nrdim) .or. &
          (s%glb_sz(2) .ne. npdim) .or. &
          (s%glb_sz(3) .ne. nzdim)) ) then
      write(*,*) '[wanring] chop: scalar dimension incompatible with the dimension defined in a kit'
    endif

    nrc = nrc + s%nrchop_offset
    npc = npc + s%npchop_offset
    nzc = nzc + s%nzchop_offset
    nzcu = nzcu - s%nzchop_offset

    if (nrc .gt. nrdim) stop 'chop: chopping in r too large'
    if (npc .gt. npdim) stop 'chop: chopping in p too large'
    if (2*s%nzchop_offset .gt. nzcu - nzc) stop 'chop: chopping in z too large'

    do mm = 1, npc
      nrcs(mm) = max(min(nrc, nrc-m(mm)), 0)
    enddo
    nrcs(npc+1:) = 0

    if (s%space(1:1) .eq. 'F') then
      do mm = 1, s%loc_sz(2)
        s%e(nrcs(s%loc_st(2)+mm)+1:, mm, :) = 0.D0
      enddo
    endif

    if (s%space(2:2) .eq. 'F') then
      do mm = 1, s%loc_sz(2)
        if (s%loc_st(2)+mm .ge. npc+1) s%e(:, mm, :) = 0.D0
      enddo
    endif

    if (s%space(3:3) .eq. 'F') then
      if (nzc .lt. nzcu) then
        do kk = 1, s%loc_sz(3)
          if ((s%loc_st(3)+kk .ge. nzc+1) .and. (s%loc_st(3)+kk .le. nzcu-1)) s%e(:, :, kk) = 0.D0
        enddo
      endif
    endif

    return
  end procedure

  module procedure trans
    integer(i4) :: curr_sp, new_sp
    integer(i4) :: kk

    if (s%space(1:3) .eq. 'PPP') then
      curr_sp = 0
    elseif (s%space(1:3) .eq. 'PFP') then
      curr_sp = 1
    elseif (s%space(1:3) .eq. 'FFP') then
      curr_sp = 2
    elseif (s%space(1:3) .eq. 'FFF') then
      curr_sp = 3
    else
      stop 'trans: scalar space info corrupted (only accepting PPP, PFP, FFP and FFF)'
    endif

    if (space(1:3) .eq. 'PPP') then
      new_sp = 0
    elseif (space(1:3) .eq. 'PFP') then
      new_sp = 1
    elseif (space(1:3) .eq. 'FFP') then
      new_sp = 2
    elseif (space(1:3) .eq. 'FFF') then
      new_sp = 3
    else
      stop 'trans: only taking PPP, PFP, FFP and FFF for spectral transformation'
    endif

    if (curr_sp .lt. new_sp) then ! forward transforamtion
      do while(.true.)
      if (curr_sp .eq. 0) then
        call horizontal_fft_forward(s, tfm)
        s%space = 'PFP'
        call s%exchange(2,1)
      endif
      if (curr_sp .eq. 1) then
        do kk = 1, nz
          s%e(:nr, 1, kk) = s%e(:nr, 1, kk) - s%ln*tfm%ln(:nr)
        enddo
        call rtrans_forward(s, tfm)
        s%space = 'FFP'
        call s%exchange(1,3)
      endif
      if (curr_sp .eq. 2) then
        call vertical_fft_forward(s, tfm)
        s%space = 'FFF'
      endif
      if (curr_sp .eq. 3) exit
      curr_sp = curr_sp + 1
      if (curr_sp .eq. new_sp) exit
      enddo
    endif

    if (curr_sp .gt. new_sp) then ! backward transforamtion
      do while(.true.)
      if (curr_sp .eq. 3) then
        call vertical_fft_backward(s, tfm)
        s%space = 'FFP'
      endif
      if (curr_sp .eq. 2) then
        call s%exchange(3,1)
        call rtrans_backward(s, tfm)
        do kk = 1, nz
          s%e(:nr, 1, kk) = s%e(:nr, 1, kk) + s%ln*tfm%ln(:nr)
        enddo
        s%space = 'PFP'
      endif
      if (curr_sp .eq. 1) then
        call s%exchange(1,2)
        call horizontal_fft_backward(s, tfm)
        s%space = 'PPP'
      endif
      if (curr_sp .eq. 0) exit
      curr_sp = curr_sp - 1
      if (curr_sp .eq. new_sp) exit
      enddo
    endif

  end procedure

! ======================================================================================================== !
! VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV INTERNAL (PRIVATE) SUBROUTINES/FUNCTIONS VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV !
! ======================================================================================================== !

  subroutine horizontal_fft_forward(s, tfm)
    implicit none
    class(scalar), intent(inout) :: s
    class(tfm_kit), intent(in) :: tfm

    complex(p8), dimension(:,:,:), allocatable :: se
    integer(i4) :: npdim, npc 

    complex(p8), dimension(:,:), allocatable :: b
    real(p8), dimension(:,:), allocatable :: c
    integer(i4) :: ii, kk, mm, nph

    external :: dzfft2d

    select type(tfm)
      type is (tfm_kit_1d)
        npdim = 1        ; npc = 1
      type is (tfm_kit_2d)
        npdim = tfm%npdim; npc = tfm%chopp
      type is (tfm_kit_3d)
        npdim = tfm%npdim; npc = tfm%chopp
    end select

    if ((is_warning) .and. (s%npchop_offset .ne. 0)) then
      write(*,*) '[warning] horizontal_fft_forward: chopping offset exists in p (!= global npchop)'
    endif

    npc = npc + s%npchop_offset
    if (npc .gt. npdim) stop 'horizontal_fft_forward: chopping in p too large'

    se = s%e
    nph = np/2

    allocate( b(npdim, 1) ) ! npdim = np/2+1 = nph + 1
    allocate( c(2*npdim, 1) )

    b(:,1) = se(1, :npdim, 1)
    c = 0.D0
    do ii = 1, nph
      c(2*ii - 1, 1) = real(b(ii, 1))
      c(2*ii, 1)     = aimag(b(ii, 1))
    enddo
    call dzfft2d(c, 2*nph, 1, 0, b)

    do ii = 1, size(se, 1)
      do kk = 1, size(se, 3)
        b(:, 1) = se(ii, :npdim, kk)
        c = 0.D0
        do mm = 1, nph
          c(2*mm - 1, 1) = real(b(mm, 1))
          c(2*mm, 1)     = aimag(b(mm, 1))
        enddo
        call dzfft2d(c, 2*nph, 1, -1, b)
        do mm = 1, nph + 1
          b(mm, 1) = c(2*mm-1, 1) + iu * c(2*mm, 1)
        enddo
        se(ii, :npdim, kk) = b(:npdim, 1)/(2*nph)
      enddo
    enddo

    s%e = se

    if (.not. all(abs(se(:, npc+1, :)) .lt. 1.d-14)) then
      if (is_warning) write(*,*) '[warning] horizontal_fft_forward: residues in the chopped region'
    endif

    deallocate( b )
    deallocate( c )

  end subroutine

! ======================================================================================================== !

  subroutine horizontal_fft_backward(s, tfm)
    implicit none
    class(scalar), intent(inout) :: s
    class(tfm_kit), intent(in) :: tfm

    complex(p8), dimension(:,:,:), allocatable :: se
    integer(i4) :: npdim, npc 

    complex(p8), dimension(:,:), allocatable :: b
    complex(p8), dimension(:,:), allocatable :: c
    integer(i4) :: ii, kk, mm, nph

    external :: zdfft2d

    select type(tfm)
      type is (tfm_kit_1d)
        npdim = 1        ; npc = 1
      type is (tfm_kit_2d)
        npdim = tfm%npdim; npc = tfm%chopp
      type is (tfm_kit_3d)
        npdim = tfm%npdim; npc = tfm%chopp
    end select

    if ((is_warning) .and. (s%npchop_offset .ne. 0)) then
      write(*,*) '[warning] horizontal_fft_backward: chopping offset exists in p (!= global npchop)'
    endif

    npc = npc + s%npchop_offset
    if (npc .gt. npdim) stop 'horizontal_fft_fbackard: chopping in p too large'

    if ((is_warning) .and. (.not. all(abs(se(:, npc+1, :)) .lt. 1.D-14))) then
      write(*,*) '[warning] horizontal_fft_backward: residues in the chopped region'
    endif

    se = s%e
    nph = np/2

    allocate( b(npdim, 1) ) ! npdim = np/2+1 = nph + 1
    allocate( c(npdim, 1) )

    b(:,1) = se(1, :npdim, 1)
    c = 0.D0
    do ii = 1, nph + 1
      c(ii, 1) = b(ii, 1)
    enddo
    call zdfft2d(c, 2*nph, 1, 0, b)

    do ii = 1, size(se, 1)
      do kk = 1, size(se, 3)
        b(:, 1) = se(ii, :npdim, kk)
        c = 0.D0
        do mm = 1, nph + 1
          c(mm, 1) = b(mm, 1)
        enddo
        call zdfft2d(c, 2*nph, 1, 1, b)
        do mm = 1, nph + 1
          b(mm, 1) = c(mm, 1)
        enddo
        se(ii, :npdim, kk) = b(:npdim, 1)*(2*nph)
      enddo
    enddo

    s%e = se

    deallocate( b )
    deallocate( c )

  end subroutine

! ======================================================================================================== !

  subroutine vertical_fft_forward(s, tfm)
    implicit none
    class(scalar), intent(inout) :: s
    class(tfm_kit), intent(in) :: tfm

    complex(p8), dimension(:,:,:), allocatable :: se
    integer(i4) :: nzdim, nzc, nzcu 

    complex(p8), dimension(:), allocatable :: b
    complex(p8), dimension(:), allocatable :: c
    integer(i4) :: ii, jj, mm

    external :: zfft1d

    select type(tfm)
      type is (tfm_kit_1d)
        nzdim = 1        ; nzc = 1; nzcu = 1
      type is (tfm_kit_2d)
        nzdim = 1        ; nzc = 1; nzcu = 1
      type is (tfm_kit_3d)
        nzdim = tfm%nzdim; nzc = tfm%chopzl; nzcu = tfm%chopzu
    end select

    if ((is_warning) .and. (s%nzchop_offset .ne. 0)) then
      write(*,*) '[warning] vertical_fft_forward: chopping offset exists in z (!= global nzchop)'
    endif

    nzc = nzc + s%nzchop_offset
    nzcu = nzcu - s%nzchop_offset
    if (2*s%nzchop_offset .gt. nzcu - nzc) stop 'vertical_fft_forward: chopping in z too large'

    se = s%e

    allocate( b(nz) ) ! npdim = np/2+1 = nph + 1
    allocate( c(2*nz) )

    b(:) = se(1, 1, :nzdim)
    c = 0.D0
    call zfft1d(b, nz, 0, c)

    do ii = 1, size(se, 1)
      do jj = 1, size(se, 2)
        b(:) = se(ii, jj, :nzdim)
        call zfft1d(b, nz, -1, c)
        se(ii, jj, :nzdim) = b(:nzdim)/nz
      enddo
    enddo

    s%e = se

    if (.not. all(abs(se(:, :, nzc+1:nzcu-1)) .lt. 1.D-14)) then
      if (is_warning) write(*,*) '[warning] vertical_fft_forward: residues in the chopped region'
    endif

    deallocate( b )
    deallocate( c )

  end subroutine

! ======================================================================================================== !

  subroutine vertical_fft_backward(s, tfm)
    implicit none
    class(scalar), intent(inout) :: s
    class(tfm_kit), intent(in) :: tfm

    complex(p8), dimension(:,:,:), allocatable :: se
    integer(i4) :: nzdim, nzc, nzcu 

    complex(p8), dimension(:), allocatable :: b
    complex(p8), dimension(:), allocatable :: c
    integer(i4) :: ii, jj, mm

    external :: zfft1d

    select type(tfm)
      type is (tfm_kit_1d)
        nzdim = 1        ; nzc = 1; nzcu = 1
      type is (tfm_kit_2d)
        nzdim = 1        ; nzc = 1; nzcu = 1
      type is (tfm_kit_3d)
        nzdim = tfm%nzdim; nzc = tfm%chopzl; nzcu = tfm%chopzu
    end select

    if ((is_warning) .and. (s%nzchop_offset .ne. 0)) then
      write(*,*) '[warning] vertical_fft_backward: chopping offset exists in z (!= global nzchop)'
    endif

    nzc = nzc + s%nzchop_offset
    nzcu = nzcu - s%nzchop_offset
    if (2*s%nzchop_offset .gt. nzcu - nzc) stop 'vertical_fft_backward: chopping in z too large'

    if (.not. all(abs(se(:, :, nzc+1:nzcu-1)) .lt. 1.D-14)) then
      if (is_warning) write(*,*) '[warning] vertical_fft_backward: residues in the chopped region'
    endif

    se = s%e

    allocate( b(nz) ) ! npdim = np/2+1 = nph + 1
    allocate( c(2*nz) )

    b(:) = se(1, 1, :nzdim)
    c = 0.D0
    call zfft1d(b, nz, 0, c)

    do ii = 1, size(se, 1)
      do jj = 1, size(se, 2)
        b(:) = se(ii, jj, :nzdim)
        call zfft1d(b, nz, 1, c)
        se(ii, jj, :nzdim) = b(:nzdim)*nz
      enddo
    enddo

    s%e = se

    deallocate( b )
    deallocate( c )

  end subroutine

! ======================================================================================================== !

  subroutine rtrans_forward(s, tfm)
    implicit none
    class(scalar), intent(inout) :: s
    class(tfm_kit), intent(in) :: tfm

    complex(p8), dimension(:,:,:), allocatable :: se
    complex(p8), dimension(:,:), allocatable :: be, bo
    integer(i4), dimension(:), allocatable :: nrcs, m
    integer(i4) :: nrh, nrdim, npdim, nzdim, nrc, npc, i, nn, mm

    select type(tfm)
      type is (tfm_kit_1d)
        nrdim = tfm%nrdim; nrc = tfm%chops(1); nrcs = tfm%chops
        npdim = 1        ; npc = 1
        nzdim = 1
        m = (/ 0 /)
      type is (tfm_kit_2d)
        nrdim = tfm%nrdim; nrc = tfm%chops(1); nrcs = tfm%chops
        npdim = tfm%npdim; npc = tfm%chopp
        nzdim = 1
        m = tfm%m
      type is (tfm_kit_3d)
        nrdim = tfm%nrdim; nrc = tfm%chops(1); nrcs = tfm%chops
        npdim = tfm%npdim; npc = tfm%chopp
        nzdim = tfm%nzdim
        m = tfm%m
    end select

    if ((is_warning) .and. (s%nrchop_offset .ne. 0)) then
      write(*,*) '[warning] rtrans_forward: chopping offset exists in r (!= global nrchop)'
    endif
    if ((is_warning) .and. (s%npchop_offset .ne. 0)) then
      write(*,*) '[warning] rtrans_forward: chopping offset exists in p (!= global npchop)'
    endif

    nrc = nrc + s%nrchop_offset
    npc = npc + s%npchop_offset
    if (nrc .gt. nrdim) stop 'rtrans_forward: chopping in r too large'
    if (npc .gt. npdim) stop 'rtrans_forward: chopping in p too large'
    do mm = 1, npc
      nrcs(mm) = max(min(nrc, nrc-m(mm)), 0)
    enddo
    nrcs(npc+1:) = 0

    allocate( se(size(s%e,1), size(s%e,2), size(s%e,3)) )
    se = 0.D0

    nrh = nr/2
    allocate( be(nrh, nzdim) )
    allocate( bo(nrh, nzdim) )

    do mm = 1, npc
      be = 0.D0; bo = 0.D0
      do i = 1, nrh
        be(i,:) = (s%e(i,mm,:)+s%e(nr-i+1,mm,:))*tfm%w(i)
        bo(i,:) = (s%e(i,mm,:)-s%e(nr-i+1,mm,:))*tfm%w(i)
      enddo
      nn = nrcs(mm)
      if (nn .ge. 1) then
        se(1:nn:2,mm,:) = transpose(tfm%pf(:nrh,1:nn:2,mm)) .mul. be
      endif
      if (nn .ge. 2) then
        se(2:nn:2,mm,:) = transpose(tfm%pf(:nrh,2:nn:2,mm)) .mul. bo
      endif
    enddo

    s%e = se

    deallocate( be )
    deallocate( bo )
    deallocate( se )

    where (abs(real(s%e)) .lt. 1.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    where (abs(aimag(s%e)) .lt. 1.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

  end subroutine

! ======================================================================================================== !

  subroutine rtrans_backward(s, tfm)
    implicit none
    class(scalar), intent(inout) :: s
    class(tfm_kit), intent(in) :: tfm

    complex(p8), dimension(:,:,:), allocatable :: se
    complex(p8), dimension(:,:), allocatable :: be, bo
    integer(i4), dimension(:), allocatable :: nrcs, m
    integer(i4) :: nrh, nrdim, npdim, nzdim, nrc, npc, i, nn, mm

    select type(tfm)
      type is (tfm_kit_1d)
        nrdim = tfm%nrdim; nrc = tfm%chops(1); nrcs = tfm%chops
        npdim = 1        ; npc = 1
        nzdim = 1
        m = (/ 0 /)
      type is (tfm_kit_2d)
        nrdim = tfm%nrdim; nrc = tfm%chops(1); nrcs = tfm%chops
        npdim = tfm%npdim; npc = tfm%chopp
        nzdim = 1
        m = tfm%m
      type is (tfm_kit_3d)
        nrdim = tfm%nrdim; nrc = tfm%chops(1); nrcs = tfm%chops
        npdim = tfm%npdim; npc = tfm%chopp
        nzdim = tfm%nzdim
        m = tfm%m
    end select

    if ((is_warning) .and. (s%nrchop_offset .ne. 0)) then
      write(*,*) '[warning] rtrans_backward: chopping offset exists in r (!= global nrchop)'
    endif
    if ((is_warning) .and. (s%npchop_offset .ne. 0)) then
      write(*,*) '[warning] rtrans_backward: chopping offset exists in p (!= global npchop)'
    endif

    nrc = nrc + s%nrchop_offset
    npc = npc + s%npchop_offset
    if (nrc .gt. nrdim) stop 'rtrans_backward: chopping in r too large'
    if (npc .gt. npdim) stop 'rtrans_backward: chopping in p too large'
    do mm = 1, npc
      nrcs(mm) = max(min(nrc, nrc-m(mm)), 0)
    enddo
    nrcs(npc+1:) = 0

    allocate( se(size(s%e,1), size(s%e,2), size(s%e,3)) )
    se = 0.D0

    nrh = nr/2
    allocate( be(nrh, nzdim) )
    allocate( bo(nrh, nzdim) )

    do mm = 1, npc
      be = 0.D0; bo = 0.D0
      nn = nrcs(mm)
      if (nn .ge. 1) then
        be = tfm%pf(:nrh,1:nn:2,mm) .mul. s%e(1:nn:2,mm,:)
      else
        be = 0.D0
      endif
      if (nn .ge. 2) then
        bo = tfm%pf(:nrh,2:nn:2,mm) .mul. s%e(2:nn:2,mm,:)
      else
        bo = 0.D0
      endif
      se(1:nrh,mm,:) = be+bo
      se(nr:nrh+1:-1,mm,:) = be-bo
    enddo

    s%e = se

    deallocate( be )
    deallocate( bo )
    deallocate( se )

    where (abs(real(s%e)) .lt. 1.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    where (abs(aimag(s%e)) .lt. 1.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

  end subroutine

end submodule