submodule (mlegs_scalar) mlegs_scalar_ops
  implicit none

contains

  module procedure chop
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, kk

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

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
        if (np .gt. 1) call horizontal_fft_forward(s, tfm)
        if (np .gt. 1) call s%exchange(2,1)
        s%space = 'PFP'
      endif
      if (curr_sp .eq. 1) then
        do kk = 1, s%loc_sz(3)
          s%e(:nr, 1, kk) = s%e(:nr, 1, kk) - s%ln*tfm%ln(:nr)
        enddo
        call rtrans_forward(s, tfm)
        s%space = 'FFP'
        if (nz .gt. 1) call s%exchange(1,3)
      endif
      if (curr_sp .eq. 2) then
        if (nz .gt. 1) call vertical_fft_forward(s, tfm)
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
        if (nz .gt. 1) call vertical_fft_backward(s, tfm)
        s%space = 'FFP'
      endif
      if (curr_sp .eq. 2) then
        if (nz .gt. 1) call s%exchange(3,1)
        call rtrans_backward(s, tfm)
        do kk = 1, s%loc_sz(3)
          s%e(:nr, 1, kk) = s%e(:nr, 1, kk) + s%ln*tfm%ln(:nr)
        enddo
        s%space = 'PFP'
      endif
      if (curr_sp .eq. 1) then
        if (np .gt. 1) call s%exchange(1,2)
        if (np .gt. 1) call horizontal_fft_backward(s, tfm)
        s%space = 'PPP'
      endif
      if (curr_sp .eq. 0) exit
      curr_sp = curr_sp - 1
      if (curr_sp .eq. new_sp) exit
      enddo
    endif

  end procedure

  module procedure calcat0
    type(scalar) :: st
    complex(p8), dimension(:), allocatable :: calc_

    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    call st%init(s%glb_sz, s%axis_comm)
    st = s

    if ((.not.(s%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] calcat0: an input scalar is not FFF; additional operations needed for tfm'
    endif
    call trans(st, 'FFF', tfm)

    allocate( calc(nz) )
    if (st%loc_st(2) .eq. 0) then
      calc = transpose(st%e(1:min(st%loc_sz(1),nrc-st%loc_st(1)), 1, :)) .mul. &
                       tfm%at0(st%loc_st(1)+1:min(st%loc_st(1)+st%loc_sz(1),nrc))
    else
      calc = 0.D0
    endif

    allocate( calc_(nz) )
    call MPI_allreduce(calc, calc_, nz, MPI_complex16, MPI_sum, comm_glb, MPI_err)
    calc = calc_

    call st%dealloc()
    deallocate( calc_ )

    return
  end procedure

  module procedure calcat1
    type(scalar) :: st
    complex(p8), dimension(:), allocatable :: calc_

    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    call st%init(s%glb_sz, s%axis_comm)
    st = s

    if ((.not.(s%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] calcat1: an input scalar is not FFF; additional operations needed for tfm'
    endif
    call trans(st, 'FFF', tfm)

    allocate( calc(nz) )
    if (st%loc_st(2) .eq. 0) then
      calc = transpose(st%e(1:min(st%loc_sz(1),nrc-st%loc_st(1)), 1, :)) .mul. &
                       tfm%at1(st%loc_st(1)+1:min(st%loc_st(1)+st%loc_sz(1),nrc))
    else
      calc = 0.D0
    endif

    allocate( calc_(nz) )
    call MPI_allreduce(calc, calc_, nz, MPI_complex16, MPI_sum, comm_glb, MPI_err)
    calc = calc_

    call st%dealloc()
    deallocate( calc_ )

    return
  end procedure

  module procedure zeroat1
    complex(p8), dimension(:), allocatable :: at1

    if ((.not.(s%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] zeroat1: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                   after this operation, the scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    at1 = calcat1(s, tfm)
    if((s%loc_st(1) .eq. 0) .and. (s%loc_st(2) .eq. 0)) then
      s%e(1,1,:) = s%e(1,1,:) - at1/tfm%at1(1)
    endif

  end procedure

  module procedure delsqp
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, nn, n

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] delsqp: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                  after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s

    so%ln = 0.D0

    if (m(1) .ne. 0) stop 'delsqp: 1st azimuthal wavenum is not zero'

    do mm = 1, min(so%loc_sz(2), npc - so%loc_st(2))
      do nn = 1, min(so%loc_sz(1), nrcs(so%loc_st(2) + mm) - so%loc_sz(1))
        n = m(so%loc_st(2) + mm) + (so%loc_st(1) + nn) - 1
        so%e(nn,mm,:) = -s%e(nn,mm,:)*n*(n+1.D0)/(ell**2.D0)
      enddo
    enddo

    if ((so%loc_st(1) .eq. 0) .and. (so%loc_st(2) .eq. 0)) then
      so%e(1,1,1) = s%ln/ell**2.D0/exp(tfm%lognorm(1,1))
    endif

    return
  end procedure

  module procedure idelsqp
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, nn, n

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] idelsqp: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                   after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s

    if (m(1) .ne. 0) stop 'idelsqp: 1st azimuthal wavenum is not zero'

    if ((so%loc_st(1) .eq. 0) .and. (so%loc_st(2) .eq. 0)) then
      so%ln = s%e(1,1,1)*ell**2.D0*exp(tfm%lognorm(1,1))
    endif
    call MPI_allreduce(so%ln, so%ln, 1, MPI_real8, MPI_sum, comm_glb, MPI_err)

    do mm = 1, min(so%loc_sz(2), npc - so%loc_st(2))
      do nn = 1, min(so%loc_sz(1), nrcs(so%loc_st(2) + mm) - so%loc_st(1))
        n = m(so%loc_st(2) + mm) + (so%loc_st(1) + nn) - 1
        so%e(nn,mm,:) = -s%e(nn,mm,:)/n/(n+1.D0)*(ell**2.D0)
      enddo
    enddo

    if ((so%loc_st(1) .eq. 0) .and. (so%loc_st(2) .eq. 0)) then
      so%e(1,1,:) = 0.D0
    endif

    return
  end procedure

  module procedure xxdx
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, nn, n
    type(real_bndm) :: xxdx_bnd

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] xxdx: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s
    if (nz .gt. 1) call so%exchange(3,1) ! whole r-data resides locally

    so%ln = 0.D0

    do mm = 1, min(so%loc_sz(2), npc - so%loc_st(2))
      nn = nrcs(so%loc_st(2) + mm)
      xxdx_bnd = leg_xxdx(m(so%loc_st(2) + mm), nn, nn, tfm)
      n = min(so%loc_sz(3), nzc - so%loc_st(3))
      so%e(:nn, mm, :n) = xxdx_bnd .mul. so%e(:nn, mm, :n)
      n = max(nzcu - so%loc_st(3), 1)
      so%e(:nn, mm, n:) = xxdx_bnd .mul. so%e(:nn, mm, n:)
    enddo
    deallocate( xxdx_bnd%e )

    !> r*d(P_L(r))/dr == P_L_1^0(r) + P_L_0^0(r)
    if ((so%loc_st(2).eq.0) .and. (so%loc_st(3).eq.0)) then
      so%e(1,1,1) = so%e(1,1,1) + 1.D0/exp(tfm%lognorm(1,1))*s%ln
      so%e(2,1,1) = so%e(2,1,1) + 1.D0/exp(tfm%lognorm(2,1))*s%ln
    endif
    if (nz .gt. 1) call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

    return
  end procedure

  module procedure del2h
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, nn, kk, n
    type(real_bndm) :: del2h_bnd

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] del2: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s
    if (nz .gt. 1) call so%exchange(3,1) ! whole r-data resides locally

    so%ln = 0.D0

    do mm = 1, min(so%loc_sz(2), npc - so%loc_st(2))
      nn = nrcs(so%loc_st(2) + mm)
      if (1 .le. min(so%loc_sz(3), nzc - so%loc_st(3))) then
        do kk = 1, min(so%loc_sz(3), nzc - so%loc_st(3))
          del2h_bnd = leg_del2h(m(so%loc_st(2) + mm), nn, nn, tfm)
          so%e(:nn, mm, kk) = del2h_bnd .mul. so%e(:nn, mm, kk)
        enddo
      endif
      if (max(nzcu - so%loc_st(3), 1) .le. so%loc_sz(3)) then
        do kk = max(nzcu - so%loc_st(3), 1), so%loc_sz(3)
          del2h_bnd = leg_del2h(m(so%loc_st(2) + mm), nn, nn, tfm)
          so%e(:nn, mm, kk) = del2h_bnd .mul. so%e(:nn, mm, kk)
        enddo
      endif
    enddo
    deallocate( del2h_bnd%e )

    ! 1/r*d(r*d(P_L(r))/dr)/dr = 4/3*P_L_0^0(r) - 2*P_L_1^0(r) + 2/3*P_L_2^0(r)
    if ((so%loc_st(2).eq.0) .and. (so%loc_st(3).eq.0)) then
      so%e(1,1,1) = so%e(1,1,1) + 4.D0/3.D0/ell**2.D0/exp(tfm%lognorm(1,1))*s%ln
      so%e(2,1,1) = so%e(2,1,1) - 2.D0/1.D0/ell**2.D0/exp(tfm%lognorm(2,1))*s%ln
      so%e(3,1,1) = so%e(3,1,1) + 2.D0/3.D0/ell**2.D0/exp(tfm%lognorm(3,1))*s%ln
    endif
    if (nz .gt. 1) call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

    return
  end procedure

  module procedure del2
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, nn, kk, n
    type(real_bndm) :: del2_bnd

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] del2: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s
    if (nz .gt. 1) call so%exchange(3,1) ! whole r-data resides locally

    so%ln = 0.D0

    do mm = 1, min(so%loc_sz(2), npc - so%loc_st(2))
      nn = nrcs(so%loc_st(2) + mm)
      if (1 .le. min(so%loc_sz(3), nzc - so%loc_st(3))) then
        do kk = 1, min(so%loc_sz(3), nzc - so%loc_st(3))
          del2_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
          so%e(:nn, mm, kk) = del2_bnd .mul. so%e(:nn, mm, kk)
        enddo
      endif
      if (max(nzcu - so%loc_st(3), 1) .le. so%loc_sz(3)) then
        do kk = max(nzcu - so%loc_st(3), 1), so%loc_sz(3)
          del2_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
          so%e(:nn, mm, kk) = del2_bnd .mul. so%e(:nn, mm, kk)
        enddo
      endif
    enddo
    deallocate( del2_bnd%e )

    ! 1/r*d(r*d(P_L(r))/dr)/dr = 4/3*P_L_0^0(r) - 2*P_L_1^0(r) + 2/3*P_L_2^0(r)
    if ((so%loc_st(2).eq.0) .and. (so%loc_st(3).eq.0)) then
      so%e(1,1,1) = so%e(1,1,1) + 4.D0/3.D0/ell**2.D0/exp(tfm%lognorm(1,1))*s%ln
      so%e(2,1,1) = so%e(2,1,1) - 2.D0/1.D0/ell**2.D0/exp(tfm%lognorm(2,1))*s%ln
      so%e(3,1,1) = so%e(3,1,1) + 2.D0/3.D0/ell**2.D0/exp(tfm%lognorm(3,1))*s%ln
    endif
    if (nz .gt. 1) call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

    return
  end procedure

  module procedure idel2_preln
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, nn, kk, n
    type(real_bndm) :: del2_bnd
    real(p8), dimension(:,:), allocatable :: del2_pre

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] idel2: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                 after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s
    if (nz .gt. 1) call so%exchange(3,1) ! whole r-data resides locally

    if ((is_warning) .and. (abs(so%ln) .gt. 5.D-14)) then
      write(*,*) '[warning]: idel2: an input scalar has non-zero ln. operation may be inaccurate.'
    endif

    if ((so%loc_st(2).eq.0) .and. (so%loc_st(3).eq.0)) then
      do mm = 1, min(so%loc_sz(2), npc - so%loc_st(2))
        nn = nrcs(so%loc_st(2) + mm)
        if (1 .le. min(so%loc_sz(3), nzc - so%loc_st(3))) then
          do kk = 1, min(so%loc_sz(3), nzc - so%loc_st(3))
            if ((mm .eq. 1) .and. (kk .eq. 1)) then
              del2_pre = fullmat(leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm))
              del2_pre(2:nn, :) = del2_pre(1:nn-1, :)
              del2_pre(1,:) = 0.D0

              del2_pre(1,1) = 1.D0
              del2_pre(2,1) = del2_pre(2,1) + 4.D0/3.D0/ell**2.D0
              del2_pre(3,1) = del2_pre(3,1) - 2.D0/1.D0/ell**2.D0*exp(tfm%lognorm(1,1)-tfm%lognorm(2,1))
              del2_pre(4,1) = del2_pre(4,1) + 2.D0/3.D0/ell**2.D0*exp(tfm%lognorm(1,1)-tfm%lognorm(3,1))

              del2_bnd = bandmat(del2_pre, suplen = 2, sublen = 3)
              deallocate(del2_pre)

              so%e(2:nn,mm,kk) = so%e(1:nn-1,mm,kk)
              so%e(1,mm,kk) = preln/exp(tfm%lognorm(1,1))
              call lsolve(del2_bnd, so%e(:nn,mm,kk))

              so%ln = so%e(1,mm,kk)*exp(tfm%lognorm(1,1))
            else
              del2_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
              call lsolve(del2_bnd, so%e(:nn,mm,kk))
            endif
          enddo
        endif
        if (max(nzcu - so%loc_st(3), 1) .le. so%loc_sz(3)) then
          do kk = max(nzcu - so%loc_st(3), 1), so%loc_sz(3)
            del2_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
            call lsolve(del2_bnd, so%e(:nn,mm,kk))
          enddo
        endif
      enddo
    else
      do mm = 1, min(so%loc_sz(2), npc - so%loc_st(2))
        nn = nrcs(so%loc_st(2) + mm)
        if (1 .le. min(so%loc_sz(3), nzc - so%loc_st(3))) then
          do kk = 1, min(so%loc_sz(3), nzc - so%loc_st(3))
            del2_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
            call lsolve(del2_bnd, so%e(:nn,mm,kk))
          enddo
        endif
        if (max(nzcu - so%loc_st(3), 1) .le. so%loc_sz(3)) then
          do kk = max(nzcu - so%loc_st(3), 1), so%loc_sz(3)
            del2_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
            call lsolve(del2_bnd, so%e(:nn,mm,kk))
          enddo
        endif
      enddo
    endif
    deallocate( del2_bnd%e )
    call MPI_allreduce(so%ln, so%ln, 1, MPI_real8, MPI_sum, comm_glb, MPI_err)

    if (nz .gt. 1) call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

    ! where (abs(real(so%e)) .lt. 5.D-14) so%e = cmplx(0.D0, aimag(so%e), p8)
    ! where (abs(aimag(so%e)) .lt. 5.D-14) so%e = cmplx(real(so%e), 0.D0, p8)

    return
  end procedure

  module procedure idel2_proln
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, nn, kk, n
    type(real_bndm) :: del2_bnd
    real(p8), dimension(:,:), allocatable :: del2_pre

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] idel2: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                 after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s
    if (nz .gt. 1) call so%exchange(3,1) ! whole r-data resides locally

    if ((is_warning) .and. (abs(so%ln) .gt. 5.D-14)) then
      write(*,*) '[warning]: idel2: an input scalar has non-zero ln. operation may be inaccurate.'
    endif

    if ((so%loc_st(2).eq.0) .and. (so%loc_st(3).eq.0)) then
      do mm = 1, min(so%loc_sz(2), npc - so%loc_st(2))
        nn = nrcs(so%loc_st(2) + mm)
        if (1 .le. min(so%loc_sz(3), nzc - so%loc_st(3))) then
          do kk = 1, min(so%loc_sz(3), nzc - so%loc_st(3))
            if ((mm .eq. 1) .and. (kk .eq. 1)) then
              del2_pre = fullmat(leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm))

              del2_pre(1,1) = del2_pre(1,1) + 4.D0/3.D0/ell**2.D0
              del2_pre(2,1) = del2_pre(2,1) - 2.D0/1.D0/ell**2.D0*exp(tfm%lognorm(1,1)-tfm%lognorm(2,1))
              del2_pre(3,1) = del2_pre(3,1) + 2.D0/3.D0/ell**2.D0*exp(tfm%lognorm(1,1)-tfm%lognorm(3,1))

              del2_bnd = bandmat(del2_pre, suplen = 2, sublen = 2)
              deallocate(del2_pre)

              call lsolve(del2_bnd, so%e(:nn,mm,kk))

              so%ln = so%e(1,mm,kk)*exp(tfm%lognorm(1,1))
            else
              del2_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
              call lsolve(del2_bnd, so%e(:nn,mm,kk))
            endif
          enddo
        endif
        if (max(nzcu - so%loc_st(3), 1) .le. so%loc_sz(3)) then
          do kk = max(nzcu - so%loc_st(3), 1), so%loc_sz(3)
            del2_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
            call lsolve(del2_bnd, so%e(:nn,mm,kk))
          enddo
        endif
      enddo
    else
      do mm = 1, min(so%loc_sz(2), npc - so%loc_st(2))
        nn = nrcs(so%loc_st(2) + mm)
        if (1 .le. min(so%loc_sz(3), nzc - so%loc_st(3))) then
          do kk = 1, min(so%loc_sz(3), nzc - so%loc_st(3))
            del2_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
            call lsolve(del2_bnd, so%e(:nn,mm,kk))
          enddo
        endif
        if (max(nzcu - so%loc_st(3), 1) .le. so%loc_sz(3)) then
          do kk = max(nzcu - so%loc_st(3), 1), so%loc_sz(3)
            del2_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
            call lsolve(del2_bnd, so%e(:nn,mm,kk))
          enddo
        endif
      enddo
    endif
    deallocate( del2_bnd%e )
    call MPI_allreduce(so%ln, so%ln, 1, MPI_real8, MPI_sum, comm_glb, MPI_err)

    if (nz .gt. 1) call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

    ! where (abs(real(so%e)) .lt. 5.D-14) so%e = cmplx(0.D0, aimag(so%e), p8)
    ! where (abs(aimag(so%e)) .lt. 5.D-14) so%e = cmplx(real(so%e), 0.D0, p8)

    return
  end procedure

  module procedure helm
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] helm: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    if ((abs(alpha) .lt. 5.D-14) .and. (is_warning)) then
      write(*,*) '[warning] helm: alpha equals to zero. This operation is exactly equivalent to del2'
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s
    so = del2(so, tfm)
    so%e = so%e + alpha*s%e
    so%ln = alpha*s%ln

    return
  end procedure

  module procedure ihelm
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, nn, kk, n
    type(real_bndm) :: helm_bnd

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] ihelm: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                 after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    if ((abs(alpha) .lt. 5.D-14)) then
      write(*,*) 'ihelm: alpha equals to zero. Inversion of 0*identity is impossible'
      stop
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s
    if (nz .gt. 1) call so%exchange(3,1) ! whole r-data resides locally

    so%ln = s%ln/alpha

    if ((so%loc_st(2).eq.0) .and. (so%loc_st(3).eq.0)) then
      so%e(1,1,1) = so%e(1,1,1) - 4.D0/3.D0/ell**2.D0/exp(tfm%lognorm(1,1))*so%ln
      so%e(2,1,1) = so%e(2,1,1) + 2.D0/1.D0/ell**2.D0/exp(tfm%lognorm(2,1))*so%ln
      so%e(3,1,1) = so%e(3,1,1) - 2.D0/3.D0/ell**2.D0/exp(tfm%lognorm(3,1))*so%ln
    endif

    do mm = 1, min(so%loc_sz(2), npc - so%loc_st(2))
      nn = nrcs(so%loc_st(2) + mm)
      if (1 .le. min(so%loc_sz(3), nzc - so%loc_st(3))) then
        do kk = 1, min(so%loc_sz(3), nzc - so%loc_st(3))
          helm_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
          helm_bnd%e(helm_bnd%suplen+1,:) = helm_bnd%e(helm_bnd%suplen+1,:) + alpha
          call lsolve(helm_bnd, so%e(:nn,mm,kk))
        enddo
      endif
      if (max(nzcu - so%loc_st(3), 1) .le. so%loc_sz(3)) then
        do kk = max(nzcu - so%loc_st(3), 1), so%loc_sz(3)
          helm_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
          helm_bnd%e(helm_bnd%suplen+1,:) = helm_bnd%e(helm_bnd%suplen+1,:) + alpha
          call lsolve(helm_bnd, so%e(:nn,mm,kk))
        enddo
      endif
    enddo
    deallocate( helm_bnd%e )

    if (nz .gt. 1) call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

    ! where (abs(real(so%e)) .lt. 5.D-14) so%e = cmplx(0.D0, aimag(so%e), p8)
    ! where (abs(aimag(so%e)) .lt. 5.D-14) so%e = cmplx(real(so%e), 0.D0, p8)

    return
  end procedure

  module procedure helmp
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: n

    type(scalar) :: s2, sp

    if (.not. (mod(power,2).eq.0) .and. (power.ge.4)) stop 'helmp: even power greater than or equal to 4'
    if (power .gt. 8) stop 'helmp: power must be less than or equal to 8 (supported power = 4, 6 or 8)'

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] helmp: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                 after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    if ((abs(alpha) .lt. 5.D-14) .and. (is_warning)) then
      write(*,*) '[warning] helmp: alpha equals to zero'
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s

    call s2%init(s%glb_sz, s%axis_comm)
    s2 = s
    s2 = del2(s2, tfm)

    call sp%init(s2%glb_sz, s2%axis_comm)
    sp = s2
    do n = 1, power/2-1
      sp = del2(sp, tfm)
    enddo

    so%e = sp%e + beta*s2%e + alpha*s%e
    so%ln = alpha*s%ln

    call s2%dealloc()
    call sp%dealloc()

    return
  end procedure

  module procedure ihelmp
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, nn, kk, n, np
    type(real_bndm) :: helmp_bnd, del2_bnd
    real(p8), dimension(:), allocatable :: bl, bl2
    real(p8), dimension(:,:), allocatable :: w

    if (.not. ((mod(power,2).eq.0) .and. (power.ge.4))) stop 'ihelmp: even power greater than or equal to 4'
    if (power .gt. 8) stop 'ihelmp: power must be less than or equal to 8 (supported power = 4, 6 or 8)'

    call chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] ihelmp: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                  after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    if ((abs(alpha) .lt. 5.D-14)) then
      write(*,*) 'ihelmp: alpha equals to zero. Inversion of 0*identity is impossible'
      stop
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s
    if (nz .gt. 1) call so%exchange(3,1) ! whole r-data resides locally

    so%ln = s%ln/alpha

    allocate( bl(nrc), bl2(nrc) )
    bl2 = 0.D0
    bl2(1) = 4.D0/3.D0/ell**2.D0/exp(tfm%lognorm(1,1))*so%ln
    bl2(2) = -2.D0    /ell**2.D0/exp(tfm%lognorm(2,1))*so%ln
    bl2(3) = 2.D0/3.D0/ell**2.D0/exp(tfm%lognorm(3,1))*so%ln
    bl = bl2

    if ((so%loc_st(2).eq.0) .and. (so%loc_st(3).eq.0)) then
      do n = 1, power/2-1
        bl = leg_del2(m(1), ak(1), nrc, nrc, tfm) .mul. bl
      enddo
      so%e(:nrc,1,1) = so%e(:nrc,1,1) - bl - beta*bl2
    endif

    do mm = 1, min(so%loc_sz(2), npc - so%loc_st(2))
      nn = nrcs(so%loc_st(2) + mm)
      np = nn + power
      if (1 .le. min(so%loc_sz(3), nzc - so%loc_st(3))) then
        do kk = 1, min(so%loc_sz(3), nzc - so%loc_st(3))
          helmp_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
          del2_bnd = helmp_bnd
          do n = 1, power/2-1
            helmp_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm) .mul. helmp_bnd
          enddo
          w = fullmat(helmp_bnd) + beta*fullmat(del2_bnd)
          helmp_bnd = bandmat(w(:nn,:nn), suplen=power, sublen=power)
          helmp_bnd%e(helmp_bnd%suplen+1,:) = helmp_bnd%e(helmp_bnd%suplen+1,:) + alpha
          deallocate( w )
          call lsolve(helmp_bnd, so%e(:nn,mm,kk))
        enddo
      endif
      if (max(nzcu - so%loc_st(3), 1) .le. so%loc_sz(3)) then
        do kk = max(nzcu - so%loc_st(3), 1), so%loc_sz(3)
          helmp_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm)
          del2_bnd = helmp_bnd
          do n = 1, power/2-1
            helmp_bnd = leg_del2(m(so%loc_st(2) + mm), ak(so%loc_st(3) + kk), nn, nn, tfm) .mul. helmp_bnd
          enddo
          w = fullmat(helmp_bnd) + beta*fullmat(del2_bnd)
          helmp_bnd = bandmat(w(:nn,:nn), suplen=power, sublen=power)
          helmp_bnd%e(helmp_bnd%suplen+1,:) = helmp_bnd%e(helmp_bnd%suplen+1,:) + alpha
          deallocate( w )
          call lsolve(helmp_bnd, so%e(:nn,mm,kk))
        enddo
      endif
    enddo
    call helmp_bnd%dealloc()
    call del2_bnd%dealloc()

    if (nz .gt. 1) call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

    ! where (abs(real(so%e)) .lt. 5.D-14) so%e = cmplx(0.D0, aimag(so%e), p8)
    ! where (abs(aimag(so%e)) .lt. 5.D-14) so%e = cmplx(real(so%e), 0.D0, p8)

    deallocate( bl )
    deallocate( bl2 )

    return
  end procedure

  module procedure fefe
    type(scalar) :: svis
    real(p8) :: dt_

    if (s%space(1:3).ne.'FFF') stop 'fefe: all input scalars must be in FFF for time stepping'
    if (s_rhs_nonlin%space(1:3).ne.'FFF') stop 'fefe: all input scalars must be in FFF for time stepping'

    if (curr_n - ni .eq. totaln) then
      dt_ = totaltime - (totaln-1)*dt ! final timestep's dt adjustment
    else
      dt_ = dt ! normally, globally defined dt (in base) must be called
    endif

    call svis%init(s%glb_sz, s%axis_comm)
    svis = s

    if (hyperpow .eq. 0) then ! zero hyperviscosity
      if (visc .lt. 5.D-14) then ! inviscid (zero viscosity) -- diffusive operator is null
        svis%e = 0.D0
      else ! viscous (nonzero viscosity) -- diffusive operator is visc*del^2
        svis = del2(s, tfm)
        svis%e = svis%e * visc
      endif
    else ! nonzero hyperviscosity -- diffusive operator is visc*del^2 + hypervisc*del^p
      svis = helmp(s, hyperpow, 0.D0, visc/(hypervisc*(-1.D0)**(hyperpow/2+1)), tfm)
      svis%e = svis%e * (hypervisc*(-1.D0)**(hyperpow/2+1))
    endif

    s%e = s%e + dt_ * (s_rhs_nonlin%e + svis%e)
    s%ln = s%ln + dt_ * (s_rhs_nonlin%ln + svis%ln)

    call svis%dealloc()

  end procedure

  module procedure abab
    logical :: is_2nd_svis_p_ ! 2nd arg is not s_p but svis_p (so that we skip its calculation)
    type(scalar) :: svis, svis_p
    real(p8) :: dt_

    if (s%space(1:3).ne.'FFF') stop 'fefe: all input scalars must be in FFF for time stepping'
    if (s_rhs_nonlin%space(1:3).ne.'FFF') stop 'fefe: all input scalars must be in FFF for time stepping'
    if (s_p%space(1:3).ne.'FFF') stop 'fefe: all input scalars must be in FFF for time stepping'
    if (s_rhs_nonlin_p%space(1:3).ne.'FFF') stop 'fefe: all input scalars must be in FFF for time stepping'

    if (present(is_2nd_svis_p)) then
      is_2nd_svis_p_ = is_2nd_svis_p
    else
      is_2nd_svis_p_ = .false. ! by default 2nd arg is s_p
    endif

    if (curr_n - ni .eq. totaln) then
      dt_ = totaltime - (totaln-1)*dt ! final timestep's dt adjustment
    else
      dt_ = dt ! normally, globally defined dt (in base) must be called
    endif

    call svis%init(s%glb_sz, s%axis_comm)
    svis = s

    if (hyperpow .eq. 0) then ! zero hyperviscosity
      if (visc .lt. 5.D-14) then ! inviscid (zero viscosity) -- diffusive operator is null
        svis%e = 0.D0
      else ! viscous (nonzero viscosity) -- diffusive operator is visc*del^2
        svis = del2(s, tfm)
        svis%e = svis%e * visc
      endif
    else ! nonzero hyperviscosity -- diffusive operator is visc*del^2 + hypervisc*del^p
      svis = helmp(s, hyperpow, 0.D0, visc/(hypervisc*(-1.D0)**(hyperpow/2+1)), tfm)
      svis%e = svis%e * (hypervisc*(-1.D0)**(hyperpow/2+1))
    endif

    call svis_p%init(s_p%glb_sz, s_p%axis_comm)
    svis_p = s_p

    if (.not. is_2nd_svis_p_) then ! so 2nd arg is s_p. by default this is not skipped.
      if (hyperpow .eq. 0) then ! zero hyperviscosity
        if (visc .lt. 5.D-14) then ! inviscid (zero viscosity) -- diffusive operator is null
          svis%e = 0.D0
        else ! viscous (nonzero viscosity) -- diffusive operator is visc*del^2
          svis_p = del2(s_p, tfm)
          svis_p%e = svis_p%e * visc
        endif
      else ! nonzero hyperviscosity -- diffusive operator is visc*del^2 + hypervisc*del^p
        svis_p = helmp(s_p, hyperpow, 0.D0, visc/(hypervisc*(-1.D0)**(hyperpow/2+1)), tfm)
        svis_p%e = svis_p%e * (hypervisc*(-1.D0)**(hyperpow/2+1))
      endif
    endif

    s%e = s%e + dt_ * (1.5D0*(s_rhs_nonlin%e + svis%e) - .5D0*(s_rhs_nonlin_p%e + svis_p%e))
    s%ln = s%ln + dt_ * (1.5D0*(s_rhs_nonlin%ln + svis%ln) - .5D0*(s_rhs_nonlin_p%ln + svis_p%ln))

    s_p = s
    s_rhs_nonlin_p = s_rhs_nonlin

    call svis%dealloc()
    call svis_p%dealloc()

  end procedure

  module procedure febe
    type(scalar) :: sh
    real(p8) :: dt_, a, b

    if (s%space(1:3).ne.'FFF') stop 'fefe: all input scalars must be in FFF for time stepping'
    if (s_rhs_nonlin%space(1:3).ne.'FFF') stop 'fefe: all input scalars must be in FFF for time stepping'

    if (curr_n - ni .eq. totaln) then
      dt_ = totaltime - (totaln-1)*dt ! final timestep's dt adjustment
    else
      dt_ = dt ! normally, globally defined dt (in base) must be called
    endif

    call sh%init(s%glb_sz, s%axis_comm)
    sh = s

    sh%e = s%e + dt_ * s_rhs_nonlin%e
    sh%ln = s%ln + dt_ * s_rhs_nonlin%ln

    if (hyperpow .eq. 0) then ! zero hyperviscosity
      if (visc .lt. 5.D-14) then ! inviscid (zero viscosity) -- diffusive operator is null
        stop 'febe: inviscid case and no linear term in rhs. semi-implicit time adv is impossible'
      else ! viscous (nonzero viscosity) -- diffusive operator is visc*del^2
           ! F(K+1) = (1-DT*VISC*DEL2)^(-1) F(K+1/2) 
           ! = IHELM( (-1/(DT*VISC))*F(K+1/2)), ALP=(-1/(DT*VISC)) )
        a = -1.D0 / (dt_*visc)
        sh%ln = a*sh%ln
        sh%e = a*sh%e
        s = ihelm(sh, a, tfm)
      endif
    else ! nonzero hyperviscosity -- diffusive operator is visc*del^2 + hypervisc*del^p
         ! F(K+1) = (1-DT*VISC*DEL2-DT*VISCP*DELP)^(-1) F(K+1/2)
         ! = [(-DT*VISCP)*(DELP+(VISC/VISCP)*DEL2+(-1/(DT*VISCP)))]^(-1) F(K+1/2)
         ! = IHELMP( (-1/(DT*VISCP))*F(K+1/2)), ALP=(-1/(DT*VISCP)), BET=(VISC/VISCP) )
      a = -1.D0/(dt_*(hypervisc*(-1.D0)**(hyperpow/2+1)))
      b = visc/(hypervisc*(-1.D0)**(hyperpow/2+1))
      sh%ln = a*sh%ln
      sh%e = a*sh%e
      s = ihelmp(sh, hyperpow, a, b, tfm)
    endif

    call sh%dealloc()

  end procedure

  module procedure abcn
    type(scalar) :: sh, svis
    real(p8) :: dt_, a, b

    if (s%space(1:3).ne.'FFF') stop 'fefe: all input scalars must be in FFF for time stepping'
    if (s_rhs_nonlin%space(1:3).ne.'FFF') stop 'fefe: all input scalars must be in FFF for time stepping'
    if (s_p%space(1:3).ne.'FFF') stop 'fefe: all input scalars must be in FFF for time stepping'
    if (s_rhs_nonlin_p%space(1:3).ne.'FFF') stop 'fefe: all input scalars must be in FFF for time stepping'

    if (curr_n - ni .eq. totaln) then
      dt_ = totaltime - (totaln-1)*dt ! final timestep's dt adjustment
    else
      dt_ = dt ! normally, globally defined dt (in base) must be called
    endif

    call sh%init(s%glb_sz, s%axis_comm)
    sh = s

    sh%e = s%e + dt_ * (1.5D0*s_rhs_nonlin%e - .5D0*s_rhs_nonlin_p%e)
    sh%ln = s%ln + dt_ * (1.5D0*s_rhs_nonlin%ln - .5D0*s_rhs_nonlin_p%ln)

    call svis%init(s%glb_sz, s%axis_comm)
    svis = s

    if (hyperpow .eq. 0) then ! zero hyperviscosity
      if (visc .lt. 5.D-14) then ! inviscid (zero viscosity) -- diffusive operator is null
        svis%e = 0.D0
      else ! viscous (nonzero viscosity) -- diffusive operator is visc*del^2
        svis = del2(s, tfm)
        svis%e = svis%e * visc
      endif
    else ! nonzero hyperviscosity -- diffusive operator is visc*del^2 + hypervisc*del^p
      svis = helmp(s, hyperpow, 0.D0, visc/(hypervisc*(-1.D0)**(hyperpow/2+1)), tfm)
      svis%e = svis%e * (hypervisc*(-1.D0)**(hyperpow/2+1))
    endif

    if (hyperpow .eq. 0) then ! zero hyperviscosity
      if (visc .lt. 5.D-14) then ! inviscid (zero viscosity) -- diffusive operator is null
        stop 'abcn: inviscid case and no linear term in rhs. semi-implicit time adv is impossible'
      else ! viscous (nonzero viscosity) -- diffusive operator is visc*del^2
           ! F(K+1) = (1-DT/2*VISC*DEL2)^(-1) [F(K+1/2) + DT/2*LTERM(K)] 
           ! = IHELM( (-2/(DT*VISC))*[F(K+1/2)+DT/2*LTERM(K)], ALP=(-2/(DT*VISC)) )
        a = -2.D0 / (dt_*visc)
        sh%ln = a*(sh%ln + dt_/2.D0*svis%ln)
        sh%e = a*(sh%e + dt_/2.D0*svis%e)
        s = ihelm(sh, a, tfm)
      endif
    else ! nonzero hyperviscosity -- diffusive operator is visc*del^2 + hypervisc*del^p
         ! F(K+1) = (1-DT/2*VISC*DEL2-DT/2*VISCP*DELP)^(-1) [F(K+1/2) + DT/2*LTERM(K)]
         ! = [(-DT/2*VISCP)*(DELP+(VISC/VISCP)*DEL2+(-2/(DT*VISCP)))]^(-1) [F(K+1/2) + DT/2*LTERM(K)]
         ! = IHELMP((-2/(DT*VISCP))*[F(K+1/2) + DT/2*LTERM(K)]), ALP=(-2/(DT*VISCP)), BET=(VISC/VISCP))
      a = -2.D0/(dt_*(hypervisc*(-1.D0)**(hyperpow/2+1)))
      b = visc/(hypervisc*(-1.D0)**(hyperpow/2+1))
      sh%ln = a*(sh%ln + dt_/2.D0*svis%ln)
      sh%e = a*(sh%e + dt_/2.D0*svis%e)
      s = ihelmp(sh, hyperpow, a, b, tfm)
    endif

    s_p = s
    s_rhs_nonlin_p = s_rhs_nonlin

    call sh%dealloc()
    call svis%dealloc()

  end procedure

  module procedure vector_product
    integer(i4) :: mm, nn, kk
    real(p8) :: a1, a2, a3, c1, c2, c3, b1, b2, b3, d1, d2, d3

    if (vr%space(1:3).ne.'PPP') stop 'vector_product: to compute v x u, vr must be in PPP'
    if (vp%space(1:3).ne.'PPP') stop 'vector_product: to compute v x u, vp must be in PPP'
    if (vz%space(1:3).ne.'PPP') stop 'vector_product: to compute v x u, vz must be in PPP'
    if (ur%space(1:3).ne.'PPP') stop 'vector_product: to compute v x u, ur must be in PPP'
    if (up%space(1:3).ne.'PPP') stop 'vector_product: to compute v x u, up must be in PPP'
    if (uz%space(1:3).ne.'PPP') stop 'vector_product: to compute v x u, uz must be in PPP'

    do kk = 1, min(vr%loc_sz(3), nz - vr%loc_st(3))
      do mm = 1, min(vr%loc_sz(2), np/2 - vr%loc_st(2))
        do nn = 1, min(vr%loc_sz(1), nr - vr%loc_st(1))
          a1 = real(vr%e(nn,mm,kk)); a2 = real(vp%e(nn,mm,kk)); a3 = real(vz%e(nn,mm,kk))
          b1 = real(ur%e(nn,mm,kk)); b2 = real(up%e(nn,mm,kk)); b3 = real(uz%e(nn,mm,kk))
          c1 = aimag(vr%e(nn,mm,kk)); c2 = aimag(vp%e(nn,mm,kk)); c3 = aimag(vz%e(nn,mm,kk))
          d1 = aimag(ur%e(nn,mm,kk)); d2 = aimag(up%e(nn,mm,kk)); d3 = aimag(uz%e(nn,mm,kk))
          vr%e(nn,mm,kk) = cmplx(a2*b3-a3*b2, c2*d3-c3*d2, p8)
          vp%e(nn,mm,kk) = cmplx(a3*b1-a1*b3, c3*d1-c1*d3, p8)
          vz%e(nn,mm,kk) = cmplx(a1*b2-a2*b1, c1*d2-c2*d1, p8)
        enddo
      enddo
    enddo
  end procedure

  module procedure vector_projection
    integer(i4) :: mm, nn, kk, i, n, mv
    real(p8) :: kv
    real(p8), dimension(:), allocatable :: r, fff
    real(p8), dimension(:,:), allocatable :: v, d, t, pfd
    complex(p8), dimension(:), allocatable :: inf, w1, w2

    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak

    type(scalar) :: ur, up, uz

    r = ell*sqrt((1.D0+tfm%x)/(1.D0-tfm%x))

    if (vr%space(1:3).ne.'PPP') stop 'vector_projection: vr must be in PPP'
    if (vp%space(1:3).ne.'PPP') stop 'vector_projection: vr must be in PPP'
    if (vz%space(1:3).ne.'PPP') stop 'vector_projection: vz must be in PPP'
    if (psi%space(1:3).ne.'FFF') stop 'vector_projection: psi must be in FFF'
    if (chi%space(1:3).ne.'FFF') stop 'vector_projection: chi must be in FFF'

    psi%e = 0.D0; chi%e = 0.D0
    psi%ln = 0.D0; chi%ln = 0.D0

    call ur%init(vr%glb_sz, vr%axis_comm); ur = vr
    call up%init(vr%glb_sz, vr%axis_comm); up = vp
    call uz%init(vr%glb_sz, vr%axis_comm); uz = vz

    if (ur%loc_st(1)+1 .le. nr) then
      do mm = 1, ur%loc_sz(2)
        do kk = 1, ur%loc_sz(3)
          ur%e(1:min(ur%loc_sz(1), nr-ur%loc_st(1)), mm, kk) = &
          ur%e(1:min(ur%loc_sz(1), nr-ur%loc_st(1)), mm, kk) &
          * r(ur%loc_st(1)+1:min(ur%loc_st(1)+ur%loc_sz(1), nr)) ! now ur actually stores r*ur
        enddo
      enddo
    endif

    if (up%loc_st(1)+1 .le. nr) then
      do mm = 1, up%loc_sz(2)
        do kk = 1, up%loc_sz(3)
          up%e(1:min(up%loc_sz(1), nr-up%loc_st(1)), mm, kk) = &
          up%e(1:min(up%loc_sz(1), nr-up%loc_st(1)), mm, kk) &
          * r(up%loc_st(1)+1:min(up%loc_st(1)+up%loc_sz(1), nr)) ! now up actually stores r*up
        enddo
      enddo
    endif

    call trans(ur, 'PFP', tfm)
    if (nz .gt. 1) then
      call ur%exchange(1, 3)
      call vertical_fft_forward(ur, tfm)
      call ur%exchange(3, 1)
    endif
    ur%space = 'PFF' ! this in non-standard yet neccesary for calculation in this subroutine

    call trans(up, 'FFF', tfm)
    inf = calcat1(up, tfm)
    psi%ln = -1.D0/2.D0*real(inf(1))
    if (nz .gt. 1) call up%exchange(3, 1)
    call rtrans_backward(up, tfm)
    up%space = 'PFF' ! this in non-standard yet neccesary for calculation in this subroutine
    if ((up%loc_st(2) .eq. 0) .and. (up%loc_st(3) .eq. 0)) then
      up%e(:nr, 1, 1) = up%e(:nr, 1, 1) + psi%ln*(1.D0 + tfm%x)
    endif

    call trans(uz, 'PFP', tfm)
    if (nz .gt. 1) then
      call uz%exchange(1, 3)
      call vertical_fft_forward(uz, tfm)
      call uz%exchange(3, 1)
    endif
    uz%space = 'PFF' ! this in non-standard yet neccesary for calculation in this subroutine

    call ur%chop_offset(3); call up%chop_offset(3); call uz%chop_offset(3)

    call chop_index(ur, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    allocate( fff(nr/2) )
    allocate( v(nrc, nr/2) )
    allocate( d(nrc, nr/2) )
    allocate( t(nrc, nr/2) )
    allocate( pfd(nr/2, nrc+1) )
    allocate( w1(nrc) )
    allocate( w2(nrc) )

    fff = (1.D0 - tfm%x(:nr/2)**2.D0)/tfm%w(:nr/2)

    if (nz .gt. 1) then
      call psi%exchange(3, 1)
      call chi%exchange(3, 1)
    endif

    do mm = 1, min(ur%loc_sz(2), npc-ur%loc_st(2))
      nn = nrcs(ur%loc_st(2) + mm)
      mv = m(ur%loc_st(2) + mm)
      pfd(:nr/2, :nn+1) = tfm%pf(:nr/2, :nn+1, mm+ur%loc_st(2)) .mul. leg_xxdx(mv, nn+1, nn+1, tfm)
      do i = 1, nn
        n = max(1, mv+i-1)
        v(i,:nr/2) = tfm%pf(:nr/2,i,mm+ur%loc_st(2))/n/(n+1)/fff
        d(i,:nr/2) = pfd(:nr/2,i)/n/(n+1)/fff
        t(i,:nr/2) = tfm%pf(:nr/2,i,mm+ur%loc_st(2))*tfm%w(:nr/2)
      enddo

      if (1 .le. min(psi%loc_sz(3), nzc - psi%loc_st(3))) then
        do kk = 1, min(psi%loc_sz(3), nzc - psi%loc_st(3))
          w1 = 0.D0; w2 = 0.D0; if (mv .ne. 0) call eomul(v(:nn,:), ur%e(:nr,mm,kk), w1(:nn))
          call oemul(d(:nn,:), up%e(:nr,mm,kk), w2(:nn))
          psi%e(:nn,mm,kk) = -iu*mv*w1(:nn) - w2(:nn)
          w1 = 0.D0; w2 = 0.D0; call oemul(d(:nn,:), ur%e(:nr,mm,kk), w1(:nn))
          if (mv .ne. 0) call eomul(v(:nn,:), up%e(:nr,mm,kk), w2(:nn))
          call eomul(t(:nn,:), uz%e(:nr,mm,kk), chi%e(:nn,mm,kk))
          kv = ak(psi%loc_st(3) + kk)
          chi%e(:nn,mm,kk) = iu*kv*w1(:nn)+mv*kv*w2(:nn)-chi%e(:nn,mm,kk)
        enddo
      endif
      if (max(nzcu - psi%loc_st(3), 1) .le. psi%loc_sz(3)) then
        do kk = max(nzcu - psi%loc_st(3), 1), psi%loc_sz(3)
          w1 = 0.D0; w2 = 0.D0; if (mv .ne. 0) call eomul(v(:nn,:), ur%e(:nr,mm,kk), w1(:nn))
          call oemul(d(:nn,:), up%e(:nr,mm,kk), w2(:nn))
          psi%e(:nn,mm,kk) = -iu*mv*w1(:nn) - w2(:nn)
          w1 = 0.D0; w2 = 0.D0; call oemul(d(:nn,:), ur%e(:nr,mm,kk), w1(:nn))
          if (mv .ne. 0) call eomul(v(:nn,:), up%e(:nr,mm,kk), w2(:nn))
          call eomul(t(:nn,:), uz%e(:nr,mm,kk), chi%e(:nn,mm,kk))
          kv = ak(psi%loc_st(3) + kk)
          chi%e(:nn,mm,kk) = iu*kv*w1(:nn)+mv*kv*w2(:nn)-chi%e(:nn,mm,kk)
        enddo
      endif

    enddo

    if (nz .gt. 1) then
      call psi%exchange(1, 3)
      call chi%exchange(1, 3)
    endif

    chi = idel2_proln(chi, tfm)
    deallocate( fff, v, d, t, pfd, w1, w2 )

    call psi%chop_offset(0); call chi%chop_offset(0)
    call chop(psi, tfm); call chop(chi, tfm)

    call ur%dealloc(); call up%dealloc(); call uz%dealloc()
  end procedure

  module procedure vector_reconstruction
    integer(i4) :: mm, nn, kk, mv
    real(p8) :: kv
    type(scalar) :: ur, up, uz

    real(p8), dimension(:), allocatable :: r

    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak

    r = ell*sqrt((1.D0+tfm%x)/(1.D0-tfm%x))

    if (psi%space(1:3).ne.'FFF') stop 'vector_projection: psi must be in FFF'
    if (chi%space(1:3).ne.'FFF') stop 'vector_projection: chi must be in FFF'
    if (vr%space(1:3).ne.'PPP') stop 'vector_projection: vr must be in PPP'
    if (vp%space(1:3).ne.'PPP') stop 'vector_projection: vp must be in PPP'
    if (vz%space(1:3).ne.'PPP') stop 'vector_projection: vz must be in PPP'

    call ur%init(chi%glb_sz, chi%axis_comm)
    call up%init(psi%glb_sz, psi%axis_comm)
    call uz%init(chi%glb_sz, chi%axis_comm)

    ur = chi; call ur%chop_offset(3) 
    up = psi; call up%chop_offset(3) 
    uz = chi; call uz%chop_offset(3) ! uz now stores chi

    call chop_index(ur, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    ur = xxdx(uz, tfm) ! ur now actually stores r*d(chi)/dr
    up = xxdx(up, tfm) ! up now actually stores r*d(psi)/dr

    do mm = 1, min(ur%loc_sz(2), npc - ur%loc_st(2))
      nn = nrcs(ur%loc_st(2) + mm)
      mv = m(ur%loc_st(2) + mm)
      do kk = 1, nzdim
        if ((kk .le. nzc) .or. (kk .ge. nzcu)) then
          kv = ak(kk)
          ur%e(:min(ur%loc_sz(1), nn-ur%loc_st(1)),mm,kk) = &
          iu*mv*psi%e(:min(psi%loc_sz(1), nn-psi%loc_st(1)),mm,kk) + &
          iu*kv*ur%e(:min(ur%loc_sz(1), nn-ur%loc_st(1)),mm,kk) ! ur now stores d(psi)/dp + r*d/dr(d/dz(chi))
          up%e(:min(up%loc_sz(1), nn-up%loc_st(1)),mm,kk) = &
          -up%e(:min(up%loc_sz(1), nn-up%loc_st(1)),mm,kk) - &
          mv*kv*uz%e(:min(uz%loc_sz(1), nn-uz%loc_st(1)),mm,kk) ! up now stores -r*d(psi)/dr + d/dp(d/dz(chi))
        endif
      enddo
    enddo

    call ur%chop_offset(0); call up%chop_offset(0)
    call chop(ur, tfm); call chop(up, tfm)

    call trans(ur, 'PPP', tfm); call trans(up, 'PPP', tfm)
    
    if (ur%loc_st(1)+1 .le. nr) then
      do mm = 1, ur%loc_sz(2)
        do kk = 1, ur%loc_sz(3)
          ur%e(1:min(ur%loc_sz(1), nr-ur%loc_st(1)), mm, kk) = &
          ur%e(1:min(ur%loc_sz(1), nr-ur%loc_st(1)), mm, kk) &
          / r(ur%loc_st(1)+1:min(ur%loc_st(1)+ur%loc_sz(1), nr)) ! now ur stores 1/r*d(psi)/dp + d/dr(d/dz(chi))
        enddo
      enddo
    endif

    if (up%loc_st(1)+1 .le. nr) then
      do mm = 1, up%loc_sz(2)
        do kk = 1, up%loc_sz(3)
          up%e(1:min(up%loc_sz(1), nr-up%loc_st(1)), mm, kk) = &
          up%e(1:min(up%loc_sz(1), nr-up%loc_st(1)), mm, kk) &
          / r(up%loc_st(1)+1:min(up%loc_st(1)+up%loc_sz(1), nr)) ! now up stores -d(psi)/dr + 1/r*d/dp(d/dz(chi))
        enddo
      enddo
    endif

    uz = del2h(uz, tfm) ! now uz stores del2h(chi)
    uz%e = -uz%e ! now uz stores -del2h(chi)

    call uz%chop_offset(0)
    call chop(uz, tfm)

    call trans(uz, 'PPP', tfm)

    vr = ur
    vp = up
    vz = uz

    call ur%dealloc()
    call up%dealloc()
    call uz%dealloc()

  end procedure

  module procedure curl_vector_reconstruction
    integer(i4) :: mm, nn, kk, mv
    real(p8) :: kv
    type(scalar) :: or, op, oz

    real(p8), dimension(:), allocatable :: r

    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak

    r = ell*sqrt((1.D0+tfm%x)/(1.D0-tfm%x))

    if (psi%space(1:3).ne.'FFF') stop 'curl_vector_projection: psi must be in FFF'
    if (chi%space(1:3).ne.'FFF') stop 'curl_vector_projection: chi must be in FFF'
    if (wr%space(1:3).ne.'PPP') stop 'curl_vector_projection: wr must be in PPP'
    if (wp%space(1:3).ne.'PPP') stop 'curl_vector_projection: wp must be in PPP'
    if (wz%space(1:3).ne.'PPP') stop 'curl_vector_projection: wz must be in PPP'

    call or%init(psi%glb_sz, psi%axis_comm)
    call op%init(chi%glb_sz, chi%axis_comm)
    call oz%init(psi%glb_sz, psi%axis_comm)

    or = psi; call or%chop_offset(3)
    op = chi; call op%chop_offset(3)
    oz = psi; call oz%chop_offset(3) ! oz now stores psi

    call chop_index(or, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    or = xxdx(oz, tfm) ! or now actually stores r*d(psi)/dr

    op = del2(op, tfm)
    op = xxdx(op, tfm) ! op now actually stores r*d(del2(chi))/dr

    do mm = 1, min(or%loc_sz(2), npc - or%loc_st(2))
      nn = nrcs(or%loc_st(2) + mm)
      mv = m(or%loc_st(2) + mm)
      do kk = 1, nzdim
        if ((kk .le. nzc) .or. (kk .ge. nzcu)) then
          kv = ak(kk)
          or%e(:min(or%loc_sz(1), nn-or%loc_st(1)),mm,kk) = &
          -iu*mv*chi%e(:min(chi%loc_sz(1), nn-chi%loc_st(1)),mm,kk) + &
          iu*kv*or%e(:min(or%loc_sz(1), nn-or%loc_st(1)),mm,kk) ! or now stores -d(del2(chi))/dp + r*d/dr(d/dz(psi))
          op%e(:min(op%loc_sz(1), nn-op%loc_st(1)),mm,kk) = &
          op%e(:min(op%loc_sz(1), nn-op%loc_st(1)),mm,kk) - &
          mv*kv*oz%e(:min(oz%loc_sz(1), nn-oz%loc_st(1)),mm,kk) ! op now stores r*d(del2(chi))/dr + d/dp(d/dz(psi))
        endif
      enddo
    enddo

    call or%chop_offset(0); call op%chop_offset(0)
    call chop(or, tfm); call chop(op, tfm)

    call trans(or, 'PPP', tfm); call trans(op, 'PPP', tfm)
    
    if (or%loc_st(1)+1 .le. nr) then
      do mm = 1, or%loc_sz(2)
        do kk = 1, or%loc_sz(3)
          or%e(1:min(or%loc_sz(1), nr-or%loc_st(1)), mm, kk) = &
          or%e(1:min(or%loc_sz(1), nr-or%loc_st(1)), mm, kk) &
          / r(or%loc_st(1)+1:min(or%loc_st(1)+or%loc_sz(1), nr)) ! now or stores -1/r*d(del2(chi))/dp + d/dr(d/dz(psi))
        enddo
      enddo
    endif

    if (op%loc_st(1)+1 .le. nr) then
      do mm = 1, op%loc_sz(2)
        do kk = 1, op%loc_sz(3)
          op%e(1:min(op%loc_sz(1), nr-op%loc_st(1)), mm, kk) = &
          op%e(1:min(op%loc_sz(1), nr-op%loc_st(1)), mm, kk) &
          / r(op%loc_st(1)+1:min(op%loc_st(1)+op%loc_sz(1), nr)) ! now op stores d(del2(chi))/dr + 1/r*d/dp(d/dz(psi))
        enddo
      enddo
    endif

    oz = del2h(oz, tfm) ! now oz stores del2h(psi)
    oz%e = -oz%e ! now oz stores -del2h(psi)

    call oz%chop_offset(0)
    call chop(oz, tfm)

    call trans(oz, 'PPP', tfm)

    wr = or
    wp = op
    wz = oz

    call or%dealloc()
    call op%dealloc()
    call oz%dealloc()

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

    if (.not. all(abs(se(:, npc+1, :)) .lt. 5.D-14)) then
      if (is_warning) write(*,*) '[warning] horizontal_fft_forward: residues in the chopped region'
    endif

    deallocate( b )
    deallocate( c )

    ! where (abs(real(s%e)) .lt. 5.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    ! where (abs(aimag(s%e)) .lt. 5.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

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

    if ((is_warning) .and. (.not. all(abs(se(:, npc+1, :)) .lt. 5.D-14))) then
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

    ! where (abs(real(s%e)) .lt. 5.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    ! where (abs(aimag(s%e)) .lt. 5.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

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

    if (.not. all(abs(se(:, :, nzc+1:nzcu-1)) .lt. 5.D-14)) then
      if (is_warning) write(*,*) '[warning] vertical_fft_forward: residues in the chopped region'
    endif

    deallocate( b )
    deallocate( c )

    ! where (abs(real(s%e)) .lt. 5.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    ! where (abs(aimag(s%e)) .lt. 5.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

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

    if (.not. all(abs(se(:, :, nzc+1:nzcu-1)) .lt. 5.D-14)) then
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

    ! where (abs(real(s%e)) .lt. 5.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    ! where (abs(aimag(s%e)) .lt. 5.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

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
    allocate( be(nrh, size(s%e,3)) )
    allocate( bo(nrh, size(s%e,3)) )

    do mm = 1, min(s%loc_sz(2), npc - s%loc_st(2))
      be = 0.D0; bo = 0.D0
      do i = 1, nrh
        be(i,:) = (s%e(i,mm,:)+s%e(nr-i+1,mm,:))*tfm%w(i)
        bo(i,:) = (s%e(i,mm,:)-s%e(nr-i+1,mm,:))*tfm%w(i)
      enddo
      nn = nrcs(s%loc_st(2) + mm)
      if (nn .ge. 1) then
        se(1:nn:2,mm,:) = transpose(tfm%pf(:nrh,1:nn:2,s%loc_st(2) + mm)) .mul. be
      endif
      if (nn .ge. 2) then
        se(2:nn:2,mm,:) = transpose(tfm%pf(:nrh,2:nn:2,s%loc_st(2) + mm)) .mul. bo
      endif
    enddo

    s%e = se

    deallocate( be )
    deallocate( bo )
    deallocate( se )

    ! where (abs(real(s%e)) .lt. 5.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    ! where (abs(aimag(s%e)) .lt. 5.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

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
    allocate( be(nrh, size(s%e,3)) )
    allocate( bo(nrh, size(s%e,3)) )

    do mm = 1, min(s%loc_sz(2), npc - s%loc_st(2))
      be = 0.D0; bo = 0.D0
      nn = nrcs(s%loc_st(2)+mm)
      if (nn .ge. 1) then
        be = tfm%pf(:nrh,1:nn:2,s%loc_st(2)+mm) .mul. s%e(1:nn:2,mm,:)
      else
        be = 0.D0
      endif
      if (nn .ge. 2) then
        bo = tfm%pf(:nrh,2:nn:2,s%loc_st(2)+mm) .mul. s%e(2:nn:2,mm,:)
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

    ! where (abs(real(s%e)) .lt. 5.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    ! where (abs(aimag(s%e)) .lt. 5.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

  end subroutine

! ======================================================================================================== !

  subroutine chop_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)
    implicit none
    class(scalar) :: s
    class(tfm_kit) :: tfm
    integer(i4), intent(inout) :: nrdim, npdim, nzdim
    integer(i4), intent(inout) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable, intent(inout) :: nrcs, m
    real(p8), dimension(:), allocatable, intent(inout) :: ak

    integer(i4) :: mm

    select type(tfm)
      type is (tfm_kit_1d)
        nrdim = tfm%nrdim; nrc = tfm%chops(1); nrcs = tfm%chops
        npdim = 1        ; npc = 1;
        nzdim = 1        ; nzc = 1; nzcu = 1
        m = (/ 0 /); ak = (/ 0 /)
      type is (tfm_kit_2d)
        nrdim = tfm%nrdim; nrc = tfm%chops(1); nrcs = tfm%chops
        npdim = tfm%npdim; npc = tfm%chopp
        nzdim = 1        ; nzc = 1; nzcu = 1
        m = tfm%m; ak = (/ 0 /)
      type is (tfm_kit_3d)
        nrdim = tfm%nrdim; nrc = tfm%chops(1); nrcs = tfm%chops
        npdim = tfm%npdim; npc = tfm%chopp
        nzdim = tfm%nzdim; nzc = tfm%chopzl; nzcu = tfm%chopzu
        m = tfm%m; ak = tfm%ak
    end select

    if ( (is_warning) .and. &
         ((s%glb_sz(1) .ne. nrdim) .or. &
          (s%glb_sz(2) .ne. npdim) .or. &
          (s%glb_sz(3) .ne. nzdim)) ) then
      write(*,*) '[wanring] chop_index: scalar dimension incompatible with the dimension defined in a kit'
    endif

    nrc = nrc + s%nrchop_offset
    npc = npc + s%npchop_offset
    nzc = nzc + s%nzchop_offset
    nzcu = nzcu - s%nzchop_offset

    if (nrc .gt. nrdim) stop 'chop_index: chopping in r too large'
    if (npc .gt. npdim) stop 'chop_index: chopping in p too large'
    if (2*s%nzchop_offset .gt. nzcu - nzc) stop 'chop_index: chopping in z too large'

    do mm = 1, npc
      nrcs(mm) = max(min(nrc, nrc-m(mm)), 0)
    enddo
    nrcs(npc+1:) = 0

    return
  end subroutine

! ======================================================================================================== !

  subroutine eomul(a,b,c)
! [usage]:
! compute c(:ni) = a(:ni, :nj) .mul. b(:nj) where a has pattern as follows:
! a(2*i  :j) =  a(2*i  :nj-j+1) and
! a(2*i-1:j) = -a(2*i-1:nj-j+1)
! only half of a should be given on input.
! [variables]:
! a >> a special real-valued matrix with the pattern above. dim: ni x nj/2
! b >> a matrix (vector) of the dimension nj x 1
! c >> a .mul. b
  implicit none

  real(p8), dimension(:,:), intent(in)   :: a
  complex(p8), dimension(:), intent(in)  :: b
  complex(p8), dimension(:), intent(out) :: c

  complex(p8), dimension(:), allocatable :: be, bo
  integer                                :: ni, nj, njh

  ni = size(a,1)
  njh = size(a,2)
  nj = size(b,1)

  if (njh*2 .ne. nj) stop 'eomul: size mismatch.'

  allocate( be(njh) )
  allocate( bo(njh) )

  be = b(1:njh) + b(nj:njh+1:-1)
  bo = b(1:njh) - b(nj:njh+1:-1)

  c(1::2) = a(1::2,:) .mul. be
  if (ni .gt. 1) c(2::2) = a(2::2,:) .mul. bo

  deallocate( be, bo )

  return
  end subroutine

! ======================================================================================================== !

  subroutine oemul(a,b,c)
! [usage]:
! compute c(:ni) = a(:ni, :nj) .mul. b(:nj) where a has pattern as follows:
! a(2*i  :j) = -a(2*i  :nj-j+1) and
! a(2*i-1:j) =  a(2*i-1:nj-j+1)
! only half of a should be given on input.
! [variables]:
! a >> a special real-valued matrix with the pattern above. dim: ni x nj/2
! b >> a matrix (vector) of the dimension nj x 1
! c >> a .mul. b
  implicit none

  real(p8), dimension(:,:), intent(in)   :: a
  complex(p8), dimension(:), intent(in)  :: b
  complex(p8), dimension(:), intent(out) :: c

  complex(p8), dimension(:), allocatable :: be, bo
  integer                                :: ni, nj, njh

  ni = size(a,1)
  njh = size(a,2)
  nj = size(b,1)

  if (njh*2 .ne. nj) stop 'oemul: size mismatch.'

  allocate( be(njh) )
  allocate( bo(njh) )

  be = b(1:njh) + b(nj:njh+1:-1)
  bo = b(1:njh) - b(nj:njh+1:-1)

  c(1::2) = a(1::2,:) .mul. bo
  if (ni .gt. 1) c(2::2) = a(2::2,:) .mul. be

  deallocate( be, bo )

  return
  end subroutine

end submodule