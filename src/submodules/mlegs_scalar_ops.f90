submodule (mlegs_scalar) mlegs_scalar_ops
  implicit none

contains

  module procedure chop
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, kk

    call chop_glb_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

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
        do kk = 1, s%loc_sz(3)
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
        do kk = 1, s%loc_sz(3)
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

  module procedure calcat0
    type(scalar) :: st
    complex(p8), dimension(:), allocatable :: calc_

    call st%init(s%glb_sz, s%axis_comm)
    st = s

    call trans(st, 'FFF', tfm)

    allocate( calc(nz) )
    if (s%loc_st(2) .eq. 0) then
      calc = transpose(st%e(1:st%loc_sz(1),1,:)) .mul. tfm%at0(st%loc_st(1)+1:st%loc_st(1)+st%loc_sz(1))
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

    call st%init(s%glb_sz, s%axis_comm)
    st = s

    if ((.not.(st%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] calcat1: an input scalar is not FFF; additional operations needed for tfm'
      call trans(st, 'FFF', tfm)
    endif

    allocate( calc(nz) )
    if (s%loc_st(2) .eq. 0) then
      calc = transpose(st%e(1:st%loc_sz(1),1,:)) .mul. tfm%at1(st%loc_st(1)+1:st%loc_st(1)+st%loc_sz(1))
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

    call chop_glb_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

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
        so%e(nn,mm,:) = -s%e(nn,mm,:)*n*(n+1)/ell**2.D0
      enddo
    enddo

    if ((so%loc_st(1) .eq. 0) .and. (so%loc_st(2) .eq. 0)) then
      so%e(1,1,1) = 1.D0/ell**2.D0/exp(tfm%lognorm(1,1))*s%ln
    endif

    return
  end procedure

  module procedure idelsqp
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, nn, n

    call chop_glb_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

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
    call MPI_allreduce(so%ln, so%ln, 1, MPI_complex16, MPI_sum, comm_glb, MPI_err)

    if (so%loc_st(2) .eq. 0) then
      if (so%loc_st(1) .eq. 0) then
        so%e(1,1,:) = 0.D0
        do nn = 2, min(so%loc_sz(1), nrcs(1) - so%loc_sz(1))
          n = nn - 1
          so%e(nn, 1, :) = -s%e(nn, 1, :)/n/(n+1)*ell**2.D0
        enddo
      else
        do nn = 1, min(so%loc_sz(1), nrcs(so%loc_st(2) + 1) - so%loc_sz(1))
          n = (so%loc_st(1) + nn) - 1
          so%e(nn, 1, :) = -s%e(nn, 1, :)/n/(n+1)*ell**2.D0
        enddo
      endif
      do mm = 2, min(so%loc_sz(2), npc)
        do nn = 1, min(so%loc_sz(1), nrcs(mm) - so%loc_sz(1))
          n = m(mm) + (so%loc_st(1) + nn) - 1
          so%e(nn,mm,:) = -s%e(nn,mm,:)/n/(n+1)*ell**2.D0
        enddo
      enddo
    else
      do mm = 1, min(so%loc_sz(2), npc - so%loc_st(2))
        do nn = 1, min(so%loc_sz(1), nrcs(so%loc_st(2) + mm) - so%loc_st(1))
          n = m(so%loc_st(2) + mm) + (so%loc_st(1) + nn) - 1
          so%e(nn,mm,:) = -s%e(nn,mm,:)/n/(n+1)*ell**2.D0
        enddo
      enddo
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

    call chop_glb_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] xxdx: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s
    call so%exchange(3,1) ! whole r-data resides locally

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
    call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

    return
  end procedure

  module procedure del2
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak
    integer(i4) :: mm, nn, kk, n
    type(real_bndm) :: del2_bnd

    call chop_glb_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] del2: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s
    call so%exchange(3,1) ! whole r-data resides locally

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
    call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

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

    call chop_glb_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] idel2: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                 after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s
    call so%exchange(3,1) ! whole r-data resides locally

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
    call MPI_allreduce(so%ln, so%ln, 1, MPI_complex16, MPI_sum, comm_glb, MPI_err)

    call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

    where (abs(real(so%e)) .lt. 5.D-14) so%e = cmplx(0.D0, aimag(so%e), p8)
    where (abs(aimag(so%e)) .lt. 5.D-14) so%e = cmplx(real(so%e), 0.D0, p8)

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

    call chop_glb_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

    if ((.not.(so%space(1:3) .eq. 'FFF')) .and. (is_warning)) then
      write(*,*) '[warning] idel2: an input scalar is not FFF; additional operations needed for tfm'
      write(*,*) '                 after this operation, both input and output scalar resides in FFF'
      call trans(s, 'FFF', tfm)
    endif

    call so%init(s%glb_sz, s%axis_comm)
    so = s
    call so%exchange(3,1) ! whole r-data resides locally

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
    call MPI_allreduce(so%ln, so%ln, 1, MPI_complex16, MPI_sum, comm_glb, MPI_err)

    call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

    where (abs(real(so%e)) .lt. 5.D-14) so%e = cmplx(0.D0, aimag(so%e), p8)
    where (abs(aimag(so%e)) .lt. 5.D-14) so%e = cmplx(real(so%e), 0.D0, p8)

    return
  end procedure

  module procedure helm
    integer(i4) :: nrdim, npdim, nzdim
    integer(i4) :: nrc, npc, nzc, nzcu
    integer(i4), dimension(:), allocatable :: nrcs, m
    real(p8), dimension(:), allocatable :: ak

    call chop_glb_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

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

    call chop_glb_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

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
    call so%exchange(3,1) ! whole r-data resides locally

    so%ln = s%ln/alpha

    if ((so%loc_st(2).eq.0) .and. (so%loc_st(3).eq.0)) then
      so%e(1,1,1) = so%e(1,1,1) - 4.D0/3.D0/ell**2.D0/exp(tfm%lognorm(1,1))*so%LN
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

    call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

    where (abs(real(so%e)) .lt. 5.D-14) so%e = cmplx(0.D0, aimag(so%e), p8)
    where (abs(aimag(so%e)) .lt. 5.D-14) so%e = cmplx(real(so%e), 0.D0, p8)

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

    call chop_glb_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

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

    if (.not. (mod(power,2).eq.0) .and. (power.ge.4)) stop 'ihelmp: even power greater than or equal to 4'
    if (power .gt. 8) stop 'ihelmp: power must be less than or equal to 8 (supported power = 4, 6 or 8)'

    call chop_glb_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)

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
    call so%exchange(3,1) ! whole r-data resides locally

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
          do n = 1, np
            w(n,n) = w(n,n) + alpha
          enddo
          helmp_bnd = bandmat(w(:nn,:nn), suplen=power, sublen=power)
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
          do n = 1, np
            w(n,n) = w(n,n) + alpha
          enddo
          helmp_bnd = bandmat(w(:nn,:nn), suplen=power, sublen=power)
          deallocate( w )
          call lsolve(helmp_bnd, so%e(:nn,mm,kk))
        enddo
      endif
    enddo
    deallocate( helmp_bnd%e )
    deallocate( del2_bnd%e )

    call so%exchange(1, 3) ! resuming to typical FFF configuration (z-data resides locally)

    where (abs(real(so%e)) .lt. 5.D-14) so%e = cmplx(0.D0, aimag(so%e), p8)
    where (abs(aimag(so%e)) .lt. 5.D-14) so%e = cmplx(real(so%e), 0.D0, p8)

    deallocate( bl )
    deallocate( bl2 )

    return
  end procedure

  module procedure fefe
    type(scalar) :: svis
    real(p8) :: dt_

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
    else ! nonzero hypeviscosity -- diffusive operator is visc*del^2 + hypervisc*del^p
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
    else ! nonzero hypeviscosity -- diffusive operator is visc*del^2 + hypervisc*del^p
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
      else ! nonzero hypeviscosity -- diffusive operator is visc*del^2 + hypervisc*del^p
        svis_p = helmp(s_p, hyperpow, 0.D0, visc/(hypervisc*(-1.D0)**(hyperpow/2+1)), tfm)
        svis_p%e = svis_p%e * (hypervisc*(-1.D0)**(hyperpow/2+1))
      endif
    endif

    s%e = s%e + dt_ * (1.5D0*(s_rhs_nonlin%e + svis%e) - .5D0*(s_rhs_nonlin_p%e + svis_p%e))
    s%ln = s%ln + dt_ * (1.5D0*(s_rhs_nonlin%ln + svis%ln) - .5D0*(s_rhs_nonlin_p%ln + svis_p%ln))

    call svis%dealloc()
    call svis_p%dealloc()

  end procedure

  module procedure febe
    type(scalar) :: sh
    real(p8) :: dt_, a, b

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
    else ! nonzero hypeviscosity -- diffusive operator is visc*del^2 + hypervisc*del^p
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
    else ! nonzero hypeviscosity -- diffusive operator is visc*del^2 + hypervisc*del^p
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
        s = ihelm(s, a, tfm)
      endif
    else ! nonzero hypeviscosity -- diffusive operator is visc*del^2 + hypervisc*del^p
         ! F(K+1) = (1-DT/2*VISC*DEL2-DT/2*VISCP*DELP)^(-1) [F(K+1/2) + DT/2*LTERM(K)]
         ! = [(-DT/2*VISCP)*(DELP+(VISC/VISCP)*DEL2+(-2/(DT*VISCP)))]^(-1) [F(K+1/2) + DT/2*LTERM(K)]
         ! = IHELMP((-2/(DT*VISCP))*[F(K+1/2) + DT/2*LTERM(K)]), ALP=(-2/(DT*VISCP)), BET=(VISC/VISCP))
      a = -2.D0/(dt_*(hypervisc*(-1.D0)**(hyperpow/2+1)))
      b = visc/(hypervisc*(-1.D0)**(hyperpow/2+1))
      sh%ln = a*(sh%ln + dt_/2.D0*svis%ln)
      sh%e = a*(sh%e + dt_/2.D0*svis%e)
      s = ihelmp(sh, hyperpow, a, b, tfm)
    endif

    call sh%dealloc()
    call svis%dealloc()

  end procedure

  ! !> richardson extrapolation (2nd order) for initial bootstrapping
  ! interface richardson
  !   module subroutine richardson(s, s_rhs_linear, s_rhs_nonlin)
  !     implicit none
  !     class(scalar), intent(inout) :: s
  !   end subroutine
  ! end interface

  ! module procedure richardson
  !   type(scalar) :: s1, s2
  !   integer(i4) :: i
  !   real(p8) :: dto, ddt, ddth

  !   if (ni .ne. curr_n) stop 'richardson: bootstrapping unnecessary except for the initial step'
  !   dto = dt

  !   ddt = dt/2
  !   ddth = ddt/2

  !   dt = ddt ! subtime stepping
  !   call s1%init(s%glb_sz, s%axis_comm)
  !   s1 = s
  !   do i = 1,2
  !     call febe(s1, s_rhs_linear, s_rhs_nonlin)
  !     ! rhs must be updated that is problem-specific
  !   enddo

  !   dt = ddth ! half-subtime stepping
  !   call s2%init(s%glb_sz, s%axis_comm)
  !   s2 = s
  !   do i = 1,4
  !     call febe(s2, s_rhs_linear, s_rhs_nonlin)
  !     ! rhs must be updated that is problem-specific
  !   enddo

  !   c%e = 2.D0*c2%e - c1%e
  !   dt = dto ! recover original dt

  !   call s1%dealloc()
  !   call s2%dealloc()

  !   return
  ! end procedure

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

    where (abs(real(s%e)) .lt. 5.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    where (abs(aimag(s%e)) .lt. 5.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

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

    where (abs(real(s%e)) .lt. 5.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    where (abs(aimag(s%e)) .lt. 5.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

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

    where (abs(real(s%e)) .lt. 5.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    where (abs(aimag(s%e)) .lt. 5.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

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

    where (abs(real(s%e)) .lt. 5.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    where (abs(aimag(s%e)) .lt. 5.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

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

    where (abs(real(s%e)) .lt. 5.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    where (abs(aimag(s%e)) .lt. 5.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

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

    where (abs(real(s%e)) .lt. 5.D-14) s%e = cmplx(0.D0, aimag(s%e), p8)
    where (abs(aimag(s%e)) .lt. 5.D-14) s%e = cmplx(real(s%e), 0.D0, p8)

  end subroutine

! ======================================================================================================== !

  subroutine chop_glb_index(s, tfm, nrdim, npdim, nzdim, nrc, npc, nzc, nzcu, nrcs, m, ak)
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

    return
  end subroutine

end submodule