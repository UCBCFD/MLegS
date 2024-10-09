submodule (mlegs_spectfm) mlegs_spectfm_tbp
  implicit none

contains

  module procedure tfm_kit_set
    integer(i4) :: i, n

    if (this%is_set .eqv. .true.) then
      write(*,*) 'tfm_kit_set: [warning] this kit has been initialized already;' &
                 // ' are you really intentionally executing it?'
      this%is_set = .false.
      call this%dealloc()
      tfm_kit_counter = tfm_kit_counter - 1
    endif
    if (tfm_kit_counter .ge. 1) stop 'tfm_kit_set: only one MLegS transformation kit' &  
                                     // ' (either 1D, 2D or 3D) is permitted while an app run'

    select type(this)
      type is (tfm_kit_1d)
        if ((np .ne. 1) .or. (npchop .ne. 1) .or. (nz .ne. 1) .or. (nzchop .ne. 1)) then
          write(*,*) '1D kit is set... azimuthal, axial discretizations (np, nz) are forcefully set to 1'
        endif
        np = 1
        npchop = 1
        nz = 1
        nzchop = 1
      type is (tfm_kit_2d)
        if ((nz .ne. 1) .or. (nzchop .ne. 1)) then
          write(*,*) '2D kit is set... axial discretization (nz) are forcefully set to 1'
        endif
        nz = 1
        nzchop = 1
    end select

    if (.not.((nr .gt. 0) .and. (mod(nr, 2) .eq. 0))) stop 'tfm_kit_set: nr must be even'
    if (.not.((np .gt. 0) .and. ((np .eq. 1) .or. (mod(np, 2)) .eq. 0))) stop 'tfm_kit_set: np must be even'
    if (np .ne. 1) then
      n = np
      do while (mod(n, 2) == 0 .or. mod(n, 3) == 0 .or. mod(n, 5) == 0)
          if (mod(n, 2) == 0) n = n / 2
          if (mod(n, 3) == 0) n = n / 3
          if (mod(n, 5) == 0) n = n / 5
      enddo
      if (n .ne. 1) stop 'tfm_kit_set: np must only have factors of 2, 3 and 5'
    endif
    if (.not.((nz .gt. 0) .and. ((nz .eq. 1) .or. (mod(nz, 2)) .eq. 0))) stop 'tfm_kit_set: nz must be even'
    if (nz .ne. 1) then
      n = nz
      do while (mod(n, 2) == 0 .or. mod(n, 3) == 0 .or. mod(n, 5) == 0)
          if (mod(n, 2) == 0) n = n / 2
          if (mod(n, 3) == 0) n = n / 3
          if (mod(n, 5) == 0) n = n / 5
      enddo
      if (n .ne. 1) stop 'tfm_kit_set: nz must only have factors of 2, 3 and 5'
    endif
    if (nrchop .gt. nr) stop 'tfm_kit_set: nrchop must be smaller than or equal to nr'
    if (np .ne. 1) then
      if (npchop*2 .gt. np + 2) stop 'tfm_kit_set: npchop <= np/2 + 1 must be satisfied'
    endif
    if (nz .ne. 1) then
      if (nzchop*2 .gt. nz + 2) stop 'tfm_kit_set: nzchop <= nz/2 + 1 must be satisfied'
    endif

    allocate( this%x(nr), this%w(nr) )
    call find_gauleg_zeros( this%x, this%w)

    allocate( this%ln(nr) )
    this%ln = -log(1.D0 - this%x)

    select type(this)
      type is (tfm_kit_1d)
        allocate( this%r(nr) )
        this%r = ell * sqrt((1.D0+this%x)/(1.D0-this%x))

        allocate( this%chops(npchop) )
        this%chops = (/ (max(nrchop - i, 0), i = 0, npchop-1) /)

        this%nrdim = nr + max(3, hyperpow)
            
        allocate( this%lognorm(nrchop + 14, npchop) )
        this%lognorm = leg_lognorm( nrchop + 14, (/ 0 /) )
        
        allocate( this%pf(nr/2, nrchop+14, npchop) )
        this%pf = leg_tbl(this%x(:nr/2), nrchop+14, (/ 0 /), this%lognorm)
      type is (tfm_kit_2d)
        allocate( this%r(nr) )
        this%r = ell * sqrt((1.D0+this%x)/(1.D0-this%x))

        allocate( this%chops(npchop) )
        this%chops = (/ (max(nrchop - i, 0), i = 0, npchop-1) /)

        this%nrdim = nr + max(3, hyperpow)

        allocate( this%p(np+1) )
        this%p = (/ (2*pi/np*i, i = 0, np) /)
        
        allocate( this%m(npchop) )
        this%m = (/ (i, i = 0, npchop-1) /)
        
        this%npdim = np/2 + 1

        this%chopp = npchop
    
        allocate( this%lognorm(nrchop + 14, npchop) )
        this%lognorm = leg_lognorm( nrchop + 14, this%m )
        
        allocate( this%pf(nr/2, nrchop+14, npchop) )
        this%pf = leg_tbl(this%x(:nr/2), nrchop+14, this%m, this%lognorm)

      type is (tfm_kit_3d)
        allocate( this%r(nr) )
        this%r = ell * sqrt((1.D0+this%x)/(1.D0-this%x))

        allocate( this%chops(npchop) )
        this%chops = (/ (max(nrchop - i, 0), i = 0, npchop-1) /)

        this%nrdim = nr + max(3, hyperpow)

        allocate( this%p(np+1) )
        this%p = (/ (2*pi/np*i, i = 0, np) /)
        
        allocate( this%m(npchop) )
        this%m = (/ (i, i = 0, npchop-1) /)
        
        this%npdim = np/2 + 1

        allocate( this%z(nz) )
        this%z = (/ (zlen/i, i = 0, nz-1) /)

        allocate( this%ak(nz) )
        this%ak = (/ (2.D0*pi/zlen*i, i = -nz, -1) /)
        this%ak(1:nz/2+1) = (/ (2.D0*pi/zlen*i, i = 0, nz/2) /) 
        
        this%nzdim = nz

        this%chopzl = nzchop
        this%chopzu = nz - nzchop + 2
    
        allocate( this%lognorm(nrchop + 14, npchop) )
        this%lognorm = leg_lognorm( nrchop + 14, this%m )
        
        allocate( this%pf(nr/2, nrchop+14, npchop) )
        this%pf = leg_tbl(this%x(:nr/2), nrchop+14, this%m, this%lognorm)
    end select
    allocate( this%at0(nrchop) )
    this%at0 = reshape(leg_tbl((/ -1.D0 /), nrchop, (/ 0 /), this%lognorm(1:nrchop,1:1)), (/nrchop/))
    allocate( this%at1(nrchop) )
    this%at1 = reshape(leg_tbl((/ 1.D0 /), nrchop, (/ 0 /), this%lognorm(1:nrchop,1:1)), (/nrchop/))

    this%is_set = .true.
    tfm_kit_counter = tfm_kit_counter + 1
  end procedure

  module procedure tfm_kit_dealloc
      deallocate(this%x, this%w, this%ln, this%lognorm, this%pf, this%at0, this%at1)
      select type(this)
        type is (tfm_kit_1d)
          deallocate(this%r, this%chops)
        type is (tfm_kit_2d)
          deallocate(this%r, this%chops, this%p, this%m)
        type is (tfm_kit_3d)
          deallocate(this%r, this%chops, this%p, this%m, this%z, this%ak)
      end select
  end procedure

! ======================================================================================================== !
! VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV INTERNAL (PRIVATE) SUBROUTINES/FUNCTIONS VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV !
! ======================================================================================================== !

  subroutine find_gauleg_zeros(x, w)
    implicit none
    real(p8), dimension(:), intent(inout) :: w, x

    call gauss_legendre(-1.D0, 1.D0, x, w, size(x))
  end subroutine

! ======================================================================================================== !

  subroutine gauss_legendre(x1, x2, x, w, n)
    implicit none
    real(p8), intent(in) :: x1, x2 ! lower and upper bounds
    real(p8), dimension(:) :: x, w ! values and weights
    integer(i4) :: n ! total number of quadratures

    real(p8), parameter :: eps=1.D-15 ! machine-precision tolerance
    real(p8) :: xm, xl, p1, p2, p3, z, z1, pp
    integer(i4) :: m, i, j

    m = ceiling((n+1)/2.D0)
    xm = (x2+x1)/2.D0
    xl = (x2-x1)/2.D0
    do i = 1, m ! refinement using Newton's method
      z = cos(pi*(i-.25D0)/(n+.5D0))
      do while (abs(z-z1) .gt. eps)
        p1 = 1.d0
        p2 = 0.d0
        do j = 1,n
          p3 = p2
          p2 = p1
          p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
        enddo
        pp = n*(z*p1-p2)/(z*z-1)
        z1 = z
        z  = z1-p1/pp
      enddo
      x(i)     = xm-xl*z
      x(n+1-i) = xm+xl*z
      w(i)     = 2*xl/((1-z*z)*pp*pp)
      w(n+1-i) = w(i)
    enddo
  end subroutine

! ======================================================================================================== !

  function leg_lognorm(ndim, ms) result(lnrm)
    implicit none
    integer(i4), intent(in) :: ndim
    integer(i4), dimension(:), intent(in) :: ms
    real(p8), dimension(1:ndim, 1:size(ms)) :: lnrm

    integer(i4) :: nn, mm, n, m, me
    real(p8), dimension(:), allocatable :: wk_log

    me = maxval(abs(ms))
    allocate(wk_log(0:me))
    wk_log(0) = log(.5D0)
    do m = 1, me
      ! eqn. (13) in Matsushima & Marcus, 1997, without sqrt
      wk_log(m) = wk_log(m-1)+log(2.D0*m+1.D0)-log(2.D0*m*(2.D0*m-1.D0)**2)  
    enddo
    do mm = 1, size(ms)
      m = abs(ms(mm))
      lnrm(1,mm) = wk_log(m)
    enddo
    do mm = 1, size(ms)
      m = abs(ms(mm))
      do nn = 2, ndim
        n = m + (nn-1)
        lnrm(nn,mm) = lnrm(nn-1,mm) + log(2.D0*n+1.D0)-log(2.D0*n-1.D0) &
                                    + log(1.D0* (n-m))-log(1.D0*(n+m) )
      enddo
    enddo
    do mm = 1, size(ms)
      lnrm(:,mm) = .5D0*lnrm(:,mm)
    enddo
    deallocate(wk_log)
  end function

! ======================================================================================================== !

  function leg_tbl(x, ne, ms, lnrm) result(tbl)
    use MPI
    use fmvals_parallel; use fmzm_parallel
    implicit none
    real(p8), dimension(:), intent(in) :: x
    integer(i4), intent(in) :: ne
    integer(i4), dimension(:), intent(in) :: ms
    real(p8), dimension(:,:), intent(in) :: lnrm
    real(p8), dimension(size(x), ne, size(ms)) :: tbl

    integer(i4) :: m, me, mm, nn, n, xx
    type(fm), allocatable, dimension(:,:) :: fm_tbl
    integer(i4) :: m_chunk(2), m_sidx(2), m_eidx(2)
    real(p8), dimension(:,:,:), allocatable :: m_tbl

    tbl = 0.D0
    me = size(ms)
    if ((size(lnrm,1) .lt. ne) .or. (size(lnrm,2) .lt. me)) then
      write(*,*) 'leg_tbl: norm table size incompatible -- too small'
      stop
    endif
    allocate( fm_tbl(size(x), ne) )
    call m_decompose( size(x), me, m_chunk, m_sidx, m_eidx)
    if (m_chunk(1)*m_chunk(2) .gt. 0) then
      do mm = m_sidx(2), m_eidx(2)
        m = abs(ms(mm))
        fm_tbl(m_sidx(1):m_eidx(1),1) = (/ ((to_fm(-1.D0)**m)*sqrt(to_fm(1.D0)-to_fm(x(xx))**2)**m, &
                                        xx = m_sidx(1), m_eidx(1)) /)
        fm_tbl(m_sidx(1):m_eidx(1),2) = to_fm(x(m_sidx(1):m_eidx(1)))*fm_tbl(m_sidx(1):m_eidx(1),1)
        do nn = 3, ne
          n = m + nn - 1
          fm_tbl(m_sidx(1):m_eidx(1),nn) = to_fm(1.D0)/(n-m)*(fm_tbl(m_sidx(1):m_eidx(1),nn-1)* &
                                           to_fm(x(m_sidx(1):m_eidx(1)))-fm_tbl(m_sidx(1):m_eidx(1),nn-2)* &
                                           (n+m-1)/(2*n-1)/(2*n-3))
        enddo
        do nn = 1, ne
          n = m + nn - 1
          fm_tbl(m_sidx(1):m_eidx(1),nn) = fm_tbl(m_sidx(1):m_eidx(1),nn)*exp(to_fm(log_fact(n))+to_fm(lnrm(nn,mm)))
        enddo
        tbl(m_sidx(1):m_eidx(1),:,mm) = to_dp(fm_tbl(m_sidx(1):m_eidx(1),:))
      enddo
    endif
    allocate( m_tbl(size(x), ne, me) )
    call MPI_Allreduce(tbl, m_tbl, size(x)*ne*me, MPI_Real8, MPI_Sum, comm_glb, mpi_ierr)
    tbl = m_tbl
    deallocate( fm_tbl, m_tbl )
  end function

! ======================================================================================================== !

  function log_fact(m) result(lf)
    implicit none
    integer(i4), intent(in) :: m
    real(p8) :: lf

    lf = log_gamma(2*m+1.D0) - m*log(2.D0) - log_gamma(m+1.D0)
  end function

end submodule