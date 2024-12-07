program vortical_flow_3d
!> This program simulates a temporal evoluation of a 3-dimensional vortex with axial jet.
!> Incompressible N-S equations are solved using the Toroidal-Poloidal decomposition
!> V(r,p,z,t=0) = (1-exp(-r^2))/r e_r + exp(-r^2)/q e_z (embedded in a uniform z-flow; U = b e_z)
!> Here q = 1 and b = -0.5 (fixed).

  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_spectfm; use mlegs_scalar 

  implicit none

  integer(i4) :: i, j, k
  integer(i4), dimension(3) :: glb_sz
  integer(i4) :: stepping_notice = 1 ! in order to suppress time stepping print-outs
  integer(i4) :: farfield_treat = 100 ! in order to remove the advection term's aliasing error 
  real(p8) :: origin_val, tsum = 0.D0 ! origin_val stores the vorticity value at origin in terms of time
  character(len=256) :: input_params_file = './input.params' ! must be created in advance. Refer to the tutorial ipynb.
  logical :: exists

  integer(i4) :: rand_seed = 0

!> uz: background flow (z-direction only)
  type(scalar) :: uz
!> psi, chi: T-P scalars of the vortical flow field, V
  type(scalar) :: psi, nlpsi, chi, nlchi, psi_prev, nlpsi_prev, chi_prev, nlchi_prev
!> dummy scalar for Richard Extrapolation
  type(scalar) :: psi_rich, chi_rich

!!!............ Pre-initialization for MPI parallelism
!> Initialize MPI
  call MPI_init(MPI_err)
  if (MPI_err .ne. MPI_success) then
    write(*,*) 'MPI initialization failed.'
    call MPI_abort(MPI_comm_world, error_flag_comm, MPI_err)
  endif

!> Set the MPI global communicator
  comm_glb = MPI_comm_world
  call MPI_comm_rank(comm_glb, rank_glb, MPI_err)
  call MPI_comm_size(comm_glb, nprocs_glb, MPI_err)

!> Set 1D processor grid (in this example, we assume a 2D scalar in r, phi -- no dependence on z)
  call set_comm_grps(comm_glb) ! 1D grid / slab decomposition

!> Start time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Program Started'
    call tic(); call print_real_time()
  endif

!!!............ Initialization for spectral transformation
  inquire( file=input_params_file, exist=exists )
  if (exists) then
    call read_input(input_params_file)
  endif

!> As the program runs transient simulations, simulation stepping and time variables are set
  call timestep_set() 

!> Initialize the global transformation kit 'tfm' (see mlegs_spectfm)
  call tfm%init()

!> Store collocation point coordinate information
  if (rank_glb .eq. 0) call colloc_pts(tfm) ! program-dependent subroutine; see below

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Initialization'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Field declaration & initialization
!> Initialize the scalar
  glb_sz = (/ tfm%nrdim, tfm%npdim, tfm%nzdim /)

  call uz%init(glb_sz, (/ 1, 0, 2 /))
  uz%space = 'PPP'

  call psi%init(glb_sz, (/ 2, 1, 0 /)); call nlpsi%init(glb_sz, (/ 2, 1, 0 /))
  psi%space = 'FFF'; nlpsi%space='FFF'
  call psi_prev%init(glb_sz, (/ 2, 1, 0 /)); call nlpsi_prev%init(glb_sz, (/ 2, 1, 0 /))
  psi_prev%space = 'FFF'; nlpsi_prev%space='FFF'

  call chi%init(glb_sz, (/ 2, 1, 0 /)); call nlchi%init(glb_sz, (/ 2, 1, 0 /))
  chi%space = 'FFF'; nlchi%space='FFF'
  call chi_prev%init(glb_sz, (/ 2, 1, 0 /)); call nlchi_prev%init(glb_sz, (/ 2, 1, 0 /))
  chi_prev%space = 'FFF'; nlchi_prev%space='FFF'

  call psi_rich%init(glb_sz, (/ 2, 1, 0 /))
  psi_rich%space = 'FFF'
  call chi_rich%init(glb_sz, (/ 2, 1, 0 /))
  chi_rich%space = 'FFF'

!> Initialize all vector components or T-P scalars (for those whose init cond is known)
  call qvort_dist_tp(psi, chi, q = 1.D0, ran_noise = 0.D0) ! program-dependent subroutine; see below
  call uniform_z_fld(uz, b = -0.D0) ! program-dependent subroutine; see below

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Swirling Flow Field Declaration & Initialization'
    write(*,*) 'V(r,p,z,t=0) = (1-exp(-r^2))/r e_r + exp(-r^2)/q e_z'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!> Store the field information (for demonstration purposes, only the 'vorticity magnitude' field will be recorded)
  call save_vort_mag(psi, chi) ! program-dependent subroutine; see below

!!!............ Time stepping

!> Richardson Extrapolation bootstrapping
  call advection_rhs(psi, chi, nlpsi, nlchi, uz)
  psi_prev = psi; nlpsi_prev = nlpsi
  chi_prev = chi; nlchi_prev = nlchi

  psi_rich = psi
  chi_rich = chi

  do i = 1, 1
    call advection_rhs(psi_rich, chi_rich, nlpsi, nlchi, uz)
    call febe(psi_rich, nlpsi, dt, tfm) ! forward euler - backward euler time integration by dt
    call febe(chi_rich, nlchi, dt, tfm) ! forward euler - backward euler time integration by dt
  enddo

  do i = 1, 2
    call advection_rhs(psi, chi, nlpsi, nlchi, uz)
    call febe(psi, nlpsi, dt/2, tfm) ! forward euler - backward euler time integration by dt/2
    call febe(chi, nlchi, dt/2, tfm) ! forward euler - backward euler time integration by dt/2
  enddo

  psi%e = 2.D0/1.D0*psi%e - 1.D0/1.D0*psi_rich%e ! Richardson extrapolation
  chi%e = 2.D0/1.D0*chi%e - 1.D0/1.D0*chi_rich%e ! Richardson extrapolation

  call psi_rich%dealloc()
  call chi_rich%dealloc()

  curr_n = 1; curr_t = dt

!> Finish time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Richardson Extrapolation (time stepping # 1)'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!> Store the first time step data
  call save_vort_mag(psi, chi) ! program-dependent subroutine; see below

!> Main time stepping begins
  do while (curr_n .lt. ni + totaln)
    curr_n = curr_n + 1
    if (curr_n - ni .eq. totaln) dt = totaltime-(totaln-1)*dt
    curr_t = curr_t + dt

!> For time integration, all scalar inputs must be in the FFF space form.
!> If they are already in the FFF space form, the below transformation will do nothing, but just ensure that they are FFF
    call trans(psi, 'FFF', tfm); call trans(chi, 'FFF', tfm)
    call trans(nlpsi, 'FFF', tfm); call trans(nlchi, 'FFF', tfm)
    call trans(psi_prev, 'FFF', tfm); call trans(chi_prev, 'FFF', tfm)
    call trans(nlpsi_prev, 'FFF', tfm); call trans(nlchi_prev, 'FFF', tfm)

!> Run the second-order time stepping
    call abcn(psi, psi_prev, nlpsi, nlpsi_prev, dt, tfm) ! adams bashforth - crank nicolson time integration by dt
    call abcn(chi, chi_prev, nlchi, nlchi_prev, dt, tfm) ! adams bashforth - crank nicolson time integration by dt
    call advection_rhs(psi, chi, nlpsi, nlchi, uz) ! compute the nonlinear term in rhs

!> Interim time record at every 100[==stepping_notice] steps
    if (rank_glb .eq. 0) then
      tsum = tsum + toc(); call tic();
      if (mod(curr_n, stepping_notice) .eq. 0) then
        write(*,102) ' time stepping # '// trim(adjustl(ntoa(curr_n,'(i10)')))//' '//'(t: ' &
                     //trim(adjustl(ntoa(curr_t,'(f11.6)')))//')'//' ' //repeat('.', 48), tsum/stepping_notice
        write(*,*) ''
        tsum = 0.D0
      endif
    endif

    if (mod(curr_n, farfield_treat) .eq. 0) then
      call fftreat(psi, tfm) ! smooth out the far field to cope with the advection term aliasing error 
      call fftreat(chi, tfm)
    endif

    if (isfldsav) then
      if (mod(curr_n, fldsavintvl) .eq. 0) then
        call save_vort_mag(psi, chi)
      endif
    endif

    if (islogsav) then ! currently islogsav == .FALSE.
!> Record the log files
      if (mod(curr_n, fldsavintvl) .eq. 0) then
        ! For customization -- create and add your own diagnositic logs here
      endif
    endif
  enddo

!!!............ Finalization 
  !> Finalize MPI
  call MPI_finalize(MPI_err)

!> Deallocate any allocated arrays
  call uz%dealloc()

  call psi%dealloc(); call nlpsi%dealloc()
  call psi_prev%dealloc(); call nlpsi_prev%dealloc()
  call chi%dealloc(); call nlchi%dealloc()
  call chi_prev%dealloc(); call nlchi_prev%dealloc()

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Finalization'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!> Finish time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Program Finished'
    write(*,101) toc(); call print_real_time()
  endif

101 format(' elapsed time: ',f15.6, 'seconds')
102 format(a48, 1x, f13.6, ' seconds per step')
103 format(i10, ' ', f10.6, ' ', 1pe23.14)

contains
!> program-dependent subroutines/functions

  subroutine colloc_pts(tfm) !> get the collocation points
    implicit none
    class(tfm_kit) :: tfm

    select type(tfm)
      type is (tfm_kit_1d)
        call msave(tfm%r, trim(adjustl(LOGDIR)) // 'coords_r.dat')
      type is (tfm_kit_2d)
        call msave(tfm%r, trim(adjustl(LOGDIR)) // 'coords_r.dat')
        if (np .gt. 1) call msave(tfm%p,  trim(adjustl(LOGDIR)) // 'coords_p.dat')
      type is (tfm_kit_3d)
        call msave(tfm%r, trim(adjustl(LOGDIR)) // 'coords_r.dat')
        if (np .gt. 1) call msave(tfm%p, trim(adjustl(LOGDIR)) // 'coords_p.dat')
        if (nz .gt. 1) call msave(tfm%z, trim(adjustl(LOGDIR)) // 'coords_z.dat')
    end select

  end subroutine

  subroutine qvort_dist_tp(psi, chi, q, ran_noise)
    implicit none
    type(scalar), intent(inout) :: psi, chi
    real(p8), intent(in) :: q, ran_noise

    complex(p8), dimension(:,:,:), allocatable :: glb
    integer(i4) :: i, j, k

    real(p8) :: xo, yo ! define the center of distribution
    real(p8) :: rr, ri ! distance from the center of dist.

    if (psi%space(1:3) .ne. 'FFF') stop
    if (chi%space(1:3) .ne. 'FFF') stop

    allocate( glb(psi%glb_sz(1), psi%glb_sz(2), psi%glb_sz(3)) )
    glb = 0.D0

    call trans(psi, 'PPP', tfm)
    call trans(chi, 'PPP', tfm)

    do i = 1, nr
      do j = 1, np/2
        do k = 1, nz
          xo = 1.5D0
          yo = 0.D0
          rr = sqrt((tfm%r(i)*cos(tfm%p(2*j-1))-xo)**2.D0 & ! real part of the j-th entry corresponds to azim. angle p(2*j-1)
                   +(tfm%r(i)*sin(tfm%p(2*j-1))-yo)**2.D0 )
          ri = sqrt((tfm%r(i)*cos(tfm%p(2*j  ))-xo)**2.D0 & ! imag part of the j-th entry corresponds to azim. angle p(2*j)
                   +(tfm%r(i)*sin(tfm%p(2*j  ))-yo)**2.D0 )
          glb(i,j,k) = cmplx( -exp(-(rr**2.D0))*2.D0/(1.D0-tfm%x(i))**2.D0, &
                              -exp(-(ri**2.D0))*2.D0/(1.D0-tfm%x(i))**2.D0, P8)

          xo = -1.5D0
          yo = 0.D0
          rr = sqrt((tfm%r(i)*cos(tfm%p(2*j-1))-xo)**2.D0 & ! real part of the j-th entry corresponds to azim. angle p(2*j-1)
                   +(tfm%r(i)*sin(tfm%p(2*j-1))-yo)**2.D0 )
          ri = sqrt((tfm%r(i)*cos(tfm%p(2*j  ))-xo)**2.D0 & ! imag part of the j-th entry corresponds to azim. angle p(2*j)
                   +(tfm%r(i)*sin(tfm%p(2*j  ))-yo)**2.D0 )
          glb(i,j,k) = glb(i,j,k) + cmplx( -exp(-(rr**2.D0))*2.D0/(1.D0-tfm%x(i))**2.D0, &
                                           -exp(-(ri**2.D0))*2.D0/(1.D0-tfm%x(i))**2.D0, P8)

          if (rr.lt.ell .or. ri.lt.ell) glb(i,j,k) = glb(i,j,k) + cmplx(rand(),rand(),P8)*ran_noise
        enddo
      enddo
    enddo
    call psi%disassemble(glb)
    
    glb = 0.D0
    call trans(psi, 'FFF', tfm);
    call idelsqp(psi, tfm)
    call zeroat1(psi, tfm)

    do i = 1, nr
      do j = 1, np/2
        do k = 1, nz
          xo = 1.5D0
          yo = 0.D0
          rr = sqrt((tfm%r(i)*cos(tfm%p(2*j-1))-xo)**2.D0 & ! real part of the j-th entry corresponds to azim. angle p(2*j-1)
                   +(tfm%r(i)*sin(tfm%p(2*j-1))-yo)**2.D0 )
          ri = sqrt((tfm%r(i)*cos(tfm%p(2*j  ))-xo)**2.D0 & ! imag part of the j-th entry corresponds to azim. angle p(2*j)
                   +(tfm%r(i)*sin(tfm%p(2*j  ))-yo)**2.D0 )
          glb(i,j,k) = cmplx( -exp(-(rr**2.D0))/q/(1.D0-tfm%x(i))**2.D0, &
                              -exp(-(ri**2.D0))/q/(1.D0-tfm%x(i))**2.D0, P8)

          xo = -1.5D0
          yo = 0.D0
          rr = sqrt((tfm%r(i)*cos(tfm%p(2*j-1))-xo)**2.D0 & ! real part of the j-th entry corresponds to azim. angle p(2*j-1)
                   +(tfm%r(i)*sin(tfm%p(2*j-1))-yo)**2.D0 )
          ri = sqrt((tfm%r(i)*cos(tfm%p(2*j  ))-xo)**2.D0 & ! imag part of the j-th entry corresponds to azim. angle p(2*j)
                   +(tfm%r(i)*sin(tfm%p(2*j  ))-yo)**2.D0 )
          glb(i,j,k) = glb(i,j,k) + cmplx( -exp(-(rr**2.D0))/q/(1.D0-tfm%x(i))**2.D0, &
                                           -exp(-(ri**2.D0))/q/(1.D0-tfm%x(i))**2.D0, P8)

          if (rr.lt.ell .or. ri.lt.ell) glb(i,j,k) = glb(i,j,k) + cmplx(rand(),rand(),P8)*ran_noise
        enddo
      enddo
    enddo
    call chi%disassemble(glb)
    
    glb = 0.D0
    call trans(chi, 'FFF', tfm)
    call idelsqp(chi, tfm)
    call zeroat1(chi, tfm)

    deallocate( glb )    

  end subroutine

  subroutine uniform_z_fld(uz, b)
    implicit none
    type(scalar), intent(inout) :: uz
    real(p8) :: b

    complex(p8), dimension(:,:,:), allocatable :: glb

    if (uz%space(1:3) .ne. 'PPP') stop

    allocate( glb(uz%glb_sz(1), uz%glb_sz(2), uz%glb_sz(3)) )
    glb = 0.D0

    do i = 1, nr
      do j = 1, np/2
        do k = 1, nz
          glb(i,j,k) = cmplx(b, b, p8)
        enddo
      enddo
    enddo
    call uz%disassemble(glb)

    deallocate( glb )

  end subroutine

  subroutine advection_rhs(psi, chi, nlpsi, nlchi, uz) !> compute the TP-projected advection term 
    implicit none
    type(scalar), intent(in) :: psi, chi, uz
    type(scalar), intent(inout) :: nlpsi, nlchi

    integer(i4), dimension(3) :: glb_sz
    type(scalar) :: vr, vp, vz, or, op, oz

    if (psi%space(1:3) .ne. 'FFF') stop
    if (chi%space(1:3) .ne. 'FFF') stop
    if (nlpsi%space(1:3) .ne. 'FFF') stop
    if (nlchi%space(1:3) .ne. 'FFF') stop
    if (uz%space(1:3) .ne. 'PPP') stop

    glb_sz = (/ tfm%nrdim, tfm%npdim, tfm%nzdim /)
    call vr%init(glb_sz, (/ 1, 0, 2 /)); call or%init(glb_sz, (/ 1, 0, 2 /))
    vr%space = 'PPP'; or%space = 'PPP'
    call vp%init(glb_sz, (/ 1, 0, 2 /)); call op%init(glb_sz, (/ 1, 0, 2 /))
    vp%space = 'PPP'; op%space = 'PPP'
    call vz%init(glb_sz, (/ 1, 0, 2 /)); call oz%init(glb_sz, (/ 1, 0, 2 /))
    vz%space = 'PPP'; oz%space = 'PPP'

    call tp2vec(psi, chi, vr, vp, vz, tfm) ! (psi, chi) -> (vr, vp, vz)
    vz%e = vz%e + uz%e ! adding uz to vz (vr, vp, vz) + (0, 0, uz)
    
    call tp2curlvec(psi, chi, or, op, oz, tfm) ! (psi, chi) -> (curl(V)r, curl(V)p, curl(V)z)
    call vecprod(vr, vp, vz, or, op, oz, tfm) ! -[curl(V) times (U+V)]
    call vec2tp(vr, vp, vz, nlpsi, nlchi, tfm) ! TP projection of -[curl(V) times (U+V)]

    call vr%dealloc(); call or%dealloc()
    call vp%dealloc(); call op%dealloc()
    call vz%dealloc(); call oz%dealloc()

    return
  end subroutine

  subroutine save_vort_mag(psi, chi)
    implicit none
    type(scalar), intent(in) :: psi, chi

    integer(i4), dimension(3) :: glb_sz
    type(scalar) :: or, op, oz, vormag

    if (psi%space(1:3) .ne. 'FFF') stop
    if (chi%space(1:3) .ne. 'FFF') stop

    glb_sz = (/ tfm%nrdim, tfm%npdim, tfm%nzdim /)
    call or%init(glb_sz, (/ 1, 0, 2 /))
    or%space = 'PPP'
    call op%init(glb_sz, (/ 1, 0, 2 /))
    op%space = 'PPP'
    call oz%init(glb_sz, (/ 1, 0, 2 /))
    oz%space = 'PPP'
    call vormag%init(glb_sz, (/ 1, 0, 2 /))
    vormag%space = 'PPP'

    call tp2curlvec(psi, chi, or, op, oz, tfm)

    vormag%e = 0.D0
    do i = 1, vormag%loc_sz(1)
      do j = 1, vormag%loc_sz(2)
        do k = 1, vormag%loc_sz(3)
          vormag%e(i, j, k) = cmplx( (real(or%e(i,j,k))**2.D0 + & ! real part of the j-th entry corresponds to azim. angle p(2*j-1)
                                      real(op%e(i,j,k))**2.D0 + &
                                      real(oz%e(i,j,k))**2.D0)**.5D0, & 
                                      (aimag(or%e(i,j,k))**2.D0 + & ! imag part of the j-th entry corresponds to azim. angle p(2*j  )
                                      aimag(op%e(i,j,k))**2.D0 + &
                                      aimag(oz%e(i,j,k))**2.D0)**.5D0, p8)
        enddo
      enddo
    enddo

    call msave(vormag, trim(adjustl(FLDDIR)) // 'vormagfld' // trim(adjustl(ntoa(curr_n,'(i0.6)'))) // &
         '_' // vormag%space, is_binary = .false., is_global = .true.)

    !> Interim time record
      if (rank_glb .eq. 0) then
        write(*,*) 'Field Data Saved at Time Step #' // trim(adjustl(ntoa(curr_n,'(i0.6)')))
        write(*,101) toc(); call tic()
        write(*,*) ''
      endif

    call or%dealloc()
    call op%dealloc()
    call oz%dealloc()
    call vormag%dealloc()

101 format(' elapsed time: ',f15.6, 'seconds')
    return
  end subroutine

  function rand() result(r) ! rand from -1 ~ 1 (pseudorandom seed)
    implicit none
    real(p8) :: r

    integer(i4) :: i
    integer(i4), dimension(:), allocatable :: seed

    call random_seed(size = i)

    allocate( seed(i) )

    rand_seed = rand_seed + 1632291
    seed = rand_seed

    call random_seed(put = seed)
    call random_number(r) ! ~ uniformRand [0, 1]
    r = r * 2.D0 - 1.D0 ! uniformRand [-1, 1]

    deallocate( seed )
    return
  end function


end program