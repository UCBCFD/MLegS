program time_integration_first
!> This is a tutorial program to demonstrate the first-order built-in time integration scheme.
!> It solves for a simple 2-dimensional axissymmetric vorticity equation, governed by dw/dt = VISC*del2(w).
!> Initial condition is: w(r,t=0) = 1/(2*pi) * exp(-r^2 / 2), whose analytic solution in time is w(r,t) = 1/(2*pi+(1+2*VISC*t)) * exp( -r^2 / (2*(1+2*VISC*t)) ).
!> Further details are presented in Tutorial 05: Time Integration in [root_dir]/tutorials/05_time_integration.ipynb

  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_spectfm; use mlegs_scalar 

  implicit none

  integer(i4), dimension(3) :: glb_sz
  integer(i4) :: stepping_notice = 10 ! in order to suppress time stepping print-outs
  real(p8) :: origin_val, tsum = 0.D0 ! origin_val stores the vorticity value at origin in terms of time
  character(len=256) :: input_params_file = './input_tutorial.params' ! must be created in advance. Refer to the tutorial ipynb.
  logical :: exists

!> w stores the vorticity scalar field, while nlw stores the advection term, which in this tutorial is just zero.
  type(scalar) :: w, nlw

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
  call set_comm_grps(comm_glb, (/ nprocs_glb, 1 /)) ! 1D grid / slab decomposition

!> Start time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Program Started'
    call tic(); call print_real_time()
  endif

!!!............ Initialization for spectral transformation
  inquire( file=input_params_file, exist=exists )
  if (exists) then
    call read_input(input_params_file)
  else
    stop 'time_integration_first: This tutorial requires the preset input params. ' // &
         'Follow the instructions in [root_dir]/tutorials/05_time_integration.ipynb'
  endif

!> As the program runs transient simulations, simulation stepping and time variables are set
  call timestep_set()

!> Initialize the global transformation kit 'tfm' (see mlegs_spectfm)
  call tfm%init()

!> Store collocation point coordinate information
  if (rank_glb .eq. 0) call colloc_pts(tfm)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Initialization'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Declare a scalar field and initialize its distribution
!> Initialize the scalar
  glb_sz = (/ tfm%nrdim, tfm%npdim, tfm%nzdim /)
  call w%init(glb_sz, (/ 1, 0, 2 /)); call nlw%init(glb_sz, (/1, 0, 2/))
  w%space = 'PPP'; nlw%space='PPP'

!> Load an initial distribution of scalar, w(r,p,z,t=0) = 1/(2*pi) * exp(-r^2 / 2)
  call init_dist(w) ! this is a program-dependent subroutine; see below

!> Store the initial distribution data
  call msave(w, trim(adjustl(FLDDIR)) // 'sfld' // trim(adjustl(ntoa(curr_n,'(i0.6)'))) // '_' // w%space, &
             is_binary = .false., is_global = .true.)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Scalar Loaded with its Initial Distribution'
    write(*,*) 'w(r,p,z,t=0) = 1/(2*pi) * exp(-r^2 / 2)'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!> Open log files in case some data are logged:
  open(12, file=trim(adjustl(LOGDIR)) // 'sval_origin_dt2e-1.dat', status="replace", action="write")
  open(13, file=trim(adjustl(LOGDIR)) // 'sval_origin_dt4e-1.dat', status="replace", action="write")
  open(14, file=trim(adjustl(LOGDIR)) // 'sval_origin_dt1e-0.dat', status="replace", action="write")

!!!............ Time stepping
  do while (curr_n .lt. ni + totaln)
    curr_n = curr_n + 1
    if (curr_n - ni .eq. totaln) dt = totaltime-(totaln-1)*dt
    curr_t = curr_t + dt

    call nonlinear_rhs(w, nlw) ! compute the nonlinear term in rhs (in the current case nls is simply zero)
!> For time integration, all scalar inputs must be in the FFF space form.
!> If they are already in the FFF space form, the below transformation will do nothing, but just ensure that they are FFF
    call trans(w, 'FFF', tfm); call trans(nlw, 'FFF', tfm)

!> Run the first-order time stepping
    call febe(w, nlw, dt, tfm) ! forward euler - backward euler time integration by dt

    if (isfldsav) then
!> Make the field stored in the PPP format (for more straightforward visualization)
      call trans(w, 'PPP', tfm)
      if (mod(curr_n, fldsavintvl) .eq. 0) then
        call msave(w, trim(adjustl(FLDDIR)) // 'sfld' // trim(adjustl(ntoa(curr_n,'(i0.6)'))) // '_' // w%space, &
             is_binary = .false., is_global = .true.)
      endif
    endif

    if (islogsav) then
      if (mod(curr_n, logsavintvl) .eq. 0) then
        origin_val = origin_eval(w) ! this is a program-dependent subroutine to calculate volume integral of s in the entire domain.
        if (rank_glb .eq. 0) then
          write(12, 103) curr_n, curr_t, origin_val
        endif
      endif
    endif

!> Interim time record at every 10[==stepping_notice] steps
    if (rank_glb .eq. 0) then
      tsum = tsum + toc(); call tic();
      if (mod(curr_n, stepping_notice) .eq. 0) then
        write(*,102) ' time stepping # '// trim(adjustl(ntoa(curr_n,'(i10)')))//' '//'(t: ' &
                     //trim(adjustl(ntoa(curr_t,'(f11.6)')))//')'//' ' //repeat('.', 48), tsum/stepping_notice
        write(*,*) ''
        tsum = 0.D0
      endif
    endif
  enddo

!!!............ Re-do the simulation, with increasing dt of 2e-1 and 1e-0. Only the origin scalar data will be logged, for the sake of error comparison

!> Reset the time variables, and change the time stepping variables to adapt dt = 2.D-1
  dt = 4.D-1; totaltime = 1.D1; ni = 0; totaln = 25
  call timestep_set()
  logsavintvl = 1
  call trans(w, 'PPP', tfm)
  call mload(trim(adjustl(FLDDIR)) // 'sfld' // trim(adjustl(ntoa(curr_n,'(i0.6)'))) // '_' // w%space, w, &
             is_binary = .false., is_global = .true.)
  call trans(nlw, 'PPP', tfm); nlw%e = 0.D0

!> Do the time stepping w/ dt = 1.D-3
  do while (curr_n .lt. ni + totaln)
    curr_n = curr_n + 1
    if (curr_n - ni .eq. totaln) dt = totaltime-(totaln-1)*dt
    curr_t = curr_t + dt

    call nonlinear_rhs(w, nlw) ! compute the nonlinear term in rhs (in the current case nls is simply zero)
    call trans(w, 'FFF', tfm); call trans(nlw, 'FFF', tfm)
    call febe(w, nlw, dt, tfm) ! forward euler - backward euler time integration by dt

    if (islogsav) then
      if (mod(curr_n, logsavintvl) .eq. 0) then
        origin_val = origin_eval(w) ! this is a program-dependent subroutine to calculate volume integral of s in the entire domain.
        if (rank_glb .eq. 0) then
          write(13, 103) curr_n, curr_t, origin_val
        endif
      endif
    endif
  enddo

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Re-do the simulation with dt=0.4'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!> Reset the time varaibles, and change the time stepping variables to adapt dt = 1.D-0
  dt = 1.D-0; totaltime = 1.D1; ni = 0; totaln = 10
  call timestep_set()
  logsavintvl = 1
  call trans(w, 'PPP', tfm)
  call mload(trim(adjustl(FLDDIR)) // 'sfld' // trim(adjustl(ntoa(curr_n,'(i0.6)'))) // '_' // w%space, w, &
             is_binary = .false., is_global = .true.)
  call trans(nlw, 'PPP', tfm); nlw%e = 0.D0

!> Do the time stepping w/ dt = 1.D-2
  do while (curr_n .lt. ni + totaln)
    curr_n = curr_n + 1
    if (curr_n - ni .eq. totaln) dt = totaltime-(totaln-1)*dt
    curr_t = curr_t + dt

    call nonlinear_rhs(w, nlw) ! compute the nonlinear term in rhs (in the current case nls is simply zero)
    call trans(w, 'FFF', tfm); call trans(nlw, 'FFF', tfm)
    call febe(w, nlw, dt, tfm) ! forward euler - backward euler time integration by dt

    if (islogsav) then
      if (mod(curr_n, logsavintvl) .eq. 0) then
        origin_val = origin_eval(w) ! this is a program-dependent subroutine to calculate volume integral of s in the entire domain.
        if (rank_glb .eq. 0) then
          write(14, 103) curr_n, curr_t, origin_val
        endif
      endif
    endif
  enddo

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Re-do the simulation with dt=1.0'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Finalization 
  !> Finalize MPI
  call MPI_finalize(MPI_err)

!> Deallocate any allocated arrays
  call w%dealloc()
  call nlw%dealloc()

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Finalization'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!> Close log files
  close(12)
  close(13)
  close(14)

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

  subroutine init_dist(s) !> set an initial distribution of s.
    implicit none
    type(scalar), intent(inout) :: s
    complex(p8), dimension(:,:,:), allocatable :: glb
    integer(i4) :: i, j, k

    real(p8) :: xo, yo ! define the center of distribution
    real(p8) :: rr, ri ! distance from the center of dist.

    allocate( glb(s%glb_sz(1), s%glb_sz(2), s%glb_sz(3)) )
    glb = 0.D0

    if (s%space .ne. 'PPP') stop

    xo = .0D0 ! origin is the center
    yo = .0D0

    do i = 1, nr
      do j = 1, np/2
        do k = 1, nz
          rr = sqrt((tfm%r(i)*cos(tfm%p(2*j-1))-xo)**2.D0 & ! real part of the j-th entry corresponds to azim. angle p(2*j-1)
                   +(tfm%r(i)*sin(tfm%p(2*j-1))-yo)**2.D0 )
          ri = sqrt((tfm%r(i)*cos(tfm%p(2*j  ))-xo)**2.D0 & ! imag part of the j-th entry corresponds to azim. angle p(2*j)
                   +(tfm%r(i)*sin(tfm%p(2*j  ))-yo)**2.D0 )
          glb(i,j,k) = cmplx( exp(-(rr**2.D0)/2.D0) , &
                              exp(-(rr**2.D0)/2.D0) , P8) / (2.D0*pi)
        enddo ! w(r,p,z,t=0) = 1/(2*pi) * exp(-r^2 / 2)
      enddo
    enddo

    call s%disassemble(glb)
    deallocate( glb )

    return
  end subroutine

  subroutine nonlinear_rhs(s, nls) !> compute the nonlinear term in the equation's RHS. In the current example, it is purely zero.
    implicit none
    type(scalar), intent(in) :: s
    type(scalar), intent(inout) :: nls

    nls = s
    nls%e = 0.D0

    return
  end subroutine

  function origin_eval(s) result(itg)
    implicit none
    type(scalar) :: s
    complex(p8), dimension(:), allocatable :: calc
    real(p8) :: itg

    calc = calcat0(s, tfm)
    itg = calc(1)

  end function

end program