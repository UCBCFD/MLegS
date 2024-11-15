program scalar_transport_2d
!> This program is to demonstrate how to simulate the scalar transport process in a transient manner, given a frozen swirling flow field.
!> Governing equation: For a scalar field s(r,p,z,t), ds/dt = - u_phi*[(1/r)*ds/dphi] + VISC*del2(s) + source(r,p,z,t).
!> Initial condition: s(r(x,y),p(x,y),z,t=0) = 0
!> Convective Velocity: u = (0, u_phi(r), 0) = (e_phi)*(1- exp(-r^2))/r (Lamb-Oseen) 
!> External source term: source(r,p,z,t) = translate( max(1-(4*r)^8, 0), xc = 0.5, yc = 0.5) * (1 - cos(pi*t)) / 2

  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_spectfm; use mlegs_scalar 

  implicit none

  integer(i4), dimension(3) :: glb_sz
  integer(i4) :: stepping_notice = 1000
  real(p8) :: origin_val, tsum = 0.D0
  character(len=256) :: input_params_file = './input_2d.params'
  logical :: exists

!> A distributed scalar 's' is declared, we will load the data generated from the previous tutorial, in '/output/scalar_FFF.dat'
  type(scalar) :: s, nls

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
!> Read input paramters. If no input file exists, proceed with default values (see mlegs_base)

!> The default input parameter setup, for animations on the documentation site, is as follows:

!> ./input_2d.params

! !!! COMPUTATIONAL DOMAIN INFO !!!
! # ---------- NR ----------- NP ----------- NZ ---------------------------------- 
!              32             48              1    
! # ------ NRCHOP ------- NPCHOP ------- NZCHOP ---------------------------------- 
!              32             25              1    
! # --------- ELL --------- ZLEN ------ ZLENxPI ---(IF ZLENxPI==T, ZLEN=ZLEN*PI)-- 
!            1.D0           1.D0              F
! #
! !!! TIME STEPPING INFO !!!
! # ---------- DT ----------- TI --------- TOTT ----------- NI --------- TOTN ----
!           1.D-3           0.D0          1.D+1              0          10000
! #
! !!! FIELD PROPERTY INFO !!!
! # -------- VISC ----- HYPERPOW ---- HYPERVISC ----------------------------------
!           5.D-3              0           0.D0
! #
! !!! FIELD SAVE INFO !!!
! # --------------------- FLDDIR -------------------------------------------------
!                   ./output/fld
! # ---- ISFLDSAV -- FLDSAVINTVL ---(IF ISFLDSAV!=T, FIELDS ARE NOT SAVED)--------
!               T            100
! #
! !!! DATA LOGGING INFO !!!
! # --------------------- LOGDIR -------------------------------------------------
!                   ./output/dat
! # ---- ISLOGSAV -- LOGSAVINTVL ---(IF ISLOGSAV!=T, LOGS ARE NOT GENERATED)------
!               T            100
! /* ------------------------------ END OF INPUT ------------------------------ */

  inquire( file=input_params_file, exist=exists )
  if (exists) then
    call read_input(input_params_file)
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
  call s%init(glb_sz, (/ 1, 0, 2 /)); call nls%init(glb_sz, (/1, 0, 2/))
  s%space = 'PPP'; nls%space='PPP'

!> Load an initial distribution of scalar, s(r,p,z,t=0) = max(1-r^2, 0)
  call init_dist(s) ! this is a program-dependent subroutine; see below

!> Store the initial distribution data
  call msave(s, trim(adjustl(FLDDIR)) // 'sfld' // trim(adjustl(ntoa(curr_n,'(i0.6)'))) // '_' // s%space, &
             is_binary = .false., is_global = .true.)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Scalar Loaded with its Initial Distribution'
    write(*,*) 's(r,p,z,t=0) = max(1-r^2, 0)'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!> Open log files in case some data are logged:
  open(12, file=trim(adjustl(LOGDIR)) // 'sval_origin.dat', status="replace", action="write")

!!!............ Time stepping
  do while (curr_n .lt. ni + totaln)
    curr_n = curr_n + 1
    if (curr_n - ni .eq. totaln) dt = totaltime-(totaln-1)*dt
    curr_t = curr_t + dt

    call nonlinear_rhs(s, nls) ! compute the nonlinear term in rhs (in the current case nls is simply zero)
!> For time integration, all scalar inputs must be in the FFF space form.
!> If they are already in the FFF space form, the below transformation will do nothing, but just ensure that they are FFF
    call trans(s, 'FFF', tfm); call trans(nls, 'FFF', tfm)
    call febe(s, nls, dt, tfm) ! forward euler - backward euler time integration by dt

    if (isfldsav) then
!> Make the field stored in the PPP format (for more straightforward visualization)
      call trans(s, 'PPP', tfm)
      if (mod(curr_n, fldsavintvl) .eq. 0) then
        call msave(s, trim(adjustl(FLDDIR)) // 'sfld' // trim(adjustl(ntoa(curr_n,'(i0.6)'))) // '_' // s%space, &
             is_binary = .false., is_global = .true.)
      endif
    endif

    if (islogsav) then
      if (mod(curr_n, logsavintvl) .eq. 0) then
        origin_val = origin_eval(s) ! this is a program-dependent subroutine to calculate volume integral of s in the entire domain.
        if (rank_glb .eq. 0) then
          write(12, 103) curr_n, curr_t, origin_val
        endif
      endif
    endif

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
  enddo

!!!............ Finalization 
  !> Finalize MPI
  call MPI_finalize(MPI_err)

!> Deallocate any allocated arrays
  call s%dealloc()
  call nls%dealloc()

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Finalization'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!> Close log files
  close(12)

!> Finish time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Program Finished'
    write(*,101) toc(); call print_real_time()
  endif

101 format(' elapsed time: ',f15.6, 'seconds')
102 format(a48, 1x, f13.6, ' seconds per step')
103 format(i10, ' ', f10.6, ' ', 1pe13.4)

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
        if (size(tfm%p) .gt. 1) call msave(tfm%p,  trim(adjustl(LOGDIR)) // 'coords_p.dat')
      type is (tfm_kit_3d)
        call msave(tfm%r, trim(adjustl(LOGDIR)) // 'coords_r.dat')
        if (size(tfm%p) .gt. 1) call msave(tfm%p, trim(adjustl(LOGDIR)) // 'coords_p.dat')
        if (size(tfm%z) .gt. 1) call msave(tfm%z, trim(adjustl(LOGDIR)) // 'coords_z.dat')
    end select

  end subroutine

  subroutine init_dist(s) !> set an initial distribution of s. In the current example, s is initially zero.
    implicit none
    type(scalar), intent(inout) :: s

    s%e = 0.D0

    return
  end subroutine

  subroutine nonlinear_rhs(s, nls) !> compute the nonlinear term in the equation's RHS.
    implicit none
    type(scalar), intent(in) :: s
    type(scalar), intent(inout) :: nls

    complex(p8), dimension(:,:,:), allocatable :: src_glb
    integer(i4) :: mm, mv, i, j, k
    real(p8) :: xo, yo ! define the center of distribution
    real(p8) :: rr, ri ! distance from the center of dist.
    type(scalar) :: src

    nls = s ! nls takes s
    call trans(nls, 'FFF', tfm)
    do mm = 1, min(nls%loc_sz(2), tfm%chopp + nls%npchop_offset - nls%loc_st(2))
      mv = tfm%m(nls%loc_st(2) + mm)
      nls%e(:, mm, :) = nls%e(:, mm, :) * iu * mv ! now nls takes ds/dphi
    enddo
    call trans(nls, 'PPP', tfm)

    call nls%exchange(2, 1)
    do j = 1, nls%loc_sz(2)
      do k = 1, nls%loc_sz(3)
        nls%e(:nr,j,k) = - nls%e(:nr,j,k) * (1.D0 - exp(-tfm%r)**2.D0) / tfm%r**2.D0
                          ! now nls takes -u_phi * (1/r)*ds/phi
      enddo
    enddo
    call nls%exchange(1, 2)

    allocate( src_glb(s%glb_sz(1), s%glb_sz(2), s%glb_sz(3)) ) ! source term
    src_glb = 0.D0

    xo = .5D0 ! (.5, .5) is the center
    yo = .5D0

    do i = 1, nr
      do j = 1, np/2
        do k = 1, nz
          ! real part of the j-th entry corresponds to azim. angle p(2*j-1)
          rr = sqrt((tfm%r(i)*cos(tfm%p(2*j-1))-xo)**2.D0 & 
                   +(tfm%r(i)*sin(tfm%p(2*j-1))-yo)**2.D0 )
          ! imag part of the j-th entry corresponds to azim. angle p(2*j)
          ri = sqrt((tfm%r(i)*cos(tfm%p(2*j  ))-xo)**2.D0 & 
                   +(tfm%r(i)*sin(tfm%p(2*j  ))-yo)**2.D0 )
          src_glb(i,j,k) = cmplx(max(1.D0 - (4*rr)**8.D0, 0.D0) , &
                             max(1.D0 - (4*ri)**8.D0, 0.D0) , P8) / 2.D0 * (1.D0 - cos(pi*curr_t))
        enddo
      enddo
    enddo

    call src%init(nls%glb_sz, nls%axis_comm); src%space = 'PPP'
    call src%disassemble(src_glb)
    
    nls%e = - nls%e + src%e ! now mls takes - u_phi * ds/dphi [convection term] + source [external term]

    deallocate( src_glb )
    call src%dealloc() 

    return
  end subroutine

  function origin_eval(s) result(itg)
    implicit none
    type(scalar) :: s
    real(p8) :: itg

    itg = 0.D0
    if ((s%loc_st(1).eq.0) .and. (s%loc_st(2).eq.0) .and. (s%loc_st(3).eq.0)) itg = s%e(1,1,1)
    call MPI_allreduce(itg, itg, 1, MPI_real8, MPI_sum, comm_glb, MPI_err)

  end function

end program