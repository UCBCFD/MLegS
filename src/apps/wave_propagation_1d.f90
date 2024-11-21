program wave_propagation_1d
!> This program is to demonstrate how to simulate 1-dimensional radial wave propagation.
!> Solved is the axisymmetric system (no phi-dependent): d^2s/dt^2 = 1/r(d(r*ds/dr)/dr)
!> Initial condition: s(r) = exp(-(r/4)^2), ds/dt(r) = 0

  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_spectfm; use mlegs_scalar 

  implicit none

  integer(i4), dimension(3) :: glb_sz
  integer(i4) :: stepping_notice = 100000
  real(p8) :: origin_val, tsum = 0.D0
  character(len=256) :: input_params_file = './input_1d.params'
  logical :: exists

!> A distributed scalar 's' is declared, we will load the data generated from the previous tutorial, in '/output/scalar_FFF.dat'
  type(scalar) :: s, dsdt, del2s

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

!> When solving 1D problem, there's no reason for using more than 1 proc as this proc always has to contain all radial info.
!> Therefore, if nprocs_glb >= 1 the program terminates
  if (nprocs_glb .ne. 1) stop 'For running a 1D program, only a serial run using a single processor is allowed.'

!> Set 1D processor grid (in this example, nprocs_glb ==1, so this processor grid setup is void)
  call set_comm_grps(comm_glb, (/ nprocs_glb, 1 /))

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
! calculations, set NZ = NP = 1
!             192              1              1    
! # ------ NRCHOP ------- NPCHOP ------- NZCHOP ----------------------------------
! calculations, set NZCHOP = NPCHOP = 1
!             192              1              1    
! # --------- ELL --------- ZLEN ------ ZLENxPI ---(IF ZLENxPI==T, ZLEN=ZLEN*PI)--
! is irrelevant. Simply set it to any positive value.
!            1.D1           1.D0              F
! #
! !!! TIME STEPPING INFO !!!
! # ---------- DT ----------- TI --------- TOTT ----------- NI --------- TOTN ----
!           1.D-4           0.D0          5.D+1              0         500000
! #
! !!! FIELD PROPERTY INFO !!!
! # -------- VISC ----- HYPERPOW ---- HYPERVISC ----------------------------------
!            0.D0              0           0.D0
! #
! !!! FIELD SAVE INFO !!!
! # --------------------- FLDDIR -------------------------------------------------
!                   ./output/fld
! # ---- ISFLDSAV -- FLDSAVINTVL ---(IF ISFLDSAV!=T, FIELDS ARE NOT SAVED)--------
!               T           1000
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
  call s%init(glb_sz, (/0, 1, 2/)); call dsdt%init(glb_sz, (/0, 1, 2/)); call del2s%init(glb_sz, (/0, 1, 2/))
  s%space = 'PPP'; dsdt%space='PPP'; del2s%space = 'PPP' ! the space for phi and z has no meaning here.

!> Load an initial distribution of scalar, s = exp(-r^2/4), dsdt = 0
  call init_dist(s, dsdt) ! this is a program-dependent subroutine; see below

!> Store the initial distribution data
  call msave(s, trim(adjustl(FLDDIR)) // 'sfld' // trim(adjustl(ntoa(curr_n,'(i0.6)'))) // '_' // s%space, &
             is_binary = .false., is_global = .true.)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Scalar Loaded with its Initial Distribution'
    write(*,*) 's(r,t=0) = exp(-(r/4)^2)'
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

    call trans(s, 'FFF', tfm)
    del2s = del2(s, tfm) ! compute del^2(s)
!> For time integration
    call trans(s, 'FFF', tfm); call trans(dsdt, 'FFF', tfm); call trans(del2s, 'FFF', tfm)
    call fefe(s, dsdt, dt, tfm) ! forward euler
    call fefe(dsdt, del2s, dt, tfm) ! forward euler

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

!> Interim time record at every 100000[==stepping_notice] steps
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
  call dsdt%dealloc()
  call del2s%dealloc()

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

  subroutine init_dist(s, dsdt) !> set an initial distribution of s and dsdt. Here, s(r,t=0) = exp(-(r/4)^2), dsdt(r,t=0) = 0
    implicit none
    type(scalar), intent(inout) :: s, dsdt
    integer(i4) :: i

    if (s%space .ne. 'PPP') stop

    do i = 1, nr
      s%e(i,1,1) = exp(-((tfm%r(i)/4.D0)**2.D0))
    enddo

    dsdt%e = 0.D0

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
    real(p8) :: itg

    itg = 0.D0
    if ((s%loc_st(1).eq.0) .and. (s%loc_st(2).eq.0) .and. (s%loc_st(3).eq.0)) itg = s%e(1,1,1)
    call MPI_allreduce(itg, itg, 1, MPI_real8, MPI_sum, comm_glb, MPI_err)

  end function

end program