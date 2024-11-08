program barebone_template
!> This is a 'barebone' template to run MLegS.
!> The program will call all MLegS modules and initialize the MPI system 
!> and the spectral transformation kit on a predefined kit class 'tfm'.
!> That is it. See further comments if one wants to understand what each code block does.

!> First of all, the program calls all MLegS modules (and MPI).
!> Unless necessary, it is recommended to use all these modules to be hassle-free.
  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_spectfm; use mlegs_scalar 

  implicit none

  character(len=256) :: input_params_file = './input.params'
  logical :: exists

!> This program requires no more program-bound variables as the MPI and kit varaibles
!> are stored in the module files (see mlegs_envir/mlegs_spectfm).
!> Below are some suggestions of program-bound variables for actual programs.
!> Iterative variables
  ! integer(i4) :: ii, jj, kk
!> 'Global' array where a locally distribued scalar data are collected
  ! complex(p8), dimension(:,:,:),allocatable :: array_glb
!> 'Local' scalar, usually storing physical states of a system to be simulated (e.g., temp, vel, pres)
  ! type(scalar) :: s

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

!> Set 2D processor grid
  call set_comm_grps(comm_glb) ! default: 2D grid / pencil (automatically balanced)
!> if one wants to do 1D grid / slab decompositions, comment above and uncomment the following.
! call set_comm_grps(comm_glb, (/ nprocs_glb, 1 /)) ! 1D grid / slab

!> Start time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Program Started'
    call tic(); call print_real_time()
  endif

!!!............ Initialization for spectral transformation
!> Read input paramters. If no input file exists, proceed with default values (see mlegs_base) 
  inquire( file=input_params_file, exist=exists )
  if (exists) then
    call read_input('./input.params')
  endif

!> Initialize the global transformation kit 'tfm' (see mlegs_spectfm)
  call tfm%init()

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Initialization'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

! ========================================================= !
! Now this is the part where customized MLegS code snippets !
! are placed. Define scalars, initialize the distributions, !
! and conduct some spatial and temporal operations.         !
! ========================================================= !
  if (rank_glb .eq. 0) then
    write(*,*) '*******************************************'
    write(*,*) 'This program is solely a barebone template.'
    write(*,*) 'No further computations shall be done.'
    write(*,*) '*******************************************'
    write(*,*) ''
  endif

!!!............ Finalization 
  !> Finalize MPI
  call MPI_finalize(MPI_err)

!> Deallocate any allocated arrays
!> In order to deallocate a global array,
  ! deallocate( A_glb )
!> In order to deallocate a scalar-type variable,
  ! call s%dealloc()

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

101 format(' elapsed time: ',F15.6, 'seconds')
end program