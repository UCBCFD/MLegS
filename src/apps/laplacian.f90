program laplacian
!> This is a tutorial program to demonstrate how to obtain Laplacian of a scalar.
!> Further details are presented in Tutorial 03: Spectral Operation in [root_dir]/tutorials/03_operation.ipynb

  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_spectfm; use mlegs_scalar 

  implicit none

  integer, dimension(3) :: glb_sz
  character(len=256) :: input_params_file = './input_2d.params'
  logical :: exists

!> A distributed scalar 's' is declared, we will load the data generated from the previous tutorial, in '/output/scalar_FFF.dat'
  type(scalar) :: s

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
  inquire( file=input_params_file, exist=exists )
  if (exists) then
    call read_input(input_params_file)
  endif

!> Initialize the global transformation kit 'tfm' (see mlegs_spectfm)
  call tfm%init()

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Initialization'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Load a scalar
!> Initialize the scalar and load the data
  glb_sz = (/ tfm%nrdim, tfm%npdim, tfm%nzdim /)
  call s%init(glb_sz, (/ 0, 1, 2 /))
  s%space = 'FFF'

  inquire (file=trim(adjustl(LOGDIR)) // 'scalar_FFF.dat', exist=exists )
  if (exists) then
    call mload(trim(adjustl(LOGDIR)) // 'scalar_FFF.dat', s, is_binary = .false., is_global = .true.)
  else
    ! just in case somebody attempts to run this tutorial independently ...
    s%e = 0.D0
    if ((s%loc_st(2)+1 .le. 2).and.(2 .le. s%loc_st(2)+s%loc_sz(2))) s%e(2, 2-s%loc_st(2), 1) = 1.D0
    if ((s%loc_st(2)+1 .le. 3).and.(3 .le. s%loc_st(2)+s%loc_sz(2))) s%e(3, 3-s%loc_st(2), 1) = 1.D0
    call msave(s, trim(adjustl(LOGDIR)) // 'scalar_FFF.dat', is_binary = .false., is_global = .true.)
  endif

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Scalar Loading'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Laplacian
!> Obtain Laplacian of s and replace the data
  s = del2(s, tfm)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Laplacian Operation'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!> Save the Laplacian for comparison.
  call msave(s, trim(adjustl(LOGDIR)) // 'scalar_FFF_Laplacian.dat', is_binary = .false., is_global = .true.)

!> This step is for instructional purposes only -- get the physical form of the scalar to verify that it is 'actually' Laplacian of s.
  call trans(s, 'PPP' ,tfm)
  call msave(s, trim(adjustl(LOGDIR)) // 'scalar_PPP_Laplacian.dat', is_binary = .false., is_global = .true.)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Laplacian Data Record'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Finalization 
  !> Finalize MPI
  call MPI_finalize(MPI_err)

!> Deallocate any allocated arrays
  call s%dealloc()

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