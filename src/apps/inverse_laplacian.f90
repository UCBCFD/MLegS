program inverse_laplacian
!> This is a tutorial program to demonstrate how to obtain a scalar from its Laplacian of a scalar.
!> Running this program must be followed by its sequel program 'laplacian'
!> Further details are presented in Tutorial 03: Spectral Operation in [root_dir]/tutorials/03_operation.ipynb

  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_spectfm; use mlegs_scalar 

  implicit none

  integer, dimension(3) :: glb_sz
  character(len=256) :: input_params_file = './input_2d.params'
  logical :: exists

!> A distributed scalar 's' is declared, we will load the data generated from the sequel program, in '/output/scalar_FFF_Laplacian.dat'
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

  call mload(trim(adjustl(LOGDIR)) // 'scalar_FFF_Laplacian.dat', s, is_binary = .false., is_global = .true.)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Scalar Loading'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Inverse Laplacian
!> Obtain Laplacian of s and replace the data
  call idel2(s, tfm)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Laplacian Operation'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!> Save the Laplacian for comparison.
  call msave(s, trim(adjustl(LOGDIR)) // 'scalar_FFF_from_inv.dat', is_binary = .false., is_global = .true.)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Inverse Laplacian Data Record'
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