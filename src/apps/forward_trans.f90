program forward_trans
!> This is a tutorial program to demonstrate how to obtain the spectral coefficient information of a scalar when its physical information is given.
!> Running this program must be followed by its sequel program 'backward_trans'
!> Further details are presented in Tutorial 02: Spectral Transformation in [root_dir]/tutorials/02_spectral_trans.ipynb

  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_spectfm; use mlegs_scalar 

  implicit none

  integer, dimension(3) :: glb_sz
  character(len=256) :: input_params_file = './input_2d.params'
  logical :: exists

!> A distributed scalar 's' is declared, but now we load the 'backward-transformed' scalar output from 'backward_trans'
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
!> First, initialize a scalar. 
  glb_sz = (/ tfm%nrdim, tfm%npdim, tfm%nzdim /)
!> If a scalar is in the 'PPP' space, typically the azimuthal information resides locally without distribution.
!> Thus, (/ 1, 0, 2 /) is the right choice for the second parameter of init()...
  call s%init(glb_sz, (/ 1, 0, 2 /))
  s%space = 'PPP'

!> Load the 'PPP' data generated from 'backward_trans'.
!> As is_binary is false, the data is stored in a ASCII-based readable format
!> As is_global is true, the data is globally collected before save 
!> If is_global is false, each processor stores locally distributed scalar data (proc 0 saves scalar.dat_0, proc 1 scalar.dat_1, ...) 
  call mload(trim(adjustl(LOGDIR)) // 'scalar_PPP.dat', s, is_binary = .false., is_global = .true.)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Scalar Loading'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Forward spectral transformation
!> Perform the forward transformation of s, from the 'PPP' space to the 'FFF' space
  call trans(s, 'FFF', tfm)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Forward Transformation'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!> Save the 'FFF' data for comparison purposes.
!> For direct comparison with the original FFF information, is_global is set to true.
  call msave(s, trim(adjustl(LOGDIR)) // 'scalar_FFF.dat', is_binary = .false., is_global = .true.)

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