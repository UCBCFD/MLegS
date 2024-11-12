program backward_trans
!> This is a tutorial program to demonstrate how to obtain the physical information of a scalar when its spectral coefficients are given.
!> Further details are presented in Tutorial 02: Spectral Transformation in [root_dir]/tutorials/02_spectral_trans.ipynb

  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_spectfm; use mlegs_scalar 

  implicit none

  integer, dimension(3) :: glb_sz
  character(len=256) :: input_params_file = './input_2d.params'
  logical :: exists

!> First, we create a 'global' array, assign values, and distribute it into a local scalar
  complex(p8), dimension(:,:,:),allocatable :: array_glb
!> This 'local' scalar is distributed; except a single-core execution, it contains partial information of a global scalar field
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

!> Store collocation point coordinate information
  if (rank_glb .eq. 0) call colloc_pts(tfm)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Initialization'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Scalar setup
!> Allocate global scalar. Globally, a scalar must take storage of a dimension of (tfm%nrdim, tfm%npdim, tfm%nzdim). 
  glb_sz = (/ tfm%nrdim, tfm%npdim, tfm%nzdim /)
  allocate(array_glb(glb_sz(1), glb_sz(2), glb_sz(3)))
!> Now that the storage is assigned, let's initialize all entries and assign some values for testing...  
  array_glb = 0.D0
  array_glb(2, 2, 1) = 1.D0
  array_glb(2, 3, 1) = 1.D0

!> Let's distribute this global array into a local scalar variable defined in spectral space
!> By 'spectral space', it means that the scalar entries represent 'spectral' coefficients 
!> after two Fourier transformations in phi, z and one mapped Legendre transformation in r.
!> What does (/ 0, 1, 2 /) in the second parameter of init() indicate? 
!> First, s is initializaed while dim 1 (r) is complete (arg 0), meaning that all radial information resides in each local processor,
!> dim 2 (phi) is distributed by comm_grp(1) (arg 1), which in this case consists of the entire MPI workers (1D slab decomposition),
!> and dim 3 is distribued by comm_grp(2) (arg 2), which in this case means that no data distribution has been made in this dimension.
  call s%init(glb_sz, (/ 0, 1, 2 /))

!> 'FFF' space == scalar entries are fully spectral (after the 'triple' transformations)
!> 'PPP' space == scalar entries possess the fully physical information in the physical domain (r_i, phi_j, z_k)
!> Intermediate spaces, like 'FFP' and 'PFP', are also present, but are less meaningful.
  s%space = 'FFF'

!> Now the global array information is loaded to s via 'disassemblying' the array
  call s%disassemble(array_glb)

!> Before transforming the scalar, save the 'FFF' data for comparison purposes.
!> As is_binary is false, the data is stored in a ASCII-based readable format
!> As is_global is true, the data is globally collected before save 
!> If is_global is false, each processor stores locally distributed scalar data (proc 0 saves scalar.dat_0, proc 1 scalar.dat_1, ...) 
  call msave(s, trim(adjustl(LOGDIR)) // 'scalar_FFF_original.dat', is_binary = .false., is_global = .true.)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Scalar Setup & Disassembly'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Backward spectral transformation
!> Perform the backward transformation of s, from the 'FFF' space to the 'PPP' space
  call trans(s, 'PPP', tfm)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Backward Transformation'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!> Save the 'PPP' data for comparison purposes.
!> As is_binary is false, the data is stored in a ASCII-based readable format
  call msave(s, trim(adjustl(LOGDIR)) // 'scalar_PPP.dat', is_binary = .false., is_global = .true.)

!!!............ Finalization 
  !> Finalize MPI
  call MPI_finalize(MPI_err)

!> Deallocate any allocated arrays
  deallocate( array_glb )
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

contains

subroutine colloc_pts(tfm)
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

end program