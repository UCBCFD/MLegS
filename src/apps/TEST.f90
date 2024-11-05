program test

  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_scalar; use mlegs_spectfm

  implicit none
  integer :: ii, jj, kk
  integer :: glb_sz(3), axis_comm(3)
  complex(p8), dimension(:,:,:),allocatable :: A_glb
  type(tfm_kit_3d) :: tfm
  logical :: tf

  ! Initialize MPI
  call MPI_init(MPI_err)
  if (MPI_err.ne.MPI_success) then
    write(*,*) 'MPI initialization failed.'
    call MPI_abort(MPI_comm_world, error_flag_comm, MPI_err)
  endif

  ! Read predefined input paramters (if no input file exists, proceed with default envir. vals)
  inquire( file='./input.params', exist=tf )
  if (tf) then
    call read_input('./input.params')
  endif

  ! Get the global communicator and start timing
  comm_glb = MPI_comm_world
  call MPI_comm_rank(comm_glb, rank_glb, MPI_err)

  if (rank_glb.eq.0) then
    write(*,*) 'program started'
    call tic(); call print_real_time()
  endif

  ! Initialize the transformation kit
  call tfm%alloc()

  if (rank_glb.eq.0) then
    write(*,*) 'tfm_kit initialized'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

  call MPI_finalize(MPI_err)

  if (rank_glb.eq.0) then
    write(*,*) 'program finished'
    write(*,101) toc(); call print_real_time()
  endif

  101 format(' elapsed time: ',F15.6, 'seconds')
end program