program test

  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_scalar; use mlegs_spectfm

  implicit none
  integer :: ii, jj, kk
  integer :: glb_sz(3), axis_comm(3)
  complex(p8), dimension(:,:,:),allocatable :: A_glb
  type(scalar) :: s
  logical :: tf

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

  !> Set 2D processor grid
  call set_comm_grps(comm_glb) ! default: 2D grid / pencil

  !!!............ Beginning of the program
  !> Start time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Program Started'
    call tic(); call print_real_time()
  endif

  !!!............ Initialization
  !> Read input paramters (if no input file exists, proceed with default values in the base file)
  inquire( file='./input.params', exist=tf )
  if (tf) then
    call read_input('./input.params')
  endif

  !> Initialize the global transformation kit
  call tfm%init()

  !> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Initialization'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

  !!!............ Scalar setup
  !> Allocate global scalar
  allocate(A_glb(tfm%nrdim, tfm%npdim, tfm%nzdim))
  A_glb = 0.D0
  A_glb(2, 3, 2) = 1.D0 + iu * 1.D0
  ! A_glb(:nr, 1, :nz) = 0

  !> Make a local scalar
  call s%init((/ tfm%nrdim, tfm%npdim, tfm%nzdim /), (/ 1, 0, 2 /))
  call s%disassemble(A_glb)
  call s%exchange(2, 3)
  s%space = 'FFF'
  if (rank_glb .eq. 0) call mcat(s%e, width=5)
  call trans(s, 'PPP', tfm)
  if (rank_glb .eq. 0) call mcat(s%e, width=5)
  call trans(s, 'FFF', tfm)
  if (rank_glb .eq. 0) call mcat(s%e, width=5)

  !!!............ Finalization 
  !> Finalize MPI
  call MPI_finalize(MPI_err)

  !> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Finalization'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

  !!!............ End of the program
  !> Finish time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Program Finished'
    write(*,101) toc(); call print_real_time()
  endif

  101 format(' elapsed time: ',F15.6, 'seconds')
end program