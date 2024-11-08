program test

  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_scalar; use mlegs_spectfm

  implicit none
  integer :: ii, jj, kk
  integer :: glb_sz(3), axis_comm(3)
  complex(p8), dimension(:,:,:),allocatable :: A_glb
  complex(p8), dimension(:), allocatable :: v
  type(scalar) :: s, st
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
  A_glb(10:11, 2:, :) = 1.D0

  !> Make a local scalar in FFF space
  call s%init((/ tfm%nrdim, tfm%npdim, tfm%nzdim /), (/ 0, 1, 2 /))
  call s%disassemble(A_glb)
  call s%exchange(1, 3)
  s%space = 'FFF'
  s%ln = -.5D0

  !> Perform forward-backward spectral transformations
  ! if (rank_glb .eq. 0) call mcat(s%e, width=5)
  call trans(s, 'PPP', tfm)
  ! if (rank_glb .eq. 0) call mcat(s%e, width=5)
  call trans(s, 'FFF', tfm)
  ! if (rank_glb .eq. 0) call mcat(s%e, width=5)
  
  ! !> Assemble local arrays into a global one at proc 0
  ! call s%exchange(3, 1)
  ! A_glb = s%assemble()
  ! if (rank_glb.eq.0)call mcat(A_glb(:,:10,1:2), width=5)
  ! call s%exchange(1, 3)

  ! !> Save and load scalar data files
  ! call msave(s, './scalar.dat', is_binary = .true., is_global = .true.)
  ! call mload('./scalar.dat', s, is_binary = .true., is_global = .true.)

  st = s
  do ii = 1, 10
  st = del2(st, tfm)
  st = idel2(st, tfm, preln=-.5D0)
  enddo

  !> Assemble local arrays into a global one at proc 0
  call st%exchange(3, 1)
  A_glb = st%assemble()
  if (rank_glb.eq.0)call mcat(A_glb(:,:,1:2), width=5)

  !> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Scalar Transformation'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

  !!!............ Finalization 
  !> Finalize MPI
  call MPI_finalize(MPI_err)

  !> deallocate any allocated arrays
  deallocate( A_glb )
  call s%dealloc()
  call st%dealloc()

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