program vecfld_reconstruction
!> This is a tutorial program to demonstrate how to reconstruct vector components from given toroidal and poloidal scalars (streamfunctions). 
!> A vector field, in a component-wise form, is described as: (vr, vp, vz) = ( 0, (1-exp(-r^2))/r, exp(-r^2)/q ), where q = 1
!> The corresponding toroidal and poloidal scalars are calculated from the formula of: psi = -(idel2h)[curl(V)_z], chi = -(idel2h)[V_z].
!> psi and chi are calculated in a program-dependent subroutine, qvort_dist() (see below)
!> Further details are presented in Tutorial 06: Vector Field Manipulation in [root_dir]/tutorials/06_vector_field.ipynb

  use MPI
  use mlegs_envir; use mlegs_base; use mlegs_misc
  use mlegs_bndmat; use mlegs_genmat
  use mlegs_spectfm; use mlegs_scalar 

  implicit none

  integer(i4), dimension(3) :: glb_sz
  character(len=256) :: input_params_file = './input_tutorial.params'
  logical :: exists

  integer(i4) :: ii, jj, kk
  complex(p8), dimension(:,:,:),allocatable :: A_glb
  type(scalar) :: psi, chi, vr, vp, vz

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

! !> Set 2D processor grid
  call set_comm_grps(comm_glb) ! default: 2D grid / pencil (automatically balanced)
! !> if one wants to do 1D grid / slab decompositions, comment above and uncomment the following.
!   call set_comm_grps(comm_glb, (/ nprocs_glb, 1 /)) ! 1D grid / slab

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

!> Initialize the global transformation kit
  call tfm%init()
  if ((nr .ne. 200) .or. (np .ne. 128) .or. (nz .ne. 8)) stop 'vecfld_reconstruction: invalid input params.'

!> Store collocation point coordinate information
  if (rank_glb .eq. 0) call colloc_pts(tfm)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Initialization'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Declare a scalar field and initialize its distribution
!> Initialize the vector components
  glb_sz = (/ tfm%nrdim, tfm%npdim, tfm%nzdim /)
  call vr%init(glb_sz, (/ 1, 0, 2 /)); call vp%init(glb_sz, (/1, 0, 2/)); call vz%init(glb_sz, (/1, 0, 2/))
  vr%space = 'PPP'; vp%space = 'PPP'; vz%space = 'PPP'

!> Initialize the TP scalars
  call psi%init(glb_sz, (/ 2, 1, 0 /)); call chi%init(glb_sz, (/2, 1, 0/))
  psi%space = 'FFF'; chi%space = 'FFF'

!> Load vector field information to psi and chi
  call qvort_dist(psi, chi, q=1.D0) ! this is a program-dependent subroutine; see below

!> Store psi and chi
  call msave(psi, trim(adjustl(FLDDIR)) // 'psi' // trim(adjustl(ntoa(ni,'(i0.6)'))) // '_' // &
             psi%space, is_binary = .false., is_global = .true.)
  call msave(chi, trim(adjustl(FLDDIR)) // 'chi' // trim(adjustl(ntoa(ni,'(i0.6)'))) // '_' // &
             chi%space, is_binary = .false., is_global = .true.)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'TP Scalar Calculation & Loading'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Perform the vector field component reconstruction
!> Run tp2vec()
  call tp2vec(psi, chi, vr, vp, vz, tfm)

!> Store 3 vector components
  call msave(vr, trim(adjustl(FLDDIR)) // 'vr' // trim(adjustl(ntoa(ni,'(i0.6)'))) // '_' // &
             vr%space, is_binary = .false., is_global = .true.)
  call msave(vp, trim(adjustl(FLDDIR)) // 'vp' // trim(adjustl(ntoa(ni,'(i0.6)'))) // '_' // &
             vp%space, is_binary = .false., is_global = .true.)
  call msave(vz, trim(adjustl(FLDDIR)) // 'vz' // trim(adjustl(ntoa(ni,'(i0.6)'))) // '_' // &
             vz%space, is_binary = .false., is_global = .true.)

!> Interim time record
  if (rank_glb .eq. 0) then
    write(*,*) 'Vector Field Reconstruction'
    write(*,101) toc(); call tic()
    write(*,*) ''
  endif

!!!............ Finalization 
  !> Finalize MPI
  call MPI_finalize(MPI_err)

!> Deallocate any allocated arrays
  call psi%dealloc()
  call chi%dealloc()
  call vr%dealloc()
  call vp%dealloc()
  call vz%dealloc()

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

  subroutine colloc_pts(tfm) !> get the collocation points
    implicit none
    class(tfm_kit) :: tfm

    select type(tfm)
      type is (tfm_kit_1d)
        call msave(tfm%r, trim(adjustl(LOGDIR)) // 'coords_r.dat')
      type is (tfm_kit_2d)
        call msave(tfm%r, trim(adjustl(LOGDIR)) // 'coords_r.dat')
        if (np .gt. 1) call msave(tfm%p,  trim(adjustl(LOGDIR)) // 'coords_p.dat')
      type is (tfm_kit_3d)
        call msave(tfm%r, trim(adjustl(LOGDIR)) // 'coords_r.dat')
        if (np .gt. 1) call msave(tfm%p, trim(adjustl(LOGDIR)) // 'coords_p.dat')
        if (nz .gt. 1) call msave(tfm%z, trim(adjustl(LOGDIR)) // 'coords_z.dat')
    end select

  end subroutine

  subroutine qvort_dist(psi, chi, q) !> set an initial distribution of s. In the current example, s(r,p,z,t=0) = translate( max(1-r^8, 0), xc = 0.5, yc =0)
    implicit none
    type(scalar), intent(inout) :: psi, chi
    real(p8) :: q
    complex(p8), dimension(:,:,:), allocatable :: glb
    integer(i4) :: i, j, k

    real(p8) :: xo, yo ! define the center of distribution
    real(p8) :: rr, ri ! distance from the center of dist.

    allocate( glb(psi%glb_sz(1), psi%glb_sz(2), psi%glb_sz(3)) )
    glb = 0.D0

    xo = 1.D0 ! (1., 1.) is the center
    yo = 1.D0

    call trans(psi, 'PPP', tfm)
    call trans(chi, 'PPP', tfm)

    do i = 1, nr
      do j = 1, np/2
        do k = 1, nz
          rr = sqrt((tfm%r(i)*cos(tfm%p(2*j-1))-xo)**2.D0 & ! real part of the j-th entry corresponds to azim. angle p(2*j-1)
                   +(tfm%r(i)*sin(tfm%p(2*j-1))-yo)**2.D0 )
          ri = sqrt((tfm%r(i)*cos(tfm%p(2*j  ))-xo)**2.D0 & ! imag part of the j-th entry corresponds to azim. angle p(2*j)
                   +(tfm%r(i)*sin(tfm%p(2*j  ))-yo)**2.D0 )
          glb(i,j,k) = cmplx( -exp(-(rr**2.D0))*2.D0/(1.D0-tfm%x(i))**2.D0, &
                              -exp(-(ri**2.D0))*2.D0/(1.D0-tfm%x(i))**2.D0, P8)
        enddo
      enddo
    enddo
    call psi%disassemble(glb)
    
    glb = 0.D0
    call trans(psi, 'FFF', tfm);
    call idelsqp(psi, tfm)
    call zeroat1(psi, tfm)

    do i = 1, nr
      do j = 1, np/2
        do k = 1, nz
          rr = sqrt((tfm%r(i)*cos(tfm%p(2*j-1))-xo)**2.D0 & ! real part of the j-th entry corresponds to azim. angle p(2*j-1)
                   +(tfm%r(i)*sin(tfm%p(2*j-1))-yo)**2.D0 )
          ri = sqrt((tfm%r(i)*cos(tfm%p(2*j  ))-xo)**2.D0 & ! imag part of the j-th entry corresponds to azim. angle p(2*j)
                   +(tfm%r(i)*sin(tfm%p(2*j  ))-yo)**2.D0 )
          glb(i,j,k) = cmplx( -exp(-(rr**2.D0))/q/(1.D0-tfm%x(i))**2.D0, &
                              -exp(-(ri**2.D0))/q/(1.D0-tfm%x(i))**2.D0, P8)
        enddo
      enddo
    enddo
    call chi%disassemble(glb)
    
    glb = 0.D0
    call trans(chi, 'FFF', tfm)
    call idelsqp(chi, tfm)
    call zeroat1(chi, tfm)

    deallocate( glb )

    return
  end subroutine

end program