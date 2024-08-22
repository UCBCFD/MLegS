module mlegs_mpi
  use mlegs_base
  implicit none
  private

  !> MPI-related global variables across multiprocessors (m_*)
  !> error code return, current processor's rank, total no. of processors
  integer(i4), public :: m_err, m_comm, m_rank, m_nprocs
  !> communicator with 2-dimensional cartesian coordinate information attached
  integer(i4), public :: m_cart_comm, m_cart_coords(2), m_cart_nprocs(2)
  integer(i4), public :: m_row_comm, m_row_rank, m_row_nprocs
  integer(i4), public :: m_col_comm, m_col_rank, m_col_nprocs
  !> cartesian coordinate map. m_map_*_(i) gives the info of * for processor of rank i.  
  integer(i4), public, dimension(:), allocatable :: m_map_row, m_map_col
  !> rank map. m_map_* (i,j) gives the processor rank of * for cart. coord (i, j)
  integer(i4), public, dimension(:,:), allocatable :: m_map_rank
  !> global array size to be decomposed & transposed
  integer(i4), public :: m_garr_nx, m_garr_ny, m_garr_nz
  !> local (sub-) MPI data type for data exchange
  integer(i4), public, dimension(:), allocatable :: m_slab_x, m_slab_y, m_pencil_x, m_pencil_y, m_pencil_z

  !> MPI initilaizer
  interface m_initialize
    module subroutine m_initialize(active_nprocs)
      implicit none
      integer(i4), optional :: active_nprocs
    end subroutine 
  end interface
  public :: m_initialize

  !> MPI finalizer
  interface m_finalize
    module subroutine m_finalize()
      implicit none
    end subroutine
  end interface
  public :: m_finalize

  !> 1-d/2-d/3-d global data index decomposition (1-d/2-d) for MPI parallelism
  interface m_decompose
    module subroutine m_3dsize_1ddecomp(size_1, size_2, size_3, decomp_dir, & 
                                        chunk, start_idx, final_idx, nprocs, rank)
      implicit none
      integer(i4), intent(in) :: size_1, size_2, size_3
      character(len=1), intent(in) :: decomp_dir
      integer(i4), intent(inout) :: chunk(3), start_idx(3), final_idx(3)
      integer(i4), optional :: nprocs, rank
    end subroutine
    module subroutine m_2dsize_1ddecomp(size_1, size_2, decomp_dir, &
                                        chunk, start_idx, final_idx, nprocs, rank)
      implicit none
      integer(i4), intent(in) :: size_1, size_2
      character(len=1), intent(in) :: decomp_dir
      integer(i4), intent(inout) :: chunk(2), start_idx(2), final_idx(2)
      integer(i4), optional :: nprocs, rank
    end subroutine
    module subroutine m_1dsize_1ddecomp(size_1, &
                                        chunk, start_idx, final_idx, nprocs, rank)
      implicit none
      integer(i4), intent(in) :: size_1
      integer(i4), intent(inout) :: chunk, start_idx, final_idx
      integer(i4), optional :: nprocs, rank
    end subroutine
    module subroutine m_3dsize_2ddecomp(size_1, size_2, size_3, decomp_dir_1, decomp_dir_2, &
                                        chunk, start_idx, final_idx, ncartprocs, coords)
      implicit none
      integer(i4), intent(in) :: size_1, size_2, size_3 
      character(len=1), intent(in) :: decomp_dir_1, decomp_dir_2
      integer(i4), intent(inout) :: chunk(3), start_idx(3), final_idx(3)
      integer(i4), optional :: ncartprocs(2), coords(2)
    end subroutine
    module subroutine m_2dsize_2ddecomp(size_1, size_2, &
                                        chunk, start_idx, final_idx, ncartprocs, coords)
      implicit none
      integer(i4), intent(in) :: size_1, size_2
      integer(i4), intent(inout) :: chunk(2), start_idx(2), final_idx(2)
      integer(i4), optional :: ncartprocs(2), coords(2)
    end subroutine
  end interface
  public :: m_decompose

  !> initialization for MPI data exchange
  interface m_exchange_init
    module subroutine m_exchange_init_2d(nx, ny)
      implicit none
      integer(i4), intent(in) :: nx, ny
    end subroutine
    module subroutine m_exchange_init_3d(nx, ny, nz)
      implicit none
      integer(i4), intent(in) :: nx, ny, nz
    end subroutine
  end interface
  public :: m_exchange_init

  !> make an exchange between the decomposed dimension(s) and the locally residing dimension
  ! interface m_transpose_3d_x_y
  !   module subroutine m_transpose_3d_x_y(input, output, nx_global, ny_global, nz_global)
  !     implicit none
  !     complex(p8), dimension(:,:,:), intent(in) :: input
  !     complex(p8), dimension(:,:,:), intent(inout) :: output
  !     integer(i4), intent(in) :: nx_global, ny_global, nz_global
  !   end subroutine
  ! end interface
  ! public :: m_transpose_3d_x_y

  ! !> save an MPI-decomposed local matrix (array) data  
  ! interface m_msave
  !   module subroutine m_msave()
  !     implicit none
  !   end subroutine
  ! end interface
  ! public :: m_msave

  ! !> load an MPI-decomposed local matrix (array) data  
  ! interface m_msave
  !   module subroutine m_mload()
  !     implicit none
  !   end subroutine
  ! end interface
  ! public :: m_mload

end module