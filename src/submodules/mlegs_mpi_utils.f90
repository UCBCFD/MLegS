submodule (mlegs_mpi) mlegs_mpi_utils
  implicit none

contains

  module procedure m_initialize
    use MPI
    integer(i4) :: i, world_group, world_nprocs, active_group
    integer(i4), dimension(:), allocatable :: rks, tmp

    call MPI_Init(m_err)
    call MPI_Comm_group(MPI_Comm_world, world_group, m_err)
    call MPI_Comm_size(MPI_Comm_world, world_nprocs, m_err)
    if (present(active_nprocs)) then
      if (active_nprocs .gt. world_nprocs) stop 'm_initialize: invalid number of processors requested'
      rks = (/ (i, i=0, active_nprocs-1) /)
      m_nprocs = active_nprocs
    else
      rks = (/ (i, i=0, world_nprocs-1) /)
      m_nprocs = world_nprocs
    endif
    call MPI_Group_incl(world_group, m_nprocs, rks, active_group, m_err)
    call MPI_Comm_create(MPI_Comm_world, active_group, m_comm, m_err)
    if (m_comm .ne. MPI_Comm_null) then
      call MPI_Dims_create(m_nprocs, 2, m_cart_nprocs, m_err)
      call MPI_Comm_rank(m_comm, m_rank, m_err)
      call MPI_Cart_create(m_comm, 2, m_cart_nprocs, (/.false.,.false./), .false., m_cart_comm, m_err)
      allocate( m_map_row(0:m_nprocs-1), m_map_col(0:m_nprocs-1), tmp(2) )
      allocate( m_map_rank(0:m_cart_nprocs(1)-1, 0:m_cart_nprocs(2)-1))
      do i = 0, m_nprocs-1
        call MPI_Cart_coords(m_cart_comm, i, 2, tmp, m_err)
        m_map_row(i) = tmp(1)
        m_map_col(i) = tmp(2)
        m_map_rank(tmp(1),tmp(2)) = i
      enddo
      deallocate( tmp )
      call MPI_Cart_coords(m_cart_comm, m_rank, 2, m_cart_coords, m_err)
      call MPI_Comm_split(m_cart_comm, m_cart_coords(2), m_cart_coords(1), m_row_comm, m_err)
      call MPI_Comm_split(m_cart_comm, m_cart_coords(1), m_cart_coords(2), m_col_comm, m_err)
      call MPI_Comm_rank(m_row_comm, m_row_rank, m_err)
      call MPI_Comm_rank(m_col_comm, m_col_rank, m_err)
      call MPI_Comm_size(m_row_comm, m_row_nprocs, m_err)
      call MPI_Comm_size(m_col_comm, m_col_nprocs, m_err)
    else
      m_rank = -1
      m_nprocs = 0
    endif
  end procedure

  module procedure m_finalize
    use MPI

    call MPI_Finalize(m_err)
  end procedure

  module procedure m_3dsize_1ddecomp
    integer(i4) :: i, sze(3)
    integer(i4) :: nprocs_, rank_
    character(len=6) :: xyz = 'XYZxyz'

    nprocs_ = merge(nprocs, m_nprocs, present(nprocs))
    rank_ = merge(rank, m_rank, present(rank))
    sze = (/ size_1, size_2, size_3 /)
    if (index(xyz, decomp_dir(1:1)) .eq. 0) then
      stop 'm_3dsize_1ddecomp: invalid decomposition direction'
    endif
    do i = 1, 3
      if (mod(index(xyz, decomp_dir(1:1))-1,3)+1 .eq. i) then
        if (sze(i) .ge. nprocs_) then
          chunk(i) = (sze(i)+nprocs_/2)/nprocs_
          start_idx(i) = rank_ * chunk(i) + 1
          final_idx(i) = (rank_ + 1) * chunk(i)
          if (rank_ .eq. nprocs_ - 1) final_idx(i) = sze(i)
        else
          if (rank_ .lt. sze(i)) then
            chunk(i) = 1
            start_idx(i) = rank_ + 1
            final_idx(i) = start_idx(i)
          else
            chunk(i) = 0
            start_idx(i) = 0
            final_idx(i) = -1
          endif
        endif
        chunk(i) = final_idx(i) - start_idx(i) + 1
      else
        start_idx(i) = 1
        final_idx(i) = sze(i)
        chunk(i) = sze(i)
      endif
    enddo
  end procedure

  module procedure m_2dsize_1ddecomp
    integer(i4) :: chunk_(3), start_idx_(3), final_idx_(3)
    integer(i4) :: nprocs_, rank_
    character(len=4) :: xy = 'XYxy'

    nprocs_ = merge(nprocs, m_nprocs, present(nprocs))
    rank_ = merge(rank, m_rank, present(rank))
    if (index(xy, decomp_dir(1:1)) .eq. 0) then
      stop 'm_2dsize_1ddecomp: invalid decomposition direction'
    endif
    call m_3dsize_1ddecomp(size_1, size_2, 1, decomp_dir, &
                           chunk_, start_idx_, final_idx_, nprocs_, rank_)
    chunk = chunk_(1:2)
    start_idx = start_idx_(1:2)
    final_idx = final_idx_(1:2)
  end procedure

  module procedure m_1dsize_1ddecomp
    integer(i4) :: chunk_(3), start_idx_(3), final_idx_(3)
    integer(i4) :: nprocs_, rank_
    character(len=1) :: decomp_dir = 'X'
    nprocs_ = merge(nprocs, m_nprocs, present(nprocs))
    rank_ = merge(rank, m_rank, present(rank))
    call m_3dsize_1ddecomp(size_1, 1, 1, decomp_dir, &
                           chunk_, start_idx_, final_idx_, nprocs_, rank_)
    chunk = chunk_(1)
    start_idx = start_idx_(1)
    final_idx = final_idx_(1)
  end procedure

  module procedure m_3dsize_2ddecomp
    integer(i4) :: sze(3), cart_comm
    integer(i4) :: decomp_dir(2), cart_procs(2)
    integer(i4) :: ncartprocs_(2), coords_(2)
    character(len=6) :: xyz = 'XYZxyz'

    ncartprocs_ = merge(ncartprocs, m_cart_nprocs, present(ncartprocs))
    coords_ = merge(coords, m_cart_coords, present(coords))
    sze = (/ size_1, size_2, size_3 /)
    if (decomp_dir_1(1:1) .eq. decomp_dir_2(1:1)) then
      stop 'm_3dsize_2ddecomp: invalid 2d decompition direction inputs'
    endif
    if (index(xyz, decomp_dir_1(1:1)) .eq. 0) then
      stop 'm_3dsize_2ddecomp: invalid decomposition direction (1st dir.)'
    endif
    if (index(xyz, decomp_dir_2(1:1)) .eq. 0) then
      stop 'm_3dsize_2ddecomp: invalid decomposition direction (2nd dir.)'
    endif
    cart_procs = (/ 0, 0 /)
    decomp_dir = (/ min(mod(index(xyz,decomp_dir_1(1:1))-1,3)+1, mod(index(xyz,decomp_dir_2(1:1))-1,3)+1), &
                    max(mod(index(xyz,decomp_dir_1(1:1))-1,3)+1, mod(index(xyz,decomp_dir_2(1:1))-1,3)+1) /)
    if ((decomp_dir(1) .eq. 1) .and. (decomp_dir(2) .eq. 2)) then ! r- and p- decompositions
      call m_1dsize_1ddecomp(sze(1), &
                             chunk(1), start_idx(1), final_idx(1), ncartprocs_(1), coords_(1))     
      call m_2dsize_1ddecomp(sze(2), sze(3), 'X', &
                             chunk(2:3), start_idx(2:3), final_idx(2:3), ncartprocs_(2), coords_(2))
    elseif ((decomp_dir(1) .eq. 1) .and. (decomp_dir(2) .eq. 3)) then ! r- and z- decompositions
      call m_1dsize_1ddecomp(sze(1), &
                             chunk(1), start_idx(1), final_idx(1), ncartprocs_(1), coords_(1))     
      call m_2dsize_1ddecomp(sze(2), sze(3), 'Y', &
                             chunk(2:3), start_idx(2:3), final_idx(2:3), ncartprocs_(2), coords_(2))
    else ! p- and z- decompositions
      call m_2dsize_1ddecomp(sze(1), sze(2), 'Y', &
                             chunk(1:2), start_idx(1:2), final_idx(1:2), ncartprocs_(1), coords_(1))
      call m_1dsize_1ddecomp(sze(3), &
                             chunk(3), start_idx(3), final_idx(3), ncartprocs_(2), coords_(2))
    endif
  end procedure

  module procedure m_2dsize_2ddecomp
    integer(i4) :: chunk_(3), start_idx_(3), final_idx_(3)
    integer(i4) :: ncartprocs_(2), coords_(2)
    character(len=2) :: decomp_dir = 'XY'

    ncartprocs_ = merge(ncartprocs, m_cart_nprocs, present(ncartprocs))
    coords_ = merge(coords, m_cart_coords, present(coords))
    call m_3dsize_2ddecomp(size_1, size_2, 1, decomp_dir(1:1), decomp_dir(2:2), &
                           chunk_, start_idx_, final_idx_, ncartprocs_, coords_)
    chunk = chunk_(1:2)
    start_idx = start_idx_(1:2)
    final_idx = final_idx_(1:2)
  end procedure

  module procedure m_exchange_init_2d
    use MPI
    integer(i4) :: chunk(2), start_idx(2), final_idx(2)

    m_garr_nx = nx
    m_garr_ny = ny
    call m_decompose(m_garr_nx, m_garr_ny, 'Y', chunk, start_idx, final_idx)
    m_slab_x = create_new_type(2, m_row_comm, chunk, 1, MPI_Double_complex)
    call m_decompose(m_garr_nx, m_garr_ny, 'X', chunk, start_idx, final_idx)
    m_slab_y = create_new_type(2, m_row_comm, chunk, 2, MPI_Double_complex)

  end procedure

  module procedure m_exchange_init_3d
    use MPI
    integer(i4) :: chunk(3), start_idx(3), final_idx(3)

    m_garr_nx = nx
    m_garr_ny = ny
    m_garr_nz = nz
    call m_decompose(m_garr_nx, m_garr_ny, m_garr_nz, 'Y', 'Z', chunk, start_idx, final_idx)
    m_pencil_x = create_new_type(3, m_row_comm, chunk, 1, MPI_Double_complex)
    call m_decompose(m_garr_nx, m_garr_ny, m_garr_nz, 'X', 'Z', chunk, start_idx, final_idx)
    m_pencil_y = create_new_type(3, m_row_comm, chunk, 2, MPI_Double_complex)
    call m_decompose(m_garr_nx, m_garr_ny, m_garr_nz, 'X', 'Y', chunk, start_idx, final_idx)
    m_pencil_z = create_new_type(3, m_row_comm, chunk, 3, MPI_Double_complex)
  end procedure

! ======================================================================================================== !
! VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV INTERNAL (PRIVATE) SUBROUTINES/FUNCTIONS VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV !
! ======================================================================================================== !

    function create_new_type(ndim, comm, data_size_proc_old, dim_old, data_type) result(dtype)
      use MPI
      integer(i4) :: ndim
      integer(i4) :: comm, dim_old, data_type
      integer(i4), dimension(ndim) :: data_size_proc_old, subdsize_proc, substarts_proc
      integer(i4), dimension(:), allocatable :: dtype

      integer(i4) :: i, nproc_comm

      call MPI_Comm_size(comm, nproc_comm, m_err)
      allocate( dtype(0:m_nprocs-1) )
      ! old local array: array_old of size_old is decomposed into subarray along dim_old
      call subarray(data_type, ndim, data_size_proc_old, dim_old, nproc_comm, dtype)
    end function

    subroutine subarray(element_data_type, ndim, datasize_proc, dim, nprocs, new_data_type)
      use MPI
      integer:: element_data_type, ndim, dim, nprocs
      integer,dimension(ndim):: datasize_proc, subdsize_proc, substarts_proc, subends_proc
      integer,dimension(0:nprocs-1):: new_data_type

      integer:: i

      do i = 1, ndim
          subdsize_proc(i) = datasize_proc(i)
          substarts_proc(i) = 0
          subends_proc(i) = 0
      enddo

      ! along that dim
      do i = 0,nprocs-1
          call m_decompose(datasize_proc(dim), subdsize_proc(dim), substarts_proc(dim), &
                           subends_proc(dim), nprocs, i)
          call MPI_Type_create_subarray(ndim, datasize_proc, subdsize_proc, substarts_proc, & 
                                        MPI_Order_fortran, element_data_type, new_data_type(i), m_err)
          call MPI_Type_commit(new_data_type(i), m_err)
      enddo

      return
    end subroutine subarray

end submodule