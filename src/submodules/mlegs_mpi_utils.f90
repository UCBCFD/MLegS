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

  module procedure m_exchange_init_2dsize
    use decomp_2d_constants; use decomp_2d_mpi; use decomp_2d
    integer(i4) :: p_row, p_col

    decomp_log = D2D_LOG_QUIET ! Change to D2D_LOG_TOFILE or D2D_LOG_STDOUT for debugging
    p_row = m_nprocs
    p_col = 1
    call decomp_2d_init(nx, ny, 1, p_row, p_col, &
                            (/ .false., .false., .false. /), m_comm)
  end procedure

  module procedure m_exchange_init_3dsize
    use decomp_2d_constants; use decomp_2d_mpi; use decomp_2d

    decomp_log = D2D_LOG_QUIET ! Change to D2D_LOG_TOFILE or D2D_LOG_STDOUT for debugging
    call decomp_2d_init(nx, ny, nz, m_cart_nprocs(1), m_cart_nprocs(2), &
                        (/ .false., .false., .false. /), m_comm)
  end procedure

  module procedure m_exchange_finalize
    use decomp_2d_constants; use decomp_2d_mpi; use decomp_2d

    call decomp_2d_finalize()
  end procedure

  module procedure m_exchange_x2y_2dsize_c
    use decomp_2d_constants; use decomp_2d_mpi; use decomp_2d
    complex(p8), dimension(:,:,:), allocatable :: wk_in, wk_out

    if (size(slab_x,1) .ne. xsize(1)) stop 'm_exchange_x2y_2dsize_c: input x-slab 1st dim. incorrect'
    if (size(slab_x,2) .ne. xsize(2)) stop 'm_exchange_x2y_2dsize_c: input x-slab 2nd dim. incorrect'
    if (size(slab_y,1) .ne. ysize(1)) stop 'm_exchange_x2y_2dsize_c: output y-slab 1st dim. incorrect'
    if (size(slab_y,2) .ne. ysize(2)) stop 'm_exchange_x2y_2dsize_c: output y-slab 2nd dim. incorrect'
    call alloc_x(wk_in)
    call alloc_y(wk_out)
    wk_in = reshape(slab_x, (/ size(slab_x,1), size(slab_x,2), 1 /))
    call transpose_x_to_y(wk_in, wk_out)
    slab_y = reshape(wk_out, (/ size(wk_out,1), size(wk_out,2) /))
    deallocate(wk_in, wk_out)
  end procedure

  module procedure m_exchange_x2y_3dsize_c
    use decomp_2d_constants; use decomp_2d_mpi; use decomp_2d
    complex(p8), dimension(:,:,:), allocatable :: wk_in, wk_out

    if (size(pencil_x,1) .ne. xsize(1)) stop 'm_exchange_x2y_3dsize_c: input x-pencil 1st dim. incorrect'
    if (size(pencil_x,2) .ne. xsize(2)) stop 'm_exchange_x2y_3dsize_c: input x-pencil 2nd dim. incorrect'
    if (size(pencil_x,3) .ne. xsize(3)) stop 'm_exchange_x2y_3dsize_c: input x-pencil 3rd dim. incorrect'
    if (size(pencil_y,1) .ne. ysize(1)) stop 'm_exchange_x2y_3dsize_c: output y-pencil 1st dim. incorrect'
    if (size(pencil_y,2) .ne. ysize(2)) stop 'm_exchange_x2y_3dsize_c: output y-pencil 2nd dim. incorrect'
    if (size(pencil_y,3) .ne. ysize(3)) stop 'm_exchange_x2y_3dsize_c: output y-pencil 3rd dim. incorrect'
    call alloc_x(wk_in)
    call alloc_y(wk_out)
    wk_in = pencil_x
    call transpose_x_to_y(wk_in, wk_out)
    pencil_y = wk_out
    deallocate(wk_in, wk_out)
  end procedure

  module procedure m_exchange_y2x_2dsize_c
    use decomp_2d_constants; use decomp_2d_mpi; use decomp_2d
    complex(p8), dimension(:,:,:), allocatable :: wk_in, wk_out

    if (size(slab_y,1) .ne. ysize(1)) stop 'm_exchange_y2x_2dsize_c: input y-slab 1st dim. incorrect'
    if (size(slab_y,2) .ne. ysize(2)) stop 'm_exchange_y2x_2dsize_c: input y-slab 2nd dim. incorrect'
    if (size(slab_x,1) .ne. xsize(1)) stop 'm_exchange_y2x_2dsize_c: output x-slab 1st dim. incorrect'
    if (size(slab_x,2) .ne. xsize(2)) stop 'm_exchange_y2x_2dsize_c: output x-slab 2nd dim. incorrect'
    call alloc_y(wk_in)
    call alloc_x(wk_out)
    wk_in = reshape(slab_y, (/ size(slab_y,1), size(slab_y,2), 1 /))
    call transpose_y_to_x(wk_in, wk_out)
    slab_x = reshape(wk_out, (/ size(wk_out,1), size(wk_out,2) /))
    deallocate(wk_in, wk_out)
  end procedure

  module procedure m_exchange_y2x_3dsize_c
    use decomp_2d_constants; use decomp_2d_mpi; use decomp_2d
    complex(p8), dimension(:,:,:), allocatable :: wk_in, wk_out

    if (size(pencil_y,1) .ne. ysize(1)) stop 'm_exchange_y2x_3dsize_c: input y-pencil 1st dim. incorrect'
    if (size(pencil_y,2) .ne. ysize(2)) stop 'm_exchange_y2x_3dsize_c: input y-pencil 2nd dim. incorrect'
    if (size(pencil_y,3) .ne. ysize(3)) stop 'm_exchange_y2x_3dsize_c: input y-pencil 3rd dim. incorrect'
    if (size(pencil_x,1) .ne. xsize(1)) stop 'm_exchange_y2x_3dsize_c: output x-pencil 1st dim. incorrect'
    if (size(pencil_x,2) .ne. xsize(2)) stop 'm_exchange_y2x_3dsize_c: output x-pencil 2nd dim. incorrect'
    if (size(pencil_x,3) .ne. xsize(3)) stop 'm_exchange_y2x_3dsize_c: output x-pencil 3rd dim. incorrect'
    call alloc_y(wk_in)
    call alloc_x(wk_out)
    wk_in = pencil_y
    call transpose_y_to_x(wk_in, wk_out)
    pencil_x = wk_out
    deallocate(wk_in, wk_out)
  end procedure

  module procedure m_exchange_y2z_3dsize_c
    use decomp_2d_constants; use decomp_2d_mpi; use decomp_2d
    complex(p8), dimension(:,:,:), allocatable :: wk_in, wk_out

    if (size(pencil_y,1) .ne. ysize(1)) stop 'm_exchange_y2z_3dsize_c: input y-pencil 1st dim. incorrect'
    if (size(pencil_y,2) .ne. ysize(2)) stop 'm_exchange_y2z_3dsize_c: input y-pencil 2nd dim. incorrect'
    if (size(pencil_y,3) .ne. ysize(3)) stop 'm_exchange_y2z_3dsize_c: input y-pencil 3rd dim. incorrect'
    if (size(pencil_z,1) .ne. zsize(1)) stop 'm_exchange_y2z_3dsize_c: output z-pencil 1st dim. incorrect'
    if (size(pencil_z,2) .ne. zsize(2)) stop 'm_exchange_y2z_3dsize_c: output z-pencil 2nd dim. incorrect'
    if (size(pencil_z,3) .ne. zsize(3)) stop 'm_exchange_y2z_3dsize_c: output z-pencil 3rd dim. incorrect'
    call alloc_y(wk_in)
    call alloc_z(wk_out)
    wk_in = pencil_y
    call transpose_y_to_z(wk_in, wk_out)
    pencil_z = wk_out
    deallocate(wk_in, wk_out)
  end procedure

  module procedure m_exchange_z2y_3dsize_c
    use decomp_2d_constants; use decomp_2d_mpi; use decomp_2d
    complex(p8), dimension(:,:,:), allocatable :: wk_in, wk_out

    if (size(pencil_z,1) .ne. zsize(1)) stop 'm_exchange_z2y_3dsize_c: input z-pencil 1st dim. incorrect'
    if (size(pencil_z,2) .ne. zsize(2)) stop 'm_exchange_z2y_3dsize_c: input z-pencil 2nd dim. incorrect'
    if (size(pencil_z,3) .ne. zsize(3)) stop 'm_exchange_z2y_3dsize_c: input z-pencil 3rd dim. incorrect'
    if (size(pencil_y,1) .ne. ysize(1)) stop 'm_exchange_z2y_3dsize_c: output y-pencil 1st dim. incorrect'
    if (size(pencil_y,2) .ne. ysize(2)) stop 'm_exchange_z2y_3dsize_c: output y-pencil 2nd dim. incorrect'
    if (size(pencil_y,3) .ne. ysize(3)) stop 'm_exchange_z2y_3dsize_c: output y-pencil 3rd dim. incorrect'
    call alloc_z(wk_in)
    call alloc_y(wk_out)
    wk_in = pencil_z
    call transpose_z_to_y(wk_in, wk_out)
    pencil_y = wk_out
    deallocate(wk_in, wk_out)
  end procedure

end submodule