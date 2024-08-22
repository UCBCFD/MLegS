submodule (mlegs_mpi) mlegs_mpi_utils
  implicit none

contains

  module procedure m_initialize
    use MPI
    integer(i4) :: i, world_group, world_nprocs, active_group, cart_rank
    integer(i4), dimension(:), allocatable :: rks

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
      call MPI_Cart_create(m_comm, 2, m_cart_nprocs, (/.false.,.false./), .true., m_cart_comm, m_err)
      call MPI_Comm_rank(m_cart_comm, cart_rank, m_err)
      call MPI_Cart_coords(m_cart_comm, cart_rank, 2, m_cart_ranks, m_err)
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
    integer(i4) :: decomp_dir(2), cart_procs(2), subrank, coords(2)
    integer(i4) :: nprocs_(2), ranks_(2)
    character(len=6) :: xyz = 'XYZxyz'

    nprocs_ = merge(nprocs, m_cart_nprocs, present(nprocs))
    ranks_ = merge(ranks, m_cart_ranks, present(ranks))
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
                             chunk(1), start_idx(1), final_idx(1), nprocs(1), ranks(1))     
      call m_2dsize_1ddecomp(sze(2), sze(3), 'X', &
                             chunk(2:3), start_idx(2:3), final_idx(2:3), nprocs(2), ranks(2))
    elseif ((decomp_dir(1) .eq. 1) .and. (decomp_dir(2) .eq. 3)) then ! r- and z- decompositions
      call m_1dsize_1ddecomp(sze(1), &
                             chunk(1), start_idx(1), final_idx(1), nprocs(1), ranks(1))     
      call m_2dsize_1ddecomp(sze(2), sze(3), 'Y', &
                             chunk(2:3), start_idx(2:3), final_idx(2:3), nprocs(2), ranks(2))
    else ! p- and z- decompositions
      call m_2dsize_1ddecomp(sze(1), sze(2), 'Y', &
                             chunk(1:2), start_idx(1:2), final_idx(1:2), nprocs(1), ranks(1))
      call m_1dsize_1ddecomp(sze(3), &
                             chunk(3), start_idx(3), final_idx(3), nprocs(2), ranks(2))
    endif
  end procedure

  module procedure m_2dsize_2ddecomp
    integer(i4) :: chunk_(3), start_idx_(3), final_idx_(3)
    integer(i4) :: nprocs_(2), ranks_(2)
    character(len=2) :: decomp_dir = 'XY'

    nprocs_ = merge(nprocs, m_cart_nprocs, present(nprocs))
    ranks_ = merge(ranks, m_cart_ranks, present(ranks))
    call m_3dsize_2ddecomp(size_1, size_2, 1, decomp_dir(1:1), decomp_dir(2:2), &
                           chunk_, start_idx_, final_idx_, nprocs_, ranks_)
    chunk = chunk_(1:2)
    start_idx = start_idx_(1:2)
    final_idx = final_idx_(1:2)
  end procedure


! ======================================================================================================== !
! VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV INTERNAL (PRIVATE) SUBROUTINES/FUNCTIONS VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV !
! ======================================================================================================== !

  subroutine find_best_proc_dims_2d(nx, ny, nprocs, proc_dims)
    implicit none
    integer(i4), intent(in) :: nx, ny        ! Global grid dimensions
    integer(i4), intent(in) :: nprocs        ! Total number of processors
    integer(i4), intent(out) :: proc_dims(2) ! Best processor grid dimensions (px, py)

    integer(i4) :: i, px, py, best_px, best_py
    real(p8) :: aspect_ratio, current_ratio, best_ratio

    ! Initialize
    best_px = 1
    best_py = nprocs
    best_ratio = real(best_px) / real(best_py)

    aspect_ratio = real(nx) / real(ny) ! Aspect ratio of the grid (nx/ny)
    do px = 1, nprocs ! Loop over all possible decompositions
      if (mod(nprocs, px) == 0) then
        py = nprocs / px
        current_ratio = real(px) / real(py)
        if (abs(current_ratio - aspect_ratio) .lt. abs(best_ratio - aspect_ratio)) then
          best_ratio = current_ratio
          best_px = px
          best_py = py
        endif
      endif
    enddo
    proc_dims(1) = best_px
    proc_dims(2) = best_py
  end subroutine

end submodule