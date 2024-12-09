submodule (mlegs_scalar) mlegs_scalar_init
  implicit none

contains
  !> initialize a scalar
  module procedure scalar_init
    integer(i4) :: i, comm_grp_idx, nproc_i, rank_i, loc_sz_i, loc_st_i

    this%glb_sz = glb_sz
    ! multi-dimension support
    select case (size(axis_comm))
      ! 1D: r
      case (1)
        if ((glb_sz(2).ne.1).or.(glb_sz(3).ne.1)) then
          if (rank_glb.eq.0) write(*,*) &
          "ERROR: scalar_initialize is expecting 1D data"
          call MPI_abort(comm_glb,error_flag_comm,MPI_err)
        endif
        if ((is_warning).and.(rank_glb.eq.0)) write(*,*) &
        "WARNING: scalar_initialize detects 1D (radial) configuration"
        this%axis_comm = 0
        this%axis_comm(1) = axis_comm(1)
      ! 2D: r-theta or r-z
      case (2)
        ! find the dimension that is not present: theta or z
        i = findloc(glb_sz,1,dim=1)
        if (i.eq.1) then
          if (rank_glb.eq.0) write(*,*) &
          "ERROR: scalar_initialize does not support 2D (theta-z) configuration"
          call MPI_abort(comm_glb,error_flag_comm,MPI_err)
        endif
        ! set axis_comm
        this%axis_comm(i) = 0
        this%axis_comm(3-i) = axis_comm(2)
        this%axis_comm(1) = axis_comm(1)
      ! 3D: r-theta-z
      case (3)
        this%axis_comm = axis_comm
    end select

    !> along each dimension
    do i = 1,3
      comm_grp_idx = this%axis_comm(i)
      !> dimension that is not distributed
      if (comm_grp_idx.eq.0) then
        this%loc_sz(i) = glb_sz(i)
        this%loc_st(i) = 0
      !> dimension that gets distributed
      else
        if (comm_grp_idx.gt.count(comm_grps.ne.0)) then
          if (rank_glb.eq.0) write(*,*) &
          "ERROR: scalar_initialize does not have enough communicator groups available"
          call MPI_abort(comm_glb,error_flag_comm,MPI_err)
        endif
        rank_i = rank_grps(comm_grp_idx)
        nproc_i = nprocs_grps(comm_grp_idx)
        call decompose(glb_sz(i),nproc_i,rank_i,loc_sz_i,loc_st_i)
        this%loc_sz(i) = loc_sz_i
        this%loc_st(i) = loc_st_i
      endif
    enddo

    !> allocate
    allocate(this%e(this%loc_sz(1),this%loc_sz(2),this%loc_sz(3)))
    this%space = 'PPP'

  end procedure

  !> deallocate a scalar
  module procedure scalar_dealloc
    if (.not.associated(this%e)) then
      if ((is_warning).and.(rank_glb.eq.0)) write(*,*) &
      "WARNING: scalar_dealloc tries to deallocate memory that has not been allocated"
      return
    endif
    deallocate(this%e)
    nullify(this%e)
    this%ln = 0.D0
    this%nrchop_offset = 0; this%npchop_offset = 0; this%nzchop_offset = 0
    this%space = ''
  end procedure

  !> set offset values for chopping
  module procedure scalar_chop_offset
    implicit none
    integer(i4) :: iof2_, iof3_

    if (present(iof2)) then
      iof2_ = iof2
    else
      iof2_ = 0
    endif
    
    if (present(iof3)) then
      iof3_ = iof3
    else
      iof3_ = 0
    endif

    this%nrchop_offset = iof1
    this%npchop_offset = iof2_
    this%nzchop_offset = iof3_
  end procedure

  !> copy a scalar
  module procedure scalar_copy
    if (associated(that%e) .eqv. .false.) then
      if ((is_warning).and.(rank_glb.eq.0)) write(*,*) &
      "WARNING: scalar_copy tries to copy from a scalar that has not been allocated"
      return
    endif
    if (associated(this%e) .eqv. .false.) then
      call this%init(that%glb_sz,that%axis_comm)
    elseif ((any(this%loc_sz.ne.that%loc_sz)).or.(any(this%axis_comm.ne.that%axis_comm))) then
      if ((is_warning).and.(rank_glb.eq.0)) write(*,*) &
      "WARNING: scalar_copy has inconsistent communicator layout and/or data size - reset target scalar ", &
      "with axis_comm configuration: OLD(", &
      trim(adjustl(ntoa(this%axis_comm(1)))), ",", trim(adjustl(ntoa(this%axis_comm(2)))), ",", &
      trim(adjustl(ntoa(this%axis_comm(3)))), ") to NEW(", &
      trim(adjustl(ntoa(that%axis_comm(1)))), ",", trim(adjustl(ntoa(that%axis_comm(2)))), ",", &
      trim(adjustl(ntoa(that%axis_comm(3)))), ")"
      call this%dealloc()
      call this%init(that%glb_sz,that%axis_comm)
    endif

    this%e = that%e

    ! below might be unnecessary, but it ensusres all data are copied exactly
    this%loc_st = that%loc_st
    this%loc_sz = that%loc_sz
    this%glb_sz = that%glb_sz
    this%axis_comm = that%axis_comm
    
    this%ln = that%ln
    this%nrchop_offset = that%nrchop_offset
    this%npchop_offset = that%npchop_offset
    this%nzchop_offset = that%nzchop_offset
  
    this%space = that%space
  end procedure

end submodule