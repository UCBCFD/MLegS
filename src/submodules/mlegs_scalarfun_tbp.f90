submodule (mlegs_scalarfun) mlegs_scalarfun_tbp
  implicit none

contains

  module procedure scalar_set
    use decomp_2d_constants; use decomp_2d_mpi; use decomp_2d

    if (kit%is_set .eqv. .false.) stop 'scalar_set: initialize the kit first -- run kit%set()'
    this%ln = 0.D0
    select type(this)
      type is (scalar_1d)
        select type(kit)
          type is (tfm_kit_2d)
            stop 'scalar_set: incompatible dimensionality between the kit and field(s) to be set'
          type is (tfm_kit_3d)
            stop 'scalar_set: incompatible dimensionality between the kit and field(s) to be set'
          type is (tfm_kit_1d)
            this%local_chops = kit%chops
            if (m_rank .eq. 0) then
              allocate( this%e(1:kit%nrdim) )
              this%e = 0.D0 + 0.D0*iu
            else
              allocate( this%e(0) )
            endif
            this%space = 'P'
        end select
      type is (scalar_2d)
        select type(kit)
          type is (tfm_kit_1d)
            stop 'scalar_set: incompatible dimensionality between the kit and field(s) to be set'
          type is (tfm_kit_3d)
            stop 'scalar_set: incompatible dimensionality between the kit and field(s) to be set'
          type is (tfm_kit_2d)
            this%local_chops = kit%chops
            allocate( this%e(xstart(1):xend(1), xstart(2):xend(2)) )
            this%e = 0.D0 + 0.D0*iu
            this%local_chopp = kit%chopp
            this%space = 'PP'
            this%slab_dir = 1
        end select
      type is (scalar_3d)
        select type (kit)
          type is (tfm_kit_1d)
            stop 'scalar_set: incompatible dimensionality between the kit and field(s) to be set'
          type is (tfm_kit_2d)
            stop 'scalar_set: incompatible dimensionality between the kit and field(s) to be set'
          type is (tfm_kit_3d)
            this%local_chops = kit%chops
            allocate( this%e(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)) )
            this%e = 0.D0 + 0.D0*iu
            this%local_chopzl = kit%chopzl
            this%local_chopzu = kit%chopzu
            this%space = 'PPP'
            this%pencil_dir = 1
        end select
    end select
  end procedure

  module procedure scalar_chopset
    integer(i4) :: i, j

    select type(this)
      type is (scalar_1d)
        this%local_chops = (/ max(min(this%local_chops + iof, size(this%e, 1)), 0) /)
      type is (scalar_2d)
        if (present(iofp)) then
          deallocate( this%local_chops )
          j = max(min(this%local_chopp + iofp, size(this%e, 2)), 0)
          allocate( this%local_chops(j) )
          this%local_chops = (/ (max(min(this%local_chops(i) + iof, size(this%e, 1)), 0), i = 1, j) /)
        else
          this%local_chops = (/ (max(min(this%local_chops(i) + iof, size(this%e, 1)), 0), i = 1, npchop) /)
        endif
        this%local_chopp = j
      type is (scalar_3d)
        if (present(iofp)) then
          deallocate( this%local_chops )
          j = max(min(this%local_chopp + iofp, size(this%e, 2)), 0)
          allocate( this%local_chops(j) )
          this%local_chops = (/ (max(min(this%local_chops(i) + iof, size(this%e, 1)), 0), i = 1, j) /)
        else
          this%local_chops = (/ (max(min(this%local_chops(i) + iof, size(this%e, 1)), 0), i = 1, npchop) /)
        endif
        this%local_chopp = j
        if (present(iofz)) then
          this%local_chopzl = max(min(this%local_chopzl + iofz, size(this%e, 3)/2+1), 0)
          this%local_chopzu = max(min(this%local_chopzu - iofz, size(this%e, 3)/2+1), 0)
        endif
    end select
  end procedure

  module procedure scalar_chopdo
    integer(i4) :: mm

    select type(this)
      type is (scalar_1d)
        if (this%space(1:1) .eq. 'F') then
          if (lbound(this%e, 1) .le. this%local_chops(1)+1) then
            this%e(this%local_chops(1)+1:) = 0.D0 + 0.D0*iu
          else
            this%e(:) = 0.D0 + 0.D0*iu
          endif
        endif
      type is (scalar_2d)
        if (this%space(1:1) .eq. 'F') then
          if (lbound(this%e, 1) .le. this%local_chopp+1) then
            this%e(this%local_chopp+1:,:) = 0.D0 + 0.D0*iu
          else
            this%e(:,:) = 0.D0 + 0.D0*iu
          endif
        endif
        if (this%space(1:1) .eq. 'F') then
          do mm = lbound(this%e, 1), min(this%local_chopp, ubound(this%e, 1))
            if (lbound(this%e, 2) .le. this%local_chops(mm)+1) then
              this%e(mm,this%local_chops(mm)+1:) = 0.D0 + 0.D0*iu
            else
              this%e(mm,:) = 0.D0 + 0.D0*iu
            endif
          enddo
        endif
      type is (scalar_3d)
        if (this%space(1:1) .eq. 'F') then
          if (lbound(this%e, 1) .le. this%local_chopp+1) then
            this%e(this%local_chopp+1:,:,:) = 0.D0 + 0.D0*iu
          else
            this%e(:,:,:) = 0.D0 + 0.D0*iu
          endif
        endif
        if (this%space(2:2) .eq. 'F') then
          do mm = lbound(this%e, 1), min(this%local_chopp, ubound(this%e, 1))
            if (lbound(this%e, 2) .le. this%local_chops(mm)+1) then
              this%e(mm,this%local_chops(mm)+1:,:) = 0.D0 + 0.D0*iu
            else
              this%e(mm,:,:) = 0.D0 + 0.D0*iu
            endif
          enddo
        endif
        if (this%space(3:3) .eq. 'F') then
          if (this%local_chopzl+1 .le. this%local_chopzu-1) then
            if (lbound(this%e, 3) .le. this%local_chopzl+1) then
              if (ubound(this%e, 3) .ge. this%local_chopzu+1) then
                this%e(:,:,this%local_chopzl+1:this%local_chopzu-1) = 0.D0 + 0.D0*iu
              else
                this%e(:,:,this%local_chopzl+1:) = 0.D0 + 0.D0*iu
              endif
            else
              if (ubound(this%e, 3) .ge. this%local_chopzu+1) then
                this%e(:,:,:this%local_chopzu-1) = 0.D0 + 0.D0*iu
              else
                this%e(:,:,:) = 0.D0 + 0.D0*iu
              endif
            endif
          endif
        endif
    end select
  end procedure

  module procedure scalar_dealloc
    select type(this)
      type is (scalar_1d)
        deallocate( this%e )
      type is (scalar_2d)
        deallocate( this%e, this%local_chops )
      type is (scalar_3d)
        deallocate( this%e, this%local_chops )
    end select
  end procedure

  module procedure scalar_1d_tfm
    integer(i4) :: i

    if (.not. ((sp .eq. 'P') .or. (sp .eq. 'F'))) then
      stop 'scalar_1d_tfm: (P) or (F) only allowed (rad-)'
    endif
    if ((this%space(1:1) .eq. 'P') .and. (sp(1:1) .eq. 'F')) then ! P to F
      ! 
      this%space(1:1) = 'F'
      call this%chopdo()
    endif
    if ((this%space(1:1) .eq. 'F') .and. (sp(1:1) .eq. 'P')) then ! F to P
      !
      this%space(1:1) = 'P'
      return
    endif
  end procedure

  module procedure scalar_2d_tfm
    integer(i4) :: i, j

    if (.not. ((sp .eq. 'PP') .or. (sp .eq. 'FP') .or. (sp .eq. 'FF'))) then
      stop 'scalar_2d_tfm: (P, P), (F, P) or (F, F) only allowed (azim-, rad-)'
    endif
    if ((this%space(1:1) .eq. 'P') .and. (sp(1:1) .eq. 'F')) then ! azim-P to azim-F
      if (this%slab_dir.eq. 2) then
        !
      endif
      !
      this%space(1:1) = 'F'
      call this%chopdo()
      this%slab_dir = 1
    endif
    if ((this%space(2:2) .eq. 'P') .and. (sp(2:2) .eq. 'F')) then ! rad-P to rad-F
      if (this%slab_dir .eq. 1) then ! under regular transform this if block should be run
        !
      endif
      !
      this%space(2:2) = 'F'
      call this%chopdo()
      this%slab_dir = 2
    endif
    if ((this%space(2:2) .eq. 'F') .and. (sp(2:2) .eq. 'P')) then ! rad-F to rad-P
      if (this%slab_dir.eq. 1) then
        !
      endif
      !
      this%space(2:2) = 'P'
      this%slab_dir = 2
    endif
    if ((this%space(1:1) .eq. 'F') .and. (sp(1:1) .eq. 'P')) then ! azim-P to azim-F
      if (this%slab_dir.eq. 2) then ! under regular transform this if block should be run
        !
      endif
      !
      this%space(1:1) = 'P'
      this%slab_dir = 1
    endif
  end procedure
  
  module procedure scalar_3d_tfm
    integer(i4) :: i, j, k

    if (.not. ((sp .eq. 'PPP') .or. (sp .eq. 'FPP') .or. (sp .eq. 'FFP') .or. (sp .eq. 'FFF'))) then
      stop 'scalar_3d_tfm: (P, P, P), (F, P, P), (F, F, P) or (F, F, F) only allowed (azim-, rad-, axi-)'
    endif

    if ((this%space(1:1) .eq. 'P') .and. (sp(1:1) .eq. 'F')) then ! azim-P to azim-F
      if (this%pencil_dir .eq. 2) then
        !
      endif
      if (this%pencil_dir .eq. 3) then
        !
      endif
      !
      this%space(1:1) = 'F'
      call this%chopdo()
      this%pencil_dir = 1
    endif
    if ((this%space(2:2) .eq. 'P') .and. (sp(2:2) .eq. 'F')) then ! rad-P to rad-F
      if (this%pencil_dir .eq. 1) then ! under regular transform this if block should be run
        !
      endif
      if (this%pencil_dir .eq. 3) then
        !
      endif
      !
      this%space(2:2) = 'F'
      call this%chopdo()
      this%pencil_dir = 2
    endif
    if ((this%space(3:3) .eq. 'P') .and. (sp(3:3) .eq. 'F')) then ! axi-P to axi-F
      if (this%pencil_dir .eq. 1) then
        !
      endif
      if (this%pencil_dir .eq. 2) then ! under regular transform this if block should be run
        !
      endif
      !
      this%space(3:3) = 'F'
      call this%chopdo()
      this%pencil_dir = 3
    endif
    if ((this%space(3:3) .eq. 'F') .and. (sp(3:3) .eq. 'P')) then ! axi-F to axi-P
      if (this%pencil_dir .eq. 1) then
        !
      endif
      if (this%pencil_dir .eq. 2) then
        !
      endif
      !
      this%space(3:3) = 'P'
      this%pencil_dir = 3
    endif
    if ((this%space(2:2) .eq. 'F') .and. (sp(2:2) .eq. 'P')) then ! rad-F to rad-P
      if (this%pencil_dir .eq. 1) then
        !
      endif
      if (this%pencil_dir .eq. 3) then ! under regular transform this if block should be run
        !
      endif
      !
      this%space(2:2) = 'P'
      this%pencil_dir = 2
    endif
    if ((this%space(1:1) .eq. 'F') .and. (sp(1:1) .eq. 'P')) then ! azim-P to azim-F
      if (this%pencil_dir .eq. 2) then ! under regular transform this if block should be run
        !
      endif
      if (this%pencil_dir .eq. 3) then
        !
      endif
      !
      this%space(1:1) = 'P'
      this%pencil_dir = 1
    endif
  end procedure

end submodule