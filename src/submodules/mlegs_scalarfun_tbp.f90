submodule (mlegs_scalarfun) mlegs_scalarfun_tbp
  implicit none

contains

  module procedure scalar_set
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
            allocate( this%e(kit%nrdim) )
            this%e = 0.D0 + 0.D0*iu
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
            allocate( this%e(kit%nrdim, kit%npdim) )
            this%e = 0.D0 + 0.D0*iu
            this%local_chopp = kit%chopp
            this%space = 'PP'
            this%slab_dir = 2
        end select
      type is (scalar_3d)
        select type (kit)
          type is (tfm_kit_1d)
            stop 'scalar_set: incompatible dimensionality between the kit and field(s) to be set'
          type is (tfm_kit_2d)
            stop 'scalar_set: incompatible dimensionality between the kit and field(s) to be set'
          type is (tfm_kit_3d)
            this%local_chops = kit%chops
            allocate( this%e(kit%nrdim, kit%npdim, kit%nzdim) )
            this%e = 0.D0 + 0.D0*iu
            this%local_chopzl = kit%chopzl
            this%local_chopzu = kit%chopzu
            this%space = 'PPP'
            this%pencil_dir = 2
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
        if (this%space(2:2) .eq. 'F') then
          if (lbound(this%e, 2) .le. this%local_chopp+1) then
            this%e(:, this%local_chopp+1:) = 0.D0 + 0.D0*iu
          else
            this%e(:,:) = 0.D0 + 0.D0*iu
          endif
        endif
        if (this%space(1:1) .eq. 'F') then
          do mm = lbound(this%e, 2), min(this%local_chopp, ubound(this%e, 2))
            if (lbound(this%e, 1) .le. this%local_chops(mm)+1) then
              this%e(this%local_chops(mm)+1:,mm) = 0.D0 + 0.D0*iu
            else
              this%e(:, mm) = 0.D0 + 0.D0*iu
            endif
          enddo
        endif
      type is (scalar_3d)
        if (this%space(2:2) .eq. 'F') then
          if (lbound(this%e, 2) .le. this%local_chopp+1) then
            this%e(:,this%local_chopp+1:,:) = 0.D0 + 0.D0*iu
          else
            this%e(:,:,:) = 0.D0 + 0.D0*iu
          endif
        endif
        if (this%space(1:1) .eq. 'F') then
          do mm = lbound(this%e, 2), min(this%local_chopp, ubound(this%e, 2))
            if (lbound(this%e, 1) .le. this%local_chops(mm)+1) then
              this%e(this%local_chops(mm)+1:,mm,:) = 0.D0 + 0.D0*iu
            else
              this%e(:,mm,:) = 0.D0 + 0.D0*iu
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
  end procedure

  module procedure scalar_2d_tfm
  end procedure
  
  module procedure scalar_3d_tfm
  end procedure

end submodule