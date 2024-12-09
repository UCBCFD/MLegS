submodule (mlegs_bndmat) mlegs_bndmat_tbp
  implicit none

contains

  module procedure real_bndm_alloc
    integer(i4) :: sublen_, suplen_, i, j
    if (present(sublen)) then
      sublen_ = max(sublen, 0)
    else
      sublen_ = 0
    endif
    if (present(suplen)) then
      suplen_ = max(suplen, 0)
    else
      suplen_ = 0
    endif
    this%sublen = sublen_
    this%suplen = suplen_
    allocate( this%e(sublen_+suplen_+1, nj) )
    this%e = huge(1.D0) ! initialized as machine inf to recognize 'unused' entries
    do j = 1, nj
      do i = max(1, j-suplen_), min(nj, j+sublen_)
        this%e(suplen_+1+i-j, j) = 0.D0 ! initialized as unity to recognize 'used' entries
      enddo
    enddo
  end procedure

  module procedure real_bndm_dealloc
    deallocate( this%e )
    this%sublen = 0
    this%suplen = 0
  end procedure

  module procedure real_bndm_transpose
    integer(i4) :: ku, kl, ni, nj, i, j

    ni = size(this%e, 1)
    nj = size(this%e, 2)
    kl = this%sublen
    ku = this%suplen
    i = count(this%e(ni,:) .eq. huge(1.D0))
    j = count(this%e(:,1) .eq. huge(1.D0))

    call tp%alloc(nj=(ni-i)+(nj-j)-1, sublen=ku, suplen=kl)
    do i = 1, ni
      j = min(size(this%e, 2), size(tp%e, 2))
      tp%e(i, 1:j) = this%e(ni+1-i,:)
      tp%e(i,:) = cshift(tp%e(i,:), -ni+i+ku )
    enddo
  end procedure

  module procedure cmpx_bndm_alloc
    integer(i4) :: sublen_, suplen_, i, j
    if (present(sublen)) then
      sublen_ = max(sublen, 0)
    else
      sublen_ = 0
    endif
    if (present(suplen)) then
      suplen_ = max(suplen, 0)
    else
      suplen_ = 0
    endif
    this%sublen = sublen_
    this%suplen = suplen_
    allocate( this%e(sublen_+suplen_+1, nj) )
    this%e = huge(1.D0) + huge(1.D0)*iu ! initialized as machine inf to recognize 'unused' entries
    do j = 1, nj
      do i = max(1, j-suplen_), min(nj, j+sublen_)
        this%e(suplen_+1+i-j, j) = 0.D0 + 0.D0*iu ! initialized as zero to recognize 'used' entries
      enddo
    enddo
  end procedure

  module procedure cmpx_bndm_dealloc
    deallocate( this%e )
    this%sublen = 0
    this%suplen = 0
  end procedure

  module procedure cmpx_bndm_transpose
    integer(i4) :: ku, kl, ni, nj, i, j

    ni = size(this%e, 1)
    nj = size(this%e, 2)
    kl = this%sublen
    ku = this%suplen
    i = count(this%e(ni,:) .eq. huge(1.D0) + huge(1.D0)*iu)
    j = count(this%e(:,1) .eq. huge(1.D0) + huge(1.D0)*iu)
    call tp%alloc(nj=(ni-i)+(nj-j)-1, sublen=ku, suplen=kl)
    do i = 1, ni
      j = min(size(this%e, 2), size(tp%e, 2))
      tp%e(i, 1:j) = this%e(ni+1-i,:)
      tp%e(i,:) = cshift(tp%e(i,:), -ni+i+ku )
    enddo
  end procedure

  module procedure cmpx_real
    call real_part%alloc(nj=size(this%e,2), sublen=this%sublen, suplen=this%suplen)
    real_part%e = real(this%e)
  end procedure

  module procedure cmpx_aimag
    call imag_part%alloc(nj=size(this%e,2), sublen=this%sublen, suplen=this%suplen)
    imag_part%e = aimag(this%e)
  end procedure

  module procedure realbnd_eq_cmpxbnd
    a%e = real(b%e)
    a%sublen = b%sublen
    a%suplen = b%suplen
  end procedure

  module procedure cmpxbnd_eq_realbnd
    a%e = b%e
    where( a%e .eq. huge(1.D0) ) a%e = huge(1.D0) + huge(1.D0)*iu
    a%sublen = b%sublen
    a%suplen = b%suplen
  end procedure

end submodule