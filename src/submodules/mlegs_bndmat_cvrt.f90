submodule (mlegs_bndmat) mlegs_bndmat_cvrt
  implicit none

contains

  module procedure fullmatr
    real(p8), dimension(:,:), allocatable :: tmp
    integer(i4) :: ku, kl, ni, nj, i, j
    ni = size(bndm%e, 1)
    nj = size(bndm%e, 2)
    ku = bndm%suplen
    kl = bndm%sublen
    i = count(bndm%e(ni,:) .eq. huge(1.D0))
    j = count(bndm%e(:,1) .eq. huge(1.D0))
    allocate(genm((ni-i)+(nj-j)-1,nj))
    genm = 0.D0
    ni = size(genm, 1)
    nj = size(genm, 2)
    do j = 1, nj
      do i = max(1, j-ku), min(ni, j+kl)
        genm(i,j) = bndm%e(ku+1+i-j,j) 
      enddo
    enddo
  end procedure

  module procedure fullmatc
    complex(p8), dimension(:,:), allocatable :: tmp
    integer(i4) :: ku, kl, ni, nj, i, j
    ni = size(bndm%e, 1)
    nj = size(bndm%e, 2)
    ku = bndm%suplen
    kl = bndm%sublen
    i = count(bndm%e(ni,:) .eq. huge(1.D0) + huge(1.D0)*iu)
    j = count(bndm%e(:,1) .eq. huge(1.D0) + huge(1.D0)*iu)
    allocate(genm((ni-i)+(nj-j)-1,nj))
    genm = 0.D0
    ni = size(genm, 1)
    nj = size(genm, 2)
    do j = 1, nj
      do i = max(1, j-ku), min(ni, j+kl)
        genm(i,j) = bndm%e(ku+1+i-j,j) 
      enddo
    enddo
  end procedure

  module procedure bandmatr
    integer(i4) :: ni, nj, i, j
    if (sublen .lt. 0) sublen = 0
    if (suplen .lt. 0) suplen = 0
    ni = size(genm, 1)
    nj = size(genm, 2)
    if (ni .lt. sublen) then
      write(*,*) 'bandmatr: incompatible subdiagonal number'
      stop
    endif
    if (nj .lt. suplen) then
      write(*,*) 'bandmatr: incompatible subdiagonal number'
      stop
    endif
    bndm%sublen = sublen
    bndm%suplen = suplen
    call bndm%alloc(nj, sublen, suplen)
    do j = 1, nj
      do i = max(1, j-suplen), min(ni, j+sublen)
        bndm%e(suplen+1+i-j, j) = genm(i, j)
      enddo
    enddo  
  end procedure

  module procedure bandmatc
    integer(i4) :: ni, nj, i, j
    if (sublen .lt. 0) sublen = 0
    if (suplen .lt. 0) suplen = 0
    ni = size(genm, 1)
    nj = size(genm, 2)
    if (ni .lt. sublen) then
      write(*,*) 'bandmatr: incompatible subdiagonal number'
      stop
    endif
    if (nj .lt. suplen) then
      write(*,*) 'bandmatr: incompatible subdiagonal number'
      stop
    endif
    bndm%sublen = sublen
    bndm%suplen = suplen
    call bndm%alloc(nj, sublen, suplen)
    do j = 1, nj
      do i = max(1, j-suplen), min(ni, j+sublen)
        bndm%e(suplen+1+i-j, j) = genm(i, j)
      enddo
    enddo
  end procedure

end submodule