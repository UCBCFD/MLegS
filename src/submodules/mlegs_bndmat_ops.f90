submodule (mlegs_bndmat) mlegs_bndmat_ops
  implicit none

contains

  module procedure mulccb
    integer(i4) :: n, ni, nj, i, j, a0, a1, b0, b1

    n = size(a, 2)
    ni = size(b%e, 1)
    nj = size(b%e, 2)
    i = count(b%e(ni,:) .eq. huge(1.D0) + huge(1.D0)*iu)
    j = count(b%e(:,1) .eq. huge(1.D0) + huge(1.D0)*iu)
    if (n .lt. (ni-i)+(nj-j)-1) then
      write(*,*) 'mul(band): incompatible matrix size for multiplication'
      stop
    else
      n = (ni-i) + (nj-j) - 1
    endif
    allocate( c(size(a,1), size(b%e,2)) )
    c = 0.D0
    do i = 1, size(c, 1)
      do j = 1, size(c, 2)
        a0 = max(1, j-b%suplen)
        a1 = min(n, j+b%sublen)
        b0 = max(2+b%suplen-j, 1)
        b1 = min(n+1+b%suplen-j, n)
        c(i,j) = sum(a(i,a0:a1)*b%e(b0:b1,j))
      enddo
    enddo
  end procedure

  module procedure mulc1cb
    c = mulccb(reshape(a, (/size(a), 1/)), b)
  end procedure

  module procedure mulcbc
    c = transpose(mulccb(transpose(b),a%t()))
  end procedure

  module procedure mulcbc1
    complex(p8), dimension(:,:), allocatable :: c2
    
    c2 = mulcbc(a, reshape(b, (/size(b),1/)))
    c = reshape(c2, (/size(c2)/))
    deallocate( c2 )
  end procedure

  module procedure mulrcb
    integer(i4) :: n, ni, nj, i, j, a0, a1, b0, b1

    n = size(a, 2)
    ni = size(b%e, 1)
    nj = size(b%e, 2)
    i = count(b%e(ni,:) .eq. huge(1.D0) + huge(1.D0)*iu)
    j = count(b%e(:,1) .eq. huge(1.D0) + huge(1.D0)*iu)
    if (n .lt. (ni-i)+(nj-j)-1) then
      write(*,*) 'mul(band): incompatible matrix size for multiplication'
      stop
    else
      n = (ni-i) + (nj-j) - 1
    endif
    allocate( c(size(a,1), size(b%e,2)) )
    c = 0.D0
    do i = 1, size(c, 1)
      do j = 1, size(c, 2)
        a0 = max(1, j-b%suplen)
        a1 = min(n, j+b%sublen)
        b0 = max(2+b%suplen-j, 1)
        b1 = min(n+1+b%suplen-j, n)
        c(i,j) = sum(a(i,a0:a1)*b%e(b0:b1,j))
      enddo
    enddo
  end procedure

  module procedure mulr1cb
    c = mulrcb(reshape(a, (/size(a), 1/)), b)
  end procedure

  module procedure mulcbr
    c = transpose(mulrcb(transpose(b),a%t()))
  end procedure

  module procedure mulcbr1
    complex(p8), dimension(:,:), allocatable :: c2
    
    c2 = mulcbr(a, reshape(b, (/size(b),1/)))
    c = reshape(c2, (/size(c2)/))
    deallocate( c2 )
  end procedure

  module procedure mulcrb
    integer(i4) :: n, ni, nj, i, j, a0, a1, b0, b1

    n = size(a, 2)
    ni = size(b%e, 1)
    nj = size(b%e, 2)
    i = count(b%e(ni,:) .eq. huge(1.D0))
    j = count(b%e(:,1) .eq. huge(1.D0))
    if (n .lt. (ni-i)+(nj-j)-1) then
      write(*,*) 'mul(band): incompatible matrix size for multiplication'
      stop
    else
      n = (ni-i) + (nj-j) - 1
    endif
    allocate( c(size(a,1), size(b%e,2)) )
    c = 0.D0
    do i = 1, size(c, 1)
      do j = 1, size(c, 2)
        a0 = max(1, j-b%suplen)
        a1 = min(n, j+b%sublen)
        b0 = max(2+b%suplen-j, 1)
        b1 = min(n+1+b%suplen-j, n)
        c(i,j) = sum(a(i,a0:a1)*b%e(b0:b1,j))
      enddo
    enddo
  end procedure

  module procedure mulc1rb
    c = mulcrb(reshape(a, (/size(a), 1/)), b)
  end procedure

  module procedure mulrbc
    c = transpose(mulcrb(transpose(b),a%t()))
  end procedure

  module procedure mulrbc1
    complex(p8), dimension(:,:), allocatable :: c2
    
    c2 = mulrbc(a, reshape(b, (/size(b),1/)))
    c = reshape(c2, (/size(c2)/))
    deallocate( c2 )
  end procedure

  module procedure mulrrb
    integer(i4) :: n, ni, nj, i, j, a0, a1, b0, b1

    n = size(a, 2)
    ni = size(b%e, 1)
    nj = size(b%e, 2)
    i = count(b%e(ni,:) .eq. huge(1.D0))
    j = count(b%e(:,1) .eq. huge(1.D0))
    if (n .lt. (ni-i)+(nj-j)-1) then
      write(*,*) 'mul(band): incompatible matrix size for multiplication'
      stop
    else
      n = (ni-i) + (nj-j) - 1
    endif
    allocate( c(size(a,1), size(b%e,2)) )
    c = 0.D0
    do i = 1, size(c, 1)
      do j = 1, size(c, 2)
        a0 = max(1, j-b%suplen)
        a1 = min(n, j+b%sublen)
        b0 = max(2+b%suplen-j, 1)
        b1 = min(n+1+b%suplen-j, n)
        c(i,j) = sum(a(i,a0:a1)*b%e(b0:b1,j))
      enddo
    enddo
  end procedure

  module procedure mulr1rb
    c = mulrrb(reshape(a, (/size(a), 1/)), b)
  end procedure

  module procedure mulrbr
    c = transpose(mulrrb(transpose(b),a%t()))
  end procedure

  module procedure mulrbr1
    real(p8), dimension(:,:), allocatable :: c2
    
    c2 = mulrbr(a, reshape(b, (/size(b),1/)))
    c = reshape(c2, (/size(c2)/))
    deallocate( c2 )
  end procedure

  module procedure mulcbcb
    integer :: n, ni, nj, i, j, k
    integer :: a0, a1, ct
    complex(p8) :: c0
    complex(p8), dimension(:,:), allocatable :: ce

    c%suplen = a%suplen + b%suplen
    c%sublen = a%sublen + b%sublen
    n = size(a%e, 2)
    ni = size(b%e, 1)
    nj = size(b%e, 2)
    i = count(b%e(ni,:) .eq. huge(1.D0) + huge(1.D0)*iu)
    j = count(b%e(:,1) .eq. huge(1.D0) + huge(1.D0)*iu)
    if (n .lt. (ni-i)+(nj-j)-1) then
      write(*,*) 'mul(band): incompatible matrix size for multiplication'
      stop
    else
      n = (ni-i)+(nj-j)-1
    endif
    allocate( ce(c%suplen+c%sublen+1, nj) )
    ce = huge(1.D0) + huge(1.D0)*iu
    do i = 1, size(ce, 1)
      do j = 1, size(ce, 2)
        a0 = max(1, j-b%suplen)
        a1 = min(n, j+b%sublen)
        c0 = 0.D0 + 0.D0*iu
        ct = 0
        do k = a0, a1
          if ((i+j-b%suplen-k .ge. 1) .and. (i+j-b%suplen-k .le. size(a%e, 1)) .and. &
              (b%suplen+1+k-j .ge. 1) .and. (b%suplen+1+k-j .le. size(b%e, 1))) then
            if ((a%e(i+j-b%suplen-k,k) .ne. huge(1.D0) + huge(1.D0)*iu) .and. &
                (b%e(b%suplen+1+k-j,j) .ne. huge(1.D0) + huge(1.D0)*iu)) then
              c0 = c0 + a%e(i+j-b%suplen-k,k) * b%e(b%suplen+1+k-j,j)
              ct = ct + 1
            endif
          endif
        enddo
        if (ct .ne. 0) then
          ce(i,j) = c0
        endif
      enddo
    enddo
    i = 1 ! truncate unnecessary superdiagonals
    do while (.true.)
      if (all( (ce(i,:) .eq. huge(1.D0) + 0.D0*iu) .or. &
               (ce(i,:) .eq. 0.D0 + huge(1.D0)*iu) .or. &
               (ce(i,:) .eq. huge(1.D0) + huge(1.D0)*iu) .or. &
               (ce(i,:) .eq. 0.D0 + 0.D0*iu))) then
        if (c%suplen .eq. 0) exit
        c%suplen = c%suplen - 1
        i = i+1
      else
        exit
      endif
    enddo
    j = size(ce, 1) ! truncate unnecessary subdiagonals
    do while (.true.)
      if (all( (ce(i,:) .eq. huge(1.D0) + 0.D0*iu) .or. &
               (ce(i,:) .eq. 0.D0 + huge(1.D0)*iu) .or. &
               (ce(i,:) .eq. huge(1.D0) + huge(1.D0)*iu) .or. &
               (ce(i,:) .eq. 0.D0 + 0.D0*iu))) then
        if (c%sublen .eq. 0) exit
        c%sublen = c%sublen - 1
        j = j-1
      else
        exit
      endif
    enddo
    allocate( c%e(j-i+1, nj) )
    c%e = ce(i:j,:)
    deallocate( ce )
  end procedure

  module procedure mulrbrb
    integer :: n, ni, nj, i, j, k
    integer :: a0, a1, ct
    real(p8) :: c0
    real(p8), dimension(:,:), allocatable :: ce

    c%suplen = a%suplen + b%suplen
    c%sublen = a%sublen + b%sublen
    n = size(a%e, 2)
    ni = size(b%e, 1)
    nj = size(b%e, 2)
    i = count(b%e(ni,:) .eq. huge(1.D0))
    j = count(b%e(:,1) .eq. huge(1.D0))
    if (n .lt. (ni-i)+(nj-j)-1) then
      write(*,*) 'mul(band): incompatible matrix size for multiplication'
      stop
    else
      n = (ni-i)+(nj-j)-1
    endif
    allocate( ce(c%suplen+c%sublen+1, nj) )
    ce = huge(1.D0)
    do i = 1, size(ce, 1)
      do j = 1, size(ce, 2)
        a0 = max(1, j-b%suplen)
        a1 = min(n, j+b%sublen)
        c0 = 0.D0
        ct = 0
        do k = a0, a1
          if ((i+j-b%suplen-k .ge. 1) .and. (i+j-b%suplen-k .le. size(a%e, 1)) .and. &
              (b%suplen+1+k-j .ge. 1) .and. (b%suplen+1+k-j .le. size(b%e, 1))) then
            if ((a%e(i+j-b%suplen-k,k) .ne. huge(1.D0)) .and. &
                (b%e(b%suplen+1+k-j,j) .ne. huge(1.D0))) then
              c0 = c0 + a%e(i+j-b%suplen-k,k) * b%e(b%suplen+1+k-j,j)
              ct = ct + 1
            endif
          endif
        enddo
        if (ct .ne. 0) then
          ce(i,j) = c0
        endif
      enddo
    enddo
    i = 1 ! truncate unnecessary superdiagonals
    do while (.true.)
      if (all( (ce(i,:) .eq. huge(1.D0)) .or. &
               (ce(i,:) .eq. 0.D0))) then
        if (c%suplen .eq. 0) exit
        c%suplen = c%suplen - 1
        i = i+1
      else
        exit
      endif
    enddo
    j = size(ce, 1) ! truncate unnecessary subdiagonals
    do while (.true.)
      if (all( (ce(i,:) .eq. huge(1.D0)) .or. &
               (ce(i,:) .eq. 0.D0))) then
        if (c%sublen .eq. 0) exit
        c%sublen = c%sublen - 1
        j = j-1
      else
        exit
      endif
    enddo
    allocate( c%e(j-i+1, nj) )
    c%e = ce(i:j,:)
    deallocate( ce )
  end procedure

  module procedure mulrbcb
    type(cmpx_bndm) :: a_
    a_ = a
    c = mulcbcb(a_, b)
  end procedure

  module procedure mulcbrb
    type(cmpx_bndm) :: b_
    b_ = b
    c = mulcbcb(a, b_)
  end procedure

  module procedure lurb
    external :: dgbtrf
    integer(i4), dimension(:), allocatable :: ipiv_
    integer(i4) :: m, n, ku, kl, info, ni, nj, i, j
    real(p8), dimension(:,:), allocatable :: ae

    ku = a%suplen
    kl = a%sublen
    ni = size(a%e, 1)
    nj = size(a%e, 2)
    i = count(a%e(ni,:) .eq. huge(1.D0))
    j = count(a%e(:,1) .eq. huge(1.D0))
    m = (ni-i)+(nj-j)-1
    n = nj
    if (present(ipiv)) then
      ipiv_ = ipiv
    else
      allocate( ipiv_(max(1, min(m,n))) )
    endif
    allocate( ae(2*kl+ku+1, n) )
    ae = huge(1.D0)
    where (a%e .ne. huge(1.D0)) ae(kl+1:,:) = a%e
    if (size(ipiv_) .ne. n) then
      write(*,*) 'lurc: inconsistency in the ipiv length'
      stop
    endif
    call dgbtrf(m, n, kl, ku, ae, 2*kl+ku+1, ipiv_, info)
    if (info .ne. 0) then
      write(*,*) 'lurc: lu factorization resulted in failure'
      stop
    endif
    call move_alloc(ae, a%e)
    a%suplen = kl + ku
    if (present(ipiv)) ipiv = ipiv_
    deallocate( ipiv_ )
  end procedure

  module procedure lucb
    external :: zgbtrf
    integer(i4), dimension(:), allocatable :: ipiv_
    integer(i4) :: m, n, ku, kl, info, ni, nj, i, j
    complex(p8), dimension(:,:), allocatable :: ae

    ku = a%suplen
    kl = a%sublen
    ni = size(a%e, 1)
    nj = size(a%e, 2)
    i = count(a%e(ni,:) .eq. huge(1.D0) + huge(1.D0)*iu)
    j = count(a%e(:,1) .eq. huge(1.D0) + huge(1.D0)*iu)
    m = (ni-i)+(nj-j)-1
    n = nj
    if (present(ipiv)) then
      ipiv_ = ipiv
    else
      allocate( ipiv_(max(1, min(m,n))) )
    endif
    allocate( ae(2*kl+ku+1, n) )
    ae = huge(1.D0) + huge(1.D0)*iu
    where (a%e .ne. huge(1.D0) + huge(1.D0)*iu) ae(kl+1:,:) = a%e
    if (size(ipiv_) .ne. n) then
      write(*,*) 'lurc: inconsistency in the ipiv length'
      stop
    endif
    call zgbtrf(m, n, kl, ku, ae, 2*kl+ku+1, ipiv_, info)
    if (info .ne. 0) then
      write(*,*) 'lurc: lu factorization resulted in failure'
      stop
    endif
    call move_alloc(ae, a%e)
    a%suplen = kl + ku
    if (present(ipiv)) ipiv = ipiv_
    deallocate( ipiv_ )
  end procedure

  module procedure solvecbc
    external :: zgbtrs
    integer(i4), dimension(:), allocatable :: ipiv
    integer(i4) :: n, nrhs, kl, ku, info, ni, nj, i, j

    if ( size(a%e, 2) .ne. size(b, 1) ) then
      write (*,*) 'solvecbc: matrix size inconsistent'
      stop
    endif
    n = size(b, 1)
    nrhs = size(b, 2)
    kl = a%sublen
    ku = a%suplen
    ni = size(a%e, 1)
    nj = size(a%e, 2)
    i = count(a%e(ni,:) .eq. huge(1.d0))
    j = count(a%e(:,1) .eq. huge(1.d0))
    allocate( ipiv(min((ni-i)+(nj-j)-1, nj)) )
    call lucb(a, ipiv)
    call zgbtrs('n', n, kl, ku, nrhs, a%e, 2*kl+ku+1, ipiv, b, n, info)
    if (info .ne. 0) then
      write(*,*) 'solvecbc: linear system unable to be solved'
      stop
    endif
  end procedure

  module procedure solvecbc1
    complex(p8), dimension(:,:), allocatable :: b_

    b_ = reshape(b, (/size(b),1/))
    call solvecbc(a, b_)
    b = reshape(b_,(/size(b)/))
    deallocate( b_ )
  end procedure

  module procedure solverbc
    type(cmpx_bndm) :: a_
    
    a_ = a
    call solvecbc(a_, b)
    a = a_
    call a_%dealloc()
  end procedure

  module procedure solverbc1
    type(cmpx_bndm) :: a_
    
    a_ = a
    call solvecbc1(a_, b)
    a = a_
    call a_%dealloc()
  end procedure

  module procedure solverbr
    external :: dgbtrs
    integer(i4), dimension(:), allocatable :: ipiv
    integer(i4) :: n, nrhs, kl, ku, info, ni, nj, i, j

    if ( size(a%e, 2) .ne. size(b, 1) ) then
      write (*,*) 'solverbr: matrix size inconsistent'
      stop
    endif
    n = size(b, 1)
    nrhs = size(b, 2)
    kl = a%sublen
    ku = a%suplen
    ni = size(a%e, 1)
    nj = size(a%e, 2)
    i = count(a%e(ni,:) .eq. huge(1.d0))
    j = count(a%e(:,1) .eq. huge(1.d0))
    allocate( ipiv(min((ni-i)+(nj-j)-1, nj)) )
    call lurb(a, ipiv)
    call dgbtrs('n', n, kl, ku, nrhs, a%e, 2*kl+ku+1, ipiv, b, n, info)
    if (info .ne. 0) then
      write(*,*) 'solverbr: linear system unable to be solved.'
      stop
    endif
  end procedure

  module procedure solverbr1
    real(p8), dimension(:,:), allocatable :: b_

    b_ = reshape(b, (/size(b),1/))
    call solverbr(a, b_)
    b = reshape(b_,(/size(b)/))
    deallocate( b_ )
  end procedure

end submodule