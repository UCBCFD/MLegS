submodule (mlegs_genmat) mlegs_genmat_ops
  implicit none

  !> zero
  real(p8), parameter :: zero_real = 0.D0
  complex(p8), parameter :: zero_cmpx = 0.D0 + 0.D0*iu
  !> one
  real(p8), parameter :: one_real = 1.D0
  complex(p8), parameter :: one_cmpx = 1.D0 + 0.D0*iu

contains

  module procedure mulrr
    external :: dgemm
    integer(i4) :: i, j, k
    if (size(a,2) .ne. size(b,1)) then
      write (*,*) 'mul: invalid shape.'
      write (*,*) 'argument1 = ',size(a,1),'x',size(a,2)
      write (*,*) 'argument2 = ',size(b,1),'x',size(b,2)
    endif
    i = size(a,1)
    j = size(b,2)
    k = size(a,2)
    call dgemm('n', 'n', i, j, k , one_real, a, i, b, k, zero_real, c, i)
  end procedure

  module procedure mulcc
    external :: zgemm
    integer(i4) :: i, j, k
    if (size(a,2) .ne. size(b,1)) then
      write (*,*) 'mul: invalid shape.'
      write (*,*) 'argument1 = ',size(a,1),'x',size(a,2)
      write (*,*) 'argument2 = ',size(b,1),'x',size(b,2)
    endif
    i = size(a,1)
    j = size(b,2)
    k = size(a,2)
    call zgemm('n', 'n', i, j, k , one_cmpx, a, i, b, k, zero_cmpx, c, i)
  end procedure

  module procedure mulrc
    external :: zgemm
    integer(i4) :: i, j, k
    complex(p8), dimension(:,:), allocatable :: a_
    if (size(a,2) .ne. size(b,1)) then
      write (*,*) 'mul: invalid shape.'
      write (*,*) 'argument1 = ',size(a,1),'x',size(a,2)
      write (*,*) 'argument2 = ',size(b,1),'x',size(b,2)
    endif
    i = size(a,1)
    j = size(b,2)
    k = size(a,2)
    allocate(a_(i,j))
    a_ = a
    call zgemm('n', 'n', i, j, k , one_cmpx, a_, i, b, k, zero_cmpx, c, i)
    deallocate(a_)
  end procedure

  module procedure mulcr
    external :: zgemm
    integer(i4) :: i, j, k
    complex(p8), dimension(:,:), allocatable :: b_
    if (size(a,2) .ne. size(b,1)) then
      write (*,*) 'mul: invalid shape.'
      write (*,*) 'argument1 = ',size(a,1),'x',size(a,2)
      write (*,*) 'argument2 = ',size(b,1),'x',size(b,2)
    endif
    i = size(a,1)
    j = size(b,2)
    k = size(a,2)
    allocate(b_(j,k))
    b_ = b
    call zgemm('n', 'n', i, j, k , one_cmpx, a, i, b_, k, zero_cmpx, c, i)
    deallocate(b_)
  end procedure

  module procedure mulrr1
    real(p8), dimension(:,:), allocatable :: b_
    allocate(b_(size(b),1))
    b_ = reshape(b,(/ size(b), 1 /))
    c = reshape(mulrr(a, b_), (/ size(a,1) /))
    deallocate(b_)
  end procedure

  module procedure mulcc1
    complex(p8), dimension(:,:), allocatable :: b_
    allocate(b_(size(b),1))
    b_ = reshape(b,(/ size(b), 1 /))
    c = reshape(mulcc(a, b_), (/ size(a,1) /))
    deallocate(b_)
  end procedure

  module procedure mulrc1
    complex(p8), dimension(:,:), allocatable :: b_
    allocate(b_(size(b),1))
    b_ = reshape(b,(/ size(b), 1 /))
    c = reshape(mulrc(a, b_), (/ size(a,1) /))
    deallocate(b_)
  end procedure

  module procedure mulcr1
    real(p8), dimension(:,:), allocatable :: b_
    allocate(b_(size(b),1))
    b_ = reshape(b,(/ size(b), 1 /))
    c = reshape(mulcr(a, b_), (/ size(a,1) /))
    deallocate(b_)
  end procedure

  module procedure mulr1r
    real(p8), dimension(:,:), allocatable :: a_
    allocate(a_(size(a),1))
    a_ = reshape(a,(/ size(a), 1 /))
    c = mulrr(a_, b)
    deallocate(a_)
  end procedure

  module procedure mulc1c
    complex(p8), dimension(:,:), allocatable :: a_
    allocate(a_(size(a),1))
    a_ = reshape(a,(/ size(a), 1 /))
    c = mulcc(a_, b)
    deallocate(a_)
  end procedure

  module procedure mulr1c
    real(p8), dimension(:,:), allocatable :: a_
    allocate(a_(size(a),1))
    a_ = reshape(a,(/ size(a), 1 /))
    c = mulrc(a_, b)
    deallocate(a_)
  end procedure

  module procedure mulc1r
    complex(p8), dimension(:,:), allocatable :: a_
    allocate(a_(size(a),1))
    a_ = reshape(a,(/ size(a), 1 /))
    c = mulcr(a_, b)
    deallocate(a_)
  end procedure

  module procedure lur
    external :: dgetrf
    integer(i4), dimension(:), allocatable :: ipiv_
    integer(i4) :: m, n, info

    m = size(a,1)
    n = size(a,2)
    if (present(ipiv)) then
      ipiv_ = ipiv
    else
      allocate( ipiv_(max(1, min(m,n))) )
    endif
    if (size(ipiv_) .ne. n) then
      write(*,*) 'lur: inconsistency in the ipiv length'
      stop
    endif
    call dgetrf(m, n, a, m, ipiv_, info)
    if (info .ne. 0) then
      write(*,*) 'lur: lu factorization resulted in failure'
      stop
    endif
    if (present(ipiv)) ipiv = ipiv_
    deallocate( ipiv_ )
  end procedure

  module procedure luc
    external :: zgetrf
    integer(i4), dimension(:), allocatable :: ipiv_
    integer(i4) :: m, n, info

    m = size(a,1)
    n = size(a,2)
    if (present(ipiv)) then
      ipiv_ = ipiv
    else
      allocate( ipiv_(max(1, min(m,n))) )
    endif
    if (size(ipiv_) .ne. n) then
      write(*,*) 'lur: inconsistency in the ipiv length'
      stop
    endif
    call zgetrf(m, n, a, m, ipiv_, info)
    if (info .ne. 0) then
      write(*,*) 'lur: lu factorization resulted in failure'
      stop
    endif
    if (present(ipiv)) ipiv = ipiv_
    deallocate( ipiv_ )
  end procedure

  module procedure solvecc
    external :: zgetrs
    integer(i4), dimension(min(size(a,1), size(a,2))) :: ipiv
    integer(i4) :: ni, nj, info

    if (size(a, 1) .ne. size(a, 2)) then
      write (*,*) 'solvecc: matrix not square'
      write (*,*) size(a,1),'x',size(a,2)
      stop
    endif
    if ( size(a, 2) .ne. size(b, 1) ) then
      write (*,*) 'solvecc: matrix size inconsistent'
      stop
    endif
    ni = size(a, 1)
    nj = size(b, 2)
    call luc(a, ipiv)
    call zgetrs('n', ni, nj, a, ni, ipiv, b, ni, info)
    if (info .ne. 0) then
      write(*,*) 'solvecc: linear system unable to be solved'
      stop
    endif
  end procedure

  module procedure solvecc1
    complex(p8), dimension(:,:), allocatable :: b_

    b_ = reshape(b, (/size(b),1/))
    call solvecc(a, b_)
    b = reshape(b_,(/size(b)/))
    deallocate( b_ )
  end procedure

  module procedure solverc
    complex(p8), dimension(:,:), allocatable :: a_

    a_ = a
    call solvecc(a_, b)
    a_ = a
    deallocate( a_ )
  end procedure

  module procedure solverc1
    complex(p8), dimension(:,:), allocatable :: a_

    a_ = a
    call solvecc1(a_, b)
    a_ = a
    deallocate( a_ )
  end procedure


  module procedure solverr
    external :: dgetrs
    integer(i4), dimension(min(size(a,1), size(a,2))) :: ipiv
    integer(i4) :: ni, nj, info

    if (size(a, 1) .ne. size(a, 2)) then
      write (*,*) 'solverr: matrix not square'
      write (*,*) size(a,1),'x',size(a,2)
      stop
    endif
    if ( size(a, 2) .ne. size(b, 1) ) then
      write (*,*) 'solverr: matrix size inconsistent'
      stop
    endif
    ni = size(a, 1)
    nj = size(b, 2)
    call lur(a, ipiv)
    call dgetrs('n', ni, nj, a, ni, ipiv, b, ni, info)
    if (info .ne. 0) then
      write(*,*) 'solvecc: linear system unable to be solved'
      stop
    endif
  end procedure

  module procedure solverr1
    real(p8), dimension(:,:), allocatable :: b_

    b_ = reshape(b, (/size(b),1/))
    call solverr(a, b_)
    b = reshape(b_,(/size(b)/))
    deallocate( b_ )
  end procedure

  module procedure invr
    external :: dgetrf, dgetri
    real(p8), dimension(size(a,1)) :: wk
    integer(i4), dimension(size(a,1)) :: ip
    integer(i4) :: ni, nj, info

    ni = size(a,1)
    nj = size(a,2)
    if (ni .ne. nj) then
      write(*,*) 'invdr: non-square matrix'
      stop
    endif
    ainv = a
    call dgetrf(ni, ni, ainv, ni, ip, info)
    if (info .ne. 0) then
      write(*,*) 'invdr: singular matrix. cannot be inverted'
      stop
    endif
    call dgetri(ni, ainv, ni, ip, wk, ni, info)
    if (info .ne. 0) then
      write(*,*) 'invdr: matrix inversion failed due to an unknown reason'
      stop
    endif
  end procedure

  module procedure invc
    external :: zgetrf, zgetri
    complex(p8), dimension(size(a,1)) :: wk
    integer(i4), dimension(size(a,1)) :: ip
    integer(i4) :: ni, nj, info

    ni = size(a,1)
    nj = size(a,2)
    if (ni .ne. nj) then
      write(*,*) 'invdr: non-square matrix'
      stop
    endif
    ainv = a
    call zgetrf(ni, ni, ainv, ni, ip, info)
    if (info .ne. 0) then
      write(*,*) 'invdr: singular matrix. cannot be inverted'
      stop
    endif
    call zgetri(ni, ainv, ni, ip, wk, ni, info)
    if (info .ne. 0) then
      write(*,*) 'invdr: matrix inversion failed due to an unknown reason'
      stop
    endif
  end procedure

  module procedure geigrr
    external :: dggev
    real(p8), dimension(:,:), allocatable :: w1, w2, vl, vr
    real(p8), dimension(:), allocatable :: alphar, alphai, beta, wk
    integer(i4) :: ni, lwk, info, i
    character(len=1) :: r, l

    if (size(a,1) .ne. size(a,2)) then
      write (*,*) 'geigrr: first matrix not square.'
      write (*,*) size(a,1),'x',size(a,2)
      stop
    endif
    if (size(b,1) .ne. size(b,2)) then
      write (*,*) 'geigrr: second matrix not square.'
      write (*,*) size(b,1),'x',size(b,2)
      stop
    endif
    if (size(a,1) .ne. size(b,1)) then
      write (*,*) 'geigrr: input matrices incompatible in size.'
      stop
    endif
    ni = size(a,1)
    lwk = 8*ni  
    allocate( w1(ni,ni), w2(ni,ni), vl(ni,ni), vr(ni,ni) )
    allocate( alphar(ni), alphai(ni), beta(ni) )
    allocate( wk(lwk) )
    w1 = a
    w2 = b
    l = 'N'
    if (present(el)) l = 'V'
    r = 'N'
    if (present(er)) r = 'V'
    call dggev(l, r, ni, w1, ni, w2, ni, alphar, alphai, beta, vl, ni, vr, ni, wk, lwk, info)
    if (info .ne. 0) then
      write(*,*) 'geigrr: evp solver resulted in failure'
      stop
    endif

    do i = 1, ni
      if (abs(beta(i)).eq.0.D0) then
        ev(i) = huge(1.D0) + huge(1.D0)*iu
      else
        ev(i) = (alphar(i) + alphai(i)*iu)/beta(i)
      endif
    enddo
    if (present(er)) then
      i = 1
      do while (i .le. ni)
        if (aimag(ev(i)) .eq. 0.D0) then
          er(:,i) = vr(:,i)
          i = i+1
        else
          er(:,i) = vr(:,i) + vr(:,i+1)*iu
          er(:,i+1) = vr(:,i) - vr(:,i+1)*iu
          i = i+2
        endif
      enddo
    endif
    if (present(el)) then
      i = 1
      do while (i .le. ni)
        if (aimag(ev(i)) .eq. 0.D0) then
          el(:,i) = vl(:,i)
          i = i+1
        else
          el(:,i) = vl(:,i) + vl(:,i+1)*iu
          el(:,i+1) = vl(:,i) - vl(:,i+1)*iu
          i = i+2
        endif
      enddo
    endif
    deallocate( w1, w2, vl, vr, alphar, alphai, beta, wk )
  end procedure

  module procedure steigr
    external :: dgeev
    real(p8), dimension(:,:), allocatable :: w1, vl, vr
    real(p8), dimension(:), allocatable :: alphar, alphai, wk
    integer(i4) :: ni, lwk, info, i
    character(len=1) :: r, l

    if (size(a,1) .ne. size(a,2)) then
      write (*,*) 'steigr: first matrix not square.'
      write (*,*) size(a,1),'x',size(a,2)
      stop
    endif
    ni = size(a,1)
    lwk = 4*ni  
    allocate( w1(ni,ni), vl(ni,ni), vr(ni,ni) )
    allocate( alphar(ni), alphai(ni) )
    allocate( wk(lwk) )
    w1 = a
    l = 'N'
    if (present(el)) l = 'V'
    r = 'N'
    if (present(er)) r = 'V'
    call dgeev(l, r, ni, w1, ni, alphar, alphai, vl, ni, vr, ni, wk, lwk, info)
    if (info .ne. 0) then
      write(*,*) 'steigr: evp solver resulted in failure'
      stop
    endif

    do i = 1, ni
      ev(i) = alphar(i) + alphai(i)*iu
    enddo
    if (present(er)) then
      i = 1
      do while (i .le. ni)
        if (aimag(ev(i)) .eq. 0.D0) then
          er(:,i) = vr(:,i)
          i = i+1
        else
          er(:,i) = vr(:,i) + vr(:,i+1)*iu
          er(:,i+1) = vr(:,i) - vr(:,i+1)*iu
          i = i+2
        endif
      enddo
    endif
    if (present(el)) then
      i = 1
      do while (i .le. ni)
        if (aimag(ev(i)) .eq. 0.D0) then
          el(:,i) = vl(:,i)
          i = i+1
        else
          el(:,i) = vl(:,i) + vl(:,i+1)*iu
          el(:,i+1) = vl(:,i) - vl(:,i+1)*iu
          i = i+2
        endif
      enddo
    endif
    deallocate( w1, vl, vr, alphar, alphai, wk )
  end procedure

  module procedure geigcc
    external :: zggev
    complex(p8), dimension(:,:), allocatable :: w1, w2, vl, vr
    complex(p8), dimension(:), allocatable :: alpha, beta, wk
    real(p8), dimension(:), allocatable :: rwk
    integer(i4) :: ni, lwk, info, i
    character(len=1) :: r, l

    if (size(a,1) .ne. size(a,2)) then
      write (*,*) 'geigrr: first matrix not square.'
      write (*,*) size(a,1),'x',size(a,2)
      stop
    endif
    if (size(b,1) .ne. size(b,2)) then
      write (*,*) 'geigrr: second matrix not square.'
      write (*,*) size(b,1),'x',size(b,2)
      stop
    endif
    if (size(a,1) .ne. size(b,1)) then
      write (*,*) 'geigrr: input matrices incompatible in size.'
      stop
    endif
    ni = size(a,1)
    lwk = 2*ni  
    allocate( w1(ni,ni), w2(ni,ni), vl(ni,ni), vr(ni,ni) )
    allocate( alpha(ni), beta(ni) )
    allocate( wk(lwk), rwk(4*lwk) )
    w1 = a
    w2 = b
    l = 'N'
    if (present(el)) l = 'V'
    r = 'N'
    if (present(er)) r = 'V'
    call zggev(l, r, ni, w1, ni, w2, ni, alpha, beta, vl, ni, vr, ni, wk, lwk, rwk, info)
    if (info .ne. 0) then
      write(*,*) 'geigrr: evp solver resulted in failure'
      stop
    endif
    do i = 1, ni
      if (abs(beta(i)).eq.0.D0) then
        ev(i) = huge(1.D0) + huge(1.D0)*iu
      else
        ev(i) = alpha(i)/beta(i)
      endif
      if (present(er)) er(:,i) = vr(:,i)
      if (present(el)) el(:,i) = vl(:,i)
    enddo
    deallocate( w1, w2, vl, vr, alpha, beta, wk, rwk )
  end procedure

  module procedure geigcr
    complex(p8), dimension(:,:), allocatable :: b_, er_, el_
    b_ = b
    allocate( er_(size(a,1),size(a,1)), el_(size(a,1),size(a,1)) )
    call geigcc(a, b_, ev, er_, el_)
    if (present(er)) er = er_
    if (present(el)) el = el_
    deallocate(b_, er_, el_)
  end procedure

  module procedure geigrc
    complex(p8), dimension(:,:), allocatable :: a_, er_, el_
    a_ = a
    allocate( er_(size(a,1),size(a,1)), el_(size(a,1),size(a,1)) )
    call geigcc(a_, b, ev, er_, el_)
    if (present(er)) er = er_
    if (present(el)) el = el_
    deallocate(a_, er_, el_)
  end procedure

  module procedure steigc
    external :: zgeev
    complex(p8), dimension(:,:), allocatable :: w1, vl, vr
    complex(p8), dimension(:), allocatable :: alpha, wk
    real(p8), dimension(:), allocatable :: rwk
    integer(i4) :: ni, lwk, info, i
    character(len=1) :: r, l

    if (size(a,1) .ne. size(a,2)) then
      write (*,*) 'steigc: first matrix not square.'
      write (*,*) size(a,1),'x',size(a,2)
      stop
    endif
    ni = size(a,1)
    lwk = 2*ni  
    allocate( w1(ni,ni), vl(ni,ni), vr(ni,ni) )
    allocate( alpha(ni) )
    allocate( wk(lwk), rwk(lwk) )
    w1 = a
    l = 'N'
    if (present(el)) l = 'V'
    r = 'N'
    if (present(er)) r = 'V'
    call zgeev(l, r, ni, w1, ni, alpha, vl, ni, vr, ni, wk, lwk, rwk, info)
    if (info .ne. 0) then
      write(*,*) 'steigc: evp solver resulted in failure'
      stop
    endif
    do i = 1, ni
      ev(i) = alpha(i)
      if (present(er)) er(:,i) = vr(:,i)
      if (present(el)) el(:,i) = vl(:,i)
    enddo
    deallocate( w1, vl, vr, alpha, wk, rwk )
  end procedure

end submodule