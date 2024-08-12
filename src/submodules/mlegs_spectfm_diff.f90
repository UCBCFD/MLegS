submodule (mlegs_spectfm) mlegs_spectfm_diff
  implicit none

contains

  module procedure leg_xxdx_with_tfm_kit
    integer(i4) :: nn, n, i, j, am
    real(p8), dimension(:,:), allocatable :: genm_xxdx

    select type(tfm_kit)
      type is (tfm_kit_1d)
        if (mval .ne. 0) then
          stop 'leg_xxdx_with_tfm_kit: m out of the transformation kit'
        endif
      type is (tfm_kit_2d)
        if (findloc(tfm_kit%m, mval, dim=1) .eq. 0) then
         stop 'leg_xxdx_with_tfm_kit: m out of the transformation kit'
        endif
      type is (tfm_kit_3d)
        if (findloc(tfm_kit%m, mval, dim=1) .eq. 0) then
          stop 'leg_xxdx_with_tfm_kit: m out of the transformation kit'
        endif
    end select

    am = abs(mval)
    allocate( genm_xxdx(n_input, n_output) )
    genm_xxdx = 0.D0
    do nn = 1, n_input
      n = am + (nn - 1)
      if ((nn-1) .ge. 1) then
        genm_xxdx(nn,nn-1) = -(n-1.D0)*(n-am)/(2.D0*n-1.D0) ! subdiagonal
      endif
      if ((nn+1) .le. n_output) then
        genm_xxdx(nn,nn+1) = (n+2.D0)*(n+am+1.D0)/(2.D0*n+3.D0) ! superdiagonal
      endif
    enddo
    do i = 1, n_input
      do j = 1, n_output
        genm_xxdx(i,j) = genm_xxdx(i,j) * exp(tfm_kit%lognorm(j, am+1) - tfm_kit%lognorm(i, am+1))
      enddo
    enddo
    xxdx = bandmat(genm_xxdx, sublen=1, suplen=1)
  end procedure

  module procedure leg_del2h_with_tfm_kit
    integer(i4) :: nn, n, i, j, am
    real(p8) :: ell2, s
    real(p8), dimension(:,:), allocatable :: genm_del2h

    select type(tfm_kit)
      type is (tfm_kit_1d)
        if (mval .ne. 0) then
          stop 'leg_del2h_with_tfm_kit: m out of the transformation kit'
        endif
      type is (tfm_kit_2d)
        if (findloc(tfm_kit%m, mval, dim=1) .eq. 0) then
         stop 'leg_del2h_with_tfm_kit: m out of the transformation kit'
        endif
      type is (tfm_kit_3d)
        if (findloc(tfm_kit%m, mval, dim=1) .eq. 0) then
          stop 'leg_del2h_with_tfm_kit: m out of the transformation kit'
        endif
    end select

    ell2 = ell**2.D0
    s = 0.D0
    am = abs(mval)
    allocate( genm_del2h(n_input, n_output) )
    genm_del2h = 0.D0
    do nn = 1, n_input
      n = am + (nn - 1)
      if ((nn-2) .ge. 1) then
        genm_del2h(nn,nn-2) = -(n-am-1.D0)*(n-am)*(n-2.D0+s)*(n-1.D0+s)/(2.D0*n-3.D0)/(2.D0*n-1.D0)
      endif
      if ((nn-1) .ge. 1) then
        genm_del2h(nn,nn-1) = 2.D0*n*(n-am)*(n-1.D0+s)/(2.D0*n-1.D0)
      endif
      if (nn .le. n_output) then
        genm_del2h(nn,nn) = (-2.D0*n*(n+1.D0)*(3.D0*n*n+3.D0*n-am*am-2.D0)+ &
                           2.D0*s*(s-2.D0)*(n*n+n+am*am-1.D0))/(2.D0*n-1.D0)/(2.D0*n+3.D0)
      endif  
      if ((nn+1) .le. n_output) then
        genm_del2h(nn,nn+1) = 2.D0*(n+1.D0)*(n+am+1.D0)*(n+2.D0-s)/(2.D0*n+3.D0)
      endif
      if ((nn+2) .le. n_output) then
        genm_del2h(nn,nn+2) = -(n+am+1.D0)*(n+am+2.D0)*(n+3.D0-s)*(n+2.D0-s)/(2.D0*n+3.D0)/(2.D0*n+5.D0)
      endif
    enddo
    genm_del2h = genm_del2h/ell2
    do i = 1, n_input
      do j = 1, n_output
        genm_del2h(i,j) = genm_del2h(i,j) * exp(tfm_kit%lognorm(j, am+1) - tfm_kit%lognorm(i, am+1))
      enddo
    enddo
    del2h = bandmat(genm_del2h, sublen=2, suplen=2)
  end procedure

  module procedure leg_del2_with_tfm_kit
    integer(i4) :: nn, n, i, j, am
    real(p8) :: ell2, s
    real(p8), dimension(:,:), allocatable :: genm_del2

    select type(tfm_kit)
      type is (tfm_kit_1d)
        if (mval .ne. 0) then
          stop 'leg_del2_with_tfm_kit: m out of the transformation kit'
        endif
      type is (tfm_kit_2d)
        if (findloc(tfm_kit%m, mval, dim=1) .eq. 0) then
         stop 'leg_del2_with_tfm_kit: m out of the transformation kit'
        endif
      type is (tfm_kit_3d)
        if (findloc(tfm_kit%m, mval, dim=1) .eq. 0) then
          stop 'leg_del2_with_tfm_kit: m out of the transformation kit'
        endif
    end select

    ell2 = ell**2.D0
    s = 0.D0
    am = abs(mval)
    allocate( genm_del2(n_input, n_output) )
    genm_del2 = 0.D0
    do nn = 1, n_input
      n = am + (nn - 1)
      if ((nn-2) .ge. 1) then
        genm_del2(nn,nn-2) = -(n-am-1.D0)*(n-am)*(n-2.D0+s)*(n-1.D0+s)/(2.D0*n-3.D0)/(2.D0*n-1.D0)
      endif
      if ((nn-1) .ge. 1) then
        genm_del2(nn,nn-1) = 2.D0*n*(n-am)*(n-1.D0+s)/(2.D0*n-1.D0)
      endif
      if (nn .le. n_output) then
        genm_del2(nn,nn) = (-2.D0*n*(n+1.D0)*(3.D0*n*n+3.D0*n-am*am-2.D0)+ &
                           2.D0*s*(s-2.D0)*(n*n+n+am*am-1.D0))/(2.D0*n-1.D0)/(2.D0*n+3.D0)
      endif  
      if ((nn+1) .le. n_output) then
        genm_del2(nn,nn+1) = 2.D0*(n+1.D0)*(n+am+1.D0)*(n+2.D0-s)/(2.D0*n+3.D0)
      endif
      if ((nn+2) .le. n_output) then
        genm_del2(nn,nn+2) = -(n+am+1.D0)*(n+am+2.D0)*(n+3.D0-s)*(n+2.D0-s)/(2.D0*n+3.D0)/(2.D0*n+5.D0)
      endif
    enddo
    genm_del2 = genm_del2/ell2
    do nn = 1, n_input
      if (nn .le. n_output) genm_del2(nn,nn) = genm_del2(nn,nn) - akval**2.D0
    enddo
    do i = 1, n_input
      do j = 1, n_output
        genm_del2(i,j) = genm_del2(i,j) * exp(tfm_kit%lognorm(j, am+1) - tfm_kit%lognorm(i, am+1))
      enddo
    enddo
    del2 = bandmat(genm_del2, sublen=2, suplen=2)
  end procedure

end submodule