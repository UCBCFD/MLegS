program addperturb
  ! J. Barranco, 01/25/1999
  ! This program reads in perturbation eigenfunctions which
  ! will be added to the initial base flow.
  ! Required: perturb.input
  ! -------------------------
      USE MOD_MISC
      USE MOD_BANDMAT
      USE MOD_EIG
      USE MOD_FD
      USE MOD_LIN_LEGENDRE
      USE MOD_SCALAR3
      USE MOD_FFT
      USE MOD_LAYOUT
      USE MOD_LEGOPS
      USE MOD_MARCH
      USE MOD_DIAGNOSTICS
      USE MOD_INIT
  ! -------------------------
  implicit none
  character(len=72) :: chifile, psifile
  ! -------------------------
  
  type(scalar):: psi,chi
  integer:: i,j,jj,err,num,nn,flag
  real(p8):: kk, rp1, ip1, rp2, ip2
  integer, allocatable:: mp(:), kp(:)
  complex(p8), allocatable:: perturb_psi(:,:), perturb_chi(:,:)
  complex(p8), allocatable:: scale(:)
  
  call readin(5)
  
  !> read in initial values
  call allocate(psi)
  call allocate(chi)
  call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//files%psii,psi)
  call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//files%chii,chi)
  
  open(10,FILE='new_perturb.input',STATUS='OLD',ACTION='READ',IOSTAT=err)
  if (err.ne.0) then
     print *, 'ERROR: addperturb -- Could not open file perturb.dat'
     STOP
  end if
  
  read(10,*) num
  allocate(perturb_psi(NRchop,num),perturb_chi(NRchop,num),&
       & mp(num),kp(num),scale(num),STAT=err)
  perturb_psi = 0.0_p8
  perturb_chi = 0.0_p8
  do j=1,num
     read(10,*) rp1, ip1
     scale(j) = cmplx(rp1,ip1,p8)
  end do
  do j=1,num
     read(10,*)
     read(10,*)
     read(10,*) nn, mp(j), kk
     do jj=1,8
        read(10,*)
     end do
     kk = zlen/2.0/pi*kk
     kp(j) = nint(kk)
     print *, 'm, k, scale = ', mp(j), kp(j), scale(j)
     if (abs(kp(j)-kk)>1.0E-6) then
        print *, "WARNING: Axial wavenumber is not an integer."
        print *, "         Try resetting value of zlen."
        print *, "         the difference is", abs(kp(j)-kk)
        stop
     end if
     flag = 0
     if (mp(j)<0) then
        mp(j) = -mp(j)
        kp(j) = -kp(j)
        flag  = 1
     end if
     !becoz arrays start from 1, m&k starts from 0
     mp(j) = mp(j)+1
     kp(j) = kp(j)+1

     !if (kp(j)<=0) kp(j)=kp(j)+nxchopdim
     if (kp(j)<=0) kp(j)=kp(j)+2*NXCHOP-1
     do i=1,nn
        if (i<=NRchop) then
           read(10,*) rp1, ip1, rp2, ip2
	   perturb_psi(i,j) = cmplx(rp1,ip1,p8)
	   perturb_chi(i,j) = cmplx(rp2,ip2,p8)
        else
           read(10,*)
        end if
     end do
     perturb_psi(:,j) = perturb_psi(:,j)*scale(j)
     perturb_chi(:,j) = perturb_chi(:,j)*scale(j)
     if (flag==1) then
	     perturb_psi(:,j) = conjg(perturb_psi(:,j))
        perturb_chi(:,j) = conjg(perturb_chi(:,j))
     end if
  end do

  do j=1,num
     psi%e(:NRchop,mp(j),kp(j)) = psi%e(:NRchop,mp(j),kp(j)) + &
          & perturb_psi(:NRchop,j)
     chi%e(:NRchop,mp(j),kp(j)) = chi%e(:NRchop,mp(j),kp(j)) + &
          & perturb_chi(:NRchop,j)
     if ((mp(j)==1).and.(kp(j)>1)) then      ! m = 0; k > 0
        !kp(j) = nxchopdim+2-kp(j)
        kp(j) = 2*NXCHOP+1-kp(j)
        psi%e(:NRchop,mp(j),kp(j)) = psi%e(:NRchop,mp(j),kp(j)) + &
             & conjg(perturb_psi(:NRchop,j))
        chi%e(:NRchop,mp(j),kp(j)) = chi%e(:NRchop,mp(j),kp(j)) + &
             & conjg(perturb_chi(:NRchop,j))
     end if
  end do
  
  where(abs(psi%e)<1.0e-24) psi%e=0.0_p8
  where(abs(chi%e)<1.0e-24) chi%e=0.0_p8

  call msave(psi, TRIM(ADJUSTL(FILES%SAVEDIR))//files%psi0)
  call msave(chi, TRIM(ADJUSTL(FILES%SAVEDIR))//files%chi0)

  call deallocate(psi)
  call deallocate(chi)
  deallocate(perturb_psi, perturb_chi, scale, mp, kp)
  
end program addperturb



