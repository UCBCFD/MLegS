program addperturb_phase
! ======================================================================
! * For MPI use ONLY *
! This program takes a name.input file and reads in the perturbation
! field. The result files are name_psi.dat and name_chi.dat, which can
! be loaded using MLOAD (must use the same amount of MPI tasks and OMP
! threads to load the files).
!
! If Q-vortex sets to True, the program adds the perturbation field to
! psi0/chi0.dat. Otherwise, if Correction sets to True, the perturbation
! field is added to new_perturb_psi/chi.dat
!
! UPDATE HISTORY
! J. Barranco, 01/25/1999: original
! J. Wang    , 01/13/2022: modified for MPI use
! J. Wang    , 03/14/2023: modified for phase adjustment
! ======================================================================
   USE OMP_LIB
   USE MPI
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
! ======================================================================
implicit none
! VARIABLE DECLARATION
character(len=72) :: chifile, psifile    
type(scalar):: psi,chi,psi0,chi0,w
integer:: i,j,jj,err,num,nn,flag
real(p8):: kk, rp1, ip1, rp2, ip2, phase_angle
integer, allocatable:: mp(:), kp(:)
complex(p8), allocatable:: perturb_psi(:,:), perturb_chi(:,:)
complex(p8), allocatable:: scale(:)
! MPI-SPECIFIC DECLARATION
integer:: mm_local, kk_local
! FILE-NAMES
CHARACTER(LEN=200):: FILENAME
logical:: save_q = .FALSE.
! ======================================================================

! START MPI
CALL MPI_INIT(IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)

! START INITILIZATION
CALL READCOM('NOECHO')
CALL READIN(5)
call LEGINIT()
IF (MPI_RANK .EQ. 0) THEN
   WRITE(*,*) 'PROGRAM STARTED'
   CALL PRINT_REAL_TIME()
   WRITE(*,*) 'Perturbation filename: '
   READ(*,10) FILENAME
ENDIF
CALL MPI_BCAST(FILENAME,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
10 FORMAT(A200)

! ALLOCATE SCALARS
call allocate(psi); call allocate(chi)
psi%e = 0.d0; psi%LN = 0.d0
chi%e = 0.d0; chi%LN = 0.d0

! READ IN VALUES
open(10,FILE=TRIM(ADJUSTL(FILENAME))//".input",STATUS='OLD',ACTION='READ',IOSTAT=err)
! open(10,FILE='new_perturb.input',STATUS='OLD',ACTION='READ',IOSTAT=err)
if (err.ne.0) then
   print *, 'ERROR: addperturb -- Could not open file: ',TRIM(ADJUSTL(FILENAME))//".input"
   STOP
end if
read(10,*) num
allocate(perturb_psi(NRchop,num),perturb_chi(NRchop,num),&
   & mp(num),kp(num),scale(num),STAT=err)
perturb_psi = 0.0_p8; perturb_chi = 0.0_p8
! For triad input, ask for phase angle.
if (num.eq.3) then
   IF (MPI_RANK .EQ. 0) THEN
      WRITE(*,*) 'Degen pair phase angle (in deg): '
      read(*,*) phase_angle
      WRITE(*,*) 'Adding phase anlge (in rad):', phase_angle*PI/180
   ENDIF
   CALL MPI_BCAST(phase_angle,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
endif
do j=1,num
   read(10,*) rp1, ip1
   scale(j) = cmplx(rp1,ip1,p8)
   if ((num.eq.3).and.(j.gt.1)) scale(j) = scale(j)*exp(iu*phase_angle*PI/180)
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
   if (MPI_RANK.EQ.0) THEN
   ! print *, 'm, k, scale = ', mp(j), kp(j), scale(j)
   WRITE(*,20) mp(j), kp(j), scale(j)
   20 FORMAT('m, k, scale = ',I0,', ',I0,', ',(F6.3,SP,F6.3,"i"))
   endif
   if (abs(kp(j)-kk)>1.0E-6) then
      if (MPI_RANK.EQ.0) THEN
      print *, "WARNING: Axial wavenumber is not an integer."
      print *, "         Try resetting value of zlen."
      print *, "         the difference is", abs(kp(j)-kk)
      endif
      stop
   elseif (abs(kp(j)).ge.NXCHOP) then
      if (MPI_RANK.EQ.0) THEN
      print *, "WARNING: Axial wave number exceeds number of axial modes."
      print *, "         Try resetting NXCHOP (currently:", NXCHOP," )."
      endif
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

! ADDING VALUES TO SCALARS
psi%e = 0.d0; chi%e = 0.d0
do j=1,num
mm_local = mp(j)-psi%inth
kk_local = kp(j)-psi%inx
IF (0.lt.mm_local .and. size(psi%e,2).ge.mm_local) THEN
   IF (0.lt.kk_local .and. size(psi%e,3).ge.kk_local) THEN
      ! WRITE(*,*) mp(j),kp(j),psi%inth,psi%inx
      WRITE(*,*) m(mp(j)),ak(mp(j),kp(j)),'RANK#',MPI_RANK!,size(psi%e,2),mm_local,size(psi%e,3),kk_local
      psi%e(:NRchop,mm_local,kk_local) = psi%e(:NRchop,mm_local,kk_local) + &
            & perturb_psi(:NRchop,j)
      chi%e(:NRchop,mm_local,kk_local) = chi%e(:NRchop,mm_local,kk_local) + &
            & perturb_chi(:NRchop,j)
      ! if (sum(abs(chi%e(:NRchop,mm_local,kk_local))).eq.0) write(*,*) 'zero:',m(mp(j)),ak(mp(j),kp(j))
   ENDIF

   ! m = 0; k > 0
   IF ((mp(j).eq.1).and.(kp(j).gt.1)) THEN 
      kk_local = (2*NXCHOP+1-kp(j)) - psi%inx
      IF (0.lt.kk_local .and. size(psi%e,3).ge.kk_local) THEN
         WRITE(*,*) m(mp(j)),ak(mp(j),kk_local+psi%inx),'RANK#',MPI_RANK
         psi%e(:NRchop,mm_local,kk_local) = psi%e(:NRchop,mm_local,kk_local) + &
               & conjg(perturb_psi(:NRchop,j))
         chi%e(:NRchop,mm_local,kk_local) = chi%e(:NRchop,mm_local,kk_local) + &
               & conjg(perturb_chi(:NRchop,j))
      ENDIF
   ENDIF 
ENDIF
end do
where(abs(psi%e)<1.0e-24) psi%e=0.0_p8
where(abs(chi%e)<1.0e-24) chi%e=0.0_p8

! ADDING SCALARS TO EXISTING FIELD
call allocate(psi0)
call allocate(chi0)
psi0%e = 0
chi0%e = 0
psi0%LN = 0
chi0%LN = 0
IF (MPI_RANK .EQ. 0) THEN
   WRITE(*,*) 'Q-VORTEX?[T/F]:'
   READ(*,*) save_q
ENDIF
CALL MPI_BCAST(save_q,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
IF (save_q) THEN
   call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//files%psi0,psi0)
   call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//files%chi0,chi0)
ELSE
   IF (MPI_RANK .EQ. 0) THEN
      WRITE(*,*) 'CORRECTION?[T/F]:'
      READ(*,*) save_q
   ENDIF
   CALL MPI_BCAST(save_q,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
   IF (save_q) THEN
      call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//"new_perturb_psi.dat",psi0)
      call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//"new_perturb_chi.dat",chi0)
   ENDIF
ENDIF
psi0%e = psi%e + psi0%e ! add to psi0, so that psi0%ln is preserved
chi0%e = chi%e + chi0%e 
call CHOPDO(psi0)
call CHOPDO(chi0)

! SAVING SCALARS
call msave(psi0, TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FILENAME))//"_psi.dat")
call msave(chi0, TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FILENAME))//"_chi.dat")

! FINALIZATION
call deallocate(psi)
call deallocate(chi)
call deallocate(psi0)
call deallocate(chi0)
deallocate(perturb_psi, perturb_chi, scale, mp, kp)

! MPI-FINALIZATION
IF (MPI_RANK.eq.0) THEN
WRITE(*,*) 'PROGRAM ENDED'
CALL PRINT_REAL_TIME()
ENDIF
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_FINALIZE(IERR)

end program addperturb_phase



