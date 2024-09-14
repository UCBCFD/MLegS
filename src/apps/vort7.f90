program vort7
! ======================================================================
! Initial Condition:
! Vq + lambda*V_bar + alpha*u1 + beta*u2 + lambda*correction
! EOM: Euler's Eqn (no viscosity)
! ======================================================================
! -------------------------
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
! -------------------------
implicit none
! -------------------------
type(scalar):: psi,chi,psi_c,chi_c,psi2,chi2
type(scalar):: dpsi,dchi,psi_o,chi_o
integer:: iii,it,mm,kk,pp,ii,jj
REAL(P8):: time_start,time_end
real(p8),dimension(1,1):: status
real(p8):: time0, tmp

CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
   WRITE(*,*) 'The threading support is lesser than that demanded.'
   CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
ENDIF
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)

CALL READCOM('NOECHO')
CALL READIN(5)
CALL LEGINIT()

status(1,1)=0
if (MPI_RANK.eq.0) then
WRITE(*,*) 'PROGRAM STARTED'
call msave(status, 'status.dat')
CALL PRINT_REAL_TIME()
endif

call allocate(psi)
call allocate(chi)
call allocate(psi2)
call allocate(chi2)
call allocate(psi_c)
call allocate(chi_c)
call allocate(psi_o)
call allocate(chi_o)
call allocate(dpsi)
call allocate(dchi)

!> initial values
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//"new_perturb_psi.dat",psi)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//"new_perturb_chi.dat",chi)
!> load correction terms
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//"cor_perturb_psi.dat",psi_c)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//"cor_perturb_chi.dat",chi_c)
!> combine
psi%e = psi%e + psi_c%e
chi%e = chi%e + chi_c%e

time_start = mpi_wtime()

!> startup : very first time step
!dpsi and dchi are initially empty, then they are assigned 
!the nonlinear part of the first step
call diagnost(psi,chi)
call rich(psi,chi,dpsi,dchi)
call diagnost(psi,chi)

psi_o = psi
chi_o = chi
time0 = tim%t

iii = tim%limit/tim%dt
files%n = 1

do it=1,iii

   if(monitorData%start.eq.tim%n) then 
      ! use this to calculate eigenvalues. this is the psi/chi at t1
      psi_o = psi
      chi_o = chi
      time0 = tim%t
   end if
   
   !> integrate - uxw
   call adamsb(psi,chi,dpsi,dchi)
   call diagnost(psi,chi)

   !> re-adjust psi/chi copy
   psi2%e = 0.d0
   chi2%e = 0.d0
   psi2 = psi_o
   chi2 = chi_o
   
   !> monitor eigenmodes
   tmp = (real(tim%n)-real(monitorData%start))/real(monitorData%interval)
   if (tmp.ge.1.0_p8 .and. abs(tmp-int(tmp)).lt.1e-08) then
      call egrowth_modified(psi_o,chi_o,psi,chi,time0,tim%t)
      call monitor_eig(psi_o,chi_o,psi,chi,time0,tim%t)
      DO II = 1,SIZE(PSI%E,2)
         DO JJ = 1,SIZE(PSI%E,3)
               IF ((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) THEN
                  CALL save_mode(tim%t,psi%E(:,II,JJ),chi%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
                  psi2%E(:,II,JJ) = psi%E(:,II,JJ)
                  chi2%E(:,II,JJ) = chi%E(:,II,JJ)
               ENDIF
               IF ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1)) THEN
                  CALL save_mode(tim%t,psi%E(:,II,JJ),chi%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
                  psi2%E(:,II,JJ) = psi%E(:,II,JJ)
                  chi2%E(:,II,JJ) = chi%E(:,II,JJ)
               ENDIF
         ENDDO
      ENDDO
   end if

    !> freestream adjust
    if(adv%sw.eq.1 .and. mod(it,adv%int).eq.0) then
        call freeadj
        call rich(psi,chi,dpsi,dchi)
        call diagnost(psi,chi)
    endif
         
   !> output
   if(files%t(files%n).le.tim%t) then
      call msave(psi, TRIM(ADJUSTL(FILES%SAVEDIR))//files%psi(files%n))
      call msave(chi, TRIM(ADJUSTL(FILES%SAVEDIR))//files%chi(files%n))
      files%n = files%n + 1
      if(files%n > files%ne) goto 999 ! abort when the last file is saved
      
      ! 'status.dat' is currently un-used
      ! Setting status.dat's stored value to 1 stops the program after 
      ! storing the results at the previous savepoint
      if (MPI_RANK.eq.0) then 
        call mload('status.dat',status)
        if(status(1,1).eq.1.0) goto 999
      endif
   endif

   !> clean modes
   ! In the moving frame, {0,0} and perturb mode shall not change
   ! Remove all the other modes except the two degen modes
   psi = psi2
   chi = chi2

enddo

999 continue
time_end = mpi_wtime()
if (MPI_RANK.eq.0) THEN
print *,tim%n,' steps'
WRITE(*,*) 'PROGRAM STARTED'
CALL PRINT_REAL_TIME()
WRITE(*,*) 'EXECUTION TIME: ',time_end-time_start,'seconds'
ENDIF
call MPI_BARRIER(MPI_COMM_IVP,IERR)
call MPI_FINALIZE(IERR)

! ======================================================================
contains
subroutine save_mode(t,psi_new,chi_new,m_i,k_i)
complex(p8),DIMENSION(:):: psi_new,chi_new
real:: t
integer:: m_i,k_i
integer:: nn

open(UNIT=777,FILE='mode_track_'//ITOA3(m_i)//'_'//ITOA3(k_i)//'.dat',&
&STATUS='UNKNOWN',ACTION='WRITE',ACCESS='APPEND')

WRITE(777,320) t,psi_new(1:size(psi_new)),chi_new(1:size(psi_new))

close(777)
320 FORMAT(F6.3(S,E14.6E3,SP,E14.6E3,'i'),*(','S,E14.6E3,SP,E14.6E3,'i'))

end subroutine save_mode
end program vort7
