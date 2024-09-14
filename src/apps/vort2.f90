program vort2
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
type(scalar):: psi,chi,psio,chio,dpsi,dchi,psi_o,chi_o
type(scalar):: ruro, rupo, ruzo
integer:: iii,it,pp
REAL(P8):: time_start,time_end
real(p8),dimension(1,1):: status
real(p8), allocatable, dimension(:,:):: tempMonitor_mk
real(p8):: time0, tmp, vel_o
character(len=72) :: psiFileName, chiFileName, num

!    ! debug:
!   COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: GLB_ARRAY


! CALL MPI_INIT(IERR)
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
call allocate(psio)
call allocate(chio)
call allocate(psi_o)
call allocate(chi_o)
call allocate(dpsi)
call allocate(dchi)

!> initial values
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%psii,psi)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%chii,chi)

time_start = mpi_wtime()
call PRINT_ENERGY_SPECTRUM(psi,chi,1)

!> startup : very first time step
call diagnost(psi,chi)
psi_o = psi
chi_o = chi
time0 = tim%t
vel_o = MAXVELP(psi_o,chi_o)

!dpsi and dchi are initially empty, then they are assigned 
!the nonlinear part of the first step
! call rich(psi,chi,dpsi,dchi)
call RICH_RE(psi,chi,dpsi,dchi,psio,chio,psi_o,chi_o,vel_o)
call diagnost(psi,chi)

iii = tim%limit/tim%dt
files%n = 1

do it=1,iii
    call WRITE_MODE_DATA_TO_FILE(psi,chi,PFF_SPACE) 

    if(monitorData%start.eq.tim%n) then 
        ! use this to calculate eigenvalues. this is the psi/chi at t1
        psi_o = psi
        chi_o = chi
        time0 = tim%t
    end if
    
    !> integrate
    call adamsb(psi,chi,dpsi,dchi)
    call diagnost(psi,chi)
    
    !> monitor eigenmodes
    call WRITE_ENERGY_DATA_TO_FILE(psi,chi)
!     call WRITE_MODE_DATA_TO_FILE(psi,chi)
    tmp = (real(tim%n)-real(monitorData%start))/real(monitorData%interval)
    if (tmp.ge.1.0_p8 .and. abs(tmp-int(tmp)).lt.1e-08) then
        call egrowth_modified(psi_o,chi_o,psi,chi,time0,tim%t)
        call monitor_eig(psi_o,chi_o,psi,chi,time0,tim%t)
    !   call PRINT_ENERGY_SPECTRUM(psi,chi,1)
        call PRINT_ENERGY_SPECTRUM(psi,chi,0)
        call WRITE_MODE_DATA_TO_FILE(psi,chi,PFF_SPACE)   
    end if
    
    ! !> freestream adjust
    ! if(adv%sw.eq.1 .and. mod(it,adv%int).eq.0) then
    !     call freeadj
    !     call rich(psi,chi,dpsi,dchi)
    !     call diagnost(psi,chi)
    ! endif
    
    !> vorticity removal
    if(rmv%sw.ne.0 .and. mod(it,rmv%int).eq.0) then
        call remove(psi,chi)
        call rich(psi,chi,dpsi,dchi)
        call diagnost(psi,chi)
    endif
    
    !> hyperviscosity adjust
    if(visc%adjsw.ne.0 .and. mod(it,visc%adjint).eq.0) then
        call hypadj(psi,chi,visc%nup)
        call diagnost(psi,chi)
    endif
    
    !> output
    if(files%t(files%n).le.tim%t) then
        call msave(psi, TRIM(ADJUSTL(FILES%SAVEDIR))//files%psi(files%n))
        call msave(chi, TRIM(ADJUSTL(FILES%SAVEDIR))//files%chi(files%n))
        files%n = files%n + 1
        if(files%n > files%ne) goto 999
        
        ! 'status.dat' is currently un-used
        ! Setting status.dat's stored value to 1 stops the program after 
        ! storing the results at the previous savepoint
        if (MPI_RANK.eq.0) then 
        call mload('status.dat',status)
        if(status(1,1).eq.1.0) goto 999
        endif

    endif

    !> clean modes
    ! if (tmp.ge.1.0_p8 .and. abs(tmp-int(tmp)).lt.1e-08) then
        call RENORMALIZE(psi,chi,psio,chio,psi_o,chi_o,vel_o)
    ! endif

enddo

999 continue
time_end = mpi_wtime()
!  ALLOCATE(GLB_ARRAY(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM))
!  call MASSEMBLE(PSI%E,GLB_ARRAY,1)
if (MPI_RANK.eq.0) THEN
print *,tim%n,' steps'
WRITE(*,*) 'PROGRAM COMPLETED'
CALL PRINT_REAL_TIME()
WRITE(*,*) 'EXECUTION TIME: ',time_end-time_start,'seconds'
! CALL MCAT(GLB_ARRAY(:,1,1))
ENDIF
call MPI_BARRIER(MPI_COMM_IVP,IERR)
!stop
!CALL MPI_ABORT(MPI_COMM_IVP,1,IERR)
call MPI_FINALIZE(IERR)

! ======================================================================
CONTAINS

! SUBROUTINE RENORMALIZE(PSI0,CHI0,PSI1,CHI1,PSIi,CHIi,ENE_i)
! ! ======================================================================
! ! NOTE:
! ! MONITOR_MK MUST BE CONSISTENT WITH THE PERTURBATION ADDED.
! ! MONITOR_MK{1,:} -> M_BAR, K_BAR
! ! MONITOR_MK{2,:} -> M_1, K_1
! ! MONITOR_MK{3,:} -> M_2, K_2
! ! ======================================================================
!     IMPLICIT NONE
!     TYPE(SCALAR):: PSI0,CHI0,PSI1,CHI1,PSIi,CHIi
!     INTEGER:: MON_SIZE, MON_II, MM, KK
!     REAL(P8), DIMENSION(NTCHOP,NXCHOPDIM):: ENE0, ENE_i

!     ! WIPE PSI1&CHI1
!     PSI1%ln = 0.D0 
!     CHI1%ln = 0.D0
!     PSI1%e = 0.D0
!     CHI1%e = 0.D0

!     ! OBTAIN ENERGY
!     ENE0 = ENERGY_SPEC_MODIFIED(PSI0,CHI0)

!     ! RENORMALIZE
!     MON_SIZE = SIZE(MONITOR_MK,1)
!     DO MON_II = 1,MIN(MON_SIZE,3)
!         ! NOTE: THE LINE BELOW IS KINDA INCORRECT. FOR {-M,-K}, {M,K} SHOULD BE MONITORED. HOWEVER THE CODE CHANGES 
!         !       IT TO {M,-K}. ALSO, THE LINE CHANGES MONITOR_MK, WHICH IS A GLOBAL VARIABLE -> THIS IS BAD.
!         !       CURRENTLY, ALL 3 ADDED MODES ARE POSITIVE, SO THIS IS FINE FOR NOW.
!         IF (MONITOR_MK(MON_II,1)<0) MONITOR_MK(MON_II,:) = -MONITOR_MK(MON_II,:)
!         MM = MONITOR_MK(MON_II,1)+1 ! INDEX
!         KK = MONITOR_MK(MON_II,2)+1 ! INDEX
!         IF (KK<1) KK=KK+2*NXCHOP-1

!         IF ((PSI1%INTH.LT.MM).AND.(PSI1%INTH+SIZE(PSI1%E,2).GE.MM)) THEN
!             IF ((PSI1%INX.LT.KK).AND.(PSI1%INX+SIZE(PSI1%E,3).GE.KK)) THEN
!                 IF (MON_II.GT.1) THEN ! RENORMALIZE THE TESTED EIGENVECTOR: alpha*V_{m1,k1}+beta*V_{m2,k2}
!                     PSI1%E(:,MM-PSI1%INTH,KK-PSI1%INX) = PSI0%E(:,MM-PSI0%INTH,KK-PSI0%INX)*((ENE_i(MM,KK))**0.5) &
!                                                         /((ENE0(MM,KK))**0.5)
!                     CHI1%E(:,MM-CHI1%INTH,KK-CHI1%INX) = CHI0%E(:,MM-CHI0%INTH,KK-CHI0%INX)*((ENE_i(MM,KK))**0.5) &
!                                                         /((ENE0(MM,KK))**0.5)
!                 ELSE ! REFIX THE BASE FLOW: V_{m_bar,k_bar}
!                     PSI1%E(:,MM-PSI1%INTH,KK-PSI1%INX) = PSIi%E(:,MM-PSI0%INTH,KK-PSI0%INX)
!                     CHI1%E(:,MM-CHI1%INTH,KK-CHI1%INX) = CHIi%E(:,MM-CHI0%INTH,KK-CHI0%INX)
!                 ENDIF
!             ENDIF
!         ENDIF
!     ENDDO
    
!     IF ((PSI1%INX.EQ.0).AND.(PSI1%INTH.EQ.0)) THEN ! REFIX THE BASE FLOW: V_{0,0}
!         PSI1%E(:,1,1) = PSIi%e(:,1,1) 
!         CHI1%E(:,1,1) = CHIi%e(:,1,1)
!     ENDIF
    
!     PSI0%E = PSI1%E
!     CHI0%E = CHI1%E
!     PSI0%LN = PSIi%LN
!     CHI0%LN = CHIi%LN

!     IF (MPI_RANK.EQ.0) WRITE(*,*) 'RENORMALIZED'
!     call WRITE_ENERGY_DATA_TO_FILE(PSI0,CHI0)


! END SUBROUTINE RENORMALIZE

end program vort2
