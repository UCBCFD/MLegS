program vort6
! ======================================================================
! Initial Condition:
! Vq + lambda*V_bar + alpha*u1 + beta*u2 + lambda*correction
! Initialization:
! 1. addperturb: new_perturb w/ Q-vortex TRUE
! 2. addperturb: cor_perturb w/ Q-vortex FALSE
! This program calculates one time step for the Euler's Eqn. And save
! the original psi,chi as well as the dpsi,dchi to file 'mode_compare_x'
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
! -------------------------
implicit none
! -------------------------
type(scalar):: psi,chi,psi_c,chi_c
type(scalar):: dpsi,dchi,psi_o,chi_o
integer:: it,pp,ii,jj,kk
REAL(P8):: time_start,time_end
real(p8),dimension(1,1):: status
real(p8), allocatable, dimension(:,:):: tempMonitor_mk
real(p8):: time0, tmp, vel_o
complex(p8):: eig_1_psi,eig_1_chi,eig_2_psi,eig_2_chi
character(LEN=200) :: FileName_R, FileName_L, FileName_C, num
complex(p8),dimension(:,:),allocatable:: prod_result
complex(p8),dimension(:),allocatable:: eig_1_psi_r,eig_2_psi_r
complex(p8),dimension(:),allocatable:: eig_1_chi_r,eig_2_chi_r


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
    WRITE(*,*) 'psi/chi(R) filename: '
    READ(*,10) FileName_R
    WRITE(*,*) 'psi/chi(C) filename: '
    READ(*,10) FileName_C
    ! WRITE(*,*) 'Uz: ', NADD%UZ
endif
CALL MPI_BCAST(FileName_R,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(FileName_C,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
10 FORMAT(A200)

! degen triad + Q-vortex
call allocate(psi)
call allocate(chi)
! correction term
call allocate(psi_c)
call allocate(chi_c)
! before nonlin
call allocate(psi_o)
call allocate(chi_o)
! nonlin result
call allocate(dpsi)
call allocate(dchi)

!> load Q-vortex + resonant triad
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"_psi.dat",psi)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"_chi.dat",chi)
!> load correction
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_C))//"_psi.dat",psi_c)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_C))//"_chi.dat",chi_c)
!> add correction
psi%e = psi%e + psi_c%e
chi%e = chi%e + chi_c%e

time_start = mpi_wtime()

! call PRINT_ENERGY_SPECTRUM(psi,chi,1)

!> R_OLD: psi,chi
psi_o = psi
chi_o = chi

!> startup
time0 = tim%t
call diagnost(psi,chi)
CALL WRITE_ENERGY_DATA_TO_FILE(psi,chi)
call WRITE_MODE_DATA_TO_FILE(psi,chi,FFF_SPACE)

! One Euler's step
CALL NONLIN(psi,chi,dpsi,dchi)
write(*,*) 'nonlin done'
psi%E = psi%E + TIM%DT*(dpsi%e)
chi%E = chi%E + TIM%DT*(dchi%e)
TIM%T = TIM%T + TIM%DT

CALL WRITE_ENERGY_DATA_TO_FILE(psi_o,chi_o)
call WRITE_MODE_DATA_TO_FILE(psi_o,chi_o,FFF_SPACE)

! CALCULATE EIGENVALUE
! ORIGINAL: [psi_o,chi_o] <- psi,del2chi
! AFTER   : [psi  ,chi  ] <- psi,del2chi
CALL RTRAN(psi_o,1)
CALL RTRAN(chi_o,1)
CALL RTRAN(dpsi,1)
CALL RTRAN(dchi,1)
allocate(eig_1_psi_r(size(psi_o%e,1)),eig_2_psi_r(size(psi_o%e,1)))
allocate(eig_1_chi_r(size(psi_o%e,1)),eig_2_chi_r(size(psi_o%e,1)))
DO II = 1,SIZE(psi%E,2)
    DO JJ = 1,SIZE(psi%E,3)
        IF ((II+psi%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+psi%INX .EQ. MONITOR_MK(2,2)+1)) THEN
            eig_1_psi_r = 0.D0
            eig_1_chi_r = 0.D0
            DO KK = 1,NR
                eig_1_psi_r(KK) = dpsi%E(KK,II,JJ)/psi_o%E(KK,II,JJ)
                eig_1_chi_r(KK) = dchi%E(KK,II,JJ)/chi_o%E(KK,II,JJ)
            ENDDO
            eig_1_psi = sum(eig_1_psi_r)/NR
            eig_1_chi = sum(eig_1_chi_r)/NR
            WRITE(*,*) "M1K1 eig (psi) = ",eig_1_psi,", std = ", (sum((abs(eig_1_psi_r(1:NR)-eig_1_psi))**2)/NR)**0.5
            WRITE(*,*) "M1K1 eig (chi) = ",eig_1_chi,", std = ", (sum((abs(eig_1_chi_r(1:NR)-eig_1_chi))**2)/NR)**0.5

            CALL save_compare(dpsi%E(:,II,JJ),dchi%E(:,II,JJ),psi_o%E(:,II,JJ),chi_o%E(:,II,JJ),II+psi%INTH-1,JJ+psi%INX-1)
            WRITE(*,*) "M1K1 eig (psi, TOP 10) = ",sum(eig_1_psi_r(1:10))/10
            WRITE(*,*) "M1K1 eig (chi, TOP 10) = ",sum(eig_1_chi_r(1:10))/10
        ENDIF
        IF ((II+psi%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+psi%INX .EQ. MONITOR_MK(3,2)+1)) THEN
            eig_2_psi_r = 0.D0
            eig_2_chi_r = 0.D0
            DO KK = 1,NR
                eig_2_psi_r(KK) = dpsi%E(KK,II,JJ)/psi_o%E(KK,II,JJ)
                eig_2_chi_r(KK) = dchi%E(KK,II,JJ)/chi_o%E(KK,II,JJ)
            ENDDO
            eig_2_psi = sum(eig_2_psi_r)/NR
            eig_2_chi = sum(eig_2_chi_r)/NR
            WRITE(*,*) "M2K2 eig (psi) = ",eig_2_psi,", std = ", (sum((abs(eig_2_psi_r(1:NR)-eig_2_psi))**2)/NR)**0.5
            WRITE(*,*) "M2K2 eig (chi) = ",eig_2_chi,", std = ", (sum((abs(eig_2_chi_r(1:NR)-eig_2_chi))**2)/NR)**0.5

            CALL save_compare(dpsi%E(:,II,JJ),dchi%E(:,II,JJ),psi_o%E(:,II,JJ),chi_o%E(:,II,JJ),II+psi%INTH-1,JJ+psi%INX-1)
            WRITE(*,*) "M2K2 eig (psi, TOP 10)  = ",sum(eig_2_psi_r(1:20))/20
            WRITE(*,*) "M2K2 eig (chi, TOP 10)  = ",sum(eig_2_chi_r(1:20))/20
        ENDIF
    ENDDO
ENDDO

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

CALL DEALLOCATE(psi)
CALL DEALLOCATE(chi)
CALL DEALLOCATE(dpsi)
CALL DEALLOCATE(dchi)
CALL DEALLOCATE(psi_o)
CALL DEALLOCATE(chi_o)

time_end = mpi_wtime()
IF (MPI_RANK.eq.0) THEN
    WRITE(*,*) 'PROGRAM COMPLETED'
    CALL PRINT_REAL_TIME()
    WRITE(*,*) 'EXECUTION TIME: ',time_end-time_start,'seconds'
ENDIF
call MPI_BARRIER(MPI_COMM_IVP,IERR)
call MPI_FINALIZE(IERR)

! ======================================================================
contains
subroutine save_compare(psi_new,chi_new,psi_old,chi_old,m_i,k_i)
complex(p8),DIMENSION(:):: psi_new,chi_new,psi_old,chi_old
integer:: m_i,k_i
integer:: nn

open(UNIT=777,FILE='mode_compare_'//ITOA3(m_i)//'_'//ITOA3(k_i)//'.dat',&
&STATUS='UNKNOWN')

do nn = 1,size(psi_new)
    WRITE(777,320) psi_new(nn),psi_old(nn),chi_new(nn),chi_old(nn)
enddo

close(777)
320 FORMAT((S,E14.6E3,SP,E14.6E3,'i'),*(','S,E14.6E3,SP,E14.6E3,'i'))

end subroutine save_compare


end program vort6
