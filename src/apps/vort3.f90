program vort3
! ======================================================================
! Initial Condition:
! 1. Vq + lambda*V_bar + alpha*u1 + beta*u2 + lambda*correction
! 2. Vq + lambda*V_bar + lambda*correction
! Initialization:
! 1. addperturb: new_perturb w/ Q-vortex TRUE
! 2. addperturb: cor_perturb w/ Q-vortex FALSE
! 3. addperturb: lft_perturb w/ Q-vortex FALSE
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
type(scalar):: psi,chi,psi_l,chi_l,psi_c,chi_c
type(scalar):: dpsi,dchi,psi_o,chi_o,dpsi_c,dchi_c
integer:: it,pp,ii,jj,kk
REAL(P8):: time_start,time_end
real(p8),dimension(1,1):: status
real(p8), allocatable, dimension(:,:):: tempMonitor_mk
real(p8):: time0, tmp, vel_o
character(LEN=200) :: FileName_R, FileName_L, FileName_C, num
complex(p8),dimension(:,:),allocatable:: prod_result
complex(p8):: eig_1_psi,eig_1_chi,eig_2_psi,eig_2_chi
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
    WRITE(*,*) 'Uz: ', NADD%UZ
endif
CALL MPI_BCAST(FileName_R,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
10 FORMAT(A200)

if (MPI_RANK.eq.0) then
    WRITE(*,*) 'psi/chi(C) filename: '
    READ(*,10) FileName_C
endif
CALL MPI_BCAST(FileName_C,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

! degen pair
call allocate(psi)
call allocate(chi)
! correction terms
call allocate(psi_c)
call allocate(chi_c)
! left eigenvector
call allocate(psi_l)
call allocate(chi_l)
! before nonlin
call allocate(psi_o)
call allocate(chi_o)
! nonlin result
call allocate(dpsi)
call allocate(dchi)
! nonlin wrt correction
call allocate(dpsi_c)
call allocate(dchi_c)

!> initial values
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"_psi.dat",psi)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"_chi.dat",chi)

!> R_OLD: PSI,DEL2CHI
psi_o = psi
! CALL DEL2(chi,chi_o)
chi_o = chi

!> load correction
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_C))//"_psi.dat",psi_c)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_C))//"_chi.dat",chi_c)
psi%e = psi%e + psi_c%e
chi%e = chi%e + chi_c%e

time_start = mpi_wtime()

call PRINT_ENERGY_SPECTRUM(psi,chi,1)

!> startup
time0 = tim%t
call diagnost(psi,chi)
CALL WRITE_ENERGY_DATA_TO_FILE(psi,chi)
call WRITE_MODE_DATA_TO_FILE(psi,chi,FFF_SPACE)

CALL NONLIN(PSI,CHI,dpsi,dchi)

! calculate N(U_bar)U' + M(Vq)U'
DO II = 1,SIZE(PSI%E,2)
    DO JJ = 1,SIZE(PSI%E,3)
        IF ((II+PSI%INTH .EQ. MONITOR_MK(1,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(1,2)+1)) THEN
            ! ADD bar term to psi_c&chi_c
            ! MONITOR_MK(1,1) MUST MATCH THE PERTURBATION(BAR) MODE
            psi_c%E(:,II,JJ) = psi%E(:,II,JJ) + psi_c%E(:,II,JJ)
            chi_c%E(:,II,JJ) = chi%E(:,II,JJ) + chi_c%E(:,II,JJ)
        ENDIF
    ENDDO
ENDDO
CALL NONLIN(psi_c,chi_c,dpsi_c,dchi_c)
PSI%E = PSI%E + TIM%DT*(dpsi%e-dpsi_c%e)
CHI%E = CHI%E + TIM%DT*(dchi%e-dchi_c%e)
TIM%T = TIM%T + TIM%DT

CALL WRITE_ENERGY_DATA_TO_FILE(PSI,CHI)
call WRITE_MODE_DATA_TO_FILE(PSI,CHI,FFF_SPACE)

! TAKEOUT N(V)U'
dpsi%e = dpsi%e - dpsi_c%e
dchi%e = dchi%e - dchi_c%e

! V*R_NEW: PSI,DEL2CHI
PSI = dpsi
CHI%E = 0
CHI%LN = 0
CALL DEL2(dchi,CHI)
call WRITE_MODE_DATA_TO_FILE(PSI,CHI,FFF_SPACE)

! LOAD LEFT EIGENVECTOR
if (MPI_RANK.eq.0) then
    WRITE(*,*) 'psi/chi(L) filename: '
    READ(*,10) FileName_L
endif
CALL MPI_BCAST(FileName_L,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_L))//"_psi.dat",psi_l)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_L))//"_chi.dat",chi_l)

! PERFORM DOT PRODUCT
! SUM(psi_l*PSI+chi_l*CHI)
allocate(prod_result(SIZE(PSI%E,2),SIZE(PSI%E,3)))
DO II = 1,SIZE(PSI%E,2)
    DO JJ = 1,SIZE(PSI%E,3)
        ! prod_result(II,JJ) = sum(psi_l%E(:,II,JJ)*conjg(PSI%E(:,II,JJ)))+sum(chi_l%E(:,II,JJ)*conjg(chi%E(:,II,JJ)))
        prod_result(II,JJ) = DOT_PRODUCT(psi_l%E(:,II,JJ),PSI%E(:,II,JJ))+DOT_PRODUCT(chi_l%E(:,II,JJ),CHI%E(:,II,JJ))
        ! IF (abs(prod_result(II,JJ)).NE.0) WRITE(*,*) "M:",II+PSI%INTH,",K:",JJ+PSI%INX,prod_result(II,JJ),",AK:",AK(II+PSI%INTH,JJ+PSI%INX)
        IF ((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) THEN
            WRITE(*,*) "M1K1 sum = ",prod_result(II,JJ)
            ! WRITE(*,*) psi_l%e(:,ii,jj)
            ! WRITE(*,*) chi_l%e(:,ii,jj)
            ! dpsi%e(:,II,JJ) = dpsi%e(:,II,JJ) - (-1.7579848914879*IU)*psi_c%e(:,II,JJ)
            ! dchi%e(:,II,JJ) = dchi%e(:,II,JJ) - (-1.7579848914879*IU)*chi_c%e(:,II,JJ)
            dpsi%e(:,II,JJ) = dpsi%e(:,II,JJ) - (-1.33544798*IU)*psi_c%e(:,II,JJ)
            dchi%e(:,II,JJ) = dchi%e(:,II,JJ) - (-1.33544798*IU)*chi_c%e(:,II,JJ)
        ENDIF
        IF ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1)) THEN
            WRITE(*,*) "M2K2 sum = ",prod_result(II,JJ)
            ! WRITE(*,*) psi_l%e(:,ii,jj)
            ! WRITE(*,*) chi_l%e(:,ii,jj)
            ! dpsi%e(:,II,JJ) = dpsi%e(:,II,JJ) - (-3.9598393663464*IU)*psi_c%e(:,II,JJ)
            ! dchi%e(:,II,JJ) = dchi%e(:,II,JJ) - (-3.9598393663464*IU)*chi_c%e(:,II,JJ)
            dpsi%e(:,II,JJ) = dpsi%e(:,II,JJ) - (-3.251027209*IU)*psi_c%e(:,II,JJ)
            dchi%e(:,II,JJ) = dchi%e(:,II,JJ) - (-3.251027209*IU)*chi_c%e(:,II,JJ)
        ENDIF
    ENDDO
ENDDO

! CALCULATE EIGENVALUE
! ORIGINAL: [psi_o,chi_o] <- psi,del2chi
! AFTER   : [PSI  ,CHI  ] <- psi,del2chi
CALL RTRAN(psi_o,1)
CALL RTRAN(chi_o,1)
CALL RTRAN(dpsi,1)
CALL RTRAN(dchi,1)
allocate(eig_1_psi_r(size(psi_o%e,1)),eig_2_psi_r(size(psi_o%e,1)))
allocate(eig_1_chi_r(size(psi_o%e,1)),eig_2_chi_r(size(psi_o%e,1)))
DO II = 1,SIZE(PSI%E,2)
    DO JJ = 1,SIZE(PSI%E,3)
        IF ((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) THEN
            eig_1_psi_r = 0.D0
            eig_1_chi_r = 0.D0
            DO KK = 1,NR
                eig_1_psi_r(KK) = dpsi%E(KK,II,JJ)/psi_o%E(KK,II,JJ)
                eig_1_chi_r(KK) = dchi%E(KK,II,JJ)/chi_o%E(KK,II,JJ)
            ENDDO
            eig_1_psi = sum(eig_1_psi_r)/NR
            eig_1_chi = sum(eig_1_chi_r)/NR
            psi_o%E(:,II,JJ) = psi_o%E(:,II,JJ)*prod_result(II,JJ)
            chi_o%E(:,II,JJ) = chi_o%E(:,II,JJ)*prod_result(II,JJ)
            CALL save_compare(dpsi%E(:,II,JJ),dchi%E(:,II,JJ),psi_o%E(:,II,JJ),chi_o%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)

            WRITE(*,*) "M1K1 eig (psi) = ",eig_1_psi
            ! WRITE(*,*) "M1K1 eig (psi, TOP 10) = ",sum(eig_1_psi_r(1:10))/10
            WRITE(*,*) "      std = ", (sum((abs(eig_1_psi_r(1:NR)-eig_1_psi))**2)/NR)**0.5
            WRITE(*,*) "max error = ", MAXVAL(ABS(dpsi%E(:,II,JJ)-psi_o%E(:,II,JJ))/ABS(prod_result(II,JJ)))

            WRITE(*,*) "M1K1 eig (chi) = ",eig_1_chi
            ! WRITE(*,*) "M1K1 eig (chi, TOP 10) = ",sum(eig_1_chi_r(1:10))/10
            WRITE(*,*) "     std = ", (sum((abs(eig_1_chi_r(1:NR)-eig_1_chi))**2)/NR)**0.5
            WRITE(*,*) "max error = ", MAXVAL(ABS(dchi%E(:,II,JJ)-chi_o%E(:,II,JJ))/ABS(prod_result(II,JJ)))
        ENDIF
        IF ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1)) THEN
            eig_2_psi_r = 0.D0
            eig_2_chi_r = 0.D0
            DO KK = 1,NR
                eig_2_psi_r(KK) = dpsi%E(KK,II,JJ)/psi_o%E(KK,II,JJ)
                eig_2_chi_r(KK) = dchi%E(KK,II,JJ)/chi_o%E(KK,II,JJ)
            ENDDO
            eig_2_psi = sum(eig_2_psi_r)/NR
            eig_2_chi = sum(eig_2_chi_r)/NR
            psi_o%E(:,II,JJ) = psi_o%E(:,II,JJ)*prod_result(II,JJ)
            chi_o%E(:,II,JJ) = chi_o%E(:,II,JJ)*prod_result(II,JJ)
            CALL save_compare(dpsi%E(:,II,JJ),dchi%E(:,II,JJ),psi_o%E(:,II,JJ),chi_o%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)

            WRITE(*,*) "M2K2 eig (psi) = ",eig_2_psi
            ! WRITE(*,*) "M1K1 eig (psi, TOP 10) = ",sum(eig_1_psi_r(1:10))/10
            WRITE(*,*) "      std = ", (sum((abs(eig_2_psi_r(1:NR)-eig_2_psi))**2)/NR)**0.5
            WRITE(*,*) "max error = ", MAXVAL(ABS(dpsi%E(1:NR,II,JJ)-psi_o%E(1:NR,II,JJ))/ABS(prod_result(II,JJ)))

            WRITE(*,*) "M2K2 eig (chi) = ",eig_2_chi
            ! WRITE(*,*) "M1K1 eig (chi, TOP 10) = ",sum(eig_1_chi_r(1:10))/10
            WRITE(*,*) "     std = ", (sum((abs(eig_2_chi_r(1:NR)-eig_2_chi))**2)/NR)**0.5
            WRITE(*,*) "max error = ", MAXVAL(ABS(dchi%E(1:NR,II,JJ)-chi_o%E(1:NR,II,JJ))/ABS(prod_result(II,JJ)))
        ENDIF
    ENDDO
ENDDO
! DO II = 1,SIZE(PSI%E,2)
!     DO JJ = 1,SIZE(PSI%E,3)
!         IF ((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) THEN
!             eig_1 = 0.D0
!             DO KK = 1,NRCHOPS(II+PSI%INTH)
!                 eig_1(KK) = PSI%E(KK,II,JJ)/psi_o%E(KK,II,JJ)
!             ENDDO
!             WRITE(*,*) "M1K1 eig = ",sum(eig_1)/NRCHOPS(II+PSI%INTH)
!             ! (psi_o,chi_o)*sigma
!             psi_o%E(:,II,JJ) = psi_o%E(:,II,JJ)*prod_result(II,JJ)
!             chi_o%E(:,II,JJ) = chi_o%E(:,II,JJ)*prod_result(II,JJ)
!             CALL save_compare(PSI%E(:,II,JJ),CHI%E(:,II,JJ),psi_o%E(:,II,JJ),chi_o%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
!             WRITE(*,*) "M1K1 eig = ",sum(eig_1(1:10))/10
!             WRITE(*,*) "ERROR: ",SUM(ABS(PSI%E(:,II,JJ)-psi_o%E(:,II,JJ))+ABS(CHI%E(:,II,JJ)-chi_o%E(:,II,JJ)))/size(PSI%E(:,II,JJ))
!         ENDIF
!         IF ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1)) THEN
!             eig_2 = 0.D0
!             DO KK = 1,NRCHOPS(II+PSI%INTH)
!                 eig_2(KK) = PSI%E(KK,II,JJ)/psi_o%E(KK,II,JJ)
!             ENDDO
!             WRITE(*,*) "M2K2 eig = ",sum(eig_2)/NRCHOPS(II+PSI%INTH)
!             ! (psi_o,chi_o)*sigma
!             psi_o%E(:,II,JJ) = psi_o%E(:,II,JJ)*prod_result(II,JJ)
!             chi_o%E(:,II,JJ) = chi_o%E(:,II,JJ)*prod_result(II,JJ)
!             CALL save_compare(PSI%E(:,II,JJ),CHI%E(:,II,JJ),psi_o%E(:,II,JJ),chi_o%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
!             WRITE(*,*) "M2K2 eig = ",sum(eig_2(1:10))/10
!             WRITE(*,*) "ERROR: ",SUM(ABS(PSI%E(:,II,JJ)-psi_o%E(:,II,JJ))+ABS(CHI%E(:,II,JJ)-chi_o%E(:,II,JJ)))/size(PSI%E(:,II,JJ))
!         ENDIF
!     ENDDO
! ENDDO


CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

DEALLOCATE(prod_result)
CALL DEALLOCATE(psi)
CALL DEALLOCATE(chi)
CALL DEALLOCATE(psi_l)
CALL DEALLOCATE(chi_l)
CALL DEALLOCATE(dpsi)
CALL DEALLOCATE(dchi)
CALL DEALLOCATE(psi_o)
CALL DEALLOCATE(chi_o)
call DEALLOCATE(psi_c)
call DEALLOCATE(chi_c)
call DEALLOCATE(dpsi_c)
call DEALLOCATE(dchi_c)

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


end program vort3
