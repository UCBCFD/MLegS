program vort8
! ======================================================================
! RE-FINING EIGENVECTOR VIA LINEARIZED N-S
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
   USE MOD_EVP, only: SAVE_PERTURB
! -------------------------
implicit none
! -------------------------
integer:: iii,it,mm,kk,pp,II,JJ
real(p8),dimension(1,1):: status
real(P8):: time_start,time_end,time0,tmp
type(scalar):: psi,chi,psio,chio,dpsi,dchi
complex(P8),dimension(:),allocatable:: vec_1, vec_2, vec_bar
complex(P8),dimension(:,:,:),allocatable:: psi_glb, chi_glb

! REFINEMENT
character(LEN=200):: FileName_R
real(p8):: eig1, eig2, eig10, eig20
complex(p8):: anchor_psi_val, anchor_factor
real(p8), allocatable, dimension(:,:):: EMK, EMK0
integer:: converg_flag1, converg_flag2, anchor_psi_loc
integer:: eig1_psi_loc, eig1_chi_loc, eig2_psi_loc, eig2_chi_loc

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
    WRITE(*,*) 'UZ: ', NADD%UZ
    WRITE(*,*) 'psi/chi(R) filename: '
    READ(*,10) FileName_R
endif
CALL MPI_BCAST(FileName_R,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
10 FORMAT(A200)

! degen triad + Q-vortex
call allocate(psi)
call allocate(chi)
! initial value
call allocate(psio)
call allocate(chio)
! nonlin result
call allocate(dpsi)
call allocate(dchi)
! energy
! allocate(EMK(NTCHOP,NXCHOP))
! allocate(EMK0(NTCHOP,NXCHOP))
allocate(EMK(NTCHOP,NXCHOPDIM))
allocate(EMK0(NTCHOP,NXCHOPDIM))

!> load Q-vortex + resonant triad
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"_psi.dat",psi)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"_chi.dat",chi)

!> adjust monitor_mk if necessary
DO II = 1,SIZE(MONITOR_MK,1)
   IF (MONITOR_MK(II,1).LT.0) MONITOR_MK(II,:) = -MONITOR_MK(II,:)
   IF (MONITOR_MK(II,2).LT.0) MONITOR_MK(II,2) = 2*NXCHOP-1+MONITOR_MK(II,2)
ENDDO

!> debug
IF (MPI_RANK.EQ.0) THEN
   WRITE(*,*) 'TRACK MODES: (M,AK) - (#MM,#KK)'
   DO II = 1,SIZE(MONITOR_MK,1)
      WRITE(*,92) M(MONITOR_MK(II,1)+1),AK(MONITOR_MK(II,1)+1,MONITOR_MK(II,2)+1),MONITOR_MK(II,1)+1,MONITOR_MK(II,2)+1
   ENDDO
92 FORMAT('(',I2,',',F9.2,') - (#',I3,', #',I3,')')
ENDIF

!> save the resonant pairs in PFF space
call allocate(psio); psio%e = psi%e; psio%ln = psi%ln; CALL RTRAN(psio,1)
call allocate(chio); chio%e = chi%e; chio%ln = chi%ln; CALL RTRAN(chio,1)
DO II = 1,SIZE(PSI%E,2)
   DO JJ = 1,SIZE(PSI%E,3)
      IF (((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) &
         .OR. ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1))) THEN
         CALL save_mode(tim%t,psio%E(:,II,JJ),chio%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
      ENDIF
   ENDDO
ENDDO
CALL DEALLOCATE(psio); CALL DEALLOCATE(chio)

!> richardson step
call rich_lin(psi,chi,dpsi,dchi)

!> save initial energy after Richardson step
EMK0 = ENERGY_SPEC_MODIFIED(psi,chi)
! ! CALL MCAT(PSI%E(1:4,1,:))
! ! CALL MPI_BARRIER(MPI_COMM_IVP,IERR)
! IF (MPI_RANK.EQ.0) THEN
!    CALL MCAT(EMK0); WRITE(*,*)
!    ! CALL MPI_ABORT(MPI_COMM_IVP,1,IERR)
! ENDIF

!> startup
!dpsi and dchi are initially empty, then they are assigned 
!the nonlinear part of the first step
time_start = mpi_wtime()
iii = tim%limit/tim%dt
files%n = 1
converg_flag1 = 0
converg_flag2 = 0
anchor_factor = 1.d0
eig10 = 0.d0; eig20 = 0.d0

do it=1,iii

   !> admam-bashforth - degen renormalized
   !> psi(n+1),chi(n+1) = factor * |psi(n),chi(n)|_degen + [1.5*|psiN(n),chiN(n)|_degen - 0.5*|psiN(n-1),chiN(n-1)|_degen]
   call adamsb_lin(psi,chi,dpsi,dchi,anchor_factor)

   !> check convergence in PFF space 
   ! Use real growth rate as the indicator
   EMK = ENERGY_SPEC_MODIFIED(psi,chi)
   eig1 = LOG(EMK(MONITOR_MK(2,1)+1,MONITOR_MK(2,2)+1)/EMK0(MONITOR_MK(2,1)+1,MONITOR_MK(2,2)+1))/(2*TIM%DT)
   eig2 = LOG(EMK(MONITOR_MK(3,1)+1,MONITOR_MK(3,2)+1)/EMK0(MONITOR_MK(3,1)+1,MONITOR_MK(3,2)+1))/(2*TIM%DT)

   ! IF ((it.EQ.10).OR.(it.EQ.1)) THEN
   !    IF (MPI_RANK.EQ.0) THEN
   !       ! CALL MCAT(EMK0(1,1:8)); WRITE(*,*) ''
   !       ! CALL MCAT(EMK(1,1:8)); WRITE(*,*) ''
   !       CALL MCAT(EMK0); WRITE(*,*)
   !       CALL MCAT(EMK); WRITE(*,*)
   !       WRITE(*,*) eig1
   !       WRITE(*,*) eig2
   !       CALL MPI_ABORT(MPI_COMM_IVP,1,IERR)
   !    ENDIF
   ! ENDIF

   IF ((MPI_RANK.EQ.0).and.(MOD(it,100).EQ.0)) THEN
      WRITE(*,*) 'STEP: ',it,'/',iii
      WRITE(*,*) '{M1,K1}: ',eig1,' - ',EMK(MONITOR_MK(2,1)+1,MONITOR_MK(2,2)+1)
      WRITE(*,*) '{M2,K2}: ',eig2,' - ',EMK(MONITOR_MK(3,1)+1,MONITOR_MK(3,2)+1)

      open(UNIT=666,FILE='eig_eneg.dat',STATUS='UNKNOWN',ACTION='WRITE',ACCESS='APPEND')
      WRITE(666,134) tim%t,eig1,eig2
      close(666)
      134 FORMAT(F10.3,2(',',SP,E26.15E3))
   ENDIF
   if ((abs(eig1-eig10)/abs(eig10).lt.1e-8) .AND. (abs(eig2-eig20)/abs(eig20).lt.1e-8)) GOTO 999
   eig10 = eig1; eig20 = eig2
   
   !> renormalization factor - keep {m1,k1}'s energy
   anchor_factor = EXP(-eig1*TIM%DT)

   !> save renormalized energy
   !> NOTE: or I can use the saved EMK0 after richardson step
   EMK0 = EMK*REAL(EXP(-2*eig1*TIM%DT))

   ! !>> DEBUG:
   ! !> print energy
   ! IF (MPI_RANK.EQ.0) THEN
   !    CALL MCAT(EMK0(1,1:8))
   ! ENDIF
   ! !> abort
   ! IF (it.EQ.10) THEN      
   !    CALL MPI_BARRIER(MPI_COMM_IVP,IERR)
   !    CALL MPI_ABORT(MPI_COMM_IVP,1,IERR)
   ! ENDIF

   !> save log
   IF (MOD(it,1000).EQ.0) THEN
      call allocate(psio); psio%e = psi%e; CALL RTRAN(psio,1)
      call allocate(chio); chio%e = chi%e; CALL RTRAN(chio,1)
      DO II = 1,SIZE(PSI%E,2)
         DO JJ = 1,SIZE(PSI%E,3)
            IF (((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) &
               .OR. ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1))) THEN
               CALL save_mode(tim%t,psio%E(:,II,JJ),chio%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
            ENDIF
         ENDDO
      ENDDO
      CALL DEALLOCATE(psio); CALL DEALLOCATE(chio)
   ENDIF
   ! IF (it.eq.1001) CALL MPI_ABORT(MPI_COMM_IVP,1,IERR)

enddo

999 continue
time_end = mpi_wtime()

!> output
!> 1. FFF space: perturbation file
ALLOCATE(PSI_GLB(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM),CHI_GLB(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM))
CALL MASSEMBLE(PSI%E,PSI_GLB,1); CALL MASSEMBLE(CHI%E,CHI_GLB,1)
IF (MPI_RANK.EQ.0) THEN
   ALLOCATE(vec_1(2*NRCHOPDIM),vec_2(2*NRCHOPDIM),vec_bar(2*NRCHOPDIM))
   vec_1(1:NRCHOPDIM) = PSI_GLB(:,MONITOR_MK(2,1)+1,MONITOR_MK(2,2)+1)
   vec_1(NRCHOPDIM+1:) = CHI_GLB(:,MONITOR_MK(2,1)+1,MONITOR_MK(2,2)+1)
   vec_2(1:NRCHOPDIM) = PSI_GLB(:,MONITOR_MK(3,1)+1,MONITOR_MK(3,2)+1)
   vec_2(NRCHOPDIM+1:) = CHI_GLB(:,MONITOR_MK(3,1)+1,MONITOR_MK(3,2)+1) 
   vec_bar(1:NRCHOPDIM) = PSI_GLB(:,MONITOR_MK(1,1)+1,MONITOR_MK(1,2)+1)
   vec_bar(NRCHOPDIM+1:) = CHI_GLB(:,MONITOR_MK(1,1)+1,MONITOR_MK(1,2)+1)
   CALL SAVE_PERTURB('./converg/new_perturb_refined.input', &
      M(MONITOR_MK(2,1)), AK(MONITOR_MK(2,1),MONITOR_MK(2,2)),vec_1, CMPLX(eig1), &
      M(MONITOR_MK(3,1)), AK(MONITOR_MK(3,1),MONITOR_MK(3,2)),vec_2, CMPLX(eig2), &
      M(MONITOR_MK(1,1)), AK(MONITOR_MK(1,1),MONITOR_MK(1,2)),vec_bar, CMPLX(0.D0))
   DEALLOCATE(vec_1,vec_2,vec_bar)
ENDIF
DEALLOCATE(PSI_GLB,CHI_GLB)

!> 2. PFF space: PSI-CHI
CALL RTRAN(psi,1); CALL RTRAN(chi,1)
DO II = 1,SIZE(PSI%E,2)
   DO JJ = 1,SIZE(PSI%E,3)
      IF (((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) &
         .OR. ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1))) THEN
         CALL save_mode(tim%t,psi%E(:,II,JJ),chi%E(:,II,JJ),II+PSI%INTH-1,JJ+PSI%INX-1)
      ENDIF
   ENDDO
ENDDO
CALL RTRAN(psi,-1); CALL RTRAN(chi,-1)

!> final printout
IF (MPI_RANK.eq.0) THEN
    print *,tim%n,' steps'
    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()
    WRITE(*,*) 'EXECUTION TIME: ',time_end-time_start,'seconds'
ENDIF
call MPI_BARRIER(MPI_COMM_IVP,IERR)
call MPI_FINALIZE(IERR)

! ======================================================================
contains
! ======================================================================

subroutine save_mode(t,psi_new,chi_new,m_i,k_i)
! ======================================================================
   complex(p8),DIMENSION(:):: psi_new,chi_new
   real:: t
   integer:: m_i,k_i
   integer:: nn

   open(UNIT=777,FILE='mode_track_'//ITOA3(m_i)//'_'//ITOA3(k_i)//'.dat',&
   &STATUS='UNKNOWN',ACTION='WRITE',ACCESS='APPEND')

   WRITE(777,320) t,psi_new(1:size(psi_new)),chi_new(1:size(psi_new))

   close(777)
320 FORMAT(F10.3,',',(S,E14.6E3,SP,E14.6E3,'i'),*(','S,E14.6E3,SP,E14.6E3,'i'))
end subroutine save_mode
! ======================================================================

subroutine NONLIN_LIN(PSI1,CHI1,PSIN,CHIN)
! ======================================================================
   TYPE(SCALAR),INTENT(IN):: PSI1,CHI1
   TYPE(SCALAR),INTENT(INOUT)::PSIN,CHIN

   INTEGER:: I1,J1

   CALL NONLIN(PSI1,CHI1,PSIN,CHIN)

   DO I1 = 1,SIZE(PSI1%E,2); DO J1 = 1,SIZE(PSI1%E,3)
      ! IF (((I1+PSI1%INTH .EQ. MONITOR_MK(2,1)+1).AND.(J1+PSI1%INX .EQ. MONITOR_MK(2,2)+1)) &
      !    .OR. ((I1+PSI1%INTH .EQ. MONITOR_MK(3,1)+1).AND.(J1+PSI1%INX .EQ. MONITOR_MK(3,2)+1))) THEN

      !    open(UNIT=888,FILE='eig_mode.dat',STATUS='UNKNOWN',ACTION='WRITE',ACCESS='APPEND')
      !    write(888,*) tim%t,SUM(PSIN%e(1:10,I1,J1)/PSI1%e(1:10,I1,J1))/10,I1+PSI1%INTH,J1+PSI1%INX
      !    close(888)

      !    CYCLE
      ! ENDIF
      IF ((I1+PSI1%INTH .EQ. MONITOR_MK(2,1)+1).AND.(J1+PSI1%INX .EQ. MONITOR_MK(2,2)+1)) THEN

         open(UNIT=888,FILE='eig_mode_1.dat',STATUS='UNKNOWN',ACTION='WRITE',ACCESS='APPEND')
         write(888,257) tim%t,SUM(PSIN%e(1:10,I1,J1)/PSI1%e(1:10,I1,J1))/10
         ! write(888,257) tim%t,SUM(PSIN%e(2:21,I1,J1)/PSI1%e(2:21,I1,J1))/20
         ! write(888,257) tim%t,SUM(PSIN%e(2:41,I1,J1)/PSI1%e(2:41,I1,J1))/40
         close(888)

         CYCLE
      ENDIF
      IF ((I1+PSI1%INTH .EQ. MONITOR_MK(3,1)+1).AND.(J1+PSI1%INX .EQ. MONITOR_MK(3,2)+1)) THEN

         open(UNIT=999,FILE='eig_mode_2.dat',STATUS='UNKNOWN',ACTION='WRITE',ACCESS='APPEND')
         write(999,257) tim%t,SUM(PSIN%e(1:10,I1,J1)/PSI1%e(1:10,I1,J1))/10
         ! write(999,257) tim%t,SUM(PSIN%e(2:21,I1,J1)/PSI1%e(2:21,I1,J1))/20
         ! write(999,257) tim%t,SUM(PSIN%e(2:41,I1,J1)/PSI1%e(2:41,I1,J1))/40
         close(999)

         CYCLE
      ENDIF
257 FORMAT(F10.3,',',(S,E26.15E3,SP,E26.15E3,'i'))
      ! M = 0, K != 0
      IF (((MONITOR_MK(2,1).EQ.0).AND.(I1+PSI1%INTH .EQ. MONITOR_MK(2,1)+1).AND.(J1+PSI1%INX .EQ. 2*NXCHOP-MONITOR_MK(2,2))) &
         .OR. ((MONITOR_MK(3,1).EQ.0).AND.(I1+PSI1%INTH .EQ. MONITOR_MK(3,1)+1).AND.(J1+PSI1%INX .EQ. 2*NXCHOP-MONITOR_MK(3,2)))) CYCLE

      PSIN%E(:,I1,J1) = 0.D0
      CHIN%E(:,I1,J1) = 0.D0
   ENDDO; ENDDO
end subroutine
! ======================================================================

subroutine ADAMSB_LIN(PSI1,CHI1,PSINO,CHINO,factor)
! ======================================================================
   IMPLICIT NONE
   TYPE(SCALAR):: PSI1,CHI1,PSINO,CHINO
   TYPE(SCALAR):: PSIN,CHIN
   COMPLEX(P8):: factor

   TYPE(SCALAR):: PSI0,CHI0
   REAL(P8):: DT
   INTEGER:: I1, J1

   DT = TIM%DT
   TIM%T = TIM%T + DT
   TIM%N = TIM%N + 1
   ADV%X = ADV%X + ADV%UX*DT
   ADV%Y = ADV%Y + ADV%UY*DT

   CALL ALLOCATE( PSIN )
   CALL ALLOCATE( CHIN )

   !> calculate nonlinear term
   CALL NONLIN_LIN(PSI1,CHI1,PSIN,CHIN)

   !> renormalize the degen pair
   DO I1 = 1,SIZE(PSI1%E,2); DO J1 = 1,SIZE(PSI1%E,3)
      IF (((I1+PSI1%INTH .EQ. MONITOR_MK(2,1)+1).AND.(J1+PSI1%INX .EQ. MONITOR_MK(2,2)+1)) &
         .OR. ((I1+PSI1%INTH .EQ. MONITOR_MK(3,1)+1).AND.(J1+PSI1%INX .EQ. MONITOR_MK(3,2)+1))) THEN
         PSI1%e(:,I1,J1) = PSI1%e(:,I1,J1)*factor
         CHI1%e(:,I1,J1) = CHI1%e(:,I1,J1)*factor
         write(*,*)factor
      ENDIF
      ! M = 0, K != 0
      IF (((MONITOR_MK(2,1).EQ.0).AND.(I1+PSI1%INTH .EQ. MONITOR_MK(2,1)+1).AND.(J1+PSI1%INX .EQ. 2*NXCHOP-MONITOR_MK(2,2))) &
         .OR. ((MONITOR_MK(3,1).EQ.0).AND.(I1+PSI1%INTH .EQ. MONITOR_MK(3,1)+1).AND.(J1+PSI1%INX .EQ. 2*NXCHOP-MONITOR_MK(3,2)))) THEN
         PSI1%e(:,I1,J1) = PSI1%e(:,I1,J1)*conjg(factor)
         CHI1%e(:,I1,J1) = CHI1%e(:,I1,J1)*conjg(factor)
         ! write(*,*) psi%e(anchor_psi_loc,I1,J1),chi%e(anchor_psi_loc,I1,J1)
      ENDIF
   ENDDO; ENDDO

   !> update psi & chi (not normalized)
   PSI1%E = PSI1%E +DT*(1.5D0*PSIN%E -0.5D0*PSINO%E)
   CHI1%E = CHI1%E +DT*(1.5D0*CHIN%E -0.5D0*CHINO%E)

   PSINO=PSIN
   CHINO=CHIN

   CALL DEALLOCATE( PSIN )
   CALL DEALLOCATE( CHIN )

   RETURN
end subroutine
! ======================================================================

SUBROUTINE RICH_LIN(PSI1,CHI1,PSIN,CHIN)
!=======================================================================
   IMPLICIT NONE
   TYPE(SCALAR):: PSI1,CHI1,PSIN,CHIN

   TYPE(SCALAR):: PSI2,CHI2
   TYPE(SCALAR):: PN,CN
   REAL(P8):: DT, HDT

   DT = TIM%DT
   HDT = DT*0.5D0
   TIM%T = TIM%T + DT
   TIM%N = TIM%N + 1
   ADV%X = ADV%X + ADV%UX*DT
   ADV%Y = ADV%Y + ADV%UY*DT

   CALL ALLOCATE(PSI2)
   CALL ALLOCATE(CHI2)

   CALL ALLOCATE( PN )
   CALL ALLOCATE( CN )

   ! FIRST HALF-STEP
   CALL NONLIN_LIN(PSI1,CHI1,PSIN,CHIN)
   PSI2%E =PSI1%E + (HDT*PSIN%E)
   CHI2%E =CHI1%E + (HDT*CHIN%E)

   ! SECOND HALF-STEP
   CALL NONLIN_LIN(PSI2,CHI2,PN,CN)

   PSI2%E =PSI2%E + (HDT*PN%E)
   CHI2%E =CHI2%E + (HDT*CN%E)

   ! FULL STEP
   PSI1%E =PSI1%E +(DT)*PSIN%E
   CHI1%E =CHI1%E +(DT)*CHIN%E

   PSI1%E =(2.0D0*PSI2%E)-PSI1%E
   CHI1%E =(2.0D0*CHI2%E)-CHI1%E

   CALL DEALLOCATE( PSI2 )
   CALL DEALLOCATE( CHI2 )
   CALL DEALLOCATE( PN )
   CALL DEALLOCATE( CN )

   RETURN
END SUBROUTINE
!=======================================================================

! ============================== ARCHIVE ===============================
! EULER STEP ===========================================================
! ======================================================================
! TIM%T = TIM%T + TIM%DT
! TIM%N = TIM%N + 1
! ADV%X = ADV%X + ADV%UX*TIM%DT
! ADV%Y = ADV%Y + ADV%UY*TIM%DT

! CALL NONLIN(psi,chi,dpsi,dchi)

!> time integration (ONLY UPDATE DEGEN MODES!)
! DO II = 1,SIZE(PSI%E,2); DO JJ = 1,SIZE(PSI%E,3)
!    IF (((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) &
!       .OR. ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1))) THEN
!       psi%e(:,II,JJ) = psi%e(:,II,JJ) + dpsi%e(:,II,JJ)*TIM%DT
!       chi%e(:,II,JJ) = chi%e(:,II,JJ) + dchi%e(:,II,JJ)*TIM%DT
!    ENDIF
!    ! M = 0, K != 0
!    IF (((MONITOR_MK(2,1).EQ.0).AND.(II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. 2*NXCHOP-MONITOR_MK(2,2))) &
!       .OR. ((MONITOR_MK(3,1).EQ.0).AND.(II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. 2*NXCHOP-MONITOR_MK(3,2)))) THEN
!       psi%e(:,II,JJ) = psi%e(:,II,JJ) + dpsi%e(:,II,JJ)*TIM%DT
!       chi%e(:,II,JJ) = chi%e(:,II,JJ) + dchi%e(:,II,JJ)*TIM%DT
!    ENDIF
! ENDDO; ENDDO
! ======================================================================
! OLD RENORMALIZATION ==================================================
! ======================================================================
!    CALL RTRAN(psi,1)
!    DO II = 1,SIZE(psi%E,2)
!       DO JJ = 1,SIZE(psi%E,3)
!          IF ((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) THEN
!             anchor_factor = anchor_psi_val/psi%E(anchor_psi_loc,II,JJ)
!             WRITE(*,*) 'anchor_factor:', anchor_factor
!          ENDIF
!       ENDDO
!    ENDDO
!    CALL RTRAN(psi,-1)
!    CALL MPI_ALLREDUCE(MPI_IN_PLACE, anchor_factor, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
!    DO II = 1,SIZE(PSI%E,2)
!       DO JJ = 1,SIZE(PSI%E,3)
!          IF (((II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(2,2)+1)) &
!             .OR. ((II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. MONITOR_MK(3,2)+1))) THEN
!             psi%e(:,II,JJ) = psi%e(:,II,JJ)*anchor_factor
!             chi%e(:,II,JJ) = chi%e(:,II,JJ)*anchor_factor
!             ! write(*,*) psi%e(anchor_psi_loc,II,JJ),chi%e(anchor_psi_loc,II,JJ)
!          ENDIF
!          ! M = 0, K != 0
!          IF (((MONITOR_MK(2,1).EQ.0).AND.(II+PSI%INTH .EQ. MONITOR_MK(2,1)+1).AND.(JJ+PSI%INX .EQ. 2*NXCHOP-MONITOR_MK(2,2))) &
!             .OR. ((MONITOR_MK(3,1).EQ.0).AND.(II+PSI%INTH .EQ. MONITOR_MK(3,1)+1).AND.(JJ+PSI%INX .EQ. 2*NXCHOP-MONITOR_MK(3,2)))) THEN
!             psi%e(:,II,JJ) = psi%e(:,II,JJ)*conjg(anchor_factor)
!             chi%e(:,II,JJ) = chi%e(:,II,JJ)*conjg(anchor_factor)
!             ! write(*,*) psi%e(anchor_psi_loc,II,JJ),chi%e(anchor_psi_loc,II,JJ)
!          ENDIF
!       ENDDO
!    ENDDO
!    ! CALL RTRAN(psi,-1); CALL RTRAN(chi,-1)
!    ! CALL RTRAN(dpsi,-1); CALL RTRAN(dchi,-1)
! ======================================================================

end program vort8
