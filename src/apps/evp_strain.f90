PROGRAM EVP_STRAIN
!=======================================================================
! [USAGE]:
! USE VPROD_PFF RATHER THAN VPROD
!
! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
! EXPRESSED IN A POLOIDAL-TOROLIDALLY DECOMPOSED FORM
! OPERATOR H CORRESPONDS TO EIGENVECTOR: [PSI,DEL2CHI]^T
!
! [UPDATES]:
! LAST UPDATE ON DEC 3, 2020
!=======================================================================
USE omp_lib
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
USE MOD_INIT
USE MOD_EVP
USE MOD_DIAGNOSTICS

IMPLICIT NONE
INTEGER     :: II, JJ, KK
INTEGER     :: M_BAR
REAL(P8)    :: K_BAR
REAL(P8), DIMENSION(:), ALLOCATABLE:: RUR0, RUP0, UZ0, ROR0, ROP0, OZ0
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BAR, RUP_BAR, UZ_BAR, PSI_BAR, CHI_BAR

! MPI:
INTEGER:: MPI_GLB_PROCS, MPI_GLB_RANK, newcomm
! ADD PERTURB:
INTEGER:: IS, NDIM
! SCANNING:
CHARACTER(LEN=72):: FILENAME
INTEGER:: FID

! CALL MPI_INIT(IERR)
CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED, MPI_THREAD_MODE, IERR)
IF (MPI_THREAD_MODE .LT. MPI_THREAD_SERIALIZED) THEN
    WRITE (*, *) 'The threading support is lesser than that demanded.'
    CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
END IF

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE (*, *) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()   ! @ MOD_MISC
END IF
CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR) ! DIVIDE NR BY #OF MPI PROCS
CALL READCOM('NOECHO')
CALL READIN(5)

! {M_bar, AK_bar}: STRAIN FIELD
!=======================================================================
M_BAR = 2
K_BAR = 1.0

IF (MPI_GLB_RANK .EQ. 0) THEN

    WRITE(*,*) 'EVP_STRAIN: SET NEW M AND AK'
    ZLEN = 2*PI/K_BAR
    NTH = 2
    NX = 4
    NTCHOP = 2 ! M = 0, M = M_ACTUAL
    NXCHOP = 2 ! K = 0, K = AK_ACTUAL
    CALL LEGINIT(newcomm,M_BAR)

    ALLOCATE(PSI_BAR(NR),CHI_BAR(NR))
    ALLOCATE(RUR_BAR(NRCHOPDIM),RUP_BAR(NRCHOPDIM),UZ_BAR(NRCHOPDIM))
    PSI_BAR = STRAIN_FIELD(RC_IN = 8.0); WRITE(*,*) REAL(PSI_BAR)
    CHI_BAR = 0.D0
    CALL RTRAN_MK(PSI_BAR,-1)
    CALL RTRAN_MK(CHI_BAR,-1)
    ! PSI_BAR = SMOOTH(PSI_BAR,SIZE(PSI_BAR))
    WRITE(*,*) 'SPECTRAL:'
    WRITE(*,*) REAL(PSI_BAR)

    CALL CHOPSET(3)
    CALL PC2VEL_MK(PSI_BAR,CHI_BAR,RUR_BAR,RUP_BAR,UZ_BAR)
    CALL RTRAN_MK(RUR_BAR,1)
    CALL RTRAN_MK(RUP_BAR,1)
    CALL RTRAN_MK( UZ_BAR,1)
    CALL CHOPSET(-3)
    CALL SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './converg/strain_field_vel.output')

    DEALLOCATE(PSI_BAR,CHI_BAR,RUR_BAR,RUP_BAR,UZ_BAR)

    WRITE (*, *) ''
    WRITE (*, *) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
END IF

CALL MPI_FINALIZE(IERR)

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================
FUNCTION STRAIN_FIELD(R_IN, RC_IN, E_IN)
! ======================================================================
    IMPLICIT NONE
    REAL(P8),DIMENSION(:),INTENT(IN),OPTIONAL:: R_IN
    REAL(P8),OPTIONAL:: RC_IN, E_IN

    REAL(P8):: RC = 4, E = 1.D0, dR = 1.6
    REAL(P8),DIMENSION(:),ALLOCATABLE:: R
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: STRAIN_FIELD

    IF (PRESENT(R_IN)) THEN
        ALLOCATE(R(SIZE(R_IN)),STRAIN_FIELD(SIZE(R_IN)))
    ELSE
        ALLOCATE(R(SIZE(TFM%R)),STRAIN_FIELD(SIZE(TFM%R)))
        R = TFM%R
    ENDIF
    STRAIN_FIELD = 0.D0

    IF (PRESENT(RC_IN)) THEN
        RC = RC_IN
    ELSE
        RC = ELL
    ENDIF
    IF (PRESENT(E_IN)) E = E_IN

    ! WHERE (R.LT.RC)
    !     STRAIN_FIELD = 0.5*E*R**2
    ! ELSEWHERE
    !     STRAIN_FIELD = 0.5*E*R**2*EXP(-(R-RC))
    ! ENDWHERE
    ! STRAIN_FIELD = SMOOTH(STRAIN_FIELD,SIZE(STRAIN_FIELD),3)

    STRAIN_FIELD = 0.5*E*R**2*0.5*(1-TANH((R-RC)/dR))

    ! WHERE (R.LT.RC)
    !     STRAIN_FIELD = 0.5*E*R**2 - 0.5*E*RC**2
    ! ELSEWHERE
    !     STRAIN_FIELD = 0
    ! ENDWHERE

    STRAIN_FIELD = STRAIN_FIELD*0.5

END FUNCTION
! ======================================================================

FUNCTION SMOOTH(fnc,N,hwidth)
! ======================================================================
! Smoothing a cmplx set using a modified moving average filter
! (note: based on the spatial distribution of R)
! ======================================================================
COMPLEX(P8),DIMENSION(:):: fnc
COMPLEX(P8):: SMOOTH(N)
INTEGER:: N, I
INTEGER, OPTIONAL:: hwidth
REAL(P8),DIMENSION(:),ALLOCATABLE:: R

IF (N.ne.SIZE(fnc)) THEN 
    WRITE(*,*) 'DIMENSION MISMATCH'
ENDIF

SMOOTH(1) = fnc(1)
SMOOTH(N) = fnc(N)

IF (.NOT.(PRESENT(hwidth))) THEN
    IF (N.ne.SIZE(TFM%R)) THEN ! IN THE FFF-SPACE
        SMOOTH(2:N-1) = (fnc(1:N-2)+fnc(2:N-1)+fnc(3:N))/3
    ELSE
        ALLOCATE(R(SIZE(TFM%R)))
        R = TFM%R
        SMOOTH(2:N-1) = fnc(3:N)*(R(2:N-1)-R(1:N-2))/(R(3:N)-R(1:N-2))
        SMOOTH(2:N-1) = SMOOTH(2:N-1) + &
                    & fnc(1:N-2)*(R(3:N)-R(2:N-1))/(R(3:N)-R(1:N-2))
        SMOOTH(2:N-1) = (SMOOTH(2:N-1)*2+fnc(2:N-1))/3
        DEALLOCATE(R)
    ENDIF
ELSE
    SMOOTH(1+hwidth:N-hwidth) = fnc(1+hwidth:N-hwidth)
    DO I=1,hwidth
        SMOOTH(1+hwidth-I) = SUM(fnc(1:2*(hwidth-I)+1))/(2*(hwidth-I)+1)
        SMOOTH(N-hwidth+I) = SUM(fnc(N-2*(hwidth-I):N))/(2*(hwidth-I)+1)
        SMOOTH(1+hwidth:N-hwidth) = SMOOTH(1+hwidth:N-hwidth) + fnc(1+hwidth-I:N-hwidth-I)
        SMOOTH(1+hwidth:N-hwidth) = SMOOTH(1+hwidth:N-hwidth) + fnc(1+hwidth+I:N-hwidth+I)
    ENDDO
    SMOOTH(1+hwidth:N-hwidth) = SMOOTH(1+hwidth:N-hwidth)/(2*hwidth+1)
ENDIF

END FUNCTION
! ======================================================================

SUBROUTINE EIG_MATRIX(MREAD, AKREAD, H, EIG_VAL, EIG_VEC_R, EIG_VEC_L, comm_grp, switch)
! ======================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: MREAD, comm_grp
    REAL(P8), INTENT(IN)   :: AKREAD
    COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT)            :: H
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT), OPTIONAL     :: EIG_VAL
    COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT), OPTIONAL  :: EIG_VEC_R, EIG_VEC_L
    LOGICAL, OPTIONAL:: switch

    LOGICAL     :: flip_switch =.FALSE.
    INTEGER     :: NR_MK, I, M_ACTUAL
    REAL(P8)    :: REI, AK_ACTUAL
    INTEGER     :: MPI_NR_SIZE, MPI_NR_INDEX
! REAL(P8),DIMENSION(:),ALLOCATABLE:: RUR0,RUP0,UZ0,ROR0,ROP0,OZ0
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RURU, RUPU, UZU, RORU, ROPU, OZU
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: PSIU, CHIU, PSI1, CHI1, PSI2, CHI2, PSI3, CHI3

    IF (MPI_GLB_RANK.EQ.0) THEN
        ! IF MREAD IS NEGATIVE, CALCULATE ITS C.C MODE
        IF (MREAD.LT.0) THEN
            flip_switch = .TRUE.
            M_ACTUAL = -MREAD
            AK_ACTUAL = -AKREAD
        ELSE
            flip_switch = .FALSE.
            M_ACTUAL = MREAD
            AK_ACTUAL = AKREAD
        ENDIF
        ZLEN = 2*PI/AK_ACTUAL
    ENDIF
    CALL MPI_BCAST(M_ACTUAL, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_BCAST(ZLEN, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)

! INIT
    NTH = 2
    NX = 4
    NTCHOP = 2 ! M = 0, M = MREAD
    NXCHOP = 2 ! K = 0, K = AKREAD
    CALL LEGINIT(comm_grp, M_ACTUAL)
    IF (.NOT. ALLOCATED(RUR0)) THEN
        ALLOCATE (RUR0(NR), RUP0(NR), UZ0(NR), ROR0(NR), ROP0(NR), OZ0(NR))
        CALL INIT_LOOP(RUR0, RUP0, UZ0, ROR0, ROP0, OZ0) ! ALL IN PFF space
    END IF

! SPLIT NRCHOPDIM
    NR_MK = NRCHOPS(2)
    IF (ALLOCATED(H) .AND. (SIZE(H, 1) .NE. NR_MK)) THEN
        DEALLOCATE (H)
        ALLOCATE (H(2*NR_MK, 2*NR_MK))
    ELSEIF (.NOT. (ALLOCATED(H))) THEN
        ALLOCATE (H(2*NR_MK, 2*NR_MK))
    END IF
    H = CMPLX(0.D0, 0.D0, P8)
    CALL DECOMPOSE(NR_MK, MPI_GLB_PROCS, MPI_GLB_RANK, MPI_NR_SIZE, MPI_NR_INDEX)

! MAIN JOB
    ALLOCATE (PSIU(NRCHOPDIM), CHIU(NRCHOPDIM))
    ALLOCATE (PSI1(NRCHOPDIM), CHI1(NRCHOPDIM))
    ALLOCATE (PSI2(NRCHOPDIM), CHI2(NRCHOPDIM))
    ALLOCATE (PSI3(NRCHOPDIM), CHI3(NRCHOPDIM))
    ALLOCATE (RURU(NRCHOPDIM), RUPU(NRCHOPDIM), UZU(NRCHOPDIM))
    ALLOCATE (RORU(NRCHOPDIM), ROPU(NRCHOPDIM), OZU(NRCHOPDIM))

    DO I = MPI_NR_INDEX + 1, MPI_NR_INDEX + MPI_NR_SIZE
! ================================ PSI =================================
        PSIU = 0.D0
        CHIU = 0.D0
        PSIU(I) = 1.D0

! VISCOSITY/HYPERV
        CALL CHOPSET(2)
        IF (VISC%SW .EQ. 1) THEN
            CALL DEL2_MK(PSIU, PSI3)
            PSI3 = PSI3*VISC%NU
        ELSE IF (VISC%SW .EQ. 2) THEN
            CALL CHOPSET(-2 + VISC%P)
            CALL HELMP_MK(VISC%P, PSIU, PSI3, 0.D0)
            IF (MOD(VISC%P/2, 2) .EQ. 0) PSI3 = -PSI3
            PSI3 = PSI3*VISC%NUP
            CALL CHOPSET(2 - VISC%P)
        ELSE
            PSI3 = 0.D0
        END IF
        CHI3 = 0.D0
        CALL CHOPSET(-2)

        CALL CHOPSET(3)
! NONLINEAR TERM: W0 X U'
        CALL PC2VEL_MK(PSIU, CHIU, RURU, RUPU, UZU)
        CALL RTRAN_MK(RURU, 1)
        CALL RTRAN_MK(RUPU, 1)
        CALL RTRAN_MK(UZU, 1)
        CALL VPROD_MK(ROR0, ROP0, OZ0, RURU, RUPU, UZU)
        CALL PROJECT_MK(RURU, RUPU, UZU, PSI1, CHI1)

! NONLINEAR TERM: W' X U0
        CALL PC2VOR_MK(PSIU, CHIU, RORU, ROPU, OZU)
        CALL RTRAN_MK(RORU, 1)
        CALL RTRAN_MK(ROPU, 1)
        CALL RTRAN_MK(OZU, 1)
        CALL VPROD_MK(RUR0, RUP0, UZ0, RORU, ROPU, OZU)
        CALL PROJECT_MK(RORU, ROPU, OZU, PSI2, CHI2)
        CALL CHOPSET(-3)
! IF (I.EQ.2) THEN
!     CALL XXDX_MK(PSIU,RURU)
!     CALL MCAT(ROP0)
! ENDIF

        PSI1 = -PSI1 + PSI2 + PSI3
        CHI1 = -CHI1 + CHI2 + CHI3
        PSI1(NRCHOPS(2) + 1:) = 0.D0
        CHI1(NRCHOPS(2) + 1:) = 0.D0

        H(:NR_MK, I) = PSI1(:NR_MK)
        H(NR_MK + 1:, I) = CHI1(:NR_MK)

! ================================ CHI =================================
        PSIU = 0.D0
        CHIU = 0.D0
        CHIU(I) = 1.D0

        CALL IDEL2_MK(CHIU, CHI1)
        CHIU = CHI1
        CHI1 = 0.D0

        IF (VISC%SW .NE. 0) THEN
            CHI3 = PSI3
        ELSE
            CHI3 = 0.D0
        END IF
        PSI3 = 0.D0

        CALL CHOPSET(3)
! NONLINEAR TERM: W0 X U'
        CALL PC2VEL_MK(PSIU, CHIU, RURU, RUPU, UZU)
        CALL RTRAN_MK(RURU, 1)
        CALL RTRAN_MK(RUPU, 1)
        CALL RTRAN_MK(UZU, 1)
        CALL VPROD_MK(ROR0, ROP0, OZ0, RURU, RUPU, UZU)
        CALL PROJECT_MK(RURU, RUPU, UZU, PSI1, CHI1)

! NONLINEAR TERM: W' X U0
        CALL PC2VOR_MK(PSIU, CHIU, RORU, ROPU, OZU)
        CALL RTRAN_MK(RORU, 1)
        CALL RTRAN_MK(ROPU, 1)
        CALL RTRAN_MK(OZU, 1)
        CALL VPROD_MK(RUR0, RUP0, UZ0, RORU, ROPU, OZU)
        CALL PROJECT_MK(RORU, ROPU, OZU, PSI2, CHI2)
        CALL CHOPSET(-3)

        PSI1 = -PSI1 + PSI2 + PSI3
        CHI1 = -CHI1 + CHI2 + CHI3
        PSI1(NRCHOPS(2) + 1:) = 0.D0
        CHI1(NRCHOPS(2) + 1:) = 0.D0

        H(:NR_MK, I + NR_MK) = PSI1(:NR_MK)
        H(NR_MK + 1:, I + NR_MK) = CHI1(:NR_MK)
    END DO
! DEALLOCATE( RUR0,RUP0,UZ0,ROR0,ROP0,OZ0 )
    DEALLOCATE (RURU, RUPU, UZU, RORU, ROPU, OZU)
    DEALLOCATE (PSIU, CHIU, PSI1, CHI1, PSI2, CHI2, PSI3, CHI3)

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, H, SIZE(H), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                        MPI_COMM_WORLD, IERR)

! OUTPUT:
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    IF (MPI_GLB_RANK .EQ. 0) THEN  ! START OF THE SERIAL PART
!WRITE OUT THE 3D MATRIX IN THE SCALAR USING UNFORMATTED MODE

        IF (flip_switch) THEN
            H = CONJG(H)
        ENDIF

        IF (PRESENT(EIG_VAL)) THEN
            IF (ALLOCATED(EIG_VAL)) DEALLOCATE (EIG_VAL)
            ALLOCATE (EIG_VAL(2*NR_MK))

            IF (PRESENT(EIG_VEC_L)) THEN
                IF (ALLOCATED(EIG_VEC_R)) DEALLOCATE (EIG_VEC_R)
                IF (ALLOCATED(EIG_VEC_L)) DEALLOCATE (EIG_VEC_L)
                ALLOCATE (EIG_VEC_R(2*NR_MK, 2*NR_MK), EIG_VEC_L(2*NR_MK, 2*NR_MK))
                CALL EIGENDECOMPOSE(H, EIG_VAL, ER=EIG_VEC_R, EL=EIG_VEC_L)
            ELSEIF (PRESENT(EIG_VEC_R)) THEN
                IF (ALLOCATED(EIG_VEC_R)) DEALLOCATE (EIG_VEC_R)
                ALLOCATE (EIG_VEC_R(2*NR_MK, 2*NR_MK))
                CALL EIGENDECOMPOSE(H, EIG_VAL, ER=EIG_VEC_R)
            ELSE
                CALL EIGENDECOMPOSE(H, EIG_VAL)
            END IF
        END IF

        IF (PRESENT(switch)) THEN
            open (10, FILE='./converg/eigM_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
                    //'.output', STATUS='unknown', ACTION='WRITE', IOSTAT=IS)
            if (IS .ne. 0) then
                print *, 'ERROR: EIG_MATRIX -- Could not creat new file'
                RETURN
            end if
            DO I = 1, 2*NR_MK
                WRITE (10, *) H(I, 1:20)
            END DO
            close (10)
            write (*, *) 'EVP MAT: max element = ', maxval(abs(H))
            write (*, *) 'EVP MAT: maxUR = ', maxval(abs(RUR0(1:NR)/TFM%R))
            write (*, *) 'EVP MAT: maxUP = ', maxval(abs(RUP0(1:NR)/TFM%R))
            write (*, *) 'EVP MAT: maxUZ = ', maxval(abs( UZ0(1:NR)))
            write (*, *) 'EVP MAT: maxOR = ', maxval(abs(ROR0(1:NR)/TFM%R))
            write (*, *) 'EVP MAT: maxOP = ', maxval(abs(ROP0(1:NR)/TFM%R))
            write (*, *) 'EVP MAT: maxOZ = ', maxval(abs( OZ0(1:NR)))
        END IF
    END IF

    RETURN

END SUBROUTINE EIG_MATRIX
!=======================================================================

END PROGRAM EVP_STRAIN
!=======================================================================
