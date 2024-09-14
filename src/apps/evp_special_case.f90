PROGRAM EVP_SPECIAL_CASE
!=======================================================================
! [USAGE]:
! GIVEN A STARTING GUESS, TUNE FOR A RESONANT TRIAD AND PERFORM DEGEN PT
! CALCULATE THE CORRECTION TERMS AND PERFORM FULL DIAGNOSTICS
! CREATE PERTURBED SYSTEM AND CALCULATE EIGVAL&EIGVEC
!
! [UPDATES]:
! LAST UPDATE ON MAY 10, 2024
! UPDATE ON FEB 28, 2023
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
INTEGER     :: N_mat1, N_mat2, II, JJ, KK
INTEGER     :: M_BAR, M_1, M_2, EIG_BAR_IND, tol_ind(2)
REAL(P8)    :: K_BAR, K_VAL_1_0, K_VAL_1_1, K_VAL_2_0, K_VAL_2_1
REAL(P8)    :: K_shift, tol, tol_diff, K_delta, EIG_LMR_1
COMPLEX(P8) :: EIG_BAR, J1, J2, JB
LOGICAL, DIMENSION(:), ALLOCATABLE:: MASK_1, MASK_2
REAL(P8), DIMENSION(:), ALLOCATABLE:: K_tune, K_tune_JJ
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig_1_0, M_eig_1_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig_2_0, M_eig_2_1
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_diff_1, M_mat_diff_2
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_1_0, EIG_R_mat_1_0, EIG_L_mat_1_0
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_2_0, EIG_R_mat_2_0, EIG_L_mat_2_0
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BAR, RUP_BAR, UZ_BAR, ROR_BAR, ROP_BAR, OZ_BAR
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BARC, RUP_BARC, UZ_BARC, ROR_BARC, ROP_BARC, OZ_BARC
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_B, VECR_B, VECR_12
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_1, VECR_1, VECR_B1, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_1_COPY, ROP_1_COPY, OZ_1_COPY
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_2, VECR_2, VECR_B2, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_2_COPY, ROP_2_COPY, OZ_2_COPY

! SPECIAL CASES
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BARC2, RUP_BARC2, UZ_BARC2, ROR_BARC2, ROP_BARC2, OZ_BARC2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_12, RUP_12, UZ_12, ROR_12, ROP_12, OZ_12
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECR_BCBC, VECR_BC1

! PERTURBED SYSTEM
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: SYS_EIG
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: P_mat_1, P_mat_2, SYS_ORG, SYS_BAR, SYS_NEW, SYS_EIG_VEC
! MATRIX WHOSE EIGVALS ARE THE SIGMA(1)s FOR NEARLY DEGEN
COMPLEX(P8), DIMENSION(2):: EIG_NDGEN
COMPLEX(P8), DIMENSION(2, 2):: M_NDGEN, VECR_NDGEN, VECL_NDGEN
! I/O:
INTEGER:: IS, FID
CHARACTER(LEN=72):: FILENAME
CHARACTER(len=6) :: K_VAL
LOGICAL:: DIAGNOST_SWITCH = .FALSE.

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

! {M_bar, AK_bar}: BASE FLOW
!=======================================================================
IF (MPI_GLB_RANK .EQ. 0) THEN

    FILENAME = 'almostdegen_read.input'
    CALL READIN(FILENAME, M_BAR, K_BAR, EIG_BAR_IND, M_1, K_VAL_1_0)

END IF ! MPI_GLB_RANK.EQ.0

CALL EIG_MATRIX(M_BAR, K_BAR, M_mat, M_eig, EIG_R_mat, EIG_L_mat, comm_grp=newcomm, print_switch=DIAGNOST_SWITCH) ! calculate operator matrix and its corresponding eigenvalues

! Pick a target Sigma_bar
IF (MPI_GLB_RANK .EQ. 0) THEN

    ! CALL TEST_ENSTROPHY(CMPLX(ROR0),CMPLX(ROP0),CMPLX(OZ0)); CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
    EIG_BAR = M_eig(EIG_BAR_IND)
    WRITE (*, *) 'TARGET EIG_BAR:', EIG_BAR
    DEALLOCATE (M_eig)
    ! =================== CHECK IF CRITICAL LAYER MODE =====================
    IF  (.NOT.((EIGRES(EIG_R_mat(:,EIG_BAR_IND),M_BAR)))) CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    ! ======================================================================

    ! Save rur_bar, rup_bar, uz_bar, ror_bar, rop_bar, oz_bar for DEGENERATE PERTURBATION THEORY
    CALL EIG2VELVOR(M_BAR,K_BAR,EIG_BAR_IND,EIG_R_mat,EIG_L_mat,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,VECR_B,VECL_B,comm_grp=newcomm)

END IF ! MPI_GLB_RANK.EQ.0
DEALLOCATE (M_mat)
CALL MPI_BCAST(EIG_BAR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)

! ! TUNE K
! !=======================================================================
! 0. Parameters:
K_delta = 0.001
tol = 1.0E-3  ! tolerance
K_shift = 1.D0 ! Initialization

! For M1,K1:
CALL MPI_BCAST(M_1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
CALL MPI_BCAST(K_VAL_1_0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
CALL EIG_MATRIX(M_1,K_VAL_1_0,M_mat_1_0,M_eig_1_0,EIG_R_mat_1_0,EIG_L_mat_1_0,comm_grp=newcomm)
! For M2,K2:
M_2 = M_1 + M_BAR
K_VAL_2_0 = K_VAL_1_0 + K_BAR
CALL EIG_MATRIX(M_2,K_VAL_2_0,M_mat_2_0,M_eig_2_0,EIG_R_mat_2_0,EIG_L_mat_2_0,comm_grp=newcomm)

DO II = 1,6

    ! 1. Check if tolerence is reached:
    ! Match the IMAG part
    IF (MPI_GLB_RANK.EQ.0) THEN

        ! Create resolving mask
        IF(.NOT.ALLOCATED(MASK_1)) THEN
            ALLOCATE(MASK_1(SIZE(EIG_R_mat_1_0,2)),MASK_2(SIZE(EIG_R_mat_2_0,2)))
        ENDIF
        MASK_1 = EIGRES(EIG_R_mat_1_0,M_1)
        MASK_2 = EIGRES(EIG_R_mat_2_0,M_2)

        ! Avoid same eigenmode for {mbar,kbar} and {m1,k1}
        IF ((M_BAR.EQ.M_1).AND.(ABS(K_BAR-K_VAL_1_0).LE.0.1)) THEN
            MASK_1 = MASK_1 .AND. (ABS(EIG_BAR-M_eig_1_0).GT.0.01)
        ENDIF
        
        ! Find the minimum distance
        tol_diff = MINVAL(ABS(AIMAG(MESHGRID(M_eig_1_0, M_eig_2_0) - EIG_BAR)),MESHGRID(MASK_1,MASK_2))/ABS(AIMAG(EIG_BAR))

        IF ((tol_diff.LT.tol).OR.(ABS(K_shift).LE.1.0E-4)) THEN

            tol_ind = MINLOC(ABS(AIMAG(MESHGRID(M_eig_1_0, M_eig_2_0) - EIG_BAR)),MESHGRID(MASK_1,MASK_2))
            JJ = tol_ind(1); KK = tol_ind(2)

            WRITE(*,*) 'Desired tol reached w/ ', ITOA3(II-1),' steps (tol_diff = ', tol_diff,')'
            ! WRITE(*,201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
            ! WRITE(*,201) M_1, K_VAL_1_0, JJ, M_eig_1_0(JJ)
            ! WRITE(*,201) M_2, K_VAL_2_0, kk, M_eig_2_0(KK)
            201  FORMAT('Eig: M# = ',I3,'; AK = ',F20.13,'; Sig (#',I3') =',1P8E20.13)

            ! TEST
            ! CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
            DEALLOCATE(MASK_1,MASK_2)
            tol = 999.0
        ENDIF
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(tol,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    IF (tol.EQ.999.0) goto 20

    ! 2. Calculate K_shift:
    ! For M1,K1:
    K_VAL_1_1 = K_VAL_1_0 + K_delta
    CALL EIG_MATRIX(M_1,K_VAL_1_1,M_mat_1_1,comm_grp=newcomm)
    ! For M2,K2:
    K_VAL_2_1 = K_VAL_2_0 + K_delta
    CALL EIG_MATRIX(M_2,K_VAL_2_1,M_mat_2_1,comm_grp=newcomm)

    ! MASTER RANK:
    IF (MPI_GLB_RANK.EQ.0) THEN

        IF(.NOT.ALLOCATED(M_mat_diff_1)) THEN
            N_mat1 = SIZE(EIG_R_mat_1_0,1)
            ALLOCATE(M_mat_diff_1(N_mat1,N_mat1))
        ENDIF
        M_mat_diff_1 = 0.D0
        M_mat_diff_1 = MATMUL((M_mat_1_1 - M_mat_1_0), EIG_R_mat_1_0)

        IF(.NOT.ALLOCATED(M_mat_diff_2)) THEN
            N_mat2 = SIZE(EIG_R_mat_2_0,1)
            ALLOCATE(M_mat_diff_2(N_mat2,N_mat2))
        ENDIF
        M_mat_diff_2 = 0.D0
        M_mat_diff_2 = MATMUL((M_mat_2_1 - M_mat_2_0), EIG_R_mat_2_0)        

        ! Find possible K_shift
        IF(.NOT.ALLOCATED(K_tune)) ALLOCATE(K_tune(N_mat1))
        IF(.NOT.ALLOCATED(K_tune_JJ)) ALLOCATE(K_tune_JJ(N_mat2))

        !$OMP parallel do private(EIG_LMR_1,KK,K_tune_JJ)
        DO JJ = 1,N_mat1
            ! eigc_b = -(m_bar*Q+k_bar*H);
            IF (.NOT.((EIGRES(EIG_R_mat_1_0(:,JJ),M_1)))) THEN
                K_tune(JJ) = 100.D0
            ! ELSEIF ((M_BAR.EQ.M_1).AND.(K_BAR.EQ.K_1)) THEN
            !     IF (EIG_BAR_IND.EQ.JJ) K_tune(JJ) = 100.D0
            ELSEIF ((M_BAR.EQ.M_1).AND.(ABS(K_BAR-K_VAL_1_0).LE.0.1)) THEN   
                IF (ABS(EIG_BAR-M_eig_1_0(JJ)).LE.0.01) K_tune(JJ) = 100.D0
            ELSE
                EIG_LMR_1 = 0.D0
                EIG_LMR_1 = &
                & AIMAG(DOT_PRODUCT(EIG_L_mat_1_0(:,JJ),M_mat_diff_1(:,JJ))/DOT_PRODUCT(EIG_L_mat_1_0(:,JJ),EIG_R_mat_1_0(:,JJ)))
    
                K_tune_JJ = 0.D0
                DO KK = 1,N_mat2
                    IF (.NOT.((EIGRES(EIG_R_mat_2_0(:,KK),M_2)))) THEN
                        K_tune_JJ(KK) = 100.D0
                    ELSE
                        K_tune_JJ(KK) = K_delta*AIMAG(M_eig_2_0(KK) - M_eig_1_0(JJ) - EIG_BAR) / ( &
                                    & EIG_LMR_1 - &
                                    & AIMAG(DOT_PRODUCT(EIG_L_mat_2_0(:,KK),M_mat_diff_2(:,KK))/DOT_PRODUCT(EIG_L_mat_2_0(:,KK),EIG_R_mat_2_0(:,KK))))
                    ENDIF
                ENDDO
                K_tune(JJ) = K_tune_JJ(MINLOC(ABS(K_tune_JJ),1))
            ENDIF
        ENDDO
        !$OMP end parallel do

        K_shift = K_tune(MINLOC(ABS(K_tune),1))
        WRITE(*,*) 'step#',ITOA3(II),'- K_shift:', K_shift
        IF (ABS(K_shift).GT.1) CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
    ENDIF

    ! 3. Update K_VAL_1_0 & K_VAL_2_0:
    K_VAL_1_0 = K_VAL_1_0 + K_shift
    CALL EIG_MATRIX(M_1,K_VAL_1_0,M_mat_1_0,M_eig_1_0,EIG_R_mat_1_0,EIG_L_mat_1_0,comm_grp=newcomm)
    K_VAL_2_0 = K_VAL_1_0 + K_BAR
    CALL EIG_MATRIX(M_2,K_VAL_2_0,M_mat_2_0,M_eig_2_0,EIG_R_mat_2_0,EIG_L_mat_2_0,comm_grp=newcomm)

ENDDO

! Maximum number of steps allowed
!=======================================================================
IF (II .EQ. 6) THEN ! Exceeds maximum number of steps allowed
    IF (MPI_GLB_RANK .EQ. 0) WRITE (*, *) 'Seem trapped w/ tol_diff:', tol_diff
    CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
END IF

! Saving and outputing the results in FFF Space
!=======================================================================
! Find JJth Eigvector of {M1,K1} & KKth Eigvector of {M2,K2}
20  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

! NEARLY DEGEN - ROUND TO NEAREST QUARTER
K_VAL_1_1 = K_VAL_1_0
K_VAL_2_1 = K_VAL_1_1 + K_BAR
CALL EIG_MATRIX(M_1, K_VAL_1_1, M_mat_1_1, M_eig_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1, comm_grp=newcomm)
CALL EIG_MATRIX(M_2, K_VAL_2_1, M_mat_2_1, M_eig_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1, comm_grp=newcomm)

IF (MPI_GLB_RANK.EQ.0) THEN

    tol_ind = MINLOC(ABS(AIMAG(MESHGRID(M_eig_1_1, M_eig_2_1) - EIG_BAR)),MESHGRID(EIGRES(EIG_R_mat_1_1,M_1),EIGRES(EIG_R_mat_2_1,M_2)))
    JJ = tol_ind(1); KK = tol_ind(2)

    WRITE (*, *) 'RESONANT TRIAD:'
    WRITE (*, 201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
    WRITE (*, 201) M_1, K_VAL_1_1, JJ, M_eig_1_1(JJ)
    WRITE (*, 201) M_2, K_VAL_2_1, kk, M_eig_2_1(KK)
    !=======================================================================

    ! Calculate the FINAL RESULT:
    ! M_bar, K_bar: RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,VECL_BAR
    ! {MBAR,KBAR}
    ALLOCATE (RUR_BARC(NDIMR), RUP_BARC(NDIMR), UZ_BARC(NDIMR), ROR_BARC(NDIMR), ROP_BARC(NDIMR), OZ_BARC(NDIMR))
    RUR_BARC = CONJG(RUR_BAR); RUP_BARC = CONJG(RUP_BAR); UZ_BARC = CONJG(UZ_BAR)
    ROR_BARC = CONJG(ROR_BAR); ROP_BARC = CONJG(ROP_BAR); OZ_BARC = CONJG(OZ_BAR)    

    ! {M1,K1}
    CALL EIG2VELVOR(M_1, K_VAL_1_1, JJ, EIG_R_mat_1_1, EIG_L_mat_1_1, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1, VECR_1, VECL_1, comm_grp=newcomm, print_switch=DIAGNOST_SWITCH)
    ALLOCATE(ROR_1_COPY(NDIMR),ROP_1_COPY(NDIMR),OZ_1_COPY(NDIMR))
    ROR_1_COPY = ROR_1; ROP_1_COPY = ROP_1; OZ_1_COPY = OZ_1

    ! {M2,K2}
    CALL EIG2VELVOR(M_2, K_VAL_2_1, KK, EIG_R_mat_2_1, EIG_L_mat_2_1, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2, VECR_2, VECL_2, comm_grp=newcomm, print_switch=DIAGNOST_SWITCH)
    ALLOCATE(ROR_2_COPY(NDIMR),ROP_2_COPY(NDIMR),OZ_2_COPY(NDIMR))
    ROR_2_COPY = ROR_2; ROP_2_COPY = ROP_2; OZ_2_COPY = OZ_2

    ! MAKE COPIES FOR SPECIAL CASE
    ALLOCATE(RUR_BARC2(NDIMR), RUP_BARC2(NDIMR), UZ_BARC2(NDIMR), ROR_BARC2(NDIMR), ROP_BARC2(NDIMR), OZ_BARC2(NDIMR))
    RUR_BARC2 = CONJG(RUR_BAR); RUP_BARC2 = CONJG(RUP_BAR); UZ_BARC2 = CONJG(UZ_BAR)
    ROR_BARC2 = CONJG(ROR_BAR); ROP_BARC2 = CONJG(ROP_BAR); OZ_BARC2 = CONJG(OZ_BAR)  
    ALLOCATE(RUR_12(NDIMR),RUP_12(NDIMR),UZ_12(NDIMR),ROR_12(NDIMR),ROP_12(NDIMR),OZ_12(NDIMR))
    RUR_12 = RUR_1; RUP_12 = RUP_1; UZ_12 = UZ_1
    ROR_12 = ROR_1; ROP_12 = ROP_1; OZ_12 = OZ_1

    ! NONLINEAR TERMS:
    ! {M2,K2}
    CALL NONLIN_MK(M_2,K_VAL_2_1,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1,VECR_B1,comm_grp=newcomm)

    ! {M1,K1}
    CALL NONLIN_MK(M_1,K_VAL_1_1,RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC,RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_B2,comm_grp=newcomm)

    ! {M1,K1} - SPECIAL CASE
    CALL NONLIN_MK(M_1,K_VAL_1_1,RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC,RUR_BARC2,RUP_BARC2,UZ_BARC2,ROR_BARC2,ROP_BARC2,OZ_BARC2,VECR_BCBC,comm_grp=newcomm)

    ! {MB,KB}
    CALL NONLIN_MK(M_BAR,K_BAR,CONJG(RUR_1),CONJG(RUP_1),CONJG(UZ_1),CONJG(ROR_1_COPY),CONJG(ROP_1_COPY),CONJG(OZ_1_COPY),RUR_2,RUP_2,UZ_2,ROR_2_COPY,ROP_2_COPY,OZ_2_COPY,VECR_12,comm_grp=newcomm)

    ! {MB,KB} - SPECIAL CASE
    CALL NONLIN_MK(M_BAR,K_BAR,RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC,RUR_12,RUP_12,UZ_12,ROR_12,ROP_12,OZ_12,VECR_BC1,comm_grp=newcomm)
    
ENDIF ! MPI_GLB_RANK.EQ.0

RANK0: IF (MPI_GLB_RANK .EQ. 0) THEN

! ============================== RESULT ================================
WRITE (*, *) "NEARLY DEGEN:"
JB = DOT_PRODUCT(VECL_B, VECR_12); WRITE (*,*) "JB: ", JB
J1 = DOT_PRODUCT(VECL_1, VECR_B2); WRITE (*,*) "J1: ", J1
J2 = DOT_PRODUCT(VECL_2, VECR_B1); WRITE (*,*) "J2: ", J2

WRITE (*,*) "JB-S: ", DOT_PRODUCT(VECL_B, VECR_BC1)
WRITE (*,*) "J1-S: ", DOT_PRODUCT(VECL_1, VECR_BCBC)
WRITE(*,*) "J1 + JB - J2: ", J1 + JB - J2

! DEALLOCATE
!=======================================================================
1   DEALLOCATE (RUR_BAR, RUP_BAR, UZ_BAR, ROR_BAR, ROP_BAR, OZ_BAR)
    DEALLOCATE (RUR_BARC, RUP_BARC, UZ_BARC, ROR_BARC, ROP_BARC, OZ_BARC)
    DEALLOCATE (RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1)
    DEALLOCATE (RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2)
    DEALLOCATE (ROR_1_COPY, ROP_1_COPY, OZ_1_COPY)
    DEALLOCATE (ROR_2_COPY, ROP_2_COPY, OZ_2_COPY)
    DEALLOCATE (VECR_B1, VECR_B2, VECR_12, VECR_1, VECR_2, VECL_1, VECL_2, VECL_B, VECR_B)
    DEALLOCATE (EIG_R_mat_1_1, EIG_L_mat_1_1, EIG_R_mat_2_1, EIG_L_mat_2_1)

END IF RANK0

! DEALLOCATE
!=======================================================================
10  CONTINUE
DEALLOCATE (RUR0, RUP0, UZ0, ROR0, ROP0, OZ0)
DEALLOCATE (M_mat_1_0, M_mat_2_0)
IF (ALLOCATED(M_mat_1_1)) THEN
    DEALLOCATE(M_mat_1_1,M_mat_2_1)
ENDIF
IF (ALLOCATED(RUR_BAR)) THEN
    DEALLOCATE(RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR)
    DEALLOCATE(RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC)
ENDIF

IF (MPI_GLB_RANK .EQ. 0) THEN
    DEALLOCATE(EIG_R_mat, EIG_L_mat)
    DEALLOCATE(M_eig_1_0, M_eig_1_1)
    DEALLOCATE(M_eig_2_0, M_eig_2_1)
    DEALLOCATE(EIG_R_mat_1_0, EIG_L_mat_1_0)
    DEALLOCATE(EIG_R_mat_2_0, EIG_L_mat_2_0)
    IF (ALLOCATED(K_tune)) THEN
        DEALLOCATE(K_tune, K_tune_JJ)
        DEALLOCATE(M_mat_diff_1)
        DEALLOCATE(M_mat_diff_2)
    ENDIF

    WRITE (*, *) ''
    WRITE (*, *) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
END IF

CALL MPI_FINALIZE(IERR)

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================

END PROGRAM EVP_SPECIAL_CASE
!=======================================================================
