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
INTEGER     :: M_BAR, M_1, M_2
REAL(P8)    :: K_BAR, K_VAL_1_0, K_VAL_1_1, K_VAL_2_0, K_VAL_2_1
REAL(P8)    :: K_delta, tol
COMPLEX(P8) :: EIG_PRIME_R, ALPHA, BETA, ALPHA_L, BETA_L
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECR_BAR, VECL_BAR
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: PSI_BAR, CHI_BAR, EIG_R_BAR
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BAR, RUP_BAR, UZ_BAR, ROR_BAR, ROP_BAR, OZ_BAR
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BARC, RUP_BARC, UZ_BARC, ROR_BARC, ROP_BARC, OZ_BARC

INTEGER     :: N_mat1
REAL(P8)    :: K_shift, EIG_LMR_1, tol_diff
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig_1_0, M_eig_1_1, M_eig_2_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_1, VECR_1, VECR_B1, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_2, VECR_2, VECR_B2, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECR_1_CORRECT, VECR_2_CORRECT
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECR_1_CORRECT2, VECR_2_CORRECT2

COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_diff_1
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_1_0, EIG_R_mat_1_0, EIG_L_mat_1_0
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1
REAL(P8), DIMENSION(:), ALLOCATABLE:: K_tune

! ADD PERTURB:
INTEGER:: IS, NDIM
! SCANNING:
CHARACTER(len=5) :: K_VAL
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
M_BAR = 2; K_BAR = 0.0

IF (MPI_GLB_RANK .EQ. 0) THEN

    WRITE(*,*) 'EVP_STRAIN: SET NEW M AND AK'
    ZLEN = 2*PI/K_BAR
    NTH = 2; NX = 4
    NTCHOP = 2 ! M = 0, M = M_ACTUAL
    NXCHOP = 2 ! K = 0, K = AK_ACTUAL
    CALL LEGINIT(newcomm,M_BAR)

    ALLOCATE(PSI_BAR(NR),CHI_BAR(NR))
    ALLOCATE(RUR_BAR(NRCHOPDIM),RUP_BAR(NRCHOPDIM),UZ_BAR(NRCHOPDIM))
    ALLOCATE(ROR_BAR(NRCHOPDIM),ROP_BAR(NRCHOPDIM),OZ_BAR(NRCHOPDIM))
    ! PSI_BAR = CMPLX(STRAIN_FIELD(RC_IN = 10.0)) !, e_in = 1.18)
    PSI_BAR = CMPLX(STRAIN_FIELD2(newcomm)) !, e_in = 1.18)
    ! PSI_BAR = CMPLX(STRAIN_DIPOLAR(newcomm)) !, e_in = 1.18)

    CHI_BAR = 0.D0

    CALL LEGINIT(newcomm,M_BAR)
    CALL RTRAN_MK(PSI_BAR,-1)
    CALL RTRAN_MK(CHI_BAR,-1)
    
    ALLOCATE(VECR_BAR(2*NRCHOPDIM))
    VECR_BAR(1:NRCHOPDIM) = PSI_BAR(1:NRCHOPDIM)
    VECR_BAR(NRCHOPDIM+1:) = CHI_BAR(1:NRCHOPDIM)
    CALL SAVE_PSICHI(VECR_BAR, './converg/check2/PC_o_MK_'//ITOA3(M_BAR)//'_'//' 0.00' &
        //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    DEALLOCATE(VECR_BAR)

    CALL CHOPSET(3)
    CALL PC2VEL_MK(PSI_BAR,CHI_BAR,RUR_BAR,RUP_BAR,UZ_BAR)
    CALL RTRAN_MK(RUR_BAR,1)
    CALL RTRAN_MK(RUP_BAR,1)
    CALL RTRAN_MK( UZ_BAR,1)

    CALL PC2VOR_MK(PSI_BAR,CHI_BAR,ROR_BAR,ROP_BAR,OZ_BAR)
    CALL RTRAN_MK(ROR_BAR,1)
    CALL RTRAN_MK(ROP_BAR,1)
    CALL RTRAN_MK( OZ_BAR,1)

    WRITE(*,*) 'ENS of STRAIN VELOCITY: ', ENSTROPHY_MK_MOD(ROR_BAR,ROP_BAR,OZ_BAR)/ZLEN0
    CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)

    CALL CHOPSET(-3)

    ! ALLOCATE(RUR_BAR(SIZE(TFM%R)),RUP_BAR(SIZE(TFM%R)),UZ_BAR(SIZE(TFM%R)))
    ! RUR_BAR = TFM%R**2*(-0.5*IU)
    ! RUP_BAR = TFM%R**2*(0.5)
    !  UZ_BAR = 0.D0

    CALL SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './converg/check2/strain_field_vel.output')
    CALL SAVE_VEL(ROR_BAR, ROP_BAR, OZ_BAR, './converg/check2/strain_field_vor.output')

    ALLOCATE(RUR_BARC(NDIMR), RUP_BARC(NDIMR), UZ_BARC(NDIMR)) 
    RUR_BARC = CONJG(RUR_BAR); RUP_BARC = CONJG(RUP_BAR); UZ_BARC = CONJG(UZ_BAR)

    ! ALLOCATE(ROR_BAR(NDIMR), ROP_BAR(NDIMR), OZ_BAR(NDIMR))
    ! ROR_BAR = 0.D0; ROP_BAR = 0.D0; OZ_BAR = 0.D0

    ALLOCATE(ROR_BARC(NDIMR), ROP_BARC(NDIMR), OZ_BARC(NDIMR))
    ROR_BARC = CONJG(ROR_BAR); ROP_BARC = CONJG(ROP_BAR); OZ_BARC = CONJG(OZ_BAR)

END IF

! {M1, K1}: 
! M_1 = 1; K_VAL_1_0 = 2.0;
! M_1 = 1; K_VAL_1_0 = 4.0;
! M_1 = 1; K_VAL_1_0 = 5.5;
M_1 = 1; K_VAL_1_0 = 7.0;

K_delta = 0.001; tol = 1.0E-6;
DO II = 1, 6
    ! For M1,K1:
    K_VAL_1_1 = K_VAL_1_0 + K_delta
    CALL EIG_MATRIX(M_1, K_VAL_1_0, M_mat_1_0, M_eig_1_0, EIG_R_mat_1_0, EIG_L_mat_1_0, comm_grp=newcomm)
    CALL EIG_MATRIX(M_1, K_VAL_1_1, M_mat_1_1, comm_grp=newcomm)

    ! TEST
    ! IF (MPI_GLB_RANK.EQ.0) WRITE(*,*) M_eig_1_0(1:10)
    ! M_2 = -M_1
    ! K_VAL_2_1 = K_VAL_1_1
    ! CALL EIG_MATRIX(M_2,K_VAL_2_1, M_mat_2_1, M_eig_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1, comm_grp=newcomm)
    ! IF (MPI_GLB_RANK.EQ.0) WRITE(*,*) M_eig_2_1(1:10)
    ! CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    ! CALL MPI_ABORT(MPI_COMM_WORLD,2,IERR)


    ! MASTER RANK:
    IF (MPI_GLB_RANK .EQ. 0) THEN
        IF (.NOT. ALLOCATED(M_mat_diff_1)) THEN
            N_mat1 = SIZE(EIG_R_mat_1_0, 1)
            ALLOCATE (M_mat_diff_1(N_mat1, N_mat1))
        END IF
        M_mat_diff_1 = 0.D0
        M_mat_diff_1 = MATMUL((M_mat_1_1 - M_mat_1_0), EIG_R_mat_1_0)

        ! Find possible K_shift
        IF (.NOT. ALLOCATED(K_tune)) ALLOCATE (K_tune(N_mat1))

        !$OMP parallel do private(EIG_LMR_1)
        DO JJ = 1, N_mat1
            IF (.NOT.(EIGRES(EIG_R_mat_1_0(:,JJ),M_1))) THEN
                M_eig_1_0(JJ) = IU*100
                K_tune(JJ) = 100.D0
            ELSE
                EIG_LMR_1 = 0.D0
                EIG_LMR_1 = &
                & AIMAG(DOT_PRODUCT(EIG_L_mat_1_0(:, JJ), M_mat_diff_1(:, JJ))/DOT_PRODUCT(EIG_L_mat_1_0(:, JJ), EIG_R_mat_1_0(:, JJ)))
                K_tune(JJ) = AIMAG(-M_eig_1_0(JJ))*K_delta/EIG_LMR_1
            END IF
        END DO
        !$OMP end parallel do

        K_shift = K_tune(MINLOC(ABS(K_tune),1))

        WRITE (*, *) 'STEP#', ITOA3(II), ' K_shift:', K_shift
        IF (ABS(K_shift) .GT. 1) CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
    END IF

    ! Calculate the actual
    K_VAL_1_1 = K_VAL_1_0 + K_shift
    CALL EIG_MATRIX(M_1, K_VAL_1_1, M_mat_1_1, M_eig_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1, comm_grp=newcomm)
    

    ! Check if tolerance is reached
    ! Match the IMAG part
    IF (MPI_GLB_RANK .EQ. 0) THEN

        M_eig_1_0 = M_eig_1_1
        DO JJ = 1, N_mat1
            IF (.NOT.(EIGRES(EIG_R_mat_1_1(:,JJ),M_1))) M_eig_1_0(JJ) = IU*100
        END DO
        tol_diff = MINVAL(ABS(AIMAG(M_eig_1_0)))
        JJ = MINLOC(ABS(AIMAG(M_eig_1_0)),1)

        WRITE (*, *) 'tol_diff:', tol_diff
        IF (tol_diff .LT. tol) THEN
            WRITE (*, *) 'Desired tol reached'
            WRITE (*, *) 'tol_diff:', tol_diff
            WRITE (*, *) 'Number of steps:', II
            tol = 999.0
            ! GOTO 20
        END IF
        K_VAL_1_0 = K_VAL_1_1

        WRITE (*, 201) M_1, K_VAL_1_1, JJ, M_eig_1_1(JJ)
        201 FORMAT('Eig: M# = ', I3, '; AK = ', F20.13, '; Sig (#', I3') =', 1P8E20.13)
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    CALL MPI_BCAST(tol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_BCAST(K_VAL_1_1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
    IF (tol .EQ. 999.0) goto 20

END DO

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

M_2 = -M_1
K_VAL_2_1 = K_VAL_1_1
CALL EIG_MATRIX(M_2, K_VAL_2_1, M_mat_2_1, M_eig_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1, comm_grp=newcomm)

IF (MPI_GLB_RANK .EQ. 0) THEN

    KK = JJ
    WRITE(K_VAL,'(F05.2)') K_VAL_1_1

    ! Q-VORTEX
    CALL SAVE_VEL(CMPLX(RUR0), CMPLX(RUP0), CMPLX(UZ0), './converg/check2/vel_MK_'//ITOA3(0)//'_'//' 0.00' &
        //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

    ! {M1,K1} = {1,2.26}
    CALL EIG2VELVOR(M_1, K_VAL_1_1, JJ, EIG_R_mat_1_1, EIG_L_mat_1_1, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1, VECR_1, VECL_1, comm_grp=newcomm)
    CALL SAVE_VEL(RUR_1, RUP_1, UZ_1, './converg/check2/vel_MK_'//ITOA3(M_1)//'_'//K_VAL &
        //'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    CALL SAVE_VEL(ROR_1, ROP_1, OZ_1, './converg/check2/vor_MK_'//ITOA3(M_1)//'_'//K_VAL &
        //'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! write(*,*) EIG_R_mat_1_1(1:5,JJ)

    ! {M2,K2} = {-1,2.26}
    CALL EIG2VELVOR(M_2, K_VAL_2_1, JJ, EIG_R_mat_2_1, EIG_L_mat_2_1, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2, VECR_2, VECL_2, comm_grp=newcomm)
    call SAVE_VEL(RUR_2, RUP_2, UZ_2, './converg/check2/vel_MK_'//ITOA3(M_2)//'_'//K_VAL &
        //'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    call SAVE_VEL(ROR_2, ROP_2, OZ_2, './converg/check2/vor_MK_'//ITOA3(M_2)//'_'//K_VAL &
        //'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! write(*,*) EIG_R_mat_2_1(1:5,JJ)

    CALL NONLIN_MK(M_2,K_VAL_2_1,RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC,RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1,VECR_B1,comm_grp=newcomm)
    CALL SAVE_PSICHI(VECR_B1, './converg/check2/PC_X_MK_'//ITOA3(M_1)//'_'//K_VAL &
        //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    CALL SAVE_PSICHI(VECL_2, './converg/check2/PC_L_MK_'//ITOA3(M_1)//'_'//K_VAL &
        //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

    ! {M1,K1}
    CALL NONLIN_MK(M_1,K_VAL_1_1,RUR_BAR ,RUP_BAR ,UZ_BAR ,ROR_BAR ,ROP_BAR ,OZ_BAR ,RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_B2,comm_grp=newcomm)

    ! ========================= DEGEN CALCULATION ==========================
    EIG_PRIME_R = DOT_PRODUCT(VECL_1, VECR_B2)*DOT_PRODUCT(VECL_2, VECR_B1)
    WRITE(*,*) 'DEGEN:'
    WRITE (*, 201) M_1, K_VAL_1_1, JJ, M_eig_1_1(JJ)
    WRITE (*, 201) M_2, K_VAL_2_1, JJ, M_eig_2_1(JJ)
    WRITE (*, *) " SIGMA(degen)'  = ", EIG_PRIME_R**0.5
    ! Beta/Alpha: ALPHA * {M1,K1} + BETA * {M2,K2}
    ALPHA = 1.D0; BETA = EIG_PRIME_R**0.5/DOT_PRODUCT(VECL_1, VECR_B2) ! BETA = (DOT_PRODUCT(VECL_2,VECR_B1)/DOT_PRODUCT(VECL_1,VECR_B2))**0.5
    WRITE (*, *) " BETA(degen)    = ", BETA
    ALPHA_L = 1.D0; BETA_L = ALPHA_L*CONJG(EIG_PRIME_R**0.5/DOT_PRODUCT(VECL_2, VECR_B1))
    WRITE (*, *) " BETA_L(degen)    = ", BETA_L

    ! ============================ DIAGNOSTICS =============================
    CALL TEST_NONLIN(VECL_1,VECR_B2,BETA)
    WRITE(*,*) "L1VR2: ", DOT_PRODUCT(VECL_1,VECR_B2)
    WRITE(*,*) "L2VR1: ", DOT_PRODUCT(VECL_2,VECR_B1)
    WRITE (*, *) "Orthogonality: (L1_new)V(R2_new):"
    WRITE (*, *) "A(+)BETA(-)L1VR2: ", DOT_PRODUCT(ALPHA_L*VECL_1, (-BETA)*VECR_B2)
    WRITE (*, *) "B(+)ALFA(-)L2VR1: ", DOT_PRODUCT(BETA_L*VECL_2, (ALPHA)*VECR_B1)
    WRITE (*, *) "Orthogonality: (L2_new)V(R1_new):"
    WRITE (*, *) "A(-)BETA(+)L1VR2: ", DOT_PRODUCT(-ALPHA_L*VECL_1, (BETA)*VECR_B2)
    WRITE (*, *) "B(-)ALFA(+)L2VR1: ", DOT_PRODUCT(-BETA_L*VECL_2, (ALPHA)*VECR_B1)
    ! WRITE(*,*) 'VECR_B2/EIG_R_1',VECR_B2(:)/EIG_R_mat_1_1(:,JJ)
    ! WRITE(*,*) 'VECR_B1/EIG_R_2',VECR_B2(:)/EIG_R_mat_2_1(:,KK)
    WRITE (*, *) "Normalization before: ", DOT_PRODUCT(VECL_1, VECR_1)
    WRITE (*, *) "Normalization before: ", DOT_PRODUCT(VECL_2, VECR_2)

    ! ====================== CORRECTION CALCULATION ========================
    ! NORMALIZE ALL EIGVEC EXCEPT #EIG_IND <- Needs to have -m support
    ! CALL EIG2NORMED(M_1, K_VAL_1_1, JJ, EIG_R_mat_1_1, EIG_L_mat_1_1, comm_grp=newcomm)
    ! CALL EIG2NORMED(M_2, K_VAL_2_1, JJ, EIG_R_mat_2_1, EIG_L_mat_2_1, comm_grp=newcomm)

    ! CALCULATE CORRECTION TERMS
    ! VECR_1_CORRECT = EIG_CORRECT_2(JJ, M_eig_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1, VECR_B2, M_1, K_VAL_1_1, newcomm)
    ! VECR_2_CORRECT = EIG_CORRECT_2(KK, M_eig_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1, VECR_B1, M_2, K_VAL_2_1, newcomm)
    ! VECR_1_CORRECT = VECR_1_CORRECT*BETA
    ! VECR_2_CORRECT = VECR_2_CORRECT*ALPHA

    VECR_1_CORRECT2 = EIG_CORRECT_3(JJ, M_eig_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1, VECR_B2, M_1, K_VAL_1_1, newcomm)
    VECR_2_CORRECT2 = EIG_CORRECT_3(KK, M_eig_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1, VECR_B1, M_2, K_VAL_2_1, newcomm)
    CALL SAVE_PERTURB('./converg/cor_perturb_org.input', M_1, K_VAL_1_1, VECR_1_CORRECT2, M_eig_1_1(JJ), &
                        M_2, K_VAL_2_1, VECR_2_CORRECT2, M_eig_2_1(kk)) ! SAVE ORIGINAL CORRECTION TERMS
    CALL SAVE_PERTURB('./converg/new_perturb_org.input', M_1, K_VAL_1_1, EIG_R_mat_1_1(:, JJ), M_eig_1_1(JJ), &
                        M_2, K_VAL_2_1, EIG_R_mat_2_1(:, KK), M_eig_2_1(kk)) ! SAVE ORIGINAL PERTURBATION TERMS

    ! FORM NEW EIGEN-VECTORS:
    VECR_1_CORRECT2 = VECR_1_CORRECT2*BETA
    VECR_2_CORRECT2 = VECR_2_CORRECT2*ALPHA    
    EIG_R_mat_1_1(:, JJ) = EIG_R_mat_1_1(:, JJ)*ALPHA
    EIG_R_mat_2_1(:, KK) = EIG_R_mat_2_1(:, KK)*BETA
    VECR_1 = VECR_1*ALPHA; VECR_2 = VECR_2*BETA
    VECL_1 = VECL_1*ALPHA_L; VECL_2 = VECL_2*BETA_L

    WRITE (*, *) "Normalization after: ", DOT_PRODUCT(VECL_1, VECR_1)
    WRITE (*, *) "Normalization after: ", DOT_PRODUCT(VECL_2, VECR_2)
    ! WRITE(*,*) 'DIFFERENCE - BetaVECR_B2 vs sigEIG_R_1',maxval(abs(BETA*VECR_B2(:) - EIG_PRIME_R**0.5*VECR_1))
    ! WRITE(*,*) 'DIFFERENCE - BetaVECR_B1 vs sigEIG_R_2',maxval(abs(ALPHA*VECR_B1(:) - EIG_PRIME_R**0.5*VECR_2))
    WRITE (*, *) "Diagnoalization:"
    WRITE (*, *) "(L1_new)V(R1_new): ", DOT_PRODUCT(VECL_1, (BETA)*VECR_B2)
    WRITE (*, *) "(L2_new)V(R2_new): ", DOT_PRODUCT(VECL_2, (ALPHA)*VECR_B1)
    WRITE(*,*) 'DIFFERENCE - BetaVECR_B2 vs sigEIG_R_1',sum((BETA*VECR_B2(:)+EIG_CORRECT_M(JJ, M_eig_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1, VECR_B2, M_1, K_VAL_1_1, newcomm))/(VECR_1))/size(VECR_1,1)
    ! WRITE(*,*) 'DIFFERENCE - BetaVECR_B1 vs sigEIG_R_2',(ALPHA*VECR_B1(:)+EIG_CORRECT_M(KK, M_eig_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1, VECR_B1, M_2, K_VAL_2_1, newcomm))/(VECR_2)

    ! WRITE(*,*) ALPHA*VECR_B1(:)/VECR_2

    ! ====================== SAVE PERTURBATION FILES =======================
    ! ======================================================================
    ALLOCATE(EIG_R_BAR(SIZE(EIG_R_mat_1_1,1)))
    EIG_R_BAR(1:SIZE(EIG_R_BAR)/2) = PSI_BAR
    EIG_R_BAR(SIZE(EIG_R_BAR)/2+1:) = CHI_BAR

    CALL SAVE_PERTURB('./converg/new_perturb.input', M_1, K_VAL_1_1, EIG_R_mat_1_1(:, JJ), M_eig_1_1(JJ), &
                        M_2, K_VAL_2_1, EIG_R_mat_2_1(:, KK), M_eig_2_1(kk), &
                        M_BAR, K_BAR, EIG_R_BAR, 0.D0*IU)
    CALL SAVE_PERTURB('./converg/lft_perturb.input', M_1, K_VAL_1_1, VECL_1, M_eig_1_1(JJ), &
                        M_2, K_VAL_2_1, VECL_2, M_eig_2_1(kk)) ! SAVE NEW LEFT EIGENVECTORS
    ! CALL SAVE_PERTURB('./converg/cor_perturb.input', M_1, K_VAL_1_1, VECR_1_CORRECT, M_eig_1_1(JJ), &
    !                     M_2, K_VAL_2_1, VECR_2_CORRECT, M_eig_2_1(kk)) ! SAVE CORRECTION TERMS
    CALL SAVE_PERTURB('./converg/cor_perturb.input', M_1, K_VAL_1_1, VECR_1_CORRECT2, M_eig_1_1(JJ), &
                        M_2, K_VAL_2_1, VECR_2_CORRECT2, M_eig_2_1(kk)) ! SAVE CORRECTION TERMS

    ! ! LAMBDA = 0.01
    ! CALL SAVE_PSICHI(EIG_R_mat_1_1(:, JJ), './converg/check2/PC_o_MK_'//ITOA3(M_1)//'_'//K_VAL &
    !     //'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! CALL SAVE_PSICHI(EIG_R_mat_2_1(:, KK), './converg/check2/PC_o_MK_'//ITOA3(M_2)//'_'//K_VAL &
    !     //'_IND_'//ITOA3(KK)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! CALL SAVE_PSICHI(VECR_1_CORRECT2, './converg/check2/PC_c_MK_'//ITOA3(M_1)//'_'//K_VAL &
    !     //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! CALL SAVE_PSICHI(VECR_2_CORRECT2, './converg/check2/PC_c_MK_'//ITOA3(M_2)//'_'//K_VAL &
    !     //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! CALL SAVE_PSICHI(EIG_R_mat_1_1(:, JJ)+0.01*VECR_1_CORRECT2, './converg/check2/PC_oc_MK_'//ITOA3(M_1)//'_'//K_VAL &
    !     //'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! CALL SAVE_PSICHI(EIG_R_mat_2_1(:, KK)+0.01*VECR_2_CORRECT2, './converg/check2/PC_oc_MK_'//ITOA3(M_2)//'_'//K_VAL &
    !     //'_IND_'//ITOA3(KK)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

    ! ============================ DEALLOCATE ==============================
    DEALLOCATE(RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1)
    DEALLOCATE(RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2)
    DEALLOCATE(VECR_1,VECL_1,VECR_2,VECL_2,VECR_B1,VECR_B2)
    ! DEALLOCATE(VECR_1_CORRECT,VECR_2_CORRECT)
    DEALLOCATE(VECR_1_CORRECT2,VECR_2_CORRECT2)

ENDIF

! FINALIZE
IF (MPI_GLB_RANK .EQ. 0) THEN
    DEALLOCATE(PSI_BAR,CHI_BAR)
    DEALLOCATE(RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR)
    DEALLOCATE(RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC)
    DEALLOCATE(M_eig_1_0, M_eig_1_1, M_eig_2_1, K_tune)
    DEALLOCATE(M_mat_diff_1)
    DEALLOCATE(M_mat_1_0, EIG_R_mat_1_0, EIG_L_mat_1_0)
    DEALLOCATE(M_mat_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1)
    DEALLOCATE(M_mat_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1)
    WRITE (*, *) ''
    WRITE (*, *) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC    
ENDIF

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

    REAL(P8):: RC = 4, E = 1.D0, dR = 2
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

    STRAIN_FIELD = -0.5*E*R**2*0.5*(1-TANH((R-RC)/dR))

    ! WHERE (R.LT.RC)
    !     STRAIN_FIELD = 0.5*E*R**2 - 0.5*E*RC**2
    ! ELSEWHERE
    !     STRAIN_FIELD = 0
    ! ENDWHERE

    STRAIN_FIELD = STRAIN_FIELD*0.5

END FUNCTION
! ======================================================================

FUNCTION STRAIN_FIELD2(comm_grp)
! ======================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN):: comm_grp
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE  :: STRAIN_FIELD2

    COMPLEX(P8), DIMENSION(:), ALLOCATABLE  :: EIG_VAL
    COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE:: H, EIG_VEC_R
    
    INTEGER     :: NR_MK, I, M_ACTUAL
    REAL(P8)    :: REI, AK_ACTUAL
    COMPLEX(P8) :: COEFF
    INTEGER     :: MPI_NR_SIZE, MPI_NR_INDEX
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: TU, T1, T2, T3, G

    ! OPERATOR:
    REAL(P8),DIMENSION(:),ALLOCATABLE:: func
    REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XXDXOP

    ! DAMPING:
    REAL(P8):: RC = 20, E = 1.D0, dR = 2
    REAL(P8),DIMENSION(:),ALLOCATABLE:: R

    ! INIT
    CALL LEGINIT(comm_grp,0)

    ! FORM Wz'/U_theta:
    ALLOCATE(func(NR))
    func = 4*(TFM%R)**4/(EXP((TFM%R)**2)-1) ! Q-VORTEX

    ! FORM RD()/DR OPERATOR FOR M = 0:
    ALLOCATE( XXDXOP(NRCHOP,3) )
    NR_MK = NRCHOPS(1)
    XXDXOP(:NR_MK,:) = BAND_LEG_XXDX(NR_MK,M(1),TFM%NORM(:,1))
    ! B(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,A(:NR_MK))

    ALLOCATE(H(NR_MK,NR_MK))
    H = CMPLX(0.D0, 0.D0, P8)

    ! MAIN JOB
    DO I = 1, NR_MK

        ALLOCATE (TU(NRCHOPDIM),T1(NRCHOPDIM),T2(NRCHOPDIM),T3(NDIMR))
        TU = 0.D0; T1 = 0.D0; T2 = 0.D0; T3 = 0.D0

        ! TEST FUNCTION:
        TU(I) = 1.D0

        ! OPERATOR:
        T1(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,TU(:NR_MK)) ! T1 = R*D[TU]/DR
        T2(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,T1(:NR_MK)) ! T2 = R*D[R*D[TU]/DR]/DR

        ! BACK TO PHYSICAL
        CALL RTRAN_MK(TU, 1)
        T3(:NR) = FUNC(:NR)*TU(:NR)
        CALL RTRAN_MK(T3,-1)

        ! FILL IN MATRIX
        H(:, I) = T2+4*T1+T3

        ! DEALLOCATE
        DEALLOCATE (TU,T1,T2,T3)

    END DO               

    ! SOLVE EVP:
    ALLOCATE(EIG_VAL(NR_MK), EIG_VEC_R(NR_MK, NR_MK))
    CALL EIGENDECOMPOSE(H, EIG_VAL, ER=EIG_VEC_R)

    ! LOCATE EIG=0
    I = MINLOC(ABS(EIG_VAL),1)
    ALLOCATE(G(NRCHOPDIM)); G(1:SIZE(EIG_VAL)) = EIG_VEC_R(:,I)
    CALL RTRAN_MK(G,1); G = G/G(NR)
    WRITE(*,*) G,TFM%R
    DEALLOCATE(EIG_VAL,EIG_VEC_R,H)

    ! FORM STREAM FUNCTION
    ALLOCATE(R(SIZE(TFM%R)),STRAIN_FIELD2(SIZE(TFM%R)))
    R = TFM%R
    STRAIN_FIELD2 = -0.5*E*R**2*G(1:SIZE(R))
    STRAIN_FIELD2 = STRAIN_FIELD2*0.5*(1-TANH((R-RC)/dR))*0.5

    DEALLOCATE(XXDXOP,FUNC)
    RETURN

END FUNCTION STRAIN_FIELD2
!=======================================================================

FUNCTION STRAIN_FIELD3(comm_grp)
! ======================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN):: comm_grp
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE  :: STRAIN_FIELD3

    COMPLEX(P8), DIMENSION(:), ALLOCATABLE  :: EIG_VAL
    COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE:: H, EIG_VEC_R
    
    INTEGER     :: NR_MK, I, M_ACTUAL
    REAL(P8)    :: REI, AK_ACTUAL
    COMPLEX(P8) :: COEFF
    INTEGER     :: MPI_NR_SIZE, MPI_NR_INDEX
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: TU, T1, T2, T3, G

    ! OPERATOR:
    REAL(P8),DIMENSION(:),ALLOCATABLE:: func
    REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XXDXOP

    ! DAMPING:
    REAL(P8):: RC = 20, E = 1.D0, dR = 2
    REAL(P8),DIMENSION(:),ALLOCATABLE:: R

    ! INIT
    CALL LEGINIT(comm_grp,0)

    ! FORM Wz'/U_theta:
    ALLOCATE(func(NR))
    func = 4*(TFM%R)**4/(EXP((TFM%R)**2)-1) ! Q-VORTEX

    ! FORM RD()/DR OPERATOR FOR M = 0:
    ALLOCATE( XXDXOP(NRCHOP,3) )
    NR_MK = NRCHOPS(1)
    XXDXOP(:NR_MK,:) = BAND_LEG_XXDX(NR_MK,M(1),TFM%NORM(:,1))
    ! B(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,A(:NR_MK))

    ALLOCATE(H(NR_MK,NR_MK))
    H = CMPLX(0.D0, 0.D0, P8)

    ! MAIN JOB
    DO I = 1, NR_MK

        ALLOCATE (TU(NRCHOPDIM),T1(NRCHOPDIM),T2(NRCHOPDIM),T3(NDIMR))
        TU = 0.D0; T1 = 0.D0; T2 = 0.D0; T3 = 0.D0

        ! TEST FUNCTION:
        TU(I) = 1.D0

        ! OPERATOR:
        T1(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,TU(:NR_MK)) ! T1 = R*D[TU]/DR
        T2(:NR_MK)=BANMUL(XXDXOP(:NR_MK,:),2,T1(:NR_MK)) ! T2 = R*D[R*D[TU]/DR]/DR

        ! BACK TO PHYSICAL
        CALL RTRAN_MK(TU, 1)
        T3(:NR) = FUNC(:NR)*TU(:NR)-4
        CALL RTRAN_MK(T3,-1)

        ! FILL IN MATRIX
        H(:, I) = T2+T3

        ! DEALLOCATE
        DEALLOCATE (TU,T1,T2,T3)

    END DO               

    ! SOLVE EVP:
    ALLOCATE(EIG_VAL(NR_MK), EIG_VEC_R(NR_MK, NR_MK))
    CALL EIGENDECOMPOSE(H, EIG_VAL, ER=EIG_VEC_R)

    ! LOCATE EIG=0
    I = MINLOC(ABS(EIG_VAL),1)
    ALLOCATE(G(NRCHOPDIM)); G(1:SIZE(EIG_VAL)) = EIG_VEC_R(:,I)
    CALL RTRAN_MK(G,1); G = G/G(NR)
    WRITE(*,*) EIG_VAL(I)
    WRITE(*,*) G,TFM%R,G/TFM%R
    DEALLOCATE(EIG_VAL,EIG_VEC_R,H)
    !!!!!!!!!!!!!!!!!!!!!!!!!
    CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)

    ! FORM STREAM FUNCTION
    ALLOCATE(R(SIZE(TFM%R)),STRAIN_FIELD3(SIZE(TFM%R)))
    R = TFM%R
    STRAIN_FIELD3 = -0.5*E*R**2*G(1:SIZE(R))
    STRAIN_FIELD3 = STRAIN_FIELD3*0.5*(1-TANH((R-RC)/dR))*0.5

    DEALLOCATE(XXDXOP,FUNC)
    RETURN

END FUNCTION STRAIN_FIELD3
! ======================================================================

FUNCTION STRAIN_DIPOLAR(comm_grp)
! ======================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN):: comm_grp
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE  :: STRAIN_DIPOLAR
  
    ! OPERATOR:
    REAL(P8),DIMENSION(:),ALLOCATABLE:: FUNC, VTHETA

    ! DAMPING:
    REAL(P8):: RC = 20, E = 1.D0, dR = 2
    REAL(P8),DIMENSION(:),ALLOCATABLE:: R

    ! INIT
    CALL LEGINIT(comm_grp,0)

    ! FORM INTEGRANTS
    ALLOCATE(R(NR),FUNC(NR),VTHETA(NR),STRAIN_DIPOLAR(NR))
    R = TFM%R
    VTHETA = (1-EXP(-R**2))/R; VTHETA(1) = 0.D0
    FUNC = 4*R*EXP(-R**2) + VTHETA

    ! NUMERICAL INTEGRATIONS
    FUNC = CUMTRAPZ(R,FUNC*VTHETA*R)
    FUNC = FUNC / (R * VTHETA**2); FUNC(1) = 0.D0
    
    ! PSI/R
    STRAIN_DIPOLAR = VTHETA*CUMTRAPZ(R,FUNC)/R; STRAIN_DIPOLAR(1) = 0.D0

    DEALLOCATE(FUNC,VTHETA,R)

    ! TEST
    WRITE(*,*) STRAIN_DIPOLAR
    CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)

    RETURN

END FUNCTION STRAIN_DIPOLAR
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

FUNCTION EIG_CORRECT_2(EIG_IND_DEGEN,EIG_VAL,EIG_VEC_R,EIG_VEC_L,VECR_B_DEGEN,MREAD,AKREAD,comm_grp)
! ======================================================================
! NOTE: VECR_B_DEGEN AND EIG_VEC_R CORRESPOND TO THE TWO WAVENUMBERS
!       e.g. VECR_B_DEGEN = N(V_BAR)V2, EIG_VEC_R/L & EIG_IND_DEGEN & EIG_VAL ~ V1
! ======================================================================
    IMPLICIT NONE
    INTEGER:: EIG_IND_DEGEN, MREAD, comm_grp
    REAL(P8):: AKREAD
    COMPLEX(P8), DIMENSION(:):: EIG_VAL, VECR_B_DEGEN
    COMPLEX(P8), DIMENSION(:, :):: EIG_VEC_R, EIG_VEC_L
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: EIG_CORRECT_2

    INTEGER:: IND, EIG_NUM, NR_MK, M_ACTUAL
    REAL(P8):: AK_ACTUAL
    COMPLEX(P8):: EIG_VAL_DEGEN, COEFF
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: CORRECTION
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: PSI, DEL2CHI, CHI
    LOGICAL:: flip_switch

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

    IF ((M(2) .NE. M_ACTUAL) .OR. (ABS(ZLEN - 2*PI/AK_ACTUAL) .GE. 1.0E-13)) THEN
        WRITE (*, *) 'EIG_CORRECT_2: SET NEW M AND AK'
        ZLEN = 2*PI/AK_ACTUAL
        NTH = 2
        NX = 4
        NTCHOP = 2 ! M = 0, M = MREAD
        NXCHOP = 2 ! K = 0, K = AKREAD
        CALL LEGINIT(comm_grp, M_ACTUAL)
    END IF

    EIG_VAL_DEGEN = EIG_VAL(EIG_IND_DEGEN)
    ALLOCATE (EIG_CORRECT_2(SIZE(EIG_VEC_R, 1)))
    ALLOCATE (CORRECTION(SIZE(EIG_VEC_R, 1)))
    ALLOCATE (PSI(NRCHOPDIM), DEL2CHI(NRCHOPDIM), CHI(NRCHOPDIM))

    EIG_CORRECT_2 = 0.D0
    EIG_NUM = SIZE(EIG_VEC_R, 1)
    NR_MK = SIZE(EIG_VEC_R, 2)/2

    DO IND = 1, SIZE(EIG_VEC_R, 2)
        IF (IND .EQ. EIG_IND_DEGEN) CYCLE

        ! calculate coefficient
        COEFF = DOT_PRODUCT(EIG_VEC_L(:, IND), VECR_B_DEGEN(:))&
        &/DOT_PRODUCT(EIG_VEC_L(:, IND), EIG_VEC_R(:, IND))&
        &/(EIG_VAL_DEGEN - EIG_VAL(IND))

        ! save original
        PSI = 0.D0
        DEL2CHI = 0.D0
        CHI = 0.D0
        PSI(:NR_MK) = EIG_VEC_R(:NR_MK, IND) !PSI
        DEL2CHI(:NR_MK) = EIG_VEC_R(NR_MK + 1:, IND) !DEL2CHI

        IF (flip_switch) DEL2CHI = CONJG(DEL2CHI)
        CALL IDEL2_MK(DEL2CHI, CHI)

        ! times coefficient
        CORRECTION(:NR_MK) = COEFF*PSI(:NR_MK)
        IF (flip_switch) CHI = CONJG(CHI)
        CORRECTION(NR_MK + 1:) = COEFF*CHI(:NR_MK)

        ! sum eig_vec
        EIG_CORRECT_2 = EIG_CORRECT_2 + CORRECTION

    END DO

    DEALLOCATE(PSI, DEL2CHI, CHI, CORRECTION)

END FUNCTION EIG_CORRECT_2
! ======================================================================

FUNCTION EIG_CORRECT_3(EIG_IND_DEGEN,EIG_VAL,EIG_VEC_R,EIG_VEC_L,VECR_B_DEGEN,MREAD,AKREAD,comm_grp)
! ======================================================================
! NOTE: VECR_B_DEGEN AND EIG_VEC_R CORRESPOND TO THE TWO WAVENUMBERS
!       e.g. VECR_B_DEGEN = N(V_BAR)V2, EIG_VEC_R/L & EIG_IND_DEGEN & EIG_VAL ~ V1
! ======================================================================
    IMPLICIT NONE
    INTEGER:: EIG_IND_DEGEN, MREAD, comm_grp
    REAL(P8):: AKREAD
    COMPLEX(P8), DIMENSION(:):: EIG_VAL, VECR_B_DEGEN
    COMPLEX(P8), DIMENSION(:, :):: EIG_VEC_R, EIG_VEC_L
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: EIG_CORRECT_3

    INTEGER:: IND, EIG_NUM, NR_MK, M_ACTUAL
    REAL(P8):: AK_ACTUAL
    COMPLEX(P8):: EIG_VAL_DEGEN, COEFF
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: CORRECTION
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: PSI, DEL2CHI, CHI
    LOGICAL:: flip_switch

    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: RUR,RUP,UZ
    CHARACTER(len=5) :: AK_VAL

    WRITE(AK_VAL,'(F05.2)') AKREAD

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

    IF ((M(2) .NE. M_ACTUAL) .OR. (ABS(ZLEN - 2*PI/AK_ACTUAL) .GE. 1.0E-13)) THEN
        WRITE (*, *) 'EIG_CORRECT_3: SET NEW M AND AK'
        ZLEN = 2*PI/AK_ACTUAL
        NTH = 2
        NX = 4
        NTCHOP = 2 ! M = 0, M = MREAD
        NXCHOP = 2 ! K = 0, K = AKREAD
        CALL LEGINIT(comm_grp, M_ACTUAL)
    END IF

    EIG_VAL_DEGEN = EIG_VAL(EIG_IND_DEGEN)
    ALLOCATE(EIG_CORRECT_3(SIZE(EIG_VEC_R, 1)))
    ALLOCATE(CORRECTION(SIZE(EIG_VEC_R, 1)))
    ALLOCATE(PSI(NRCHOPDIM), DEL2CHI(NRCHOPDIM), CHI(NRCHOPDIM))
    ALLOCATE(RUR(NRCHOPDIM),RUP(NRCHOPDIM),UZ(NRCHOPDIM))

    EIG_CORRECT_3 = 0.D0
    EIG_NUM = SIZE(EIG_VEC_R, 1)
    NR_MK = SIZE(EIG_VEC_R, 2)/2

    DO IND = 1, SIZE(EIG_VEC_R, 2)
        IF (IND .EQ. EIG_IND_DEGEN) CYCLE

        ! calculate coefficient
        COEFF = DOT_PRODUCT(EIG_VEC_L(:, IND), VECR_B_DEGEN(:))&
        &/DOT_PRODUCT(EIG_VEC_L(:, IND), EIG_VEC_R(:, IND))&
        &/(EIG_VAL_DEGEN - EIG_VAL(IND))

        ! sum eig_vec
        EIG_CORRECT_3 = EIG_CORRECT_3 + COEFF*EIG_VEC_R(:, IND) ! [PSI,DEL2CHI]

        IF (ABS(COEFF).GT.0.1) THEN
            WRITE(*,*) IND, abs(COEFF), ITOA3(MREAD), AK_VAL
        ENDIF

    END DO

    CALL SAVE_PSICHI(EIG_CORRECT_3, './converg/check2/PC_cd_MK_'//ITOA3(MREAD)//'_'//AK_VAL &
        //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

    PSI = 0.D0; DEL2CHI = 0.D0; CHI = 0.D0
    PSI(:NR_MK) = EIG_CORRECT_3(:NR_MK)
    DEL2CHI(:NR_MK) = EIG_CORRECT_3(NR_MK+1:)

    IF (flip_switch) THEN
        PSI = CONJG(PSI)
        DEL2CHI = CONJG(DEL2CHI)
    ENDIF
    CALL IDEL2_MK(DEL2CHI, CHI) !< add negative M support

    CALL CHOPSET(3)
    CALL PC2VEL_MK(PSI,CHI,RUR,RUP,UZ)
    CALL RTRAN_MK(RUR,1); CALL RTRAN_MK(RUP,1); CALL RTRAN_MK( UZ,1)
    CALL CHOPSET(-3)
    
    IF (flip_switch) THEN
        PSI = CONJG(PSI)
        CHI = CONJG(CHI)
        RUR = CONJG(RUR); RUP = CONJG(RUP); UZ = CONJG( UZ)
    ENDIF
    EIG_CORRECT_3(:NR_MK) = PSI(:NR_MK)
    EIG_CORRECT_3(NR_MK+1:) = CHI(:NR_MK)

    CALL SAVE_VEL(RUR, RUP, UZ, './converg/check2/vec_MK_'//ITOA3(MREAD)//'_'//AK_VAL &
        //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

    DEALLOCATE(PSI,DEL2CHI,CHI,CORRECTION)
    DEALLOCATE(RUR,RUP,UZ)

END FUNCTION EIG_CORRECT_3
! ======================================================================

FUNCTION EIG_CORRECT_M(EIG_IND_DEGEN,EIG_VAL,EIG_VEC_R,EIG_VEC_L,VECR_B_DEGEN,MREAD,AKREAD,comm_grp)
! ======================================================================
! NOTE: VECR_B_DEGEN AND EIG_VEC_R CORRESPOND TO THE TWO WAVENUMBERS
!       e.g. VECR_B_DEGEN = N(V_BAR)V2, EIG_VEC_R/L & EIG_IND_DEGEN & EIG_VAL ~ V1
! ======================================================================
    IMPLICIT NONE
    INTEGER:: EIG_IND_DEGEN, MREAD, comm_grp
    REAL(P8):: AKREAD
    COMPLEX(P8), DIMENSION(:):: EIG_VAL, VECR_B_DEGEN
    COMPLEX(P8), DIMENSION(:, :):: EIG_VEC_R, EIG_VEC_L
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: EIG_CORRECT_M

    INTEGER:: IND, EIG_NUM, NR_MK, M_ACTUAL
    REAL(P8):: AK_ACTUAL
    COMPLEX(P8):: EIG_VAL_DEGEN, COEFF
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: PSI, DEL2CHI, CHI
    LOGICAL:: flip_switch

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

    IF ((M(2) .NE. M_ACTUAL) .OR. (ABS(ZLEN - 2*PI/AK_ACTUAL) .GE. 1.0E-13)) THEN
        WRITE (*, *) 'EIG_CORRECT_M: SET NEW M AND AK'
        ZLEN = 2*PI/AK_ACTUAL
        NTH = 2
        NX = 4
        NTCHOP = 2 ! M = 0, M = MREAD
        NXCHOP = 2 ! K = 0, K = AKREAD
        CALL LEGINIT(comm_grp, M_ACTUAL)
    END IF

    EIG_VAL_DEGEN = EIG_VAL(EIG_IND_DEGEN)
    ALLOCATE (EIG_CORRECT_M(SIZE(EIG_VEC_R, 1)))
    ! ALLOCATE (PSI(NRCHOPDIM), DEL2CHI(NRCHOPDIM), CHI(NRCHOPDIM))

    EIG_CORRECT_M = 0.D0
    EIG_NUM = SIZE(EIG_VEC_R, 1)
    NR_MK = SIZE(EIG_VEC_R, 2)/2

    DO IND = 1, SIZE(EIG_VEC_R, 2)
        IF (IND .EQ. EIG_IND_DEGEN) CYCLE

        ! calculate coefficient
        COEFF = DOT_PRODUCT(EIG_VEC_L(:, IND), VECR_B_DEGEN(:))&
        &/DOT_PRODUCT(EIG_VEC_L(:, IND), EIG_VEC_R(:, IND))&
        &/(EIG_VAL_DEGEN - EIG_VAL(IND))
        ! write(*,*) IND, EIG_VAL(IND)

        ! sum eig_vec
        EIG_CORRECT_M = EIG_CORRECT_M + COEFF*EIG_VEC_R(:, IND)*EIG_VAL(IND) ! [PSI,DEL2CHI]

    END DO

    ! PSI = 0.D0; DEL2CHI = 0.D0; CHI = 0.D0
    ! PSI(:NR_MK) = EIG_CORRECT_M(:NR_MK)
    ! DEL2CHI(:NR_MK) = EIG_CORRECT_M(NR_MK+1:)

    ! IF (flip_switch) DEL2CHI = CONJG(DEL2CHI)
    ! CALL IDEL2_MK(DEL2CHI, CHI) !< add negative M support

    ! EIG_CORRECT_M(:NR_MK) = PSI(:NR_MK)
    ! IF (flip_switch) CHI = CONJG(CHI)
    ! EIG_CORRECT_M(NR_MK+1:) = CHI(:NR_MK)


    ! DEALLOCATE (PSI, DEL2CHI, CHI)

END FUNCTION EIG_CORRECT_M
! ======================================================================

SUBROUTINE TEST_NONLIN(L,VR,factor)
! ======================================================================
COMPLEX(P8), DIMENSION(:), ALLOCATABLE, INTENT(IN):: L,VR ! FFF
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: L2,VR2
COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE:: WW
COMPLEX(P8),OPTIONAL:: factor

ALLOCATE(L2(SIZE(L)),VR2(SIZE(VR)),WW(NRH,1))
L2 = L; VR2 = VR

CALL RTRAN_MK(L2,1)
CALL RTRAN_MK(VR2,1)

IF (PRESENT(factor)) THEN
    WW(:,1) = TFM%W(1:NRH)
    WRITE(*,*) 'PFF DOT PROD: ', DOT_PRODUCT(L2,factor*VR2), DOT_PRODUCT(L,factor*VR)
    WRITE(*,*) 'TFM^T*TFM = ', TRANSPOSE(TFM%PF(1:NRH,1:NRCHOPS(2),2))*WW*TFM%PF(1:NRH,1:NRCHOPS(2),2)
ELSE
    WRITE(*,*) 'PFF DOT PROD: ', DOT_PRODUCT(L2,VR2), DOT_PRODUCT(L,VR)
ENDIF

END SUBROUTINE TEST_NONLIN

END PROGRAM EVP_STRAIN
!=======================================================================
