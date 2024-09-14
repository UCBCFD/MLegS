PROGRAM EVP_CORRECTION
!=======================================================================
! [USAGE]:
! GIVEN A STARTING GUESS, TUNE FOR A RESONANT TRIAD AND PERFORM DEGEN PT
! CALCULATE THE CORRECTION TERMS AND PERFORM FULL DIAGNOSTICS
! CREATE PERTURBED SYSTEM AND CALCULATE EIGVAL&EIGVEC
!
! [UPDATES]:
! LAST UPDATE ON FEB 28, 2023
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
INTEGER     :: N_mat1, N_mat2, N_matS, II, JJ, KK
INTEGER     :: M_IND1, M_IND2
INTEGER     :: M_BAR, M_1, M_2, EIG_BAR_IND, tol_ind(2)
REAL(P8)    :: K_BAR, K_VAL_1_0, K_VAL_1_1, K_VAL_2_0, K_VAL_2_1
REAL(P8)    :: K_shift, tol, tol_diff, K_delta, LAMB_CUST, EIG_LMR_1
COMPLEX(P8) :: EIG_BAR, EIG_PRIME, EIG_PRIME_R, ALPHA, BETA, ALPHA_L, BETA_L
LOGICAL, DIMENSION(:), ALLOCATABLE:: MASK_1, MASK_2
REAL(P8), DIMENSION(:), ALLOCATABLE:: K_tune, K_tune_JJ, tol_diff_JJ
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig_1_0, M_eig_1_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig_2_0, M_eig_2_1
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_diff_1, M_mat_diff_2
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_1_0, EIG_R_mat_1_0, EIG_L_mat_1_0
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_2_0, EIG_R_mat_2_0, EIG_L_mat_2_0
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECR_BAR, VECL_BAR
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BAR, RUP_BAR, UZ_BAR, ROR_BAR, ROP_BAR, OZ_BAR
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BARC, RUP_BARC, UZ_BARC, ROR_BARC, ROP_BARC, OZ_BARC
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_1, VECR_1, VECR_B1, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_2, VECR_2, VECR_B2, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECR_1_CORRECT, RUR_C_1, RUP_C_1, UZ_C_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECR_2_CORRECT, RUR_C_2, RUP_C_2, UZ_C_2
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
    CALL EIG2VELVOR(M_BAR,K_BAR,EIG_BAR_IND,EIG_R_mat,EIG_L_mat,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,VECR_BAR,VECL_BAR,comm_grp=newcomm)

END IF ! MPI_GLB_RANK.EQ.0
DEALLOCATE (M_mat)
CALL MPI_BCAST(EIG_BAR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)

! ! TUNE K
! !=======================================================================
! 0. Parameters:
K_delta = 0.001
tol = 1.0E-6  ! tolerance
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
            WRITE(*,201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
            WRITE(*,201) M_1, K_VAL_1_0, JJ, M_eig_1_0(JJ)
            WRITE(*,201) M_2, K_VAL_2_0, kk, M_eig_2_0(KK)
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
! K_VAL_1_1 = real(NINT(K_VAL_1_0*100))/100
! K_VAL_1_1 = real(NINT(K_VAL_1_0*50))/50
! K_VAL_1_1 = real(NINT(K_VAL_1_0*40))/40
! K_VAL_1_1 = real(NINT(K_VAL_1_0*20))/20
! K_VAL_1_1 = real(NINT(K_VAL_1_0*10))/10
! K_VAL_1_1 = real(FLOOR(K_VAL_1_0*10))/10
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
    ! WRITE(K_VAL,'(F06.2)') K_BAR
    ! call SAVE_VEL(ROR_BAR, ROP_BAR, OZ_BAR, './converg/vor_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !             //'_IND_'//ITOA3(EIG_BAR_IND)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

    ! {M1,K1}
    ! WRITE(K_VAL,'(F06.2)') K_VAL_1_1
    CALL EIG2VELVOR(M_1, K_VAL_1_1, JJ, EIG_R_mat_1_1, EIG_L_mat_1_1, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1, VECR_1, VECL_1, comm_grp=newcomm, print_switch=DIAGNOST_SWITCH)
    ! call SAVE_VEL(ROR_1, ROP_1, OZ_1, './converg/vor_MK_'//ITOA3(M_1)//'_'//K_VAL &
    !             //'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! call SAVE_VEL(RUR_1, RUP_1, UZ_1, './converg/vel_MK_'//ITOA3(M_1)//'_'//K_VAL &
    !             //'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

    ! {M2,K2}
    ! WRITE(K_VAL,'(F06.2)') K_VAL_2_1
    CALL EIG2VELVOR(M_2, K_VAL_2_1, KK, EIG_R_mat_2_1, EIG_L_mat_2_1, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2, VECR_2, VECL_2, comm_grp=newcomm, print_switch=DIAGNOST_SWITCH)
    ! call SAVE_VEL(ROR_2, ROP_2, OZ_2, './converg/vor_MK_'//ITOA3(M_2)//'_'//K_VAL &
    !             //'_IND_'//ITOA3(KK)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! call SAVE_VEL(RUR_2, RUP_2, UZ_2, './converg/vel_MK_'//ITOA3(M_2)//'_'//K_VAL &
    !             //'_IND_'//ITOA3(KK)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    CALL NONLIN_MK(M_2,K_VAL_2_1,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,RUR_1,RUP_1,UZ_1,ROR_1,ROP_1,OZ_1,VECR_B1,comm_grp=newcomm)

    ! {M1,K1}
    ALLOCATE (RUR_BARC(NDIMR), RUP_BARC(NDIMR), UZ_BARC(NDIMR), ROR_BARC(NDIMR), ROP_BARC(NDIMR), OZ_BARC(NDIMR))
    RUR_BARC = CONJG(RUR_BAR); RUP_BARC = CONJG(RUP_BAR); UZ_BARC = CONJG(UZ_BAR)
    ROR_BARC = CONJG(ROR_BAR); ROP_BARC = CONJG(ROP_BAR); OZ_BARC = CONJG(OZ_BAR)
    CALL NONLIN_MK(M_1,K_VAL_1_1,RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC,RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_B2,comm_grp=newcomm)

ENDIF ! MPI_GLB_RANK.EQ.0

! CHECK PERTURBATION MATRIX
IF (.NOT.(ALLOCATED(RUR_BAR))) THEN
    ALLOCATE (RUR_BAR(NDIMR), RUP_BAR(NDIMR), UZ_BAR(NDIMR), ROR_BAR(NDIMR), ROP_BAR(NDIMR), OZ_BAR(NDIMR))
    ALLOCATE (RUR_BARC(NDIMR), RUP_BARC(NDIMR), UZ_BARC(NDIMR), ROR_BARC(NDIMR), ROP_BARC(NDIMR), OZ_BARC(NDIMR))
    RUR_BAR = 0.D0; RUP_BAR = 0.D0; UZ_BAR = 0.D0; ROR_BAR = 0.D0; ROP_BAR = 0.D0; OZ_BAR = 0.D0
    RUR_BARC = 0.D0; RUP_BARC = 0.D0; UZ_BARC = 0.D0; ROR_BARC = 0.D0; ROP_BARC = 0.D0; OZ_BARC = 0.D0
ENDIF
CALL MPI_ALLREDUCE(MPI_IN_PLACE, RUR_BAR, NDIMR, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE, RUP_BAR, NDIMR, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE, UZ_BAR, NDIMR, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE, ROR_BAR, NDIMR, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE, ROP_BAR, NDIMR, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE, OZ_BAR, NDIMR, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE, RUR_BARC, NDIMR, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE, RUP_BARC, NDIMR, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE, UZ_BARC, NDIMR, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE, ROR_BARC, NDIMR, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE, ROP_BARC, NDIMR, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
CALL MPI_ALLREDUCE(MPI_IN_PLACE, OZ_BARC, NDIMR, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, IERR)
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

CALL PERTURB_MATRIX2(M_BAR,K_BAR,RUR_BAR,RUP_BAR,UZ_BAR,ROR_BAR,ROP_BAR,OZ_BAR,M_1,K_VAL_1_1,P_mat_2,comm_grp=newcomm,print_switch=DIAGNOST_SWITCH)
CALL PERTURB_MATRIX2(-M_BAR,-K_BAR,RUR_BARC,RUP_BARC,UZ_BARC,ROR_BARC,ROP_BARC,OZ_BARC,M_2,K_VAL_2_1,P_mat_1,comm_grp=newcomm,print_switch=DIAGNOST_SWITCH)

RANK0: IF (MPI_GLB_RANK .EQ. 0) THEN

! =============================== DEGEN ================================
WRITE (*, *) "NEARLY DEGEN:"
EIG_PRIME_R = DOT_PRODUCT(VECL_1, VECR_B2)*DOT_PRODUCT(VECL_2, VECR_B1)
! Beta/Alpha: ALPHA * {M1,K1} + BETA * {M2,K2}
ALPHA = 1.D0; BETA = EIG_PRIME_R**0.5/DOT_PRODUCT(VECL_1, VECR_B2) 
! BETA = (DOT_PRODUCT(VECL_2,VECR_B1)/DOT_PRODUCT(VECL_1,VECR_B2))**0.5
ALPHA_L = 1.D0; BETA_L = ALPHA_L*CONJG(EIG_PRIME_R**0.5/DOT_PRODUCT(VECL_2, VECR_B1))
WRITE (*, *) " (dSigma(0))^2/4   = ", (M_eig_1_1(JJ) + EIG_BAR - M_eig_2_1(KK))**2/4
WRITE (*, *) " <L1NR2><L2NR1>    = ", EIG_PRIME_R
WRITE (*, *) " <L1NR2><L2NR1>^.5 = ", EIG_PRIME_R**0.5
WRITE (*, *) " BETA(degen)       = ", BETA
WRITE (*, *) " BETA_L(degen)     = ", BETA_L

! ============================ DIAGNOSTICS =============================
IF (DIAGNOST_SWITCH) THEN
WRITE (*, *) "DIAGNOSTICS:"
WRITE (*, *) "L1VR2: ", DOT_PRODUCT(VECL_1, VECR_B2)
WRITE (*, *) "L2VR1: ", DOT_PRODUCT(VECL_2, VECR_B1)
WRITE (*, *) "Orthogonality: (L1_new)V(R2_new):"
WRITE (*, *) "A(+)BETA(-)L1VR2: ", DOT_PRODUCT(ALPHA_L*VECL_1, (-BETA)*VECR_B2)
WRITE (*, *) "B(+)ALFA(-)L2VR1: ", DOT_PRODUCT(BETA_L*VECL_2, (ALPHA)*VECR_B1)
WRITE (*, *) "Orthogonality: (L2_new)V(R1_new):"
WRITE (*, *) "A(-)BETA(+)L1VR2: ", DOT_PRODUCT(ALPHA_L*VECL_1, (BETA)*VECR_B2)
WRITE (*, *) "B(-)ALFA(+)L2VR1: ", DOT_PRODUCT(-BETA_L*VECL_2, (ALPHA)*VECR_B1)
! WRITE(*,*) 'VECR_B2/EIG_R_1',VECR_B2(:)/EIG_R_mat_1_1(:,JJ)
! WRITE(*,*) 'VECR_B1/EIG_R_2',VECR_B2(:)/EIG_R_mat_2_1(:,KK)
WRITE (*, *) "Normalization before: ", DOT_PRODUCT(VECL_1, VECR_1)
WRITE (*, *) "Normalization before: ", DOT_PRODUCT(VECL_2, VECR_2)
ENDIF
! =========================== NEARLY DEGEN =============================
WRITE(*,*) "ENTER LAMB_CUST: "
READ(5,'(F20.10)') LAMB_CUST
WRITE(*,*) "NEARLY DEGEN   : ", ((M_eig_1_1(JJ) + EIG_BAR - M_eig_2_1(KK))**2/4 + LAMB_CUST**2 * EIG_PRIME_R)**0.5

! NEARLY DEGEN CALCULAITON
! M_NDGEN, EIG_NDGEN, VECR_NDGEN, VECL_NDGEN
M_NDGEN(1, 1) = M_eig_1_1(JJ) - K_VAL_1_1*EIG_BAR/K_BAR
M_NDGEN(1, 2) = LAMB_CUST*DOT_PRODUCT(VECL_1, VECR_B2)
M_NDGEN(2, 1) = LAMB_CUST*DOT_PRODUCT(VECL_2, VECR_B1)
M_NDGEN(2, 2) = M_eig_2_1(KK) - K_VAL_2_1*EIG_BAR/K_BAR
CALL EIGENDECOMPOSE(M_NDGEN, EIG_NDGEN, VECR_NDGEN, VECL_NDGEN)
WRITE (*, *) 'NEARLY DEGEN(1): ', EIG_NDGEN(1)
WRITE (*, *) 'NEARLY DEGEN(2): ', EIG_NDGEN(2)

! ========================= PERTURBED SYSTEM ===========================
! {m1,k1} -> N_mat1 & {m2,k2} -> N_mat2
! The perturbed system looks like below:
! +----------+----------+   +---+
! |  N_mat1  |  N_mat1  |   |   |
! |    x     |    x     |   | 1 |
! |  N_mat1  |  N_mat2  |   |   |
! +----------+----------+ X +---+
! |  N_mat2  |  N_mat2  |   |   |
! |    x     |    x     |   | 2 |
! |  N_mat1  |  N_mat2  |   |   |
! +----------+----------+   +---+
! Note that the two off-diagonal blocks are NOT square!
! ======================================================================
N_mat1 = SIZE(P_mat_1,1); N_mat2 = SIZE(P_mat_2,1); N_matS = N_mat1+N_mat2
ALLOCATE(SYS_ORG(N_matS,N_matS),SYS_BAR(N_matS,N_matS),SYS_NEW(N_matS,N_matS))
! ORIGINAL SYS:
SYS_ORG = 0.D0
SYS_ORG(1:N_mat1,1:N_mat1) = M_mat_1_1(:,:) + CMPLX(EYE(N_mat1))*(-K_VAL_1_1*EIG_BAR/K_BAR)
SYS_ORG(N_mat1+1:,N_mat1+1:) = M_mat_2_1(:,:) + CMPLX(EYE(N_mat2))*(-K_VAL_2_1*EIG_BAR/K_BAR)
! PERTURBATION:
SYS_BAR = 0.D0
SYS_BAR(1:N_mat1,N_mat1+1:) = P_mat_1(:,:)
SYS_BAR(N_mat1+1:,1:N_mat1) = P_mat_2(:,:)
! PERTURBD SYS:
SYS_NEW = SYS_ORG + LAMB_CUST * SYS_BAR
ALLOCATE(SYS_EIG(N_matS),SYS_EIG_VEC(N_matS,N_matS))
CALL EIGENDECOMPOSE(SYS_NEW, SYS_EIG, SYS_EIG_VEC)
WRITE(*,*) "PERTURBED SYS  : "
DO II = 1,SIZE(SYS_EIG)
    IF (ABS(REAL(SYS_EIG(II))) .GT. 0.05*LAMB_CUST) THEN 
        WRITE(*,*) '#',ITOA3(II),': ',SYS_EIG(II), ' - ', &
            EIGRES(SYS_EIG_VEC(1:N_mat1,II),M_1), EIGRES(SYS_EIG_VEC(N_mat1+1:,II),M_2)
        !RTOA6(REAL(SYS_EIG(II))),'+',RTOA6(AIMAG(SYS_EIG(II))),'i'
        IF (REAL(SYS_EIG(II)) .GT. 0.D0) THEN
            ! THE TUNED TRIAD
            IF (ABS(AIMAG(SYS_EIG(II))-AIMAG(M_eig_1_0(JJ)-K_VAL_1_1*EIG_BAR/K_BAR)).LT.0.01) THEN

                CALL PERTURB_SAVE(SYS_EIG_VEC(:,II),M_1, K_VAL_1_1, M_2, K_VAL_2_1, newcomm, print_switch=.true.)
                CALL SAVE_PERTURB('./perturb/per_perturb.input', &
                    M_1, K_VAL_1_1, SYS_EIG_VEC(1:N_mat1,II), SYS_EIG(II), &
                    M_2, K_VAL_2_1, SYS_EIG_VEC(N_mat1+1:,II), SYS_EIG(II), &
                    M_BAR, K_BAR, EIG_R_mat(:, EIG_BAR_IND), EIG_BAR, LAMB_CUST)
                
                DEALLOCATE(RUR_1, RUP_1, UZ_1, RUR_2, RUP_2, UZ_2)
                CALL EIG2VEL(M_1, K_VAL_1_1, SYS_EIG_VEC(1:N_mat1,II), RUR_1, RUP_1, UZ_1, newcomm, .true.)
                CALL EIG2VEL(M_2, K_VAL_2_1, SYS_EIG_VEC(N_mat1+1:,II), RUR_2, RUP_2, UZ_2, newcomm, .true.)

            ! OTHER POSTIVELY GROWING TRIAD
            ELSEIF (EIGRES(SYS_EIG_VEC(1:N_mat1,II),M_1).AND.EIGRES(SYS_EIG_VEC(N_mat1+1:,II),M_2)) THEN
                CALL PERTURB_SAVE(SYS_EIG_VEC(:,II),M_1, K_VAL_1_1, M_2, K_VAL_2_1, newcomm, II)
            ENDIF
        ENDIF
    ENDIF
ENDDO
DEALLOCATE(SYS_ORG,SYS_BAR,SYS_NEW,SYS_EIG,SYS_EIG_VEC)
! ======================================================================

! ! NORMALIZE ALL EIGVEC EXCEPT #EIG_IND
! CALL EIG2NORMED(M_1, K_VAL_1_1, JJ, EIG_R_mat_1_1, EIG_L_mat_1_1, comm_grp=newcomm)
! CALL EIG2NORMED(M_2, K_VAL_2_1, KK, EIG_R_mat_2_1, EIG_L_mat_2_1, comm_grp=newcomm)

! CALCULATE CORRECTION TERMS
! VECR_1_CORRECT = EIG_CORRECT(JJ,M_eig_1_1,EIG_R_mat_1_1,EIG_L_mat_1_1,VECR_B2)
! VECR_2_CORRECT = EIG_CORRECT(KK,M_eig_2_1,EIG_R_mat_2_1,EIG_L_mat_2_1,VECR_B1)
! EIG_CORRECT_2(EIG_IND_DEGEN,EIG_VAL,EIG_VEC_R,EIG_VEC_L,VECR_B_DEGEN,MREAD,KREAD,factor_self,factor_other,comm_grp)
VECR_1_CORRECT = EIG_CORRECT_3(JJ, M_eig_1_1, EIG_R_mat_1_1, EIG_L_mat_1_1, VECR_B2, M_1, K_VAL_1_1, ALPHA, BETA, newcomm, DIAGNOST_SWITCH)
VECR_2_CORRECT = EIG_CORRECT_3(KK, M_eig_2_1, EIG_R_mat_2_1, EIG_L_mat_2_1, VECR_B1, M_2, K_VAL_2_1, BETA, ALPHA, newcomm, DIAGNOST_SWITCH)
VECR_1_CORRECT = VECR_1_CORRECT*BETA
VECR_2_CORRECT = VECR_2_CORRECT*ALPHA

! FORM NEW EIGEN-VECTORS:
EIG_R_mat_1_1(:, JJ) = EIG_R_mat_1_1(:, JJ)*ALPHA
EIG_R_mat_2_1(:, KK) = EIG_R_mat_2_1(:, KK)*BETA
VECR_1 = VECR_1*ALPHA  ; VECR_2 = VECR_2*BETA
VECL_1 = VECL_1*ALPHA_L; VECL_2 = VECL_2*BETA_L

! ============================ DIAGNOSTICS =============================
IF (DIAGNOST_SWITCH) THEN
WRITE (*, *) "Normalization after: ", DOT_PRODUCT(VECL_1, VECR_1)
WRITE (*, *) "Normalization after: ", DOT_PRODUCT(VECL_2, VECR_2)
! WRITE(*,*) 'DIFFERENCE - BetaVECR_B2 vs sigEIG_R_1',sum(abs(BETA*VECR_B2(:) - EIG_PRIME_R**0.5*EIG_R_mat_1_1(:,JJ)))/(2*SIZE(VECR_B2))
! WRITE(*,*) 'DIFFERENCE - BetaVECR_B1 vs sigEIG_R_2',sum(abs(ALPHA*VECR_B1(:) - EIG_PRIME_R**0.5*EIG_R_mat_2_1(:,KK)))/(2*SIZE(VECR_B1))
WRITE(*,*) 'DIFFERENCE - BetaVECR_B2 vs sigEIG_R_1',maxval(abs(BETA*VECR_B2(:) - EIG_PRIME_R**0.5*VECR_1))
WRITE(*,*) 'DIFFERENCE - BetaVECR_B1 vs sigEIG_R_2',maxval(abs(ALPHA*VECR_B1(:) - EIG_PRIME_R**0.5*VECR_2))
ENDIF
! ======================================================================

! CHECK IF EIG_VEC_R & VEC_R ARE CONSISTENT W/ EACH OTHER:
! CALL EIG2PT(M_1, K_VAL_1_1, JJ, EIG_R_mat_1_1, VECR_1, comm_grp=newcomm)
! CALL EIG2PT(M_2, K_VAL_2_1, KK, EIG_R_mat_2_1, VECR_2, comm_grp=newcomm)
WRITE(K_VAL,'(F06.2)') K_BAR 
CALL SAVE_PSICHI(VECR_BAR,'./perturb/pt_MK_'//ITOA3(M_BAR)//'_'//K_VAL//'_O.output')

WRITE(K_VAL,'(F06.2)') K_VAL_1_1 
CALL SAVE_PSICHI(VECR_1,'./perturb/pt_MK_'//ITOA3(M_1)//'_'//K_VAL//'_O.output')
! CALL SAVE_PSICHI(EIG_R_mat_1_1(:, JJ),'./perturb/pt_MK_'//ITOA3(M_1)//'_'//K_VAL//'_O.output',FFF_SWITCH=.TRUE.)
WRITE(K_VAL,'(F06.2)') K_VAL_2_1
CALL SAVE_PSICHI(VECR_2,'./perturb/pt_MK_'//ITOA3(M_2)//'_'//K_VAL//'_O.output')
! CALL SAVE_PSICHI(EIG_R_mat_2_1(:, KK),'./perturb/pt_MK_'//ITOA3(M_2)//'_'//K_VAL//'_O.output',FFF_SWITCH=.TRUE.)

! ! SAVE DEGEN PREDICTION RESULT
! ! ======================================================================
! IF (abs(REAL(EIG_PRIME_R**0.5)) .GT. 0.01) THEN
!     open (unit=FID, FILE='./converg/perturb_scan.GROW', STATUS='unknown', ACTION='WRITE', ACCESS='APPEND', IOSTAT=IS)
!     WRITE (FID, *) ''
!     WRITE (FID, *) 'NEARLY DEGEN: *GROWTH*'
! ELSE
!     open (unit=FID, FILE='./converg/perturb_scan.NEGLIGIBLE', STATUS='unknown', ACTION='WRITE', ACCESS='APPEND', IOSTAT=IS)
!     WRITE (FID, *) ''
!     WRITE (FID, *) 'NEARLY DEGEN: *NEGLIGIBLE*'
! END IF
! WRITE (FID, 201) M_BAR, K_BAR, EIG_BAR_IND, EIG_BAR
! WRITE (FID, 201) M_1, K_VAL_1_1, JJ, M_eig_1_1(JJ)
! WRITE (FID, 201) M_2, K_VAL_2_1, kk, M_eig_2_1(KK)
! WRITE (FID, *) "SIGMA(degen)'= RE:", abs(real(EIG_PRIME_R**0.5)), ',IM:', imag(EIG_PRIME_R**0.5)
! WRITE (FID, *) "BETA/ALPHA   = RE:", abs(real(BETA)), ',IM:', imag(BETA)
! WRITE (FID, *) ''
! close (FID)

! SAVE COLLOCATION INFO
! ======================================================================
CALL COLLOC_INFO()

! SAVE PERTURBATION FILES
! ======================================================================
CALL SAVE_PERTURB('./converg/new_perturb.input', M_1, K_VAL_1_1, EIG_R_mat_1_1(:, JJ), M_eig_1_1(JJ), &
                    M_2, K_VAL_2_1, EIG_R_mat_2_1(:, KK), M_eig_2_1(kk), &
                    M_BAR, K_BAR, EIG_R_mat(:, EIG_BAR_IND), EIG_BAR)
! CALL SAVE_PERTURB('./converg/lft_perturb.input', M_1, K_VAL_1_1, VECL_1, M_eig_1_1(JJ), &
!                     M_2, K_VAL_2_1, VECL_2, M_eig_2_1(kk)) ! SAVE NEW LEFT EIGENVECTORS
! CALL SAVE_PERTURB('./converg/cor_perturb.input', M_1, K_VAL_1_1, VECR_1_CORRECT, M_eig_1_1(JJ), &
!                     M_2, K_VAL_2_1, VECR_2_CORRECT, M_eig_2_1(kk)) ! SAVE CORRECTION TERMS

! ANALYSIS
! ======================================================================
! 1. SAVE UR, UP, UZ & OR, OP, OZ
! call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './converg/vel_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
!                 //'_IND_'//ITOA3(EIG_BAR_IND)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
! call SAVE_VEL(ROR_BAR, ROP_BAR, OZ_BAR, './converg/vor_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
!                 //'_IND_'//ITOA3(EIG_BAR_IND)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
! WRITE(K_VAL,'(F06.2)') K_VAL_1_1              
! call SAVE_VEL(RUR_1, RUP_1, UZ_1, './converg/vel_MK_'//ITOA3(M_1)//'_'//K_VAL &
!                 //'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
! WRITE(K_VAL,'(F06.2)') K_VAL_2_1
! call SAVE_VEL(RUR_2, RUP_2, UZ_2, './converg/vel_MK_'//ITOA3(M_2)//'_'//K_VAL &
!                 //'_IND_'//ITOA3(KK)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
WRITE(K_VAL,'(F06.2)') K_BAR              
call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './perturb/vel_MK_'//ITOA3(M_BAR)//'_'//K_VAL//'_O.output')
call SAVE_VEL(ROR_BAR, ROP_BAR, OZ_BAR, './perturb/vor_MK_'//ITOA3(M_BAR)//'_'//K_VAL//'_O.output')
WRITE(K_VAL,'(F06.2)') K_VAL_1_1              
call SAVE_VEL(RUR_1*ALPHA, RUP_1*ALPHA, UZ_1*ALPHA, './perturb/vel_MK_'//ITOA3(M_1)//'_'//K_VAL//'_O.output')
call SAVE_VEL(ROR_1*ALPHA, ROP_1*ALPHA, OZ_1*ALPHA, './perturb/vor_MK_'//ITOA3(M_1)//'_'//K_VAL//'_O.output')
WRITE(K_VAL,'(F06.2)') K_VAL_2_1
call SAVE_VEL(RUR_2*BETA , RUP_2*BETA , UZ_2*BETA , './perturb/vel_MK_'//ITOA3(M_2)//'_'//K_VAL//'_O.output')
call SAVE_VEL(ROR_2*BETA , ROP_2*BETA , OZ_2*BETA , './perturb/vor_MK_'//ITOA3(M_2)//'_'//K_VAL//'_O.output')

! 2. SAVE UR, UP, UZ OF THE CORRECTION TERMS
! VECR_1_CORRECT, VECR_2_CORRECT
CALL EIG2VEL(M_1, K_VAL_1_1, VECR_1_CORRECT, RUR_C_1, RUP_C_1, UZ_C_1, comm_grp=newcomm,chi_switch=.true.)
CALL EIG2VEL(M_2, K_VAL_2_1, VECR_2_CORRECT, RUR_C_2, RUP_C_2, UZ_C_2, comm_grp=newcomm,chi_switch=.true.)
WRITE(K_VAL,'(F06.2)') K_VAL_1_1 
call SAVE_VEL(RUR_C_1, RUP_C_1, UZ_C_1, './perturb/vel_MK_'//ITOA3(M_1)//'_'//K_VAL//'_C.output')
CALL SAVE_PSICHI(VECR_1_CORRECT,'./perturb/pt_MK_'//ITOA3(M_1)//'_'//K_VAL//'_C.output')
! CALL SAVE_PSICHI(VECR_1_CORRECT,'./converg/COR_pt_MK_'//ITOA3(M_1)//'_'//K_VAL &
!         //'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
WRITE(K_VAL,'(F06.2)') K_VAL_2_1
call SAVE_VEL(RUR_C_2, RUP_C_2, UZ_C_2, './perturb/vel_MK_'//ITOA3(M_2)//'_'//K_VAL//'_C.output')
CALL SAVE_PSICHI(VECR_2_CORRECT,'./perturb/pt_MK_'//ITOA3(M_2)//'_'//K_VAL//'_C.output')
! CALL SAVE_PSICHI(VECR_2_CORRECT,'./converg/COR_pt_MK_'//ITOA3(M_2)//'_'//K_VAL &
!         //'_IND_'//ITOA3(KK)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

! 3. SAVE UR0, UP0, UZ0 & OR0, OP0, OZ0
! call SAVE_VEL(CMPLX(RUR0), CMPLX(RUP0), CMPLX(UZ0), './converg/vel_MK_'//ITOA3(0)//'_'//ITOA4(0) &
!                 //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
! call SAVE_VEL(CMPLX(ROR0), CMPLX(ROP0), CMPLX(OZ0), './converg/vor_MK_'//ITOA3(0)//'_'//ITOA4(0) &
!                 //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

! 4. SAVE UR, UP, UZ OF THE EIGVEC OF THE PERTURBED EVP
call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './perturb/vel_MK_0.output')
call SAVE_VEL(ROR_BAR, ROP_BAR, OZ_BAR, './perturb/vor_MK_0.output')
call SAVE_VEL(RUR_1, RUP_1, UZ_1, './perturb/vel_MK_1.output')
call SAVE_VEL(RUR_2, RUP_2, UZ_2, './perturb/vel_MK_2.output')

! DEALLOCATE
!=======================================================================
1   DEALLOCATE (RUR_BAR, RUP_BAR, UZ_BAR, ROR_BAR, ROP_BAR, OZ_BAR)
    DEALLOCATE (RUR_BARC, RUP_BARC, UZ_BARC, ROR_BARC, ROP_BARC, OZ_BARC)
    DEALLOCATE (RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1)
    DEALLOCATE (RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2)
    DEALLOCATE (VECR_B1, VECR_B2, VECR_1, VECR_2, VECL_1, VECL_2, VECL_BAR)
    DEALLOCATE (EIG_R_mat_1_1, EIG_L_mat_1_1, EIG_R_mat_2_1, EIG_L_mat_2_1)
    DEALLOCATE (VECR_1_CORRECT, VECR_2_CORRECT)
    DEALLOCATE (RUR_C_1, RUP_C_1, UZ_C_1)
    DEALLOCATE (RUR_C_2, RUP_C_2, UZ_C_2)

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

FUNCTION EIG_CORRECT(EIG_IND_DEGEN, EIG_VAL, EIG_VEC_R, EIG_VEC_L, VECR_B_DEGEN)
! ======================================================================
! NOTE: VECR_B_DEGEN AND EIG_VEC_R CORRESPOND TO THE TWO WAVENUMBERS
!       e.g. VECR_B_DEGEN = N(V_BAR)V2, EIG_VEC_R/L & EIG_IND_DEGEN & EIG_VAL ~ V1
    IMPLICIT NONE
! INTEGER,INTENT(IN)    :: EIG_IND_DEGEN
! COMPLEX(P8),DIMENSION(:),ALLOCATABLE,INTENT(IN):: EIG_VAL, VECR_B_DEGEN
! COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: EIG_VEC_R, EIG_VEC_L
    INTEGER:: EIG_IND_DEGEN
    COMPLEX(P8), DIMENSION(:):: EIG_VAL, VECR_B_DEGEN
    COMPLEX(P8), DIMENSION(:, :):: EIG_VEC_R, EIG_VEC_L
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: EIG_CORRECT

    INTEGER:: IND
    COMPLEX(P8):: EIG_VAL_DEGEN, COEFF
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE::CORRECTION

    EIG_VAL_DEGEN = EIG_VAL(EIG_IND_DEGEN)
    ALLOCATE (EIG_CORRECT(SIZE(EIG_VEC_R, 1)))
    ALLOCATE (CORRECTION(SIZE(EIG_VEC_R, 1)))
    EIG_CORRECT = 0.D0

    DO IND = 1, SIZE(EIG_VEC_R, 2)
        IF (IND .EQ. EIG_IND_DEGEN) CYCLE
        COEFF = DOT_PRODUCT(EIG_VEC_L(:, IND), VECR_B_DEGEN(:))&
        &/DOT_PRODUCT(EIG_VEC_L(:, IND), EIG_VEC_R(:, IND))&
        &/(EIG_VAL_DEGEN - EIG_VAL(IND))
        CORRECTION = COEFF*EIG_VEC_R(:, IND)
        EIG_CORRECT = EIG_CORRECT + CORRECTION
        IF (MAXVAL(ABS(COEFF*EIG_VEC_R(:, IND))) .GE. 1E-1) THEN
            WRITE (*, *) IND, ',', MAXVAL(ABS(COEFF*EIG_VEC_R(:, IND))), ',', 1/(EIG_VAL_DEGEN - EIG_VAL(IND))
        END IF
    END DO

END FUNCTION EIG_CORRECT
! ======================================================================

FUNCTION EIG_CORRECT_2(EIG_IND_DEGEN,EIG_VAL,EIG_VEC_R,EIG_VEC_L,VECR_B_DEGEN,MREAD,AKREAD,factor_self,factor_other,comm_grp,save_switch_in)
! ======================================================================
! NOTE: VECR_B_DEGEN AND EIG_VEC_R CORRESPOND TO THE TWO WAVENUMBERS
!       e.g. VECR_B_DEGEN = N(V_BAR)V2, EIG_VEC_R/L & EIG_IND_DEGEN & EIG_VAL ~ V1
! ======================================================================
    IMPLICIT NONE
    ! INTEGER,INTENT(IN)    :: EIG_IND_DEGEN
    ! COMPLEX(P8),DIMENSION(:),ALLOCATABLE,INTENT(IN):: EIG_VAL, VECR_B_DEGEN
    ! COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: EIG_VEC_R, EIG_VEC_L
    INTEGER:: EIG_IND_DEGEN, MREAD, comm_grp
    REAL(P8):: AKREAD
    COMPLEX(P8), DIMENSION(:):: EIG_VAL, VECR_B_DEGEN
    COMPLEX(P8), DIMENSION(:, :):: EIG_VEC_R, EIG_VEC_L
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: EIG_CORRECT_2
    LOGICAL, OPTIONAL:: save_switch_in
    LOGICAL:: save_switch

    INTEGER:: IND, EIG_NUM, NR_MK
    COMPLEX(P8):: EIG_VAL_DEGEN, COEFF, factor_self, factor_other
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: CORRECTION, EIG_LEFTOVER
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_new, RUP_new, UZ_new
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_old, RUP_old, UZ_old
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: PSI, DEL2CHI, CHI

    IF ((PRESENT(save_switch_in)) .AND. (save_switch_in)) THEN
        save_switch = .TRUE.
    ELSE
        save_switch = .FALSE.
    END IF

    IF ((M(2) .NE. MREAD) .OR. (ABS(ZLEN - 2*PI/AKREAD) .GE. 1.0E-13)) THEN
        WRITE (*, *) 'EIG_CORRECT_2: SET NEW M AND AK'
        ZLEN = 2*PI/AKREAD
        NTH = 2
        NX = 4
        NTCHOP = 2 ! M = 0, M = MREAD
        NXCHOP = 2 ! K = 0, K = AKREAD
        CALL LEGINIT(comm_grp, MREAD)
    END IF

    EIG_VAL_DEGEN = EIG_VAL(EIG_IND_DEGEN)
    ALLOCATE (EIG_CORRECT_2(SIZE(EIG_VEC_R, 1)))
    ALLOCATE (CORRECTION(SIZE(EIG_VEC_R, 1)))
    ALLOCATE (EIG_LEFTOVER(SIZE(EIG_VEC_R, 1)))
    ALLOCATE (RUR_new(NRCHOPDIM), RUP_new(NRCHOPDIM), UZ_new(NRCHOPDIM))
    ALLOCATE (RUR_old(NRCHOPDIM), RUP_old(NRCHOPDIM), UZ_old(NRCHOPDIM))
    ALLOCATE (PSI(NRCHOPDIM), DEL2CHI(NRCHOPDIM), CHI(NRCHOPDIM))

    EIG_CORRECT_2 = 0.D0
    EIG_LEFTOVER = 0.D0
    EIG_NUM = SIZE(EIG_VEC_R, 1)
    NR_MK = SIZE(EIG_VEC_R, 2)/2

    IF (save_switch) THEN
        open (10, FILE='./converg/check/contribution_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
                //'.output', STATUS='unknown', ACTION='WRITE', IOSTAT=IS)
        if (IS .ne. 0) then
            print *, 'ERROR: createperturb -- Could not creat new file contribution.output'
            RETURN
        end if
    END IF

    ! ============================ DEGEN MODE ==============================
    ! save original
    RUR_old = 0.D0
    RUP_old = 0.D0
    UZ_old = 0.D0
    PSI = 0.D0
    DEL2CHI = 0.D0
    CHI = 0.D0
    PSI(:NR_MK) = EIG_VEC_R(:NR_MK, EIG_IND_DEGEN) !PSI
    CHI(:NR_MK) = EIG_VEC_R(NR_MK + 1:, EIG_IND_DEGEN) !CHI
    CALL CHOPSET(3)
    CALL PC2VEL_MK(PSI, CHI, RUR_old, RUP_old, UZ_old)
    CALL RTRAN_MK(RUR_old, 1)
    CALL RTRAN_MK(RUP_old, 1)
    CALL RTRAN_MK(UZ_old, 1)
    CALL CHOPSET(-3)

    ! save new
    RUR_new = 0.D0
    RUP_new = 0.D0
    UZ_new = 0.D0
    PSI = 0.D0
    DEL2CHI = 0.D0
    CHI = 0.D0
    PSI(:NR_MK) = factor_self*EIG_VEC_R(:NR_MK, EIG_IND_DEGEN) !PSI
    CHI(:NR_MK) = factor_self*EIG_VEC_R(NR_MK + 1:, EIG_IND_DEGEN) !DEL2CHI
    CALL DEL2_MK(CHI, DEL2CHI)
    CALL CHOPSET(3)
    CALL PC2VEL_MK(PSI, CHI, RUR_new, RUP_new, UZ_new)
    CALL RTRAN_MK(RUR_new, 1)
    CALL RTRAN_MK(RUP_new, 1)
    CALL RTRAN_MK(UZ_new, 1)
    CALL CHOPSET(-3)
    IF (save_switch) THEN
        WRITE (10, 1345) EIG_IND_DEGEN, ABS(EIG_VAL_DEGEN), ABS(factor_self),&
        &MAXVAL(ABS(RUP_old/TFM%R(1:NDIMR))),&
        &MAXVAL(ABS(RUP_new/TFM%R(1:NDIMR))),&
        &MAX(MAXVAL(ABS(PSI(:NR_MK)/factor_self)), MAXVAL(ABS(DEL2CHI(:NR_MK)/factor_self))),&
        &MAX(MAXVAL(ABS(PSI(:NR_MK))), MAXVAL(ABS(DEL2CHI(:NR_MK)))),&
        &MAXVAL(ABS(RUP_old)),&
        &MAXVAL(ABS(RUP_new))
    END IF
    1345    FORMAT(I4, 8(',', G20.12))
    ! WRITE(10,*) IND,',',ABS(EIG_VAL_DEGEN-EIG_VAL(IND)),',',ABS(COEFF),&
    ! &',',MAXVAL(ABS(RUP/TFM%R)),&
    ! &',',MAXVAL(ABS(EIG_VEC_R(:,IND))),&
    ! &',',MAXVAL(ABS(COEFF*EIG_VEC_R(:,IND))),&
    ! &',',MAXVAL(ABS(RUP))
    ! ======================================================================

    DO IND = 1, SIZE(EIG_VEC_R, 2)
        IF (IND .EQ. EIG_IND_DEGEN) CYCLE
        IF (SUM(ABS(EIG_VEC_L(:, IND))) .GT. 10E50) CYCLE ! TEMPORARY - FOR EIGVAL = 0 MODES

        ! calculate coefficient
        COEFF = DOT_PRODUCT(EIG_VEC_L(:, IND), VECR_B_DEGEN(:))&
        &/DOT_PRODUCT(EIG_VEC_L(:, IND), EIG_VEC_R(:, IND))&
        &/(EIG_VAL_DEGEN - EIG_VAL(IND))

        ! debug
        ! write(*,*) IND,DOT_PRODUCT(EIG_VEC_L(:, IND), EIG_VEC_R(:, IND)) !,EIG_VAL_DEGEN - EIG_VAL(IND)
        write(*,*) IND,SUM(ABS(EIG_VEC_R(2:, IND))),SUM(ABS(EIG_VEC_L(2:, IND)))

        ! save original
        RUR_old = 0.D0
        RUP_old = 0.D0
        UZ_old = 0.D0
        PSI = 0.D0
        DEL2CHI = 0.D0
        CHI = 0.D0
        PSI(:NR_MK) = EIG_VEC_R(:NR_MK, IND) !PSI
        DEL2CHI(:NR_MK) = EIG_VEC_R(NR_MK + 1:, IND) !DEL2CHI
        CALL IDEL2_MK(DEL2CHI, CHI)
        CALL CHOPSET(3)
        CALL PC2VEL_MK(PSI, CHI, RUR_old, RUP_old, UZ_old)
        CALL RTRAN_MK(RUR_old, 1)
        CALL RTRAN_MK(RUP_old, 1)
        CALL RTRAN_MK(UZ_old, 1)
        CALL CHOPSET(-3)

        IF ((IND .LE. 61).AND.(save_switch)) THEN
            CALL SAVE_VEL(RUR_old, RUP_old, UZ_old, './converg/check/before_vel_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
                            //'_IND_'//ITOA3(IND)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
        END IF

        ! times coefficient
        CORRECTION(:NR_MK) = COEFF*PSI(:NR_MK)
        CORRECTION(NR_MK + 1:) = COEFF*CHI(:NR_MK)

        ! save new
        RUR_new = 0.D0; RUP_new = 0.D0; UZ_new = 0.D0
        PSI = 0.D0; DEL2CHI = 0.D0; CHI = 0.D0
        PSI(:NR_MK) = factor_other*COEFF*EIG_VEC_R(:NR_MK, IND) !PSI
        DEL2CHI(:NR_MK) = factor_other*COEFF*EIG_VEC_R(NR_MK + 1:, IND) !DEL2CHI
        CALL IDEL2_MK(DEL2CHI, CHI)
        CALL CHOPSET(3)
        CALL PC2VEL_MK(PSI, CHI, RUR_new, RUP_new, UZ_new)
        CALL RTRAN_MK(RUR_new, 1); CALL RTRAN_MK(RUP_new, 1); CALL RTRAN_MK(UZ_new, 1)
        CALL CHOPSET(-3)
        IF (IND .LE. 61) THEN
            IF (save_switch) THEN
                CALL SAVE_VEL(RUR_new, RUP_new, UZ_new, './converg/check/after_vel_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
                                //'_IND_'//ITOA3(IND)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
            END IF
        ELSE
            EIG_LEFTOVER = EIG_LEFTOVER + factor_other*CORRECTION
        END IF

        ! sum eig_vec
        EIG_CORRECT_2 = EIG_CORRECT_2 + CORRECTION

        ! IF (MAXVAL(ABS(COEFF*EIG_VEC_R(:,IND))) .GE. 1E-1) THEN
        IF (save_switch) THEN
            WRITE (10, 1345) IND, ABS(EIG_VAL(IND)), ABS(factor_other*COEFF),&
            &MAXVAL(ABS(RUP_old/TFM%R(1:NDIMR))),&
            &MAXVAL(ABS(RUP_new/TFM%R(1:NDIMR))),&
            &MAXVAL(ABS(EIG_VEC_R(:, IND))),&
            &MAXVAL(ABS(factor_other*COEFF*EIG_VEC_R(:, IND))),&
            &MAXVAL(ABS(RUP_old)),&
            &MAXVAL(ABS(RUP_new))
        END IF
        ! ENDIF
    END DO
    IF (save_switch) THEN
        close (10)

        RUR_new = 0.D0; RUP_new = 0.D0; UZ_new = 0.D0
        PSI = 0.D0; DEL2CHI = 0.D0; CHI = 0.D0
        PSI(:NR_MK) = EIG_LEFTOVER(:NR_MK) !PSI
        DEL2CHI(:NR_MK) = EIG_LEFTOVER(NR_MK + 1:) !DEL2CHI
        CALL IDEL2_MK(DEL2CHI, CHI)
        CALL CHOPSET(3)
        CALL PC2VEL_MK(PSI, CHI, RUR_new, RUP_new, UZ_new)
        CALL RTRAN_MK(RUR_new, 1); CALL RTRAN_MK(RUP_new, 1); CALL RTRAN_MK(UZ_new, 1)
        CALL CHOPSET(-3)
    ! IF (save_switch) THEN
        CALL SAVE_VEL(RUR_new, RUP_new, UZ_new, './converg/check/lft_vel_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
                        //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    END IF

    DEALLOCATE (PSI, DEL2CHI, CHI, RUR_old, RUP_old, UZ_old, RUR_new, RUP_new, UZ_new, CORRECTION, EIG_LEFTOVER)

END FUNCTION EIG_CORRECT_2
! ======================================================================

FUNCTION EIG_CORRECT_3(EIG_IND_DEGEN,EIG_VAL,EIG_VEC_R,EIG_VEC_L,VECR_B_DEGEN,MREAD,AKREAD,factor_self,factor_other,comm_grp,print_switch)
! ======================================================================
! NOTE: VECR_B_DEGEN AND EIG_VEC_R CORRESPOND TO THE TWO WAVENUMBERS
!       e.g. VECR_B_DEGEN = N(V_BAR)V2, EIG_VEC_R/L & EIG_IND_DEGEN & EIG_VAL ~ V1
! ======================================================================
    IMPLICIT NONE
    INTEGER:: EIG_IND_DEGEN, MREAD, comm_grp
    REAL(P8):: AKREAD,AK_ACTUAL
    COMPLEX(P8), DIMENSION(:):: EIG_VAL, VECR_B_DEGEN
    COMPLEX(P8), DIMENSION(:, :):: EIG_VEC_R, EIG_VEC_L
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: EIG_CORRECT_3
    LOGICAL, OPTIONAL:: print_switch
    LOGICAL:: flip_switch

    INTEGER:: IND, EIG_NUM, NR_MK, M_ACTUAL
    COMPLEX(P8):: EIG_VAL_DEGEN, COEFF, factor_self, factor_other
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: CORRECTION, EIG_LEFTOVER
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_new, RUP_new, UZ_new
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_old, RUP_old, UZ_old
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: PSI, DEL2CHI, CHI

    IF (MREAD.LT.0) THEN
        flip_switch = .TRUE.
        M_ACTUAL = -MREAD; AK_ACTUAL = -AKREAD
    ELSE
        flip_switch = .FALSE.
        M_ACTUAL = MREAD; AK_ACTUAL = AKREAD
    ENDIF
    IF ((M(2) .NE. M_ACTUAL) .OR. (ABS(ZLEN - 2*PI/AK_ACTUAL) .GE. 1.0E-13)) THEN
        IF ((PRESENT(print_switch)).AND.(print_switch)) WRITE (*, *) 'EIG_CORRECT_3: SET NEW M AND AK'
        ZLEN = 2*PI/AK_ACTUAL
        NTH = 2
        NX = 4
        NTCHOP = 2 ! M = 0, M = MREAD
        NXCHOP = 2 ! K = 0, K = AKREAD
        CALL LEGINIT(comm_grp, M_ACTUAL)
    ENDIF

    EIG_VAL_DEGEN = EIG_VAL(EIG_IND_DEGEN)
    ALLOCATE (EIG_CORRECT_3(SIZE(EIG_VEC_R, 1)))
    ALLOCATE (CORRECTION(SIZE(EIG_VEC_R, 1)))
    ALLOCATE (EIG_LEFTOVER(SIZE(EIG_VEC_R, 1)))

    EIG_CORRECT_3 = 0.D0
    EIG_LEFTOVER = 0.D0
    EIG_NUM = SIZE(EIG_VEC_R, 1)
    NR_MK = SIZE(EIG_VEC_R, 2)/2

    DO IND = 1, SIZE(EIG_VEC_R, 2)
        IF (IND .EQ. EIG_IND_DEGEN) CYCLE
        ! IF (SUM(ABS(EIG_VEC_L(:, IND))) .GT. 10E50) CYCLE ! TEMPORARY - FOR EIGVAL = 0 MODES
        IF (EIG_VAL(IND).EQ.0.D0) THEN
            IF (ISZERO(EIG_VEC_R,IND,NR_MK,flip_switch)) THEN
                WRITE(*,*) 'ISZERO:',IND; CYCLE
            ENDIF
        ENDIF

        ! calculate coefficient
        COEFF = DOT_PRODUCT(EIG_VEC_L(:, IND), VECR_B_DEGEN(:))&
        &/DOT_PRODUCT(EIG_VEC_L(:, IND), EIG_VEC_R(:, IND))&
        &/(EIG_VAL_DEGEN - EIG_VAL(IND))

        ! debug
        IF (ABS(COEFF).GT.1E2) THEN
            ! write(*,*) IND,DOT_PRODUCT(EIG_VEC_L(:, IND), EIG_VEC_R(:, IND))
            write(*,*) IND,DOT_PRODUCT(EIG_VEC_L(:, IND), VECR_B_DEGEN(:))
            write(*,*) IND,EIG_VAL_DEGEN - EIG_VAL(IND)
            ! write(*,*) IND,SUM(ABS(EIG_VEC_R(2:, IND))),SUM(ABS(EIG_VEC_L(2:, IND)))
        ENDIF

        ! obtain correction term
        CORRECTION = 0.D0
        CORRECTION = EIG2PTVEL(EIG_VEC_R,IND,NR_MK,COEFF,factor_other,flip_switch,print_switch,MREAD,AKREAD,EIG_VAL(IND))
        EIG_CORRECT_3 = EIG_CORRECT_3 + CORRECTION

        ! save left-over
        IF (IND .GT. 61) THEN
            EIG_LEFTOVER = EIG_LEFTOVER + factor_other*CORRECTION
        END IF     

    END DO

    ! save high-order terms
    IF ((PRESENT(print_switch)).AND.(print_switch)) THEN
        ALLOCATE (RUR_new(NRCHOPDIM), RUP_new(NRCHOPDIM), UZ_new(NRCHOPDIM))
        ALLOCATE (PSI(NRCHOPDIM), DEL2CHI(NRCHOPDIM), CHI(NRCHOPDIM))
        RUR_new = 0.D0; RUP_new = 0.D0; UZ_new = 0.D0
        PSI = 0.D0; DEL2CHI = 0.D0; CHI = 0.D0
        PSI(:NR_MK) = EIG_LEFTOVER(:NR_MK) !PSI
        DEL2CHI(:NR_MK) = EIG_LEFTOVER(NR_MK + 1:) !DEL2CHI
        IF (flip_switch) THEN
            PSI = CONJG(PSI); DEL2CHI = CONJG(DEL2CHI)
        ENDIF
        CALL IDEL2_MK(DEL2CHI, CHI)
        CALL CHOPSET(3)
        CALL PC2VEL_MK(PSI, CHI, RUR_new, RUP_new, UZ_new)
        CALL RTRAN_MK(RUR_new, 1); CALL RTRAN_MK(RUP_new, 1); CALL RTRAN_MK(UZ_new, 1)
        CALL CHOPSET(-3)
        IF (flip_switch) THEN
            RUR_new = CONJG(RUR_new); RUP_new = CONJG(RUP_new); UZ_new = CONJG(UZ_new)
        ENDIF
        CALL SAVE_VEL(RUR_new, RUP_new, UZ_new, './converg/check/lft_vel_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
                        //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
        DEALLOCATE(PSI, DEL2CHI, CHI, RUR_new, RUP_new, UZ_new)
    END IF

    DEALLOCATE(CORRECTION, EIG_LEFTOVER)

END FUNCTION EIG_CORRECT_3
! ======================================================================

FUNCTION EIG2PTVEL(EIG_VEC_R,IND,NR_MK,coeff,factor,flip_switch,save_switch,MREAD,AKREAD,eig_val)
! ======================================================================
! SUBROUTINE FOR EIG_CORRECT_3
! CONVERT EIG_VEC_R(:,IND) TO PSI,CHI
! IF SAVE_SWITH == TRUE:
!   SAVE THE CORRESPONDING VELOCITY FIELDS
!   IF FACTOR IS PRESENTED, SAVE THE FACTORED VELOCITY FIELDS AS WELL
! ======================================================================
    IMPLICIT NONE
    COMPLEX(P8):: coeff, factor, eig_val
    INTEGER:: IND, NR_MK, MREAD
    REAL(P8):: AKREAD
    LOGICAL:: flip_switch, save_switch
    COMPLEX(P8), DIMENSION(:,:):: EIG_VEC_R
    COMPLEX(P8):: EIG2PTVEL(2*NR_MK)

    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: PSI,CHI,DEL2CHI
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_new, RUP_new, UZ_new
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_old, RUP_old, UZ_old
      
    ! convert to PSI,CHI
    ALLOCATE (PSI(NRCHOPDIM), DEL2CHI(NRCHOPDIM), CHI(NRCHOPDIM))  
    PSI = 0.D0; DEL2CHI = 0.D0; CHI = 0.D0
    PSI(:NR_MK) = EIG_VEC_R(:NR_MK, IND) !PSI
    DEL2CHI(:NR_MK) = EIG_VEC_R(NR_MK + 1:, IND) !DEL2CHI
    IF (flip_switch) THEN
        PSI = CONJG(PSI); DEL2CHI = CONJG(DEL2CHI)
    ENDIF
    CALL IDEL2_MK(DEL2CHI, CHI)

    ! save original
    IF (save_switch) THEN
        ALLOCATE (RUR_old(NRCHOPDIM), RUP_old(NRCHOPDIM), UZ_old(NRCHOPDIM))
        RUR_old = 0.D0; RUP_old = 0.D0; UZ_old = 0.D0
        CALL CHOPSET(3)
        CALL PC2VEL_MK(PSI, CHI, RUR_old, RUP_old, UZ_old)
        CALL RTRAN_MK(RUR_old, 1); CALL RTRAN_MK(RUP_old, 1); CALL RTRAN_MK(UZ_old, 1)
        CALL CHOPSET(-3)
        IF (flip_switch) THEN
            RUR_old = CONJG(RUR_old); RUP_old = CONJG(RUP_old); UZ_old = CONJG(UZ_old)
        ENDIF
        IF (IND .LE. 61) CALL SAVE_VEL(RUR_old, RUP_old, UZ_old, './converg/check/before_vel_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
            //'_IND_'//ITOA3(IND)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ENDIF

    ! times coefficient
    EIG2PTVEL(:NR_MK) = PSI(:NR_MK)
    EIG2PTVEL(NR_MK + 1:) = CHI(:NR_MK)
    IF (flip_switch) THEN
        EIG2PTVEL = coeff*CONJG(EIG2PTVEL)
    ELSE
        EIG2PTVEL = coeff*EIG2PTVEL
    ENDIF

    ! save new
    IF (save_switch) THEN
        ALLOCATE (RUR_new(NRCHOPDIM), RUP_new(NRCHOPDIM), UZ_new(NRCHOPDIM))
        RUR_new = 0.D0; RUP_new = 0.D0; UZ_new = 0.D0
        IF (flip_switch) THEN
            PSI(:) = CONJG(factor*coeff)*PSI(:)
            CHI(:) = CONJG(factor*coeff)*CHI(:)
        ELSE
            PSI(:) = factor*coeff*PSI(:)
            CHI(:) = factor*coeff*CHI(:)
        ENDIF
        CALL CHOPSET(3)
        CALL PC2VEL_MK(PSI, CHI, RUR_new, RUP_new, UZ_new)
        CALL RTRAN_MK(RUR_new, 1); CALL RTRAN_MK(RUP_new, 1); CALL RTRAN_MK(UZ_new, 1)
        CALL CHOPSET(-3)
        IF (flip_switch) THEN
            RUR_new = CONJG(RUR_new); RUP_new = CONJG(RUP_new); UZ_new = CONJG(UZ_new)
        ENDIF
        IF (IND .LE. 61) CALL SAVE_VEL(RUR_new, RUP_new, UZ_new, './converg/check/after_vel_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
            //'_IND_'//ITOA3(IND)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    END IF

    ! print data
    IF (save_switch) THEN
        open (10, FILE='./converg/check/contribution_MK_'//ITOA3(MREAD)//'_'//ITOA4(INT(AKREAD*100)) &
            //'.output', STATUS='unknown', ACTION='WRITE', ACCESS='APPEND', IOSTAT=IS)
        if (IS .ne. 0) then
            print *, 'ERROR: createperturb -- Could not creat new file contribution.output'
            RETURN
        end if
        WRITE (10, 1345) IND, ABS(EIG_VAL), ABS(factor*coeff),&
        &MAXVAL(ABS(RUP_old/TFM%R(1:NDIMR))),&
        &MAXVAL(ABS(RUP_new/TFM%R(1:NDIMR))),&
        &MAXVAL(ABS(EIG_VEC_R(:, IND))),&
        &MAXVAL(ABS(factor*COEFF*EIG_VEC_R(:, IND))),&
        &MAXVAL(ABS(RUP_old)),&
        &MAXVAL(ABS(RUP_new))
        close(10)
        DEALLOCATE(RUR_old, RUP_old, UZ_old)
        DEALLOCATE(RUR_new,RUP_new,UZ_new)
    END IF
    1345    FORMAT(I4, 8(',', G20.12))

    DEALLOCATE(PSI,CHI,DEL2CHI)

END FUNCTION EIG2PTVEL
! ======================================================================

FUNCTION ISZERO(EIG_VEC_R,IND,NR_MK,flip_switch)
! ======================================================================
! CHECK IF AN ERIGENVECTOR IS PRURELY ZERO
! ======================================================================
    IMPLICIT NONE
    INTEGER:: IND, NR_MK
    LOGICAL:: flip_switch, ISZERO
    COMPLEX(P8), DIMENSION(:,:):: EIG_VEC_R

    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: PSI,CHI,DEL2CHI
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_new, RUP_new, UZ_new
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_old, RUP_old, UZ_old

    ALLOCATE (PSI(NRCHOPDIM), DEL2CHI(NRCHOPDIM), CHI(NRCHOPDIM))  
    PSI = 0.D0; DEL2CHI = 0.D0; CHI = 0.D0
    PSI(:NR_MK) = EIG_VEC_R(:NR_MK, IND) !PSI
    DEL2CHI(:NR_MK) = EIG_VEC_R(NR_MK + 1:, IND) !DEL2CHI
    IF (flip_switch) THEN
        PSI = CONJG(PSI); DEL2CHI = CONJG(DEL2CHI)
    ENDIF
    CALL IDEL2_MK(DEL2CHI, CHI)

    ! CALL RTRAN_MK(PSI,1); CALL RTRAN_MK(CHI,1);
    ! IF (SUM(ABS(PSI(1:NR)))+SUM(ABS(CHI(1:NR))) .EQ. 0.D0) THEN
    !     ISZERO = .TRUE.
    ! ELSE
    !     ISZERO = .FALSE.
    !     WRITE(*,*) 'PSI:',IND,PSI(1:NR)
    !     WRITE(*,*) 'CHI:',CHI(1:NR)
    ! ENDIF

    ALLOCATE (RUR_old(NRCHOPDIM), RUP_old(NRCHOPDIM), UZ_old(NRCHOPDIM))
    RUR_old = 0.D0; RUP_old = 0.D0; UZ_old = 0.D0
    CALL CHOPSET(3)
    CALL PC2VEL_MK(PSI, CHI, RUR_old, RUP_old, UZ_old)
    CALL RTRAN_MK(RUR_old, 1); CALL RTRAN_MK(RUP_old, 1); CALL RTRAN_MK(UZ_old, 1)
    CALL CHOPSET(-3)
    IF (flip_switch) THEN
        RUR_old = CONJG(RUR_old); RUP_old = CONJG(RUP_old); UZ_old = CONJG(UZ_old)
    ENDIF
    IF (SUM(ABS(RUR_old(1:NR)))+SUM(ABS(RUP_old(1:NR)))+SUM(ABS(UZ_old(1:NR))) .EQ. 0.D0) THEN
        ISZERO = .TRUE.
    ELSE
        ISZERO = .FALSE.
        WRITE(*,*) IND,SUM(ABS(RUR_old(1:NR)))+SUM(ABS(RUP_old(1:NR))),SUM(ABS(UZ_old(1:NR)))
    ENDIF
    DEALLOCATE(RUR_old,RUP_old,UZ_old)


END FUNCTION ISZERO
! ======================================================================

SUBROUTINE PERTURB_MATRIX(M_0,AK_0,RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0,MREAD,AKREAD,H,comm_grp,switch)
! ======================================================================
! MREAD,AKREAD: mode interacted with perturbation
! M_0  ,  AK_0: mode of perturbation
! MREAD+M_0, AKREAD+AK_0: H
! ======================================================================
    IMPLICIT NONE
    INTEGER,INTENT(IN):: M_0, MREAD, comm_grp
    REAL(P8),INTENT(IN):: AK_0
    REAL(P8),INTENT(IN):: AKREAD
    COMPLEX(P8),DIMENSION(:),INTENT(IN):: RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0
    COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT):: H
    LOGICAL, OPTIONAL:: switch

    INTEGER:: M_END
    REAL(P8):: AK_END
    INTEGER:: NR_MK, I
    INTEGER:: MPI_NR_SIZE, MPI_NR_INDEX
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RURU, RUPU, UZU, RORU, ROPU, OZU
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: PSIU, CHIU, PSI1, CHI1, PSI2, CHI2
    COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE:: RUR_PSI, RUP_PSI, UZ_PSI
    COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE:: RUR_CHI, RUP_CHI, UZ_CHI
    CHARACTER(len=6) :: K_WRITTEN

    ! INIT
    NTH = 2
    NX = 4
    NTCHOP = 2 ! M = 0, M = MREAD
    NXCHOP = 2 ! K = 0, K = AKREAD
    ZLEN = 2*PI/AKREAD
    CALL MPI_BCAST(MREAD, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_BCAST(ZLEN, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
    CALL LEGINIT(comm_grp, MREAD)

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
    ALLOCATE (RURU(NRCHOPDIM), RUPU(NRCHOPDIM), UZU(NRCHOPDIM))
    ALLOCATE (RORU(NRCHOPDIM), ROPU(NRCHOPDIM), OZU(NRCHOPDIM))
    ALLOCATE (RUR_PSI(NRCHOPDIM,MPI_NR_SIZE),RUP_PSI(NRCHOPDIM,MPI_NR_SIZE),UZ_PSI(NRCHOPDIM,MPI_NR_SIZE))
    ALLOCATE (RUR_CHI(NRCHOPDIM,MPI_NR_SIZE),RUP_CHI(NRCHOPDIM,MPI_NR_SIZE),UZ_CHI(NRCHOPDIM,MPI_NR_SIZE))

    DO I = 1, MPI_NR_SIZE
        ! ================================ PSI =================================
        PSIU = 0.D0
        CHIU = 0.D0
        PSIU(MPI_NR_INDEX + I) = 1.D0

        CALL CHOPSET(3)
        ! NONLINEAR TERM: U0 X W' - W0 X U'
        CALL PC2VEL_MK(PSIU, CHIU, RURU, RUPU, UZU)
        CALL PC2VOR_MK(PSIU, CHIU, RORU, ROPU, OZU)
        CALL RTRAN_MK(RURU, 1)
        CALL RTRAN_MK(RUPU, 1)
        CALL RTRAN_MK(UZU, 1)
        CALL RTRAN_MK(RORU, 1)
        CALL RTRAN_MK(ROPU, 1)
        CALL RTRAN_MK(OZU, 1)
        CALL NONLIN_MK(MREAD,AKREAD,RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0,RURU,RUPU,UZU,RORU,ROPU,OZU,comm_grp=comm_grp)
        CALL CHOPSET(-3)

        RUR_PSI(:,I) = RORU
        RUP_PSI(:,I) = ROPU
        UZ_PSI(:,I) =  OZU

        ! ================================ CHI =================================
        PSIU = 0.D0
        CHIU = 0.D0
        CHIU(MPI_NR_INDEX + I) = 1.D0

        CALL IDEL2_MK(CHIU, CHI1)
        CHIU = CHI1
        CHI1 = 0.D0

        CALL CHOPSET(3)
        ! NONLINEAR TERM: U0 X W' - W0 X U'
        CALL PC2VEL_MK(PSIU, CHIU, RURU, RUPU, UZU)
        CALL PC2VOR_MK(PSIU, CHIU, RORU, ROPU, OZU)
        CALL RTRAN_MK(RURU, 1)
        CALL RTRAN_MK(RUPU, 1)
        CALL RTRAN_MK(UZU, 1)
        CALL RTRAN_MK(RORU, 1)
        CALL RTRAN_MK(ROPU, 1)
        CALL RTRAN_MK(OZU, 1)
        CALL NONLIN_MK(MREAD,AKREAD,RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0,RURU,RUPU,UZU,RORU,ROPU,OZU,comm_grp=comm_grp)
        CALL CHOPSET(-3)

        RUR_CHI(:,I) = RORU
        RUP_CHI(:,I) = ROPU
        UZ_CHI(:,I) =  OZU

    END DO

    ! INIT
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    NTH = 2
    NX = 4
    NTCHOP = 2 ! M = 0, M = M_END
    NXCHOP = 2 ! K = 0, K = AK_END        
    M_END = M_0 + MREAD
    AK_END = AK_0 + AKREAD
    ZLEN = 2*PI/AK_END
    CALL MPI_BCAST(M_END, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_BCAST(ZLEN, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
    CALL LEGINIT(comm_grp, M_END)

    CALL CHOPSET(3)
    DO I = 1, MPI_NR_SIZE
        PSI1 = 0.D0
        CHI1 = 0.D0
        PSI2 = 0.D0
        CHI2 = 0.D0
        CALL PROJECT_MK(RUR_PSI(:,I),RUP_PSI(:,I),UZ_PSI(:,I),PSI1,CHI1)
        CALL PROJECT_MK(RUR_CHI(:,I),RUP_CHI(:,I),UZ_CHI(:,I),PSI2,CHI2)

        H(:NR_MK,MPI_NR_INDEX+I) = PSI1(:NR_MK)
        H(NR_MK+1:,MPI_NR_INDEX+I) = CHI1(:NR_MK)
        H(:NR_MK,MPI_NR_INDEX+I+NR_MK) = PSI2(:NR_MK)
        H(NR_MK+1:,MPI_NR_INDEX+I+NR_MK) = CHI2(:NR_MK)

    END DO
    CALL CHOPSET(-3)
    DEALLOCATE (RURU, RUPU, UZU, RORU, ROPU, OZU)
    DEALLOCATE (PSIU, CHIU, PSI1, CHI1, PSI2, CHI2)
    DEALLOCATE (RUR_PSI, RUP_PSI, UZ_PSI, RUR_CHI, RUP_CHI, UZ_CHI)

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, H, SIZE(H), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                        MPI_COMM_WORLD, IERR)

    ! OUTPUT:
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    IF (MPI_GLB_RANK .EQ. 0) THEN  ! START OF THE SERIAL PART
        IF (PRESENT(switch)) THEN
            WRITE(K_WRITTEN,'(F06.2)') AK_END
            open (10, FILE='./converg/eigM_MK_'//ITOA3(M_END)//'_'//K_WRITTEN &
                    //'.output', STATUS='unknown', ACTION='WRITE', IOSTAT=IS)
            if (IS .ne. 0) then
                print *, 'ERROR: PERTURB_MATRIX -- Could not creat new file'
                RETURN
            end if
            DO I = 1, 2*NR_MK
                WRITE (10, 1252) H(I, 1:20)
                ! 1252 FORMAT((F5.3,SP,F5.3,'i'),19(',',F5.3,SP,F5.3,'i'))
                1252 FORMAT((S,E14.6E3,SP,E14.6E3,'i'),*(','S,E14.6E3,SP,E14.6E3,'i'))
            END DO
            close (10)
            write (*, *) 'PERTURB MAT (M=',M_END,', AK=',AK_END,'):  max element = ', maxval(abs(H))
            write (*, *) 'PERTURB MAT: maxUR = ', maxval(abs(RUR_0(1:NR)/TFM%R))
            write (*, *) 'PERTURB MAT: maxUP = ', maxval(abs(RUP_0(1:NR)/TFM%R))
            write (*, *) 'PERTURB MAT: maxUZ = ', maxval(abs( UZ_0(1:NR)))
            write (*, *) 'PERTURB MAT: maxOR = ', maxval(abs(ROR_0(1:NR)/TFM%R))
            write (*, *) 'PERTURB MAT: maxOP = ', maxval(abs(ROP_0(1:NR)/TFM%R))
            write (*, *) 'PERTURB MAT: maxOZ = ', maxval(abs( OZ_0(1:NR)))
        END IF
    END IF

    RETURN

END SUBROUTINE PERTURB_MATRIX
!=======================================================================

SUBROUTINE PERTURB_MATRIX2(M_0,AK_0,RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0,MREAD,AKREAD,H,comm_grp,print_switch)
! ======================================================================
! MREAD,AKREAD: mode interacted with perturbation
! M_0  ,  AK_0: mode of perturbation
! MREAD+M_0, AKREAD+AK_0: H
! ======================================================================
    IMPLICIT NONE
    INTEGER,INTENT(IN):: M_0, MREAD, comm_grp
    REAL(P8),INTENT(IN):: AK_0
    REAL(P8),INTENT(IN):: AKREAD
    COMPLEX(P8),DIMENSION(:),INTENT(IN):: RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0
    COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT):: H
    LOGICAL, OPTIONAL:: print_switch

    LOGICAL:: flip_switch
    INTEGER:: M_END, M_ACTUAL
    REAL(P8):: AK_END, AK_ACTUAL
    INTEGER:: NR_MK_0, NR_MK_1, I
    INTEGER:: MPI_NR_SIZE, MPI_NR_INDEX
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RURU, RUPU, UZU, RORU, ROPU, OZU
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: PSIU, CHIU, PSI1, CHI1, PSI2, CHI2
    COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE:: RUR_PSI, RUP_PSI, UZ_PSI
    COMPLEX(P8), DIMENSION(:,:), ALLOCATABLE:: RUR_CHI, RUP_CHI, UZ_CHI
    CHARACTER(len=6) :: K_WRITTEN

    ! INIT - {MREAD,AKREAD}
    NTH = 2; NX = 4
    NTCHOP = 2 ! M = 0, M = MREAD
    NXCHOP = 2 ! K = 0, K = AKREAD
    IF (MPI_GLB_RANK.EQ.0) THEN
        ! IF MREAD IS NEGATIVE, CALCULATE ITS C.C MODE
        IF (MREAD.LT.0) THEN
            flip_switch = .TRUE.
            M_ACTUAL = -MREAD; AK_ACTUAL = -AKREAD
        ELSE
            flip_switch = .FALSE.
            M_ACTUAL = MREAD; AK_ACTUAL = AKREAD
        ENDIF
    ENDIF
    CALL MPI_BCAST(M_ACTUAL, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_BCAST(AK_ACTUAL, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_BCAST(flip_switch, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERR)
    ZLEN = 2*PI/AK_ACTUAL
    CALL LEGINIT(comm_grp, M_ACTUAL)

    ! SPLIT NRCHOPDIM
    NR_MK_0 = NRCHOPS(2)
    CALL DECOMPOSE(NR_MK_0, MPI_GLB_PROCS, MPI_GLB_RANK, MPI_NR_SIZE, MPI_NR_INDEX)

    ! ALLOCATION
    ALLOCATE (PSIU(NRCHOPDIM), CHIU(NRCHOPDIM))
    ALLOCATE (PSI1(NRCHOPDIM), CHI1(NRCHOPDIM))
    ALLOCATE (PSI2(NRCHOPDIM), CHI2(NRCHOPDIM))
    ALLOCATE (RUR_PSI(NDIMR,MPI_NR_SIZE),RUP_PSI(NDIMR,MPI_NR_SIZE),UZ_PSI(NDIMR,MPI_NR_SIZE))
    ALLOCATE (RUR_CHI(NDIMR,MPI_NR_SIZE),RUP_CHI(NDIMR,MPI_NR_SIZE),UZ_CHI(NDIMR,MPI_NR_SIZE))

    ! GET NONLINEAR TERM FOR EACH TEST FUNC IN PFF SPACE - {MREAD,AKREAD}
    DO I = 1, MPI_NR_SIZE
        ! ================================ PSI =================================
        PSIU = 0.D0; CHIU = 0.D0
        PSIU(MPI_NR_INDEX + I) = 1.D0
        ALLOCATE (RURU(NRCHOPDIM), RUPU(NRCHOPDIM), UZU(NRCHOPDIM))
        ALLOCATE (RORU(NRCHOPDIM), ROPU(NRCHOPDIM), OZU(NRCHOPDIM))

        CALL CHOPSET(3)
        ! NONLINEAR TERM: U0 X W' - W0 X U'
        CALL PC2VEL_MK(PSIU, CHIU, RURU, RUPU, UZU)
        CALL PC2VOR_MK(PSIU, CHIU, RORU, ROPU, OZU)
        CALL RTRAN_MK(RURU, 1); CALL RTRAN_MK(RUPU, 1); CALL RTRAN_MK(UZU, 1)
        CALL RTRAN_MK(RORU, 1); CALL RTRAN_MK(ROPU, 1); CALL RTRAN_MK(OZU, 1)
        ! BACK TO {MREAD,AKREAD}
        IF (flip_switch) THEN 
            RURU = CONJG(RURU); RUPU = CONJG(RUPU); UZU = CONJG(UZU)
            RORU = CONJG(RORU); ROPU = CONJG(ROPU); OZU = CONJG(OZU)
        ENDIF
        ! {M_0,AK_0} + {MREAD,AKREAD}
        CALL NONLIN_MK(MREAD,AKREAD,RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0,RURU,RUPU,UZU,RORU,ROPU,OZU,comm_grp=comm_grp)
        CALL CHOPSET(-3)

        RUR_PSI(:,I) = RORU; RUP_PSI(:,I) = ROPU; UZ_PSI(:,I) = OZU
        DEALLOCATE (RURU, RUPU, UZU)
        DEALLOCATE (RORU, ROPU, OZU)

        ! ================================ CHI =================================
        PSIU = 0.D0; CHIU = 0.D0
        CHIU(MPI_NR_INDEX + I) = 1.D0
        ALLOCATE (RURU(NRCHOPDIM), RUPU(NRCHOPDIM), UZU(NRCHOPDIM))
        ALLOCATE (RORU(NRCHOPDIM), ROPU(NRCHOPDIM), OZU(NRCHOPDIM))

        CALL IDEL2_MK(CHIU, CHI1)
        CHIU = CHI1
        CHI1 = 0.D0

        CALL CHOPSET(3)
        ! NONLINEAR TERM: U0 X W' - W0 X U'
        CALL PC2VEL_MK(PSIU, CHIU, RURU, RUPU, UZU)
        CALL PC2VOR_MK(PSIU, CHIU, RORU, ROPU, OZU)
        CALL RTRAN_MK(RURU, 1); CALL RTRAN_MK(RUPU, 1); CALL RTRAN_MK(UZU, 1)
        CALL RTRAN_MK(RORU, 1); CALL RTRAN_MK(ROPU, 1); CALL RTRAN_MK(OZU, 1)
        ! BACK TO {MREAD,AKREAD}
        IF (flip_switch) THEN 
            RURU = CONJG(RURU); RUPU = CONJG(RUPU); UZU = CONJG(UZU)
            RORU = CONJG(RORU); ROPU = CONJG(ROPU); OZU = CONJG(OZU)
        ENDIF
        ! {M_0,AK_0} + {MREAD,AKREAD}
        CALL NONLIN_MK(MREAD,AKREAD,RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0,RURU,RUPU,UZU,RORU,ROPU,OZU,comm_grp=comm_grp)
        CALL CHOPSET(-3)

        RUR_CHI(:,I) = RORU; RUP_CHI(:,I) = ROPU; UZ_CHI(:,I) = OZU
        DEALLOCATE (RURU, RUPU, UZU)
        DEALLOCATE (RORU, ROPU, OZU)

    END DO

    ! INIT - {M_END,AK_END}
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    NTH = 2; NX = 4
    NTCHOP = 2 ! M = 0, M = M_END
    NXCHOP = 2 ! K = 0, K = AK_END        
    M_END = M_0 + MREAD; AK_END = AK_0 + AKREAD
    IF (MPI_GLB_RANK.EQ.0) THEN
        ! IF MREAD IS NEGATIVE, CALCULATE ITS C.C MODE
        IF (M_END.LT.0) THEN
            flip_switch = .TRUE.
            M_ACTUAL = -M_END; AK_ACTUAL = -AK_END;
        ELSE
            flip_switch = .FALSE.
            M_ACTUAL = M_END; AK_ACTUAL = AK_END;
        ENDIF
    ENDIF
    CALL MPI_BCAST(M_ACTUAL, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_BCAST(AK_ACTUAL, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_BCAST(flip_switch, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, IERR)
    ZLEN = 2*PI/AK_ACTUAL
    CALL LEGINIT(comm_grp, M_ACTUAL)

    IF (flip_switch) THEN
        RUR_PSI = CONJG(RUR_PSI); RUP_PSI = CONJG(RUP_PSI); UZ_PSI = CONJG(UZ_PSI)
        RUR_CHI = CONJG(RUR_CHI); RUP_CHI = CONJG(RUP_CHI); UZ_CHI = CONJG(UZ_CHI)
    ENDIF

    ! PERTURBATION MATRIX - 2*NR_MK_END X 2*NR_MK_READ
    NR_MK_1 = NRCHOPS(2)
    IF (ALLOCATED(H) .AND. (SIZE(H,1).NE.2*NR_MK_1) .AND. (SIZE(H,2).NE.2*NR_MK_0)) THEN
        DEALLOCATE (H)
        ALLOCATE (H(2*NR_MK_1, 2*NR_MK_0))
    ELSEIF (.NOT. (ALLOCATED(H))) THEN
        ALLOCATE (H(2*NR_MK_1, 2*NR_MK_0))
    END IF
    H = CMPLX(0.D0, 0.D0, P8)

    ! PROJECT NONLINEAR TERM TO PSI,DEL2CHI IN {M_END,AK_END}
    CALL CHOPSET(3)
    DO I = 1, MPI_NR_SIZE
        PSI1 = 0.D0
        CHI1 = 0.D0
        PSI2 = 0.D0
        CHI2 = 0.D0
        CALL PROJECT_MK(RUR_PSI(:,I),RUP_PSI(:,I),UZ_PSI(:,I),PSI1,CHI1)
        CALL PROJECT_MK(RUR_CHI(:,I),RUP_CHI(:,I),UZ_CHI(:,I),PSI2,CHI2)

        H(:NR_MK_1,MPI_NR_INDEX+I) = PSI1(:NR_MK_1)
        H(NR_MK_1+1:,MPI_NR_INDEX+I) = CHI1(:NR_MK_1)
        H(:NR_MK_1,MPI_NR_INDEX+I+NR_MK_0) = PSI2(:NR_MK_1)
        H(NR_MK_1+1:,MPI_NR_INDEX+I+NR_MK_0) = CHI2(:NR_MK_1)

    END DO
    CALL CHOPSET(-3)
    ! DEALLOCATE (RURU, RUPU, UZU, RORU, ROPU, OZU)
    DEALLOCATE (PSIU, CHIU, PSI1, CHI1, PSI2, CHI2)
    DEALLOCATE (RUR_PSI, RUP_PSI, UZ_PSI, RUR_CHI, RUP_CHI, UZ_CHI)

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, H, SIZE(H), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                        MPI_COMM_WORLD, IERR)
    IF (flip_switch) H = CONJG(H)

    ! OUTPUT:
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
    IF (MPI_GLB_RANK .EQ. 0) THEN  ! START OF THE SERIAL PART

        ! WRITE(K_WRITTEN,'(F06.2)') AK_END
        ! write (*, *) 'PERTURB MAT (M=',ITOA3(M_END),', AK=',ADJUSTL(RTOA6(AK_END)),'):  max element = ', ADJUSTL(RTOA6(maxval(abs(H)))),' - LOCATION:',maxloc(abs(H))

        IF ((PRESENT(print_switch)).AND.(print_switch)) THEN
            WRITE(K_WRITTEN,'(F06.2)') AK_END
            open (10, FILE='./converg/eigM_MK_'//ITOA3(M_END)//'_'//K_WRITTEN &
                    //'.output', STATUS='unknown', ACTION='WRITE', IOSTAT=IS)
            if (IS .ne. 0) then
                print *, 'ERROR: PERTURB_MATRIX -- Could not creat new file'
                RETURN
            end if
            DO I = 1, 2*NR_MK_1
                WRITE (10, 1252) H(I, 1:20)
                ! 1252 FORMAT((F5.3,SP,F5.3,'i'),19(',',F5.3,SP,F5.3,'i'))
                1252 FORMAT((S,E14.6E3,SP,E14.6E3,'i'),*(','S,E14.6E3,SP,E14.6E3,'i'))
            END DO
            close (10)
            ! write (*, *) 'PERTURB MAT (M=',ITOA3(M_END),', AK=',TRIM(RTOA6(AK_END)),'):  max element = ', TRIM(RTOA6(maxval(abs(H)))),' - LOCATION:',maxloc(abs(H))
            write (*, *) 'PERTURB MAT: maxUR = ', maxval(abs(RUR_0(1:NR)/TFM%R))
            write (*, *) 'PERTURB MAT: maxUP = ', maxval(abs(RUP_0(1:NR)/TFM%R))
            write (*, *) 'PERTURB MAT: maxUZ = ', maxval(abs( UZ_0(1:NR)))
            write (*, *) 'PERTURB MAT: maxOR = ', maxval(abs(ROR_0(1:NR)/TFM%R))
            write (*, *) 'PERTURB MAT: maxOP = ', maxval(abs(ROP_0(1:NR)/TFM%R))
            write (*, *) 'PERTURB MAT: maxOZ = ', maxval(abs( OZ_0(1:NR)))
        END IF
    END IF

    RETURN

END SUBROUTINE PERTURB_MATRIX2
!=======================================================================

SUBROUTINE TEST_ENSTROPHY(ROR,ROP,OZ)
! ======================================================================
! This enstrophy calculation is indepedent from the wavenumbers.
! The calculation for (0,0) mode needs to be divided by 2.
! ======================================================================
    IMPLICIT NONE
    COMPLEX(P8):: ROR(:), ROP(:), OZ(:)
    COMPLEX(P8):: OR(NR), OP(NR), OZ2(NR)
    REAL(P8):: ENS

    OR = ROR(1:NR)/(1-TFM%X)**2/(TFM%R)**2
    OP = ROP(1:NR)/(1-TFM%X)**2/(TFM%R)**2
    OZ2 = OZ(1:NR)/(1-TFM%X)**2

    ENS = PROD_MK(ROR,OR)+PROD_MK(ROP,OP)+PROD_MK(OZ,OZ2)

    write(*,*) 'M,K = ',M(2),AK(2,2), ', ZLEN = ',ZLEN
    write(*,*) 'ENS = ',ENS/ZLEN0,MAXVAL(ABS(OZ))
    
END SUBROUTINE
! ======================================================================

SUBROUTINE PERTURB_SAVE(EIG_VEC,M1,AK1,M2,AK2,comm_grp,IND,print_switch)
! ======================================================================
! convert the perturbed system's eigvec to velocity in PFF space and sa-
! ve both the vel and the psi,chi to files. if IND is NOT provided, 
! eigvec is converted from [psi,del2chi] to [psi,chi]
! ======================================================================
    IMPLICIT NONE
    INTEGER:: M1, M2, comm_grp
    REAL(P8):: AK1, AK2
    COMPLEX(P8),DIMENSION(:),INTENT(INOUT):: EIG_VEC
    INTEGER,OPTIONAL:: IND
    LOGICAL,OPTIONAL:: print_switch

    LOGICAL:: flip_switch
    REAL(P8):: AK_ACTUAL, ENE_MK
    CHARACTER(LEN=6) :: AK_PRINT
    INTEGER:: NR_MK_1, NR_MK_2, M_ACTUAL
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: PSI1, DEL2CHI1, CHI1, PSI2, DEL2CHI2, CHI2
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: RUR, RUP, UZ, VEC_R

    NR_MK_1 = MAX(MIN(NRCHOP,NRCHOP-ABS(M1)),0)
    NR_MK_2 = MAX(MIN(NRCHOP,NRCHOP-ABS(M2)),0)
    IF (SIZE(EIG_VEC).NE.2*(NR_MK_1+NR_MK_2)) THEN
        WRITE(*,*) 'PERTURB_SAVE: DIMENSION MISMATCH'
        WRITE(*,*) '     NR_MK_1:', ITOA4(NR_MK_1)
        WRITE(*,*) '     NR_MK_2:', ITOA4(NR_MK_2)
        WRITE(*,*) ' ACTUAL SIZE:', ITOA4(SIZE(EIG_VEC)/2)
        RETURN
    ENDIF

    ALLOCATE(PSI1(NRCHOPDIM),DEL2CHI1(NRCHOPDIM),CHI1(NRCHOPDIM))
    ALLOCATE(PSI2(NRCHOPDIM),DEL2CHI2(NRCHOPDIM),CHI2(NRCHOPDIM))

    PSI1(1:NR_MK_1)     = EIG_VEC(1:NR_MK_1)
    DEL2CHI1(1:NR_MK_1) = EIG_VEC(NR_MK_1+1:2*NR_MK_1)
    PSI2(1:NR_MK_2)     = EIG_VEC(2*NR_MK_1+1:2*NR_MK_1+NR_MK_2)
    DEL2CHI2(1:NR_MK_2) = EIG_VEC(2*NR_MK_1+NR_MK_2+1:2*NR_MK_1+2*NR_MK_2)

    ! INIT
    NTH = 2
    NX = 4
    NTCHOP = 2 ! M = 0, M = MREAD
    NXCHOP = 2 ! K = 0, K = AKREAD

    ! {m1,k1}
    ! IF (M1.LT.0) THEN
    !     flip_switch = .TRUE. ; M_ACTUAL = -M1; AK_ACTUAL = -AK1
    !     PSI1 = CONJG(PSI1); DEL2CHI1 = CONJG(DEL2CHI1)
    ! ELSE
    !     flip_switch = .FALSE.; M_ACTUAL =  M1; AK_ACTUAL =  AK1
    ! ENDIF
    ! ZLEN = 2*PI/AK_ACTUAL; CALL LEGINIT(comm_grp, M_ACTUAL)
    ! ALLOCATE(RUR(NRCHOPDIM),RUP(NRCHOPDIM),UZ(NRCHOPDIM),VEC_R(2*NR_MK_1))
    ! CALL IDEL2_MK(DEL2CHI1,CHI1)
    ! CALL CHOPSET(3); CALL PC2VEL_MK(PSI1,CHI1,RUR,RUP,UZ); 
    ! CALL RTRAN_MK(RUR,1); CALL RTRAN_MK(RUP,1); CALL RTRAN_MK(UZ,1)
    ! CALL CHOPSET(-3)
    ! VEC_R(1:NR_MK_1) = PSI1(1:NR_MK_1)
    ! VEC_R(NR_MK_1+1:) = CHI1(1:NR_MK_1)
    ! WRITE(AK_PRINT,'(F06.2)') AK1
    ! IF (flip_switch) THEN
    !     CALL SAVE_VEL(CONJG(RUR),CONJG(RUP),CONJG(UZ),'./perturb/vel_MK_'//ITOA3(M1)//'_'//AK_PRINT//'_P.output')
    !     CALL SAVE_PSICHI(CONJG(VEC_R),'./perturb/pt_MK_'//ITOA3(M1)//'_'//AK_PRINT//'_P.output')
    ! ELSE
    !     CALL SAVE_VEL(RUR,RUP,UZ,'./perturb/vel_MK_'//ITOA3(M1)//'_'//AK_PRINT//'_P.output')
    !     CALL SAVE_PSICHI(VEC_R,'./perturb/pt_MK_'//ITOA3(M1)//'_'//AK_PRINT//'_P.output')
    ! ENDIF
    ! DEALLOCATE(RUR,RUP,UZ,VEC_R)

    WRITE(AK_PRINT,'(F06.2)') AK1
    ! SAVE VELOCITY
    IF ((PRESENT(print_switch)).AND.(print_switch)) THEN
        ALLOCATE(VEC_R(2*NR_MK_1))
        VEC_R(1:2*NR_MK_1) = EIG_VEC(1:2*NR_MK_1)
        CALL EIG2VEL(M1,AK1,VEC_R,RUR,RUP,UZ,comm_grp,print_switch=.false.)
        
        IF (PRESENT(IND)) THEN
            CALL SAVE_VEL(RUR,RUP,UZ,'./perturb/vel_MK_'//ITOA3(M1)//'_'//AK_PRINT//'_P2'//ITOA3(IND)//'.output')
        ELSE
            CALL SAVE_VEL(RUR,RUP,UZ,'./perturb/vel_MK_'//ITOA3(M1)//'_'//AK_PRINT//'_P2.output')
        ENDIF
        DEALLOCATE(VEC_R,RUR,RUP,UZ)
    ENDIF
    ! SAVE PSI,CHI
    IF (M1.LT.0) THEN
        flip_switch = .TRUE. ; M_ACTUAL = -M1; AK_ACTUAL = -AK1
        DEL2CHI1 = CONJG(DEL2CHI1)
    ELSE
        flip_switch = .FALSE.; M_ACTUAL =  M1; AK_ACTUAL =  AK1
    ENDIF
    ZLEN = 2*PI/AK_ACTUAL; CALL LEGINIT(comm_grp, M_ACTUAL)
    CALL IDEL2_MK(DEL2CHI1,CHI1)
    IF (flip_switch) CHI1 = CONJG(CHI1)
    ALLOCATE(VEC_R(2*NR_MK_1))
    VEC_R(1:NR_MK_1) = PSI1(1:NR_MK_1)
    VEC_R(NR_MK_1+1:) = CHI1(1:NR_MK_1)
    IF (PRESENT(IND)) THEN
        CALL SAVE_PSICHI(VEC_R,'./perturb/pt_MK_'//ITOA3(M1)//'_'//AK_PRINT//'_P2'//ITOA3(IND)//'.output')
    ELSE
        CALL SAVE_PSICHI(VEC_R,'./perturb/pt_MK_'//ITOA3(M1)//'_'//AK_PRINT//'_P2.output')
        ! NORMALIZE & SAVE VORTICITY
        ALLOCATE(RUR(NRCHOPDIM),RUP(NRCHOPDIM),UZ(NRCHOPDIM))
        CALL CHOPSET(3)
        CALL PC2VOR_MK(PSI1,CHI1,RUR,RUP,UZ)
        CALL RTRAN_MK(RUR,1); CALL RTRAN_MK(RUP,1); CALL RTRAN_MK( UZ,1)
        ENE_MK = ENSTROPHY_MK_MOD(RUR,RUP,UZ)/ZLEN0
        CALL CHOPSET(-3)
        CALL SAVE_VEL(RUR/ENE_MK**0.5,RUP/ENE_MK**0.5,UZ/ENE_MK**0.5,'./perturb/vor_MK_1.output')
        DEALLOCATE(RUR,RUP,UZ)
        ! CONVERT TO [PSI,CHI]
        EIG_VEC(1:2*NR_MK_1) = VEC_R(1:2*NR_MK_1)/ENE_MK**0.5
    ENDIF
    DEALLOCATE(VEC_R)

    ! {m2,k2}
    ! IF (M2.LT.0) THEN
    !     flip_switch = .TRUE. ; M_ACTUAL = -M2; AK_ACTUAL = -AK2
    !     PSI2 = CONJG(PSI2); DEL2CHI2 = CONJG(DEL2CHI2)
    ! ELSE
    !     flip_switch = .FALSE.; M_ACTUAL =  M2; AK_ACTUAL =  AK2
    ! ENDIF
    ! ZLEN = 2*PI/AK_ACTUAL; CALL LEGINIT(comm_grp, M_ACTUAL)
    ! ALLOCATE(RUR(NRCHOPDIM),RUP(NRCHOPDIM),UZ(NRCHOPDIM),VEC_R(2*NR_MK_2))
    ! CALL IDEL2_MK(DEL2CHI2,CHI2)
    ! CALL CHOPSET(3); CALL PC2VEL_MK(PSI2,CHI2,RUR,RUP,UZ); 
    ! CALL RTRAN_MK(RUR,1); CALL RTRAN_MK(RUP,1); CALL RTRAN_MK(UZ,1)
    ! CALL CHOPSET(-3)
    ! VEC_R(1:NR_MK_2) = PSI2(1:NR_MK_2)
    ! VEC_R(NR_MK_2+1:) = CHI2(1:NR_MK_2)
    ! WRITE(AK_PRINT,'(F06.2)') AK2
    ! IF (flip_switch) THEN
    !     CALL SAVE_VEL(CONJG(RUR),CONJG(RUP),CONJG(UZ),'./perturb/vel_MK_'//ITOA3(M2)//'_'//AK_PRINT//'_P.output')
    !     CALL SAVE_PSICHI(CONJG(VEC_R),'./perturb/pt_MK_'//ITOA3(M2)//'_'//AK_PRINT//'_P.output')
    ! ELSE
    !     CALL SAVE_VEL(RUR,RUP,UZ,'./perturb/vel_MK_'//ITOA3(M2)//'_'//AK_PRINT//'_P.output')
    !     CALL SAVE_PSICHI(VEC_R,'./perturb/pt_MK_'//ITOA3(M2)//'_'//AK_PRINT//'_P.output')
    ! ENDIF
    ! DEALLOCATE(RUR,RUP,UZ,VEC_R)

    WRITE(AK_PRINT,'(F06.2)') AK2
    ! SAVE VELOCITY
    IF ((PRESENT(print_switch)).AND.(print_switch)) THEN
        ALLOCATE(VEC_R(2*NR_MK_2))
        VEC_R(1:2*NR_MK_2) = EIG_VEC(2*NR_MK_1+1:2*NR_MK_1+2*NR_MK_2)
        CALL EIG2VEL(M2,AK2,VEC_R,RUR,RUP,UZ,comm_grp,print_switch=.false.)
        
        IF (PRESENT(IND)) THEN
            CALL SAVE_VEL(RUR,RUP,UZ,'./perturb/vel_MK_'//ITOA3(M2)//'_'//AK_PRINT//'_P2'//ITOA3(IND)//'.output')
        ELSE
            CALL SAVE_VEL(RUR,RUP,UZ,'./perturb/vel_MK_'//ITOA3(M2)//'_'//AK_PRINT//'_P2.output')
        ENDIF
        DEALLOCATE(VEC_R,RUR,RUP,UZ)
    ENDIF
    ! CONVERT TO [PSI,CHI]
    IF (M2.LT.0) THEN
        flip_switch = .TRUE. ; M_ACTUAL = -M2; AK_ACTUAL = -AK2
        PSI2 = CONJG(PSI2); DEL2CHI2 = CONJG(DEL2CHI2)
    ELSE
        flip_switch = .FALSE.; M_ACTUAL =  M2; AK_ACTUAL =  AK2
    ENDIF
    ZLEN = 2*PI/AK_ACTUAL; CALL LEGINIT(comm_grp, M_ACTUAL)
    CALL IDEL2_MK(DEL2CHI2,CHI2)
    ALLOCATE(VEC_R(2*NR_MK_2))
    VEC_R(1:NR_MK_2) = PSI2(1:NR_MK_2)
    VEC_R(NR_MK_2+1:) = CHI2(1:NR_MK_2)
    IF (PRESENT(IND)) THEN
        CALL SAVE_PSICHI(VEC_R,'./perturb/pt_MK_'//ITOA3(M2)//'_'//AK_PRINT//'_P2'//ITOA3(IND)//'.output')
    ELSE
        CALL SAVE_PSICHI(VEC_R,'./perturb/pt_MK_'//ITOA3(M2)//'_'//AK_PRINT//'_P2.output')
        ! SAVE VORTICITY
        ALLOCATE(RUR(NRCHOPDIM),RUP(NRCHOPDIM),UZ(NRCHOPDIM))
        CALL CHOPSET(3)
        CALL PC2VOR_MK(PSI2,CHI2,RUR,RUP,UZ)
        CALL RTRAN_MK(RUR,1); CALL RTRAN_MK(RUP,1); CALL RTRAN_MK( UZ,1)
        CALL CHOPSET(-3)
        CALL SAVE_VEL(RUR/ENE_MK**0.5,RUP/ENE_MK**0.5,UZ/ENE_MK**0.5,'./perturb/vor_MK_2.output')
        DEALLOCATE(RUR,RUP,UZ)
        ! CONVERT TO [PSI,CHI]
        EIG_VEC(2*NR_MK_1+1:2*NR_MK_1+2*NR_MK_2) = VEC_R(1:2*NR_MK_2)/ENE_MK**0.5
    ENDIF
    DEALLOCATE(VEC_R)

    ! DEALLOCATE
    DEALLOCATE(PSI1, DEL2CHI1, CHI1, PSI2, DEL2CHI2, CHI2)

END SUBROUTINE
! ======================================================================

END PROGRAM EVP_CORRECTION
!=======================================================================
