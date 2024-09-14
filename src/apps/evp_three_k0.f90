PROGRAM EVP_THREE_K0
!=======================================================================
! [USAGE]:
! Calculate interaction coefficients for h = 0, k = 0 case
! All three triad memebers are CL modes. 
! Note: for k = 0 and h = 0, each linear CL mode is continuous at the CL
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
INTEGER     :: II, JJ, KK
INTEGER     :: M_0, M_1, M_2
REAL(P8)    :: K_0, K_1, K_2, SIG_0, SIG_1, SIG_2, R_CTL, OMEGA_CTL
COMPLEX(P8) :: EIG_BAR, J1, J2, J0
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig_0, M_eig_1, M_eig_2
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_0, EIG_R_mat_0, EIG_L_mat_0
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_1, EIG_R_mat_1, EIG_L_mat_1
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_2, EIG_R_mat_2, EIG_L_mat_2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_0, VECR_0, VECR_12, RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_0_COPY, ROP_0_COPY, OZ_0_COPY
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_1, VECR_1, VECR_01, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_1_COPY, ROP_1_COPY, OZ_1_COPY
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_2, VECR_2, VECR_02, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_2_COPY, ROP_2_COPY, OZ_2_COPY

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

! ! =============================== CHECK ================================
! ! SPATIAL SYMMETRY
! QPAIR%H(1) = 0.D0
! M_0 = 2; M_1 = -1; M_2 = M_0 + M_1
! K_0 = -2; K_1 = -0.29; K_2 = K_0 + K_1 ! ALL K = 0
! SIG_0 = -2.1614; SIG_1 = 1.0264; SIG_2 = -1.1348
! ! ======================================================================

! SPATIAL SYMMETRY
QPAIR%H(1) = 0.D0
M_0 = 1; M_1 = 2; M_2 = M_0 + M_1
K_0 = 0; K_1 = 0; K_2 = K_0 + K_1 ! ALL K = 0

! CRITICAL LAYER
R_CTL = 0.21
OMEGA_CTL = (1-EXP(-R_CTL**2))/(R_CTL**2)
SIG_0 = -M_0*OMEGA_CTL; SIG_1 = -M_1*OMEGA_CTL; SIG_2 = -M_2*OMEGA_CTL

! CALCULATE EVP FOR EACH TRIAD
CALL EIG_MATRIX(M_0, K_0, M_mat_0, M_eig_0, EIG_R_mat_0, EIG_L_mat_0, comm_grp=newcomm) 
CALL EIG_MATRIX(M_1, K_1, M_mat_1, M_eig_1, EIG_R_mat_1, EIG_L_mat_1, comm_grp=newcomm)
CALL EIG_MATRIX(M_2, K_2, M_mat_2, M_eig_2, EIG_R_mat_2, EIG_L_mat_2, comm_grp=newcomm)

IF (MPI_GLB_RANK.EQ.0) THEN

    II = MINLOC(ABS(AIMAG(M_eig_0)-SIG_0),DIM=1)
    JJ = MINLOC(ABS(AIMAG(M_eig_1)-SIG_1),DIM=1)
    KK = MINLOC(ABS(AIMAG(M_eig_2)-SIG_2),DIM=1)

    ! PRINT TRIAD INFO
    WRITE (*, *) 'RESONANT TRIAD:'
    WRITE (*, 201) M_0, K_0, II, M_eig_0(II)
    WRITE (*, 201) M_1, K_1, JJ, M_eig_1(JJ)
    WRITE (*, 201) M_2, K_2, kk, M_eig_2(KK)
201  FORMAT('Eig: M# = ',I3,'; AK = ',F20.13,'; Sig (#',I3') =',1P8E20.13)

    ! SAVE PERTURBATION
    ! CALL SAVE_PERTURB('./qvortex/Lamb Oseen (h=0)/triad/new_perturb.input', &
    !                     M_0, K_0, EIG_R_mat_0(:, II), M_eig_0(II), &
    !                     M_1, K_1, EIG_R_mat_1(:, JJ), M_eig_1(JJ), &
    !                     M_2, K_2, EIG_R_mat_2(:, KK), M_eig_2(kk))

    ! FINAL RESULT:
    ! {M0,K0}
    CALL EIG2VELVOR(M_0, K_0, II, EIG_R_mat_0, EIG_L_mat_0, RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0, VECR_0, VECL_0,comm_grp=newcomm, print_switch=DIAGNOST_SWITCH)
    ALLOCATE(ROR_0_COPY(NDIMR), ROP_0_COPY(NDIMR), OZ_0_COPY(NDIMR))
    ROR_0_COPY = ROR_0; ROP_0_COPY = ROP_0; OZ_0_COPY = OZ_0

    ! {M1,K1}
    CALL EIG2VELVOR(M_1, K_1, JJ, EIG_R_mat_1, EIG_L_mat_1, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1, VECR_1, VECL_1, comm_grp=newcomm, print_switch=DIAGNOST_SWITCH)
    ALLOCATE(ROR_1_COPY(NDIMR),ROP_1_COPY(NDIMR),OZ_1_COPY(NDIMR))
    ROR_1_COPY = ROR_1; ROP_1_COPY = ROP_1; OZ_1_COPY = OZ_1

    ! {M2,K2}
    CALL EIG2VELVOR(M_2, K_2, KK, EIG_R_mat_2, EIG_L_mat_2, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2, VECR_2, VECL_2, comm_grp=newcomm, print_switch=DIAGNOST_SWITCH)
    ALLOCATE(ROR_2_COPY(NDIMR),ROP_2_COPY(NDIMR),OZ_2_COPY(NDIMR))
    ROR_2_COPY = ROR_2; ROP_2_COPY = ROP_2; OZ_2_COPY = OZ_2

    ! SAVE VELOCITTY
    ! WRITE(K_VAL,'(F06.2)') K_0; 
    ! CALL SAVE_VEL(RUR_0, RUP_0, UZ_0, './qvortex/Lamb Oseen (h=0)/triad/vel_MK_'//ITOA3(M_0)//'_'//K_VAL//'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! WRITE(K_VAL,'(F06.2)') K_1; 
    ! CALL SAVE_VEL(RUR_1, RUP_1, UZ_1, './qvortex/Lamb Oseen (h=0)/triad/vel_MK_'//ITOA3(M_1)//'_'//K_VAL//'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! WRITE(K_VAL,'(F06.2)') K_2; 
    ! CALL SAVE_VEL(RUR_2, RUP_2, UZ_2, './qvortex/Lamb Oseen (h=0)/triad/vel_MK_'//ITOA3(M_2)//'_'//K_VAL//'_IND_'//ITOA3(KK)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

    ! NONLINEAR TERMS:
    ! {M2,K2}: 0  X 1 -> 2
    CALL NONLIN_MK(M_2,K_2,RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0,RUR_1,RUP_1,UZ_1,ROR_1_COPY,ROP_1_COPY,OZ_1_COPY,VECR_01,comm_grp=newcomm)

    ! {M1,K1}: 0* X 2 -> 1
    CALL NONLIN_MK(M_1,K_1,CONJG(RUR_0),CONJG(RUP_0),CONJG(UZ_0),CONJG(ROR_0),CONJG(ROP_0),CONJG(OZ_0),RUR_2,RUP_2,UZ_2,ROR_2,ROP_2,OZ_2,VECR_02,comm_grp=newcomm)

    ! {M0,K0}: 1* X 2 -> 0
    CALL NONLIN_MK(M_0,K_0,CONJG(RUR_1),CONJG(RUP_1),CONJG(UZ_1),CONJG(ROR_1),CONJG(ROP_1),CONJG(OZ_1),RUR_2,RUP_2,UZ_2,ROR_2_COPY,ROP_2_COPY,OZ_2_COPY,VECR_12,comm_grp=newcomm)
    
ENDIF ! MPI_GLB_RANK.EQ.0

RANK0: IF (MPI_GLB_RANK .EQ. 0) THEN

! ============================== RESULT ================================
WRITE (*, *) "NEARLY DEGEN:"
J0 = DOT_PRODUCT(VECL_0, VECR_12); WRITE (*,*) "J0: ", J0
J1 = DOT_PRODUCT(VECL_1, VECR_02); WRITE (*,*) "J1: ", J1
J2 = DOT_PRODUCT(VECL_2, VECR_01); WRITE (*,*) "J2: ", J2
WRITE(*,*) "J1 + J0 - J2: ", J1 + J0 - J2

! DEALLOCATE
!=======================================================================
1   DEALLOCATE(RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0)
    DEALLOCATE(RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1)
    DEALLOCATE(RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2)
    DEALLOCATE(ROR_0_COPY, ROP_0_COPY, OZ_0_COPY)
    DEALLOCATE(ROR_1_COPY, ROP_1_COPY, OZ_1_COPY)
    DEALLOCATE(ROR_2_COPY, ROP_2_COPY, OZ_2_COPY)
    DEALLOCATE(VECR_0, VECL_0, VECR_1, VECR_2, VECL_1, VECL_2, VECR_01, VECR_02, VECR_12)
    DEALLOCATE(EIG_R_mat_0, EIG_L_mat_0, EIG_R_mat_1, EIG_L_mat_1, EIG_R_mat_2, EIG_L_mat_2)
    DEALLOCATE(M_eig_0, M_eig_1, M_eig_2)

END IF RANK0

! DEALLOCATE
!=======================================================================
10  CONTINUE
DEALLOCATE (RUR0, RUP0, UZ0, ROR0, ROP0, OZ0)
IF (ALLOCATED(M_mat_1)) THEN
    DEALLOCATE(M_mat_0,M_mat_1,M_mat_2)
ENDIF

IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE (*, *) ''
    WRITE (*, *) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
END IF

CALL MPI_FINALIZE(IERR)

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================

END PROGRAM EVP_THREE_K0
!=======================================================================
