PROGRAM EVP_CLTRIAD
!=======================================================================
! [USAGE]:
! GENERATING PERTURBATION FILE FOR TRIADS INVOLVING CRITICAL LAYERS
!
! CURRENT CONFIGURATION:
! #     SIGMA       M       K
! 0     -0.153i     1       1.00
! 1     -0.144i     1       4.42
! 2     -0.305i     2       5.32
!
! [UPDATES]:
! CREATED ON JUN 11, 2024
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
REAL(P8)    :: K_0, K_1, K_2
COMPLEX(P8) :: EIG_0, EIG_1, EIG_2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig_0, M_eig_1, M_eig_2
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_0, EIG_R_mat_0, EIG_L_mat_0
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_1, EIG_R_mat_1, EIG_L_mat_1
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat_2, EIG_R_mat_2, EIG_L_mat_2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECR_0, VECL_0, RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_1, VECR_1, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: VECL_2, VECR_2, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2

! ADD PERTURB:
INTEGER:: IS
! SCANNING:
CHARACTER(LEN=72):: FILENAME
CHARACTER(len=5) :: K_VAL
INTEGER:: FID

! INITILIZATION - MPI
! ======================================================================
CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED, MPI_THREAD_MODE, IERR)
IF (MPI_THREAD_MODE .LT. MPI_THREAD_SERIALIZED) THEN
    WRITE (*, *) 'The threading support is lesser than that demanded.'
    CALL MPI_ABORT(MPI_COMM_WORLD, 1, IERR)
END IF

! INITILIZATION - MPI PROCESSORS
! ======================================================================
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR) ! DIVIDE NR BY #OF MPI PROCS

! INITILIZATION - BASIC INFO FROM READ.INPUT
! ======================================================================
IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE (*, *) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()   ! @ MOD_MISC
END IF
CALL READCOM('NOECHO')
CALL READIN(5)

! print '(l1)', .TRUE.
! CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)

! TRIAD INFO
! ======================================================================
! M_0 = 1; K_0 = 1.D0; EIG_0 = -0.153*IU
! M_1 = 1; K_1 = 4.42; EIG_1 = -0.144*IU
! M_2 = 2; K_2 = 5.42; EIG_2 = -0.305*IU
M_0 = 1; K_0 = 0.90; EIG_0 = -0.159*IU
M_1 = 1; K_1 = 4.50; EIG_1 = -0.143*IU
M_2 = 2; K_2 = 5.40; EIG_2 = -0.305*IU

! CALCULATION
!=======================================================================
! calculate operator matrix and its corresponding eigenvalues
CALL EIG_MATRIX(M_0,K_0,M_mat_0,M_eig_0,EIG_R_mat_0,EIG_L_mat_0,comm_grp=newcomm)
CALL EIG_MATRIX(M_1,K_1,M_mat_1,M_eig_1,EIG_R_mat_1,EIG_L_mat_1,comm_grp=newcomm)
CALL EIG_MATRIX(M_2,K_2,M_mat_2,M_eig_2,EIG_R_mat_2,EIG_L_mat_2,comm_grp=newcomm)

IF (MPI_GLB_RANK .EQ. 0) THEN

    WRITE (*, *) 'TRIAD:'
    201  FORMAT('Eig: M# = ',I3,'; AK = ',F20.13,'; Sig (#',I3,') =',(s1pE20.13,SP,s1pE20.13),'; Res - ',L1) 

    ! JJ = tol_ind(1); KK = tol_ind(2)

    ! {M0,K0}
    ! ======================================================================
    ! find the correct index
    II = MINLOC(abs(AIMAG(M_eig_0-EIG_0)),DIM=1,MASK=EIGRES(EIG_R_mat_0,M_0)); EIG_0 = M_eig_0(II)
    ! print mode info
    WRITE (*, 201) M_0, K_0, II, M_eig_0(II), EIGRES(EIG_R_mat_0(:,II), M_0)
    ! save velocity vectors
    CALL EIG2VELVOR(M_0, K_0, II, EIG_R_mat_0, EIG_L_mat_0, RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0, VECR_0, VECL_0,comm_grp=newcomm)
    WRITE(K_VAL,'(F05.2)') K_0; CALL SAVE_VEL(RUR_0, RUP_0, UZ_0, './converg/CriticalLayer_240605/vel_MK_'//ITOA3(M_0)//'_'//K_VAL//'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! call SAVE_VEL(ROR_0, ROP_0, OZ_0, './converg/vor_MK_'//ITOA3(M_0)//'_'//K_VAL &
    !             //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

    ! {M1,K1}
    ! ======================================================================
    ! find the correct index
    JJ = MINLOC(abs(AIMAG(M_eig_1-EIG_1)),DIM=1,MASK=EIGRES(EIG_R_mat_1,M_1)); EIG_1 = M_eig_1(JJ)
    ! print mode info
    WRITE (*, 201) M_1, K_1, JJ, M_eig_1(JJ), EIGRES(EIG_R_mat_1(:,JJ), M_1)
    ! save velocity vectors
    CALL EIG2VELVOR(M_1, K_1, JJ, EIG_R_mat_1, EIG_L_mat_1, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1, VECR_1, VECL_1, comm_grp=newcomm)
    WRITE(K_VAL,'(F05.2)') K_1; CALL SAVE_VEL(RUR_1, RUP_1, UZ_1, './converg/CriticalLayer_240605/vel_MK_'//ITOA3(M_1)//'_'//K_VAL//'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! call SAVE_VEL(ROR_1, ROP_1, OZ_1, './converg/vor_MK_'//ITOA3(M_1)//'_'//K_VAL &
    !             //'_IND_'//ITOA3(JJ)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

    ! {M2,K2}
    ! ======================================================================
    ! find the correct index
    KK = MINLOC(abs(AIMAG(M_eig_2-EIG_2)),DIM=1,MASK=EIGRES(EIG_R_mat_2,M_2)); EIG_2 = M_eig_2(KK)
    ! print mode info
    WRITE (*, 201) M_2, K_2, kk, M_eig_2(KK), EIGRES(EIG_R_mat_2(:,KK), M_2)
    ! save velocity vectors
    CALL EIG2VELVOR(M_2, K_2, KK, EIG_R_mat_2, EIG_L_mat_2, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2, VECR_2, VECL_2, comm_grp=newcomm)
    WRITE(K_VAL,'(F05.2)') K_2; CALL SAVE_VEL(RUR_2, RUP_2, UZ_2, './converg/CriticalLayer_240605/vel_MK_'//ITOA3(M_2)//'_'//K_VAL//'_IND_'//ITOA3(KK)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! call SAVE_VEL(ROR_2, ROP_2, OZ_2, './converg/vor_MK_'//ITOA3(M_2)//'_'//K_VAL &
    !             //'_IND_'//ITOA3(KK)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')

    ! BASEFLOW
    ! ======================================================================
    CALL COLLOC_INFO()
    CALL SAVE_VEL(CMPLX(RUR0), CMPLX(RUP0), CMPLX(UZ0), './converg/CriticalLayer_240605/vel_MK_'//ITOA3(0)//'_'//ITOA4(0) &
                    //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    ! CALL SAVE_VEL(CMPLX(ROR0), CMPLX(ROP0), CMPLX(OZ0), './converg/vor_MK_'//ITOA3(0)//'_'//ITOA4(0) &
    !                 //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    Write (*,*) 'MODE---0 ENSTROPHY = ', ENSTROPHY_MK_MOD(ROR_0, ROP_0, OZ_0)/ZLEN0
    Write (*,*) 'MODE---1 ENSTROPHY = ', ENSTROPHY_MK_MOD(ROR_1, ROP_1, OZ_1)/ZLEN0
    Write (*,*) 'MODE---2 ENSTROPHY = ', ENSTROPHY_MK_MOD(ROR_2, ROP_2, OZ_2)/ZLEN0
    Write (*,*) 'BASEFLOW ENSTROPHY = ', ENSTROPHY_MK_MOD(CMPLX(ROR0),CMPLX(ROP0),CMPLX(OZ0))/ZLEN0

    Write (*,*) 'MODE---0 ENERGY = ', ENSTROPHY_MK_MOD(RUR_0, RUP_0, UZ_0)/ZLEN0
    Write (*,*) 'MODE---1 ENERGY = ', ENSTROPHY_MK_MOD(RUR_1, RUP_1, UZ_1)/ZLEN0
    Write (*,*) 'MODE---2 ENERGY = ', ENSTROPHY_MK_MOD(RUR_2, RUP_2, UZ_2)/ZLEN0
    Write (*,*) 'BASEFLOW ENERGY = ', ENSTROPHY_MK_MOD(CMPLX(RUR0),CMPLX(RUP0),CMPLX(UZ0))/ZLEN0

    ! PERTURBATION FILES
    CALL SAVE_PERTURB('./converg/CriticalLayer_240605/new_perturb.input', &
                        M_0, K_0, EIG_R_mat_0(:, II), M_eig_0(II), &
                        M_1, K_1, EIG_R_mat_1(:, JJ), M_eig_1(JJ), &
                        M_2, K_2, EIG_R_mat_2(:, KK), M_eig_2(kk))

    ! DEALLOCATE
    DEALLOCATE(M_eig_0, RUR_0, RUP_0, UZ_0, ROR_0, ROP_0, OZ_0, VECR_0, VECL_0, EIG_R_mat_0, EIG_L_mat_0)
    DEALLOCATE(M_eig_1, RUR_1, RUP_1, UZ_1, ROR_1, ROP_1, OZ_1, VECR_1, VECL_1, EIG_R_mat_1, EIG_L_mat_1)
    DEALLOCATE(M_eig_2, RUR_2, RUP_2, UZ_2, ROR_2, ROP_2, OZ_2, VECR_2, VECL_2, EIG_R_mat_2, EIG_L_mat_2)
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

! DEALLOCATE
!=======================================================================
10  CONTINUE
DEALLOCATE(M_mat_0, M_mat_1, M_mat_2)
DEALLOCATE(RUR0, RUP0, UZ0, ROR0, ROP0, OZ0)

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

FUNCTION TARGET_DIST(KVAL, MIND1, MIND2, KBAR, SIG)
! ======================================================================
! Calcualte Eig1 + Sig_0 - Eig2
    IMPLICIT NONE
    INTEGER:: MIND1, MIND2, III
    REAL(P8):: KVAL, KVAL2, KBAR, TARGET_DIST
    REAL(P8), DIMENSION(:), ALLOCATABLE:: EIG_DIST_LIST, EIG_DIST_LIST_I
    COMPLEX(P8):: SIG
    COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: Meig1, Meig2
    COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: Mmat1, Mmat2

    CALL EIG_MATRIX(MIND1, KVAL, Mmat1, Meig1, comm_grp=newcomm)
    IF (MPI_GLB_RANK .EQ. 0) THEN
        WRITE (*, *) ''
        WRITE (*, *) KVAL
    END IF

    KVAL2 = KVAL + KBAR
    CALL EIG_MATRIX(MIND2, KVAL2, Mmat2, Meig2, comm_grp=newcomm)

    IF (MPI_GLB_RANK .EQ. 0) THEN
        ALLOCATE (EIG_DIST_LIST(SIZE(Meig1)))
        ALLOCATE (EIG_DIST_LIST_I(SIZE(Meig2)))

        DO III = 1, SIZE(Meig1)
            EIG_DIST_LIST_I = AIMAG(Meig1(III) + SIG - Meig2)
            EIG_DIST_LIST(III) = EIG_DIST_LIST_I(MINLOC(ABS(EIG_DIST_LIST_I), 1)) ! CLOSEST DISTANCE for III
        END DO
        TARGET_DIST = EIG_DIST_LIST(MINLOC(ABS(EIG_DIST_LIST), 1)) ! CLOSEST DISTANCE overall
        WRITE (*, *) TARGET_DIST

        DEALLOCATE (Meig1, Meig2, Mmat1, Mmat2, EIG_DIST_LIST, EIG_DIST_LIST_I)
    END IF
    CALL MPI_BCAST(TARGET_DIST, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR)

END FUNCTION TARGET_DIST
!=======================================================================

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

END PROGRAM EVP_CLTRIAD
!=======================================================================
