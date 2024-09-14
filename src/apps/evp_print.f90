PROGRAM EVP_PRINT
!=======================================================================
! [USAGE]:
! FIND EIGENVALUES AND EIGENVECTORS OF THE LINEARIZED N-S EQUATIONS
! EXPRESSED IN A POLOIDAL-TOROLIDALLY DECOMPOSED FORM
! OPERATOR H CORRESPONDS TO EIGENVECTOR: [PSI,DEL2CHI]^T
!
! [UPDATES]:
! LAST UPDATE ON SEP 20, 2022
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
INTEGER     :: II, JJ, KK, K_ind
INTEGER     :: M_BAR, EIG_BAR_IND, M_1
REAL(P8)    :: K_BAR, K_VAL_1_0
REAL(P8),DIMENSION(:),ALLOCATABLE:: K_RANGE

COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: M_eig
COMPLEX(P8), DIMENSION(:, :), ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BAR, RUP_BAR, UZ_BAR, EIG_R_BAR
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_BAR, ROP_BAR, OZ_BAR, EIG_L_BAR

COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: RUR_BAR2, RUP_BAR2, UZ_BAR2, EIG_R_BAR2
COMPLEX(P8), DIMENSION(:), ALLOCATABLE:: ROR_BAR2, ROP_BAR2, OZ_BAR2, EIG_L_BAR2

! SCANNING:
CHARACTER(len=6) :: K_VAL
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

! {M_bar, AK_bar}: BASE FLOW
!=======================================================================
! IF (MPI_GLB_RANK .EQ. 0) THEN
!     FILENAME = 'almostdegen_read.input'
!     CALL READIN(FILENAME, M_BAR, K_BAR, EIG_BAR_IND, M_1, K_VAL_1_0)
! END IF ! MPI_GLB_RANK.EQ.0

! calculate operator matrix and its corresponding eigenvalues
ALLOCATE(K_RANGE(201))
CALL LINSPACE(-10.d0,10.d0,K_RANGE)

DO M_BAR = 1,3
! DO KK = 1,SIZE(K_RANGE)

    K_BAR = 0.d0 !K_RANGE(KK)
    ! K_BAR = 1.64
    ! DO K_ind = -1,1,2
    ! K_BAR = K_ind*2.2608213943428
    
    WRITE(K_VAL,'(F06.2)') K_BAR
    IF (MPI_GLB_RANK.EQ.0) WRITE(*,*) MPI_GLB_RANK,':',M_BAR,K_BAR

    ! CALL EIG_MATRIX_SERIAL(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm, print_switch=.true.)
    CALL EIG_MATRIX(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm, serial_switch=.false.) 

    IF (MPI_GLB_RANK .NE. 0) THEN
        GOTO 129
    ELSE
        ! SAVE EIGVALS
        ! open(FID,FILE='./converg/CriticalLayer_240618/qvortex_0.1/eig_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
        !         //'_NRCHOP_'//ITOA3(NRCHOP)//'.output',STATUS='unknown',ACTION='WRITE',IOSTAT=IERR)
        open(FID,FILE='./qvortex/h0/eig_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
                //'_NRCHOP_'//ITOA3(NRCHOP)//'.output',STATUS='unknown',ACTION='WRITE',IOSTAT=IERR)
        DO II = 1,SIZE(M_eig)
            ! WRITE(FID,*) REAL(M_eig(II)),',',AIMAG(M_eig(II)),',',EIGRES(EIG_R_mat(:,II),M_BAR)
            WRITE(FID,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_mat(:,II),M_BAR)
            ! IF (REAL(M_eig(II)).GT.0.10) WRITE(*,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_mat(:,II),M_BAR)
        ENDDO
        close(FID)

        ! SAVE EIGVECS
        ALLOCATE(EIG_R_BAR(SIZE(EIG_R_mat,1)))
        DO II = 1,SIZE(M_eig)
            ! II = 688
            EIG_R_BAR = EIG_R_mat(:,II)

            ! WRITE(*,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_BAR,M_BAR)
            CALL EIG2VEL(M_BAR, K_BAR, EIG_R_BAR, RUR_BAR, RUP_BAR, UZ_BAR, comm_grp=newcomm)
            ! call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './converg/CriticalLayer_240605/vel_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
            !             //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
            call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './qvortex/h0/vel_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
                        //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
            ! ! call SAVE_VEC('./converg/CriticalLayer_240605/PC_o_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
            !             ! //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R_BAR)
            DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)

            ! IF (EIGRES(EIG_R_BAR,M_BAR)) THEN
            ! ! IF ((EIGRES(EIG_R_BAR)).OR.(II.LE.4)) THEN
            !     ! WRITE(*,*) 'II = ',II,':',M_eig(II)
            !     CALL EIG2VEL(M_BAR, K_BAR, EIG_R_BAR, RUR_BAR, RUP_BAR, UZ_BAR, comm_grp=newcomm)
            !     call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './converg/check5/vel_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
            !                 //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
            !     call SAVE_VEC('./converg/check5/PC_o_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
            !                 //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R_BAR)
            !     DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)
            ! ENDIF
        ENDDO
        DEALLOCATE(EIG_R_BAR,M_eig)
    ENDIF

    ! DO JJ = 0,MPI_GLB_PROCS-1
    !     IF (MPI_GLB_RANK .EQ. JJ) THEN
    !         IF (MPI_GLB_RANK .NE. 0) GOTO 129

    !         ! open(FID,FILE='./converg/check5/eig_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !         !     //'_NRCHOP_'//ITOA3(NRCHOP)//'.output',STATUS='unknown',ACTION='WRITE',IOSTAT=IERR)
    !         open(FID,FILE='./converg/CriticalLayer_240618/eig_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !             //'_NRCHOP_'//ITOA3(NRCHOP)//'.output',STATUS='unknown',ACTION='WRITE',IOSTAT=IERR)
    !         DO II = 1,SIZE(M_eig)
    !             ! WRITE(FID,*) REAL(M_eig(II)),',',AIMAG(M_eig(II)),',',EIGRES(EIG_R_mat(:,II),M_BAR)
    !             WRITE(FID,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_BAR,M_BAR)
    !         ENDDO
    !         close(FID)

    !         ! EIG2VEL:
    !         ALLOCATE(EIG_R_BAR(SIZE(EIG_R_mat,1)))
    !         DO II = 1,SIZE(M_eig)
    !             ! II = 688
    !             EIG_R_BAR = EIG_R_mat(:,II)

    !             WRITE(*,*) 'II = ',II,':',M_eig(II),'-',EIGRES(EIG_R_BAR,M_BAR)
    !             CALL EIG2VEL(M_BAR, K_BAR, EIG_R_BAR, RUR_BAR, RUP_BAR, UZ_BAR, comm_grp=newcomm)
    !             call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './converg/CriticalLayer_240605/vel_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !                         //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    !             ! call SAVE_VEC('./converg/CriticalLayer_240605/PC_o_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !                         ! //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R_BAR)
    !             DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)

    !             ! IF (EIGRES(EIG_R_BAR,M_BAR)) THEN
    !             ! ! IF ((EIGRES(EIG_R_BAR)).OR.(II.LE.4)) THEN
    !             !     ! WRITE(*,*) 'II = ',II,':',M_eig(II)
    !             !     CALL EIG2VEL(M_BAR, K_BAR, EIG_R_BAR, RUR_BAR, RUP_BAR, UZ_BAR, comm_grp=newcomm)
    !             !     call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './converg/check5/vel_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !             !                 //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    !             !     call SAVE_VEC('./converg/check5/PC_o_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !             !                 //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R_BAR)
    !             !     DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)
    !             ! ENDIF
    !             ! IF ((VISC%SW.EQ.1).AND.(REAL(M_eig(II)).GT.-VISC%NU)) THEN
    !             ! ! IF ((EIGRES(EIG_R_BAR)).OR.(II.LE.4)) THEN
    !             !     WRITE(*,*) 'II = ',II,':',M_eig(II)
    !             !     CALL EIG2VEL(M_BAR, K_BAR, EIG_R_BAR, RUR_BAR, RUP_BAR, UZ_BAR, comm_grp=newcomm)
    !             !     call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './converg/check5/vel_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !             !                 //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    !             !     call SAVE_VEC('./converg/check5/PC_o_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !             !                 //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R_BAR)
    !             !     DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)
    !             ! ENDIF
    !         ENDDO
    !         DEALLOCATE(EIG_R_BAR,M_eig)
            
    !         ! EIG2VELVOR:
    !         ! DO II = 1,30!size(M_eig,1)

    !         ! WRITE(*,*) II
    !         ! CALL EIG2VELVOR(M_BAR,K_BAR,II,EIG_R_mat,EIG_L_mat,RUR_BAR, RUP_BAR, UZ_BAR, ROR_BAR, ROP_BAR, OZ_BAR, EIG_R_BAR, EIG_L_BAR, comm_grp=newcomm)
    !         ! call SAVE_VEL(RUR_BAR, RUP_BAR, UZ_BAR, './converg/check2/vel_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !         !             //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    !         ! call SAVE_VEL(ROR_BAR, ROP_BAR, OZ_BAR, './converg/check2/vor_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !         !             //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
    !         ! call SAVE_VEC('./converg/check2/vec_MK_'//ITOA3(M_BAR)//'_'//K_VAL &
    !         !             //'_IND_'//ITOA3(II)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R_BAR, EIG_L_BAR)
    !         ! WRITE(*,*) 'EVP_PRINT: NORMALIZATION = ', DOT_PRODUCT(EIG_L_BAR,EIG_R_BAR)
    !         ! DEALLOCATE(EIG_R_BAR,EIG_L_BAR)
    !         ! DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)
    !         ! DEALLOCATE(ROR_BAR, ROP_BAR, OZ_BAR)
    !         ! ENDDO

    !         ! DEALLOCATE(M_eig)

    !     ENDIF
129     CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    ! ENDDO

    ! ENDDO
ENDDO
! ENDDO

! calcualte NONLIN
! II = 25;
! DO M_BAR = -1,1,2

!     K_BAR = 2.27*M_BAR
    

!     CALL EIG_MATRIX(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat, comm_grp=newcomm, print_switch=.true.) 
!     IF (MPI_GLB_RANK.EQ.0) THEN
!         CALL EIG2VELVOR(M_BAR,K_BAR,II,EIG_R_mat,EIG_L_mat,RUR_BAR, RUP_BAR, UZ_BAR, ROR_BAR, ROP_BAR, OZ_BAR, EIG_R_BAR, EIG_L_BAR, comm_grp=newcomm)

!         KK = SIZE(RUR_BAR)
!         ALLOCATE(RUR_BAR2(KK), RUP_BAR2(KK), UZ_BAR2(KK), ROR_BAR2(KK), ROP_BAR2(KK), OZ_BAR2(KK))
!         RUR_BAR2(1:NR) = RUR_BAR(1:NR); RUP_BAR2 = RUP_BAR(1:NR); UZ_BAR2 = UZ_BAR(1:NR);
!         ROR_BAR2(1:NR) = ROR_BAR(1:NR); ROP_BAR2 = ROP_BAR(1:NR); OZ_BAR2 = OZ_BAR(1:NR);

!         DEALLOCATE(EIG_R_BAR)
!         CALL NONLIN_MK(M_BAR*2,K_BAR*2,RUR_BAR, RUP_BAR, UZ_BAR, ROR_BAR, ROP_BAR, OZ_BAR, &
!             RUR_BAR2, RUP_BAR2, UZ_BAR2, ROR_BAR2, ROP_BAR2, OZ_BAR2, EIG_R_BAR, comm_grp=newcomm)
            
!         WRITE(K_VAL,'(F05.2)') K_BAR*2
!         call SAVE_VEL(ROR_BAR2, ROP_BAR2, OZ_BAR2, './converg/check2/vor_MK_'//ITOA3(M_BAR*2)//'_'//K_VAL &
!                     //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output')
!         call SAVE_VEC('./converg/check2/vec_MK_'//ITOA3(M_BAR*2)//'_'//K_VAL &
!                     //'_IND_'//ITOA3(0)//'_NRCHOP_'//ITOA3(NRCHOP)//'.output',EIG_R_BAR)

!         DEALLOCATE(EIG_R_BAR,EIG_L_BAR,M_eig)
!         DEALLOCATE(RUR_BAR, RUP_BAR, UZ_BAR)
!         DEALLOCATE(ROR_BAR, ROP_BAR, OZ_BAR)
!         DEALLOCATE(RUR_BAR2, RUP_BAR2, UZ_BAR2)
!         DEALLOCATE(ROR_BAR2, ROP_BAR2, OZ_BAR2)
!     ENDIF

! ENDDO

IF (ALLOCATED(EIG_R_mat)) DEALLOCATE(EIG_R_mat)
IF (ALLOCATED(EIG_L_mat)) DEALLOCATE(EIG_L_mat)

IF (MPI_GLB_RANK .EQ. 0) THEN
    WRITE (*, *) ''
    WRITE (*, *) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
END IF

CALL MPI_FINALIZE(IERR)

END PROGRAM EVP_PRINT
!=======================================================================
