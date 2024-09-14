PROGRAM EVP_SPARSITY
!=======================================================================
! [USAGE]:  
! SAVE EIGENVALUES FOR DIFFERENT SETS OF {M,AK}
!
! [UPDATES]:
! 4/12/2022
! 9/12/2022
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

IMPLICIT NONE
INTEGER     :: II,JJ,KK
INTEGER     :: M_BAR,K_BAR_INT
REAL(P8)    :: K_BAR
REAL(P8),DIMENSION(:),ALLOCATABLE:: K_tune, K_tune_JJ, tol_diff_JJ
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat

! SCANNING:
CHARACTER(LEN=72):: FILENAME
INTEGER:: FID,IS

CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
    WRITE(*,*) 'The threading support is lesser than that demanded.'
    CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
ENDIF

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
IF (MPI_GLB_RANK.EQ.0) THEN
    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()   ! @ MOD_MISC
ENDIF
CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR) ! DIVIDE NR BY #OF MPI PROCS
CALL READCOM('NOECHO')
CALL READIN(5)

DO M_BAR = -10,10
    ! IF (M_BAR.GT.0) GOTO 82
    K_BAR = 20.D0
    ! DO WHILE (K_BAR.LE.20.D0)

        IF ((M_BAR.EQ.0).AND.(K_BAR.EQ.0)) CYCLE

        ! calculate operator matrix and its corresponding eigenvalues
        CALL EIG_MATRIX(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat, EIG_VEC_L = EIG_L_mat,comm_grp=newcomm)
        IF (MPI_GLB_RANK .EQ. 0) CALL SAVE_MODE(M_BAR, K_BAR, M_mat, M_eig, EIG_R_mat, EIG_L_mat)
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)

        ! ! save file
        ! FILENAME = 'EIGVAL_'//ITOA3(M_BAR)//'_'//ITOA3(K_BAR_INT)
        ! IF (MPI_GLB_RANK.EQ.0) THEN
        !     open(10,FILE='./converg/eigval/'//TRIM(ADJUSTL(FILENAME))//'.input',STATUS='unknown',ACTION='WRITE',IOSTAT=IS)
        !     if (IS.ne.0) then
        !         print *, 'ERROR: Could not creat eigen-value file'
        !     end if
        !     DO II = 1,SIZE(M_eig)
        !         WRITE(10,*) real(M_eig(II)),',',aimag(M_eig(II))
        !     ENDDO
        !     M_eig = 0.D0
        !     close(10)
        ! ENDIF

        K_BAR = K_BAR + 0.1

    ! ENDDO
ENDDO

82 CONTINUE
! DEALLOCATE
IF (MPI_GLB_RANK.EQ.0) THEN
    DEALLOCATE(M_eig)
ENDIF
DEALLOCATE(M_mat)

CALL MPI_FINALIZE(IERR)

END PROGRAM EVP_SPARSITY
!=======================================================================
