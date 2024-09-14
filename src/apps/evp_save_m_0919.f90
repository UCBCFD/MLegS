PROGRAM EVP_SAVE_M_0919
!=======================================================================
! [USAGE]:  
! SAVE EIGENVALUES FOR DIFFERENT SETS OF {M,AK}
!
! [UPDATES]:
! 4/12/2022
! 9/12/2022
! 9/19/2022
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
INTEGER     :: M_BAR,K_BAR_COUNT
REAL(P8)    :: K_BAR,K_BAR_START,K_BAR_END
REAL(P8),DIMENSION(:),ALLOCATABLE:: K_BAR_RANGE
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat

! SCANNING:
CHARACTER(LEN=72):: FILENAME
INTEGER:: FID,IS
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_eig_M
LOGICAL,DIMENSION(:,:),ALLOCATABLE:: M_res_M

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

! DO M_BAR = 1,2
M_BAR = MPI_GLB_RANK !-1

    ! create a range for K_BAR
    K_BAR_START = -10.D0 !0.D0
    K_BAR_END   = 10.D0
    ! K_BAR_COUNT = (K_BAR_END-K_BAR_START)/1E-2 + 1
    K_BAR_COUNT = (K_BAR_END-K_BAR_START)/1E-1 + 1

    IF (ALLOCATED(K_BAR_RANGE)) DEALLOCATE(K_BAR_RANGE)
    ALLOCATE(K_BAR_RANGE(K_BAR_COUNT))
    CALL LINSPACE(K_BAR_START,K_BAR_END,K_BAR_RANGE)

    ! DEBUG:
    ! IF (MPI_GLB_RANK.EQ.0) WRITE(*,*) K_BAR_RANGE
    WRITE(*,*) MPI_GLB_RANK,':',M_BAR,K_BAR_START,K_BAR_END

    ! save eig_val and mode resolving info
    DO II = 1,K_BAR_COUNT

        K_BAR = K_BAR_RANGE(II)
        IF ((M_BAR.EQ.0).AND.(K_BAR.EQ.0)) THEN
            M_eig_M(:,II) = 0.D0
            M_res_M(:,II) = .FALSE.
            CYCLE
        ENDIF

        ! calculate operator matrix and its corresponding eigenvalues
        CALL EIG_MATRIX(M_BAR, K_BAR, M_mat, M_eig, EIG_VEC_R = EIG_R_mat,comm_grp=newcomm,serial_switch=.true.)

        ! IF (MPI_GLB_RANK .EQ. 0) THEN
            ! create M_eig_M and M_res_M
            IF (.NOT.(ALLOCATED(M_eig_M))) ALLOCATE(M_eig_M(SIZE(M_eig),K_BAR_COUNT))
            IF (.NOT.(ALLOCATED(M_res_M))) ALLOCATE(M_res_M(SIZE(M_eig),K_BAR_COUNT))

            ! update M_eig_M(:,II) and M_res_M(:,II)
            M_eig_M(:,II) = M_eig(:)
            M_res_M(:,II) = EIGRES(EIG_R_mat,M_BAR)

            WRITE(*,*) 'RANK#',MPI_GLB_RANK,' PROGRESS: ',ITOA4(II),'/',ITOA4(K_BAR_COUNT)
        ! ENDIF

        ! CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        ! CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
    ENDDO

    ! save M_eig_M and M_res_M (e.g. EIG_0003.dat)
    ! IF (MPI_GLB_RANK .EQ. 0) THEN

        ! FILENAME = 'EIG_'//ITOA4(M_BAR)//'.dat'
        FILENAME = RTOA6(QPAIR%H(1))//'_EIG_'//ITOA4(M_BAR)//'.dat'
        OPEN(UNIT=7,FILE='./data/'//TRIM(FILENAME),STATUS='UNKNOWN',FORM='UNFORMATTED',IOSTAT=IS)

        IF(IS.NE.0) THEN
            WRITE(*,*) 'SAVE_MODE_M: FAILED TO OPEN ',TRIM(FILENAME)
            STOP
        ENDIF

        WRITE(7) M_BAR,K_BAR_START,K_BAR_END
        WRITE(7) SIZE(M_eig_M,1),K_BAR_COUNT
        WRITE(7) M_eig_M(:,:)
        WRITE(7) M_res_M(:,:)
        CLOSE(7)

        DEALLOCATE(M_eig_M,M_res_M,K_BAR_RANGE)

        ! DEBUG:
        ! CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)

    ! ENDIF
! ENDDO

82 CONTINUE
! DEALLOCATE
! IF (MPI_GLB_RANK.EQ.0) THEN
    DEALLOCATE(M_eig)
! ENDIF
DEALLOCATE(M_mat)

CALL MPI_FINALIZE(IERR)

! ======================================================================
CONTAINS
! ======================================================================
SUBROUTINE READ_MODE_M(M_READ, K_RANGE, EIG_MAT, RES_MAT)
! ======================================================================

! ======================================================================
IMPLICIT NONE
INTEGER:: M_READ, M_FILE
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: EIG_MAT
LOGICAL,DIMENSION(:,:),ALLOCATABLE:: RES_MAT

INTEGER:: K_COUNT, EIG_SIZE
REAL(P8):: K_START,K_END
REAL(P8),DIMENSION(:),ALLOCATABLE:: K_RANGE

CHARACTER(LEN=72):: FN
INTEGER:: STATUS

! open file
FN = 'EIG_'//ITOA4(M_READ)//'.dat'
OPEN(UNIT=7,FILE='./data/'//TRIM(FN),STATUS='OLD',FORM='UNFORMATTED',ACTION='READ',IOSTAT=STATUS)
IF (STATUS.NE.0) THEN
    WRITE(*,*) 'READ_MODE_M: FAILED TO OPEN ',TRIM(FN),' - RANK#',MPI_GLB_RANK
    STOP
ENDIF

READ(7) M_FILE,K_START,K_END
READ(7) EIG_SIZE,K_COUNT

IF (M_FILE.NE.M_READ) THEN
    WRITE(*,*) 'READ_MODE_M: M_READ = ', M_READ,' .NE. M_FILE = ', M_FILE,' - RANK#',MPI_GLB_RANK
    STOP
ENDIF

! obtain K_RANGE
IF (ALLOCATED(K_RANGE)) DEALLOCATE(K_RANGE)
ALLOCATE(K_RANGE(K_COUNT))
CALL LINSPACE(K_START,K_END,K_RANGE)

! read EIG_MAT and 
IF (ALLOCATED(EIG_MAT)) DEALLOCATE(EIG_MAT)
IF (ALLOCATED(RES_MAT)) DEALLOCATE(RES_MAT)
ALLOCATE(EIG_MAT(EIG_SIZE,K_COUNT))
ALLOCATE(RES_MAT(EIG_SIZE,K_COUNT))

READ(7) EIG_MAT(:,:)
READ(7) RES_MAT(:,:)
CLOSE(7)

END SUBROUTINE
! ======================================================================

END PROGRAM EVP_SAVE_M_0919
!=======================================================================
