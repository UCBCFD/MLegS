PROGRAM EVP_SAVE_M_0_MOD
!=======================================================================
! [USAGE]:  
! FOR M = 0, MARK EIGVALS THAT ARE STACKED CLOSE TO EACH OTHER AS "UNRES
! OLVED". THE CURRENT THRESHOLD IS 1.5E-3.
!
! THIS IS A SERIAL PROGRAM!
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
LOGICAL,DIMENSION(:),ALLOCATABLE:: M_res
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat, EIG_R_mat, EIG_L_mat

! SCANNING:
CHARACTER(LEN=72):: FILENAME
INTEGER:: FID,IS
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_eig_M
LOGICAL,DIMENSION(:,:),ALLOCATABLE:: M_res_M, M_res_M_org

CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
    WRITE(*,*) 'The threading support is lesser than that demanded.'
    CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
ENDIF

IF (MPI_GLB_RANK.EQ.0) THEN
    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()   ! @ MOD_MISC
ENDIF

CALL READCOM('NOECHO')
CALL READIN(5)

IF (MPI_GLB_RANK.EQ.0) THEN
    WRITE(*,*) 'M_BAR = '
    READ(5,'(I72)') M_BAR
ENDIF
CALL MPI_BCAST(M_BAR,1,MPI_INT,0,MPI_COMM_WORLD,IERR)

CALL READ_MODE_M(M_BAR, K_BAR_RANGE, M_eig_M, M_res_M)
K_BAR_COUNT = SIZE(K_BAR_RANGE)

! save eig_val and mode resolving info
ALLOCATE(M_res_M_org(size(M_res_M,1),size(M_res_M,2)))
M_res_M_org = M_res_M
DO II = 1,K_BAR_COUNT

    K_BAR = K_BAR_RANGE(II)

    DO JJ = 1,SIZE(M_eig_M(:,II))

        ! skip k = 0 modes
        IF (K_BAR.EQ.0.D0) THEN
            M_res_M(JJ,II) = .FALSE.
            CYCLE
        ENDIF

        ! skip unresolved mode
        IF (.NOT.(M_res_M(JJ,II))) CYCLE

        ALLOCATE(M_res(SIZE(M_eig_M(:,II))))
        M_res = M_res_M_org(:,II); M_res(JJ) = .FALSE.

        ! if the mode's eigval is close to other resolved eigvals, remove that mode
        IF (MINVAL(ABS(M_eig_M(JJ,II) - M_eig_M(:,II)),M_res).LE.1E-3) M_res_M(JJ,II) = .FALSE.

        DEALLOCATE(M_res)

    ENDDO

ENDDO

! save M_eig_M and M_res_M (e.g. EIG_0003.dat)
! IF (MPI_GLB_RANK .EQ. 0) THEN

! FILENAME = 'EIG_'//ITOA4(M_BAR)//'.dat'
FILENAME = RTOA6(QPAIR%H(1))//'_EIG_'//ITOA4(M_BAR)//'_mod.dat'
OPEN(UNIT=7,FILE='./data/'//TRIM(FILENAME),STATUS='UNKNOWN',FORM='UNFORMATTED',IOSTAT=IS)

IF(IS.NE.0) THEN
    WRITE(*,*) 'SAVE_MODE_M: FAILED TO OPEN ',TRIM(FILENAME)
    STOP
ENDIF

WRITE(7) M_BAR,K_BAR_RANGE(1),K_BAR_RANGE(K_BAR_COUNT)
WRITE(7) SIZE(M_eig_M,1),K_BAR_COUNT
WRITE(7) M_eig_M(:,:)
WRITE(7) M_res_M(:,:)
CLOSE(7)

DEALLOCATE(M_eig_M,M_res_M,M_res_M_org,K_BAR_RANGE)


82 CONTINUE

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
FN = RTOA6(QPAIR%H(1))//'_EIG_'//ITOA4(M_READ)//'.dat'
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

END PROGRAM EVP_SAVE_M_0_MOD
!=======================================================================
