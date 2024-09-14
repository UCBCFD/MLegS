PROGRAM EVP_DISPERSION_M
!=======================================================================
! [USAGE]: 
! For a given M, save the IM(eig) as a function of k
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
INTEGER :: II, JJ, KK
! SCAN:
INTEGER :: M_BAR
REAL(P8),DIMENSION(:),ALLOCATABLE:: K_BAR_RANGE
LOGICAL,    DIMENSION(:,:),ALLOCATABLE:: M_res_bar_mat
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_eig_bar_mat

! LOG:
INTEGER:: IS
CHARACTER(LEN=72):: FILENAME
INTEGER:: FID

QPAIR%H = 0.4

! M_BAR = 1
WRITE(*,*) 'M_BAR = '
READ(5,'(I72)') M_BAR
CALL READ_MODE_M(M_BAR, K_BAR_RANGE, M_eig_bar_mat, M_res_bar_mat)

! IF (.NOT.(M_res_bar_mat(EIG_BAR_INDEX,K_BAR_INDEX))) CYCLE loop_ind
open(unit=FID,FILE='eigres_M_'//ITOA3(M_BAR)//'.output',STATUS='unknown',ACTION='WRITE',IOSTAT=IS)

DO II = 1,SIZE(M_eig_bar_mat,2)

    WRITE(FID,'(F20.13)',ADVANCE="NO") K_BAR_RANGE(II)
    
    DO JJ = 1,SIZE(M_eig_bar_mat,1)
        IF (M_res_bar_mat(JJ,II)) WRITE(FID,47,ADVANCE="NO") AIMAG(M_eig_bar_mat(JJ,II))
    ENDDO

    WRITE(FID,'')

ENDDO
47 FORMAT(',', 1P8E20.13)

DEALLOCATE(M_res_bar_mat,M_eig_bar_mat,K_BAR_RANGE)

! ======================================================================
CONTAINS
! ======================================================================
SUBROUTINE READ_MODE_M(M_READ, K_RANGE, EIG_MAT, RES_MAT)
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
FN = RTOA6(QPAIR%H(1))//'_EIG_'//ITOA4(M_READ)//'_mod.dat'
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

END PROGRAM EVP_DISPERSION_M
!=======================================================================
