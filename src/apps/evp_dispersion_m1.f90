PROGRAM EVP_DISPERSION_M1
!=======================================================================
! [USAGE]:  
! CREATE DISPERSION CURVES OF THE Q-VORTEX AT |M|=1
!
! [UPDATES]:
! 6/29/2022
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
INTEGER     :: II,JJ,KK, N_K
INTEGER     :: M_BAR,K_BAR_INT
REAL(P8)    :: K_BAR,K_BAR_END
! REAL(P8),DIMENSION(:),ALLOCATABLE::  K_BAR_RANGE
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: M_eig
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat

! SCANNING:
CHARACTER(len=5) :: K_VAL
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

N_K = 100
K_BAR_END = 10.D0

! DO M_BAR = -1,1,2
!     DO K_BAR_INT = 0,N_K
M_BAR = 0
    DO K_BAR_INT = 1,5

        IF ((K_BAR_INT.NE.1).AND.(K_BAR_INT.NE.2).AND.(K_BAR_INT.NE.5)) CYCLE

        ! K_BAR = (K_BAR_INT)*K_BAR_END/N_K
        K_BAR = K_BAR_INT + 2

        WRITE(K_VAL,'(F05.2)') K_BAR

        FILENAME = 'EIGVAL_'//ITOA3(M_BAR)//'_'//K_VAL

        ! calculate operator matrix and its corresponding eigenvalues
        CALL EIG_MATRIX(M_BAR,K_BAR,M_mat,M_eig,comm_grp=newcomm)
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

        ! save file
        IF (MPI_GLB_RANK.EQ.0) THEN
            open(10,FILE='./converg/eigval3/'//TRIM(ADJUSTL(FILENAME))//'.input',STATUS='unknown',ACTION='WRITE',IOSTAT=IS)
            if (IS.ne.0) then
                print *, 'ERROR: Could not creat eigen-value file'
            end if
            DO II = 1,SIZE(M_eig)
                WRITE(10,*) real(M_eig(II)),',',aimag(M_eig(II))
            ENDDO
            ! write(*,*) aimag(M_eig)
            M_eig = 0.D0
            close(10)
        ENDIF

    ENDDO
! ENDDO

82 CONTINUE
! DEALLOCATE
IF (MPI_GLB_RANK.EQ.0) THEN
    DEALLOCATE(M_eig)
ENDIF
DEALLOCATE(M_mat)

CALL MPI_FINALIZE(IERR)

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================


END PROGRAM EVP_DISPERSION_M1
!=======================================================================
