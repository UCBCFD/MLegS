program addperturb_split
! J. Barranco, 01/25/1999
! J. Wang, 03/30/2023
! This program splits the full perturbation file which has 3
! eigenfunctions into two separate files for perturbation and degen pair
! separately
! ======================================================================
   USE OMP_LIB
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
   USE MOD_DIAGNOSTICS
   USE MOD_INIT
! ======================================================================
IMPLICIT NONE
! VARIABLE DECLARATION
INTEGER:: I,J,JJ,ERR,NN
REAL(P8):: KK, RP1, IP1, RP2, IP2
INTEGER:: MM
COMPLEX(P8), ALLOCATABLE:: PERTURB_PSI(:,:), PERTURB_CHI(:,:)
COMPLEX(P8), ALLOCATABLE:: SCALE(:)
! FILE
INTEGER:: FILENUM
CHARACTER(LEN=200):: FILENAME, LINE
! ======================================================================

WRITE(*,*) 'PROGRAM STARTED'
CALL PRINT_REAL_TIME()
WRITE(*,*) 'Perturbation filename: '
READ(*,10) FILENAME
10 FORMAT(A200)

! original file
OPEN(10,FILE=TRIM(ADJUSTL(FILENAME))//".input",STATUS='OLD',ACTION='READ',IOSTAT=err)
IF (ERR.NE.0) THEN
   PRINT *, 'ERROR: addperturb -- Could not open file: ',TRIM(ADJUSTL(FILENAME))//".input"
   STOP
END IF

! copy
OPEN(20,FILE=TRIM(ADJUSTL(FILENAME))//"0.input",ACTION='WRITE',IOSTAT=err) ! perturbation file
OPEN(30,FILE=TRIM(ADJUSTL(FILENAME))//"1.input",ACTION='WRITE',IOSTAT=err) ! degen pair file

! number of modes
WRITE(20,'(I2)') 1; WRITE(30,'(I2)') 2

! scaling
READ(10,*)
DO J=1,3
    FILENUM = 30
    IF (J.EQ.1) FILENUM = 20;

    READ(10,*) RP1, IP1
    WRITE(FILENUM,'(F9.7,1X,F9.7)') RP1, IP1
END DO

! copy contents
DO J=1,3
    FILENUM = 30
    IF (J.EQ.1) FILENUM = 20;

    ! header
    READ(10,'(A)') LINE; WRITE(FILENUM,'(A)') trim(LINE)
    READ(10,'(A)') LINE; WRITE(FILENUM,'(A)') trim(LINE)
    READ(10,*) NN, MM, KK
    BACKSPACE 10
    READ(10,'(A)') LINE; WRITE(FILENUM,'(A)') trim(LINE)
    ! WRITE(FILENUM,74) NN, MM, KK
74 FORMAT(2(12X,I12),12X,F12.8)
    DO JJ=1,8
        READ(10,'(A)') LINE; WRITE(FILENUM,'(A)') trim(LINE)
    END DO

    ! coeffs
    DO I=1,NN
        READ(10,'(A)') LINE; WRITE(FILENUM,'(A)') trim(LINE)
    END DO
END DO

! finalize
close(10); close(20); close(30)
WRITE(*,*) 'PROGRAM ENDED'
CALL PRINT_REAL_TIME()

END PROGRAM ADDPERTURB_SPLIT



