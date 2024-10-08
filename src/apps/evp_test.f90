PROGRAM EVP_SPARSITY
!=======================================================================
! [USAGE]:  
! SAVE EIGENVALUES FOR DIFFERENT SETS OF {M,AK}
!
! [UPDATES]:
! 4/12/2022
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
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: M_mat
REAL(P8),DIMENSION(:),ALLOCATABLE:: RUR0,RUP0,UZ0,ROR0,ROP0,OZ0

! MPI:
INTEGER:: MPI_GLB_PROCS, MPI_GLB_RANK, newcomm
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

M_BAR = 0
K_BAR_INT = 1

FILENAME = 'EIGVAL_'//ITOA3(M_BAR)//'_'//ITOA3(K_BAR_INT)
K_BAR = K_BAR_INT

! calculate operator matrix and its corresponding eigenvalues
CALL EIG_MATRIX(M_BAR,K_BAR,M_mat,M_eig,comm_grp=newcomm)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

! save file
IF (MPI_GLB_RANK.EQ.0) THEN
    open(10,FILE='./converg/eigval/'//TRIM(ADJUSTL(FILENAME))//'.input',STATUS='unknown',ACTION='WRITE',IOSTAT=IS)
    if (IS.ne.0) then
        print *, 'ERROR: Could not creat eigen-value file'
    end if
    DO II = 1,SIZE(M_eig)
        WRITE(10,*) real(M_eig(II)),',',aimag(M_eig(II))
    ENDDO
    M_eig = 0.D0
    close(10)
ENDIF

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
SUBROUTINE EIG_MATRIX(MREAD,AKREAD,H,EIG_VAL,EIG_VEC_R,EIG_VEC_L,comm_grp)
    IMPLICIT NONE
    INTEGER,INTENT(IN)    :: MREAD, comm_grp
    REAL(P8),INTENT(INOUT)    :: AKREAD
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: H
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE,INTENT(INOUT),OPTIONAL:: EIG_VAL
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT),OPTIONAL:: EIG_VEC_R, EIG_VEC_L

    INTEGER     :: NR_MK, I
    REAL(P8)    :: REI
    INTEGER     :: MPI_NR_SIZE,MPI_NR_INDEX
    ! REAL(P8),DIMENSION(:),ALLOCATABLE:: RUR0,RUP0,UZ0,ROR0,ROP0,OZ0
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: RURU,RUPU,UZU,RORU,ROPU,OZU
    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: PSIU,CHIU,PSI1,CHI1,PSI2,CHI2,PSI3,CHI3

    ZLEN = 2*PI/AKREAD
    CALL MPI_BCAST(MREAD,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(ZLEN,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
    
    ! INIT 
    NTH = 2
    NX = 4
    NTCHOP = 2 ! M = 0, M = MREAD
    NXCHOP = 2 ! K = 0, K = AKREAD
    CALL LEGINIT(comm_grp,MREAD)
    IF (.NOT.ALLOCATED(RUR0)) THEN
        ALLOCATE(RUR0(NR),RUP0(NR),UZ0(NR),ROR0(NR),ROP0(NR),OZ0(NR))
        CALL INIT_LOOP(RUR0,RUP0,UZ0,ROR0,ROP0,OZ0) ! ALL IN PFF space
    ENDIF

    ! SPLIT NRCHOPDIM
    NR_MK = NRCHOPS(2)
    IF (ALLOCATED(H).AND.(SIZE(H,1).NE.NR_MK)) THEN
        DEALLOCATE(H)
        ALLOCATE(H(2*NR_MK,2*NR_MK))
    ELSEIF (.NOT.(ALLOCATED(H))) THEN
        ALLOCATE(H(2*NR_MK,2*NR_MK))
    ENDIF 
    H = CMPLX(0.D0,0.D0,P8)
    CALL DECOMPOSE(NR_MK, MPI_GLB_PROCS, MPI_GLB_RANK, MPI_NR_SIZE, MPI_NR_INDEX)

    ! MAIN JOB
    ALLOCATE(PSIU(NRCHOPDIM),CHIU(NRCHOPDIM))
    ALLOCATE(PSI1(NRCHOPDIM),CHI1(NRCHOPDIM))
    ALLOCATE(PSI2(NRCHOPDIM),CHI2(NRCHOPDIM))
    ALLOCATE(PSI3(NRCHOPDIM),CHI3(NRCHOPDIM))
    ALLOCATE(RURU(NRCHOPDIM),RUPU(NRCHOPDIM),UZU(NRCHOPDIM))
    ALLOCATE(RORU(NRCHOPDIM),ROPU(NRCHOPDIM),OZU(NRCHOPDIM))

    DO I = MPI_NR_INDEX+1,MPI_NR_INDEX+MPI_NR_SIZE
    ! ================================ PSI =================================
        PSIU = 0.D0
        CHIU = 0.D0
        PSIU(I) = 1.D0

        ! VISCOSITY/HYPERV
        CALL CHOPSET(2)
        IF (VISC%SW .EQ. 1) THEN
        CALL DEL2_MK(PSIU,PSI3)
        PSI3 = PSI3 * VISC%NU
        ELSE IF (VISC%SW .EQ. 2) THEN
        CALL CHOPSET(-2+VISC%P)
        CALL HELMP_MK(VISC%P,PSIU,PSI3,0.D0)
        IF (MOD(VISC%P/2,2).EQ.0) PSI3 = -PSI3
        PSI3 = PSI3 * VISC%NUP
        CALL CHOPSET( 2-VISC%P)
        ELSE
        PSI3 = 0.D0
        ENDIF
        CHI3 = 0.D0
        CALL CHOPSET(-2)

        CALL CHOPSET(3)
        ! NONLINEAR TERM: W0 X U'
        CALL PC2VEL_MK(PSIU,CHIU,RURU,RUPU,UZU)
        CALL RTRAN_MK(RURU,1)
        CALL RTRAN_MK(RUPU,1)
        CALL RTRAN_MK(UZU,1)
        CALL VPROD_MK(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        CALL PROJECT_MK(RURU,RUPU,UZU,PSI1,CHI1)

        ! NONLINEAR TERM: W' X U0
        CALL PC2VOR_MK(PSIU,CHIU,RORU,ROPU,OZU)
        CALL RTRAN_MK(RORU,1)
        CALL RTRAN_MK(ROPU,1)
        CALL RTRAN_MK(OZU,1)
        CALL VPROD_MK(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        CALL PROJECT_MK(RORU,ROPU,OZU,PSI2,CHI2)
        CALL CHOPSET(-3)
        ! IF (I.EQ.2) THEN
        !     CALL XXDX_MK(PSIU,RURU)
        !     CALL MCAT(ROP0)
        ! ENDIF

        PSI1 = -PSI1 + PSI2 + PSI3
        CHI1 = -CHI1 + CHI2 + CHI3
        PSI1(NRCHOPS(2)+1:) = 0.D0
        CHI1(NRCHOPS(2)+1:) = 0.D0

        H(:NR_MK,  I) = PSI1(:NR_MK)
        H(NR_MK+1:,I) = CHI1(:NR_MK)

    ! ================================ CHI =================================
        PSIU = 0.D0
        CHIU = 0.D0
        CHIU(I) = 1.D0

        CALL IDEL2_MK(CHIU,CHI1)
        CHIU = CHI1
        CHI1 = 0.D0

        IF (VISC%SW .NE. 0) THEN
        CHI3 = PSI3
        ELSE
        CHI3 = 0.D0
        ENDIF
        PSI3 = 0.D0

        CALL CHOPSET(3)
        ! NONLINEAR TERM: W0 X U'
        CALL PC2VEL_MK(PSIU,CHIU,RURU,RUPU,UZU)
        CALL RTRAN_MK(RURU,1)
        CALL RTRAN_MK(RUPU,1)
        CALL RTRAN_MK(UZU,1)
        CALL VPROD_MK(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        CALL PROJECT_MK(RURU,RUPU,UZU,PSI1,CHI1)

        ! NONLINEAR TERM: W' X U0
        CALL PC2VOR_MK(PSIU,CHIU,RORU,ROPU,OZU)
        CALL RTRAN_MK(RORU,1)
        CALL RTRAN_MK(ROPU,1)
        CALL RTRAN_MK(OZU,1)
        CALL VPROD_MK(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        CALL PROJECT_MK(RORU,ROPU,OZU,PSI2,CHI2)
        CALL CHOPSET(-3)

        PSI1 = -PSI1 + PSI2 + PSI3
        CHI1 = -CHI1 + CHI2 + CHI3
        PSI1(NRCHOPS(2)+1:) = 0.D0
        CHI1(NRCHOPS(2)+1:) = 0.D0

        H(:NR_MK,  I+NR_MK) = PSI1(:NR_MK)
        H(NR_MK+1:,I+NR_MK) = CHI1(:NR_MK)
    ENDDO
    ! DEALLOCATE( RUR0,RUP0,UZ0,ROR0,ROP0,OZ0 )
    DEALLOCATE( RUPU,UZU,RORU,ROPU,OZU )
    DEALLOCATE( PSIU,CHIU,PSI1,CHI1,PSI2,CHI2,PSI3,CHI3 )

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,H,SIZE(H),MPI_DOUBLE_COMPLEX,MPI_SUM,&
    MPI_COMM_WORLD,IERR)

! OUTPUT:
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    IF (MPI_GLB_RANK.EQ.0) THEN  ! START OF THE SERIAL PART
        !WRITE OUT THE 3D MATRIX IN THE SCALAR USING UNFORMATTED MODE
        
        IF (PRESENT(EIG_VAL)) THEN
            IF (ALLOCATED(EIG_VAL)) DEALLOCATE(EIG_VAL)
            ALLOCATE( EIG_VAL(2*NR_MK) )

            IF (PRESENT(EIG_VEC_L)) THEN
                IF (ALLOCATED(EIG_VEC_R)) DEALLOCATE(EIG_VEC_R) 
                IF (ALLOCATED(EIG_VEC_L)) DEALLOCATE(EIG_VEC_L) 
                ALLOCATE( EIG_VEC_R(2*NR_MK,2*NR_MK),EIG_VEC_L(2*NR_MK,2*NR_MK) )
                CALL EIGENDECOMPOSE(H,EIG_VAL,ER = EIG_VEC_R,EL = EIG_VEC_L)
            ELSEIF (PRESENT(EIG_VEC_R)) THEN
                IF (ALLOCATED(EIG_VEC_R)) DEALLOCATE(EIG_VEC_R) 
                ALLOCATE( EIG_VEC_R(2*NR_MK,2*NR_MK) )
                CALL EIGENDECOMPOSE(H,EIG_VAL,ER = EIG_VEC_R)
            ELSE
                CALL EIGENDECOMPOSE(H,EIG_VAL)
            ENDIF
        ENDIF

        ! WRITE(*,*) 'NON-HERMITIAN MATRIX EIGENVALUE SOLVER'
        ! WRITE(*,*) 'CURRENTLY RUNNING IN SERIAL...'
        ! WRITE(*,*) 'EIGVALS:'
        ! CALL MCAT(EIG_VAL)
    ENDIF

    RETURN

END SUBROUTINE EIG_MATRIX
!=======================================================================
END PROGRAM EVP_SPARSITY
!=======================================================================
