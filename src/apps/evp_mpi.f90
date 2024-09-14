PROGRAM EVP_MPI

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

implicit NONE
! ======================================================================
INTEGER     :: MREAD, NR_MK, II
REAL(P8)    :: AKREAD, REI
INTEGER     :: MPI_GLB_PROCS,MPI_GLB_RANK,newcomm,MPI_NR_SIZE,MPI_NR_INDEX

REAL(P8),DIMENSION(:),ALLOCATABLE:: RUR0,RUP0,UZ0,ROR0,ROP0,OZ0
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: RURU,RUPU,UZU,RORU,ROPU,OZU
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: PSIU,CHIU,PSI1,CHI1,PSI2,CHI2,PSI3,CHI3

COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: H ! Assembled NS operator
COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: EIGVAL
COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: EIGVEC_L,EIGVEC_R

CHARACTER(LEN=72),DIMENSION(20):: COM
CHARACTER(LEN=72):: FILENAME
INTEGER:: PRINT_SWITCH
CHARACTER*72:: STR1, STR2, STR3, STR4, STR5
COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: EIGVECS
! ======================================================================

! BEGIN
CALL MPI_INIT(IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_GLB_PROCS, IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_GLB_RANK, IERR)
IF (MPI_GLB_RANK.EQ.0) THEN
    WRITE(*,*) 'PROGRAM STARTED'
    CALL PRINT_REAL_TIME()   ! @ MOD_MISC
ENDIF
CALL READCOM('NOECHO')
CALL READIN(5)

! INPUT M & K
IF (MPI_GLB_RANK.EQ.0) THEN
! ! INPUT M AND AK FROM KEYBOARD:
! ======================================================================
!     WRITE(*,*) ''
!     WRITE(*,*) 'CHOOSE AN AZIMUTHAL WAVENUMBER (M >= 0) ; M ='
!     READ(5,'(I72)') MREAD
!     WRITE(*,*) ''
!     WRITE(*,103) ZLEN
! 103  FORMAT(' AK = 2*PI/', F0.3)    
!     WRITE(*,*) 'CHOOSE AN AXIAL WAVENUMBER (AK) INDEX; AK ='
!     READ(5,'(F5.3)') AKREAD ! INPUT MUST HAVE DIGISTS; E.G. 1.0
!     ZLEN = 2*PI/AKREAD
! ======================================================================

! INPT M AND AK FROM *.INPUT FILE:
! ======================================================================
    FILENAME = 'matlab_read.input'
    CALL READIN(FILENAME,MREAD,AKREAD,PRINT_SWITCH)

    WRITE(*,*) 'MREAD:',MREAD
    WRITE(*,*) 'AKREAD',AKREAD
    WRITE(*,*) 'PRINT',PRINT_SWITCH
    WRITE(*,*) ''
    ZLEN = 2*PI/AKREAD
! ======================================================================

ENDIF
CALL MPI_BCAST(MREAD,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(ZLEN,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

! SPLIT PROCS: EA PROC GROUP ONLY CONTAINS ONE PROC
CALL MPI_Comm_split(MPI_COMM_WORLD, MPI_GLB_RANK, 0, newcomm, IERR)

! INIT 
NTH = 2
NX = 4
NTCHOP = 2 ! M = 0, M = MREAD
NXCHOP = 2 ! K = 0, K = AKREAD
CALL LEGINIT(newcomm,MREAD)
ALLOCATE(RUR0(NR),RUP0(NR),UZ0(NR),ROR0(NR),ROP0(NR),OZ0(NR))
CALL INIT_LOOP(RUR0,RUP0,UZ0,ROR0,ROP0,OZ0) ! ALL IN PFF space

! SPLIT NRCHOPDIM
NR_MK = NRCHOPS(2)
ALLOCATE(H(2*NR_MK,2*NR_MK))
H = CMPLX(0.D0,0.D0,P8)
CALL DECOMPOSE(NR_MK, MPI_GLB_PROCS, MPI_GLB_RANK, MPI_NR_SIZE, MPI_NR_INDEX)

! MAIN JOB
ALLOCATE(PSIU(NRCHOPDIM),CHIU(NRCHOPDIM))
ALLOCATE(PSI1(NRCHOPDIM),CHI1(NRCHOPDIM))
ALLOCATE(PSI2(NRCHOPDIM),CHI2(NRCHOPDIM))
ALLOCATE(PSI3(NRCHOPDIM),CHI3(NRCHOPDIM))
ALLOCATE(RURU(NRCHOPDIM),RUPU(NRCHOPDIM),UZU(NRCHOPDIM))
ALLOCATE(RORU(NRCHOPDIM),ROPU(NRCHOPDIM),OZU(NRCHOPDIM))

DO II = MPI_NR_INDEX+1,MPI_NR_INDEX+MPI_NR_SIZE
! ================================ PSI =================================
    PSIU = 0.D0
    CHIU = 0.D0
    PSIU(II) = 1.D0

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
    ! IF (II.EQ.2) THEN
    !     CALL XXDX_MK(PSIU,RURU)
    !     CALL MCAT(ROP0)
    ! ENDIF

    PSI1 = -PSI1 + PSI2 + PSI3
    CHI1 = -CHI1 + CHI2 + CHI3
    PSI1(NRCHOPS(2)+1:) = 0.D0
    CHI1(NRCHOPS(2)+1:) = 0.D0

    H(:NR_MK,  II) = PSI1(:NR_MK)
    H(NR_MK+1:,II) = CHI1(:NR_MK)

! ================================ CHI =================================
    PSIU = 0.D0
    CHIU = 0.D0
    CHIU(II) = 1.D0

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

    H(:NR_MK,  II+NR_MK) = PSI1(:NR_MK)
    H(NR_MK+1:,II+NR_MK) = CHI1(:NR_MK)
ENDDO

CALL MPI_ALLREDUCE(MPI_IN_PLACE,H,SIZE(H),MPI_DOUBLE_COMPLEX,MPI_SUM,&
MPI_COMM_WORLD,IERR)

! END
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
IF (MPI_GLB_RANK.EQ.0) THEN  ! START OF THE SERIAL PART
  !WRITE OUT THE 3D MATRIX IN THE SCALAR USING UNFORMATTED MODE
    ALLOCATE( EIGVAL(2*NR_MK) )
    ALLOCATE( EIGVEC_R(2*NR_MK,2*NR_MK) )
    ALLOCATE( EIGVEC_L(2*NR_MK,2*NR_MK) )
    WRITE(*,*) 'NON-HERMITIAN MATRIX EIGENVALUE SOLVER'
    WRITE(*,*) 'CURRENTLY RUNNING IN SERIAL...'

    CALL EIGENDECOMPOSE(H,EIGVAL,EIGVEC_R,EIGVEC_L)
    WRITE(*,*) 'EIGVALS:'
    CALL MCAT(EIGVAL)

    IF (PRINT_SWITCH.EQ.1) THEN ! SAVE FILES
        WRITE( STR1, '(f10.3)' )  AKREAD
        WRITE( STR2, '(f10.3)' )  ELL
        WRITE( STR3, '(f10.3)' )  QPAIR%H(1)

      ALLOCATE(EIGVECS(NDIMR,2*NR_MK,3))
      
      DO II = 1,2*NR_MK

        PSIU = 0.D0
        CHIU = 0.D0
        RURU = 0.D0
        RUPU = 0.D0
        UZU = 0.D0

        PSIU(:NR_MK) = EIGVEC_R(:NR_MK,II)
        CHIU(:NR_MK) = EIGVEC_R(NR_MK+1:,II)

        CALL CHOPSET(2)
        CALL IDEL2_MK(CHIU,RURU)
        CHIU = RURU
        RURU = 0.D0
        CALL CHOPSET(1)
        CALL PC2VEL_MK(PSIU,CHIU,RURU,RUPU,UZU)
        CALL CHOPSET(-3)
        CALL RTRAN_MK(RURU,1)
        CALL RTRAN_MK(RUPU,1)
        CALL RTRAN_MK( UZU,1)

        EIGVECS(:NDIMR,II,1) = RURU(:NDIMR)
        EIGVECS(:NDIMR,II,2) = RUPU(:NDIMR)
        EIGVECS(:NDIMR,II,3) =  UZU(:NDIMR)

      ENDDO

      CALL MSAVE(EIGVECS(:,:,1), TRIM(ADJUSTL(FILES%SAVEDIR)) //&
                          'efun_phys_rur_m_'// ITOA4(MREAD) //&
                          '_k_'       // TRIM(ADJUSTL(STR1))  //&
                          '_L_'       // TRIM(ADJUSTL(STR2))  //&
                          '_QI_'      // TRIM(ADJUSTL(STR3))  //&
                          '_N_'       // ITOA4(NR_MK) // '.dat')
      CALL MSAVE(EIGVECS(:,:,2), TRIM(ADJUSTL(FILES%SAVEDIR)) //&
                          'efun_phys_rup_m_'// ITOA4(MREAD) //&
                          '_k_'       // TRIM(ADJUSTL(STR1))  //&
                          '_L_'       // TRIM(ADJUSTL(STR2))  //&
                          '_QI_'      // TRIM(ADJUSTL(STR3))  //&
                          '_N_'       // ITOA4(NR_MK) // '.dat')
      CALL MSAVE(EIGVECS(:,:,3), TRIM(ADJUSTL(FILES%SAVEDIR)) //&
                          'efun_phys__uz_m_'// ITOA4(MREAD) //&
                          '_k_'       // TRIM(ADJUSTL(STR1))  //&
                          '_L_'       // TRIM(ADJUSTL(STR2))  //&
                          '_QI_'      // TRIM(ADJUSTL(STR3))  //&
                          '_N_'       // ITOA4(NR_MK) // '.dat')

      DEALLOCATE(EIGVECS)

    ENDIF

ENDIF

! END
DEALLOCATE( RUR0,RUP0,UZ0,ROR0,ROP0,OZ0 )
DEALLOCATE( RURU,RUPU,UZU,RORU,ROPU,OZU )
DEALLOCATE( PSIU,CHIU,PSI1,CHI1,PSI2,CHI2,PSI3,CHI3 )
IF (MPI_GLB_RANK.EQ.0) THEN
    WRITE(*,*) ''
    WRITE(*,*) 'PROGRAM FINISHED'
    CALL PRINT_REAL_TIME()  ! @ MOD_MISC
ENDIF
CALL MPI_FINALIZE(IERR)


! ======================================================================
! ========================= FAST SUBROUTINES ===========================
! ======================================================================
! CONTAINS
! ! XXDX
! ! ======================================================================
! SUBROUTINE XXDX_MK(A,B)
! !=======================================================================
! COMPLEX(P8),DIMENSION(:),INTENT(IN):: A
! COMPLEX(P8),DIMENSION(:),INTENT(INOUT):: B

! INTEGER:: NN
! REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XXDXOP

! IF (M(2).EQ.0) THEN
!     IF (MPI_RANK.EQ.0) WRITE(*,*) 'XXDX_MK: M CANNOT EQUAL 0'
!     RETURN
! ENDIF
! B = 0.D0

! ALLOCATE( XXDXOP(NRCHOP,3) )

! NN = NRCHOPS(2)
! XXDXOP(:NN,:) = BAND_LEG_XXDX(NN,M(2),TFM%NORM(:,2))
! B(:NN)=BANMUL(XXDXOP(:NN,:),2,A(:NN))

! DEALLOCATE( XXDXOP )

! RETURN
! END SUBROUTINE XXDX_MK
! ! ======================================================================

! ! MULXM
! SUBROUTINE MULXM_MK(A,B)
! !=======================================================================
! COMPLEX(P8),DIMENSION(:),INTENT(IN):: A
! COMPLEX(P8),DIMENSION(:),INTENT(INOUT):: B

! INTEGER:: NN
! REAL(P8),DIMENSION(:,:),ALLOCATABLE:: XM

! IF (M(2).EQ.0) THEN
!     IF (MPI_RANK.EQ.0) WRITE(*,*) 'MULXM_MK: M CANNOT EQUAL 0'
!     RETURN
! ENDIF

! B = 0.D0

! ALLOCATE( XM(NRCHOP,3) )

! NN = NRCHOPS(2)
! XM(:NN,:) = BAND_LEG_XM(NN,M(2),TFM%NORM(:,2))
! B(:NN)=BANMUL(XM(:NN,:),2,A(:NN))

! DEALLOCATE( XM )

! RETURN
! END SUBROUTINE MULXM_MK
! ! ======================================================================

! ! DESLQH
! SUBROUTINE DELSQH_MK(A,B)
! !=======================================================================
! IMPLICIT NONE
! COMPLEX(P8),DIMENSION(:),INTENT(IN):: A
! COMPLEX(P8),DIMENSION(:),INTENT(INOUT):: B

! INTEGER:: NN,N

! IF (M(2).EQ.0) THEN
!     IF (MPI_RANK.EQ.0) WRITE(*,*) 'DELSQH_MK: M CANNOT EQUAL 0'
!     RETURN
! ENDIF

! B = 0.D0

! DO NN=1,NRCHOPS(2)
!     N = M(2) + (NN-1)
!     B(NN) = -A(NN)*(N*(N+1)/ELL2)
! ENDDO

! RETURN
! END SUBROUTINE DELSQH_MK
! ! ======================================================================

! ! DEL2
! SUBROUTINE DEL2_MK(A,B)
! !=======================================================================
! IMPLICIT NONE
! COMPLEX(P8),DIMENSION(:),INTENT(IN):: A
! COMPLEX(P8),DIMENSION(:),INTENT(INOUT):: B

! INTEGER:: NN
! REAL(P8),DIMENSION(:,:),ALLOCATABLE:: DEL2OP

! IF (M(2).EQ.0) THEN
!     IF (MPI_RANK.EQ.0) WRITE(*,*) 'DEL2_MK: M CANNOT EQUAL 0'
!     RETURN
! ENDIF
! IF (AK(2,2).EQ.0) THEN
!     IF (MPI_RANK.EQ.0) WRITE(*,*) 'DEL2_MK: AK CANNOT EQUAL 0'
!     RETURN
! ENDIF

! B = 0.D0

! ALLOCATE( DEL2OP(NRCHOP,5) )

! NN = NRCHOPS(2)
! DEL2OP(:NN,:) = BAND_LEG_RAT_DEL2H(NN,M(2),ELL,TFM%NORM(:,2))
! B(:NN)=BANMUL(DEL2OP(:NN,:),3,A(:NN))
! B(:NN)=B(:NN)-(AK(2,2)**2)*A(:NN)

! DEALLOCATE( DEL2OP )

! RETURN
! END SUBROUTINE DEL2_MK
! ! ======================================================================

! ! IDEL2
! SUBROUTINE IDEL2_MK(B,A)
! !=======================================================================
! IMPLICIT NONE
! COMPLEX(P8),DIMENSION(:),INTENT(IN):: B
! COMPLEX(P8),DIMENSION(:),INTENT(INOUT):: A
! REAL(P8),DIMENSION(:,:,:),ALLOCATABLE :: DEL2OP
! REAL(P8),DIMENSION(:,:),ALLOCATABLE :: DEL2OP_0

! INTEGER::KK,NN,MM,KV,NP,NM

! IF (M(2).EQ.0) THEN
!     IF (MPI_RANK.EQ.0) WRITE(*,*) 'IDEL2_MK: M CANNOT EQUAL 0'
!     RETURN
! ENDIF
! IF (AK(2,2).EQ.0) THEN
!     IF (MPI_RANK.EQ.0) WRITE(*,*) 'IDEL2_MK: AK CANNOT EQUAL 0'
!     RETURN
! ENDIF

! A = 0.D0

! A = B

! !> ALL M .NE. 0
! ALLOCATE( DEL2OP_0(NRCHOP,5))
! ALLOCATE( DEL2OP(NRCHOP,5,1) )  !NXCHOPDIM) )

! ! each M, all K
! NN=NRCHOPS(2)
! DEL2OP_0(:NN,1:5) = BAND_LEG_RAT_DEL2H(NN,M(2),ELL,TFM%NORM(:,2))
! DEL2OP(:NN,1:2,1)=DEL2OP_0(:NN,1:2)
! DEL2OP(:NN,3  ,1)=DEL2OP_0(:NN,3  )-AK(2,2)**2
! DEL2OP(:NN,4:5,1)=DEL2OP_0(:NN,4:5)

! CALL LUB(DEL2OP(:NN,1:5,1),3)
! CALL SOLVEB(DEL2OP(:NN,1:5,1),3,A(:NN))

! DEALLOCATE( DEL2OP, DEL2OP_0 )
! A(NRCHOPS(2)+1:) = CMPLX(0.D0,0.D0)

! RETURN
! END SUBROUTINE IDEL2_MK
! ! ======================================================================
        
! ! HELMP
! SUBROUTINE HELMP_MK(P,A,B,ALP,NU)
! !=======================================================================
! IMPLICIT NONE
! INTEGER:: P
! COMPLEX(P8),DIMENSION(:),INTENT(IN):: A
! COMPLEX(P8),DIMENSION(:),INTENT(INOUT):: B
! REAL(P8),INTENT(IN):: ALP
! REAL(P8),OPTIONAL:: NU

! COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: W,D2
! INTEGER:: IP

! ALLOCATE(W(SIZE(A)),D2(SIZE(A)))
! CALL DEL2_MK(A,D2)

! B=D2

! DO IP=4,P,2
! CALL DEL2_MK(B,W)
! B=W
! ENDDO

! B = B-ALP*A

! IF(PRESENT(NU)) THEN
! B = B+NU*D2
! ENDIF

! DEALLOCATE(W,D2)

! RETURN
! END SUBROUTINE HELMP_MK
! ! ======================================================================

! ! PC2VOR
! SUBROUTINE PC2VOR_MK(PSI,CHI,ROR,ROP,OZ,C)
! !=======================================================================
! IMPLICIT NONE
! COMPLEX(P8),DIMENSION(:):: PSI,CHI,ROR,ROP,OZ
! INTEGER:: NN,MV
! REAL(P8):: KV
! COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: W
! CHARACTER(LEN=1),OPTIONAL:: C

! CALL DEL2_MK(CHI,OZ)                                                  ! OZ NOW HAS (DEL^2)(CHI)
! CALL XXDX_MK(PSI,ROR)                                                 ! ROR NOW HAS (1-X^2)D/DX(PSI) = R*D/DR(PSI)
! CALL XXDX_MK(OZ,ROP)                                                  ! ROP NOW HAS (1-X^2)D/DX((DEL^2)(CHI)) = R*D/DR((DEL^2)(CHI))

! NN=NRCHOPS(2)
! MV=M(2)
! KV=AK(2,2)
! ROR(:NN)=-IU*MV*OZ(:NN)+IU*KV*ROR(:NN)                                ! ROR = -D/DPHI(DEL^2(CHI))+ R*D/DR(D/DZ(PSI))
! ROP(:NN)=    ROP(:NN)-(KV*MV)*PSI(:NN)                                ! ROP = R*D/DR(DEL^2(CHI))+ D/DPHI(D/DZ(PSI))

! CALL DELSQH_MK(PSI,OZ)                                                ! OZ/(1-X)**2 = - DELSQH(PSI) (MINUS SIGN ASSIGNED AT THE END)

! IF(.NOT.PRESENT(C)) THEN
! ALLOCATE(W(SIZE(PSI)))
! CALL MULXM_MK(OZ,W)                                                   ! W <- (1-X)*[OZ/(1-X)**2] = OZ/(1-X)
! CALL MULXM_MK(W,OZ)                                                   ! OZ <- (1-X)*[W] = (1-X)*[OZ/(1-X)]= OZ
! DEALLOCATE(W)
! ENDIF
! OZ = -OZ

! RETURN
! END SUBROUTINE PC2VOR_MK
! ! ======================================================================

! ! PC2VEL
! SUBROUTINE PC2VEL_MK(PSI,CHI,RUR,RUP,UZ,C)
! !=======================================================================
! IMPLICIT NONE
! COMPLEX(P8),DIMENSION(:):: PSI,CHI,RUR,RUP,UZ
! INTEGER:: NN,MV
! REAL(P8):: KV
! COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: W
! CHARACTER(LEN=1),OPTIONAL:: C

! CALL XXDX_MK(CHI,RUR)                                                 ! RUR NOW HAS (1-X^2)D/DX(CHI) = R*D/DR(CHI)
! CALL XXDX_MK(PSI,RUP)                                                 ! RUP NOW HAS (1-X^2)D/DX(PSI) = R*D/DR(PSI)

! NN=NRCHOPS(2)
! MV=M(2)
! KV=AK(2,2)
! RUR(:NN)=IU*MV*PSI(:NN)+(IU*KV)*RUR(:NN)                              ! RUR = D/DPHI(PSI)+ R*D/DR(D/DZ(CHI))
! RUP(:NN)=     -RUP(:NN)-(MV*KV)*CHI(:NN)                              ! RUP =-R*D/DR(PSI)+ D/DPHI(D/DZ(CHI))

! CALL DELSQH_MK(CHI,UZ)                                                ! UZ/(1-X)**2 = - DELSQH(CHI) (MINUS SIGN ASSIGNED AT THE END)
! IF(.NOT.PRESENT(C)) THEN
!     ALLOCATE(W(SIZE(PSI)))
!     CALL MULXM_MK(UZ,W)                                               ! W <- (1-X)*[UZ/(1-X)**2] = UZ/(1-X)
!     CALL MULXM_MK(W,UZ)                                               ! UZ <- (1-X)*[W] = (1-X)*[UZ/(1-X)]= UZ
!     DEALLOCATE(W)
! ENDIF
! UZ = -UZ

! RETURN
! END SUBROUTINE PC2VEL_MK
! ! ======================================================================

! ! VPROD
! SUBROUTINE VPROD_MK(RUR_0,RUP_0,UZ_0,ROR_U,ROP_U,OZ_U)
! !=======================================================================
! IMPLICIT NONE
! REAL(P8),DIMENSION(:),INTENT(IN):: RUR_0,RUP_0,UZ_0
! COMPLEX(P8),DIMENSION(:),INTENT(INOUT):: ROR_U,ROP_U,OZ_U

! INTEGER:: NN
! REAL(P8):: A1,A2,A3,C1,C2,C3,B1,B2,B3,D1,D2,D3

! !> DIMENSION CHECK
! IF(SIZE(RUR_0).NE.NR) STOP 'VPROD_MK: RUR_0 HAS IMPROPER DIMENSION'
! IF(SIZE(RUP_0).NE.NR) STOP 'VPROD_MK: RUP_0 HAS IMPROPER DIMENSION'
! IF(SIZE(UZ_0 ).NE.NR) STOP 'VPROD_MK: UZ_0 HAS IMPROPER DIMENSION'
! IF(SIZE(ROR_U ).NE.NDIMR) STOP 'VPROD_MK: UZ_0 HAS IMPROPER DIMENSION'
! IF(SIZE(ROP_U ).NE.NDIMR) STOP 'VPROD_MK: UZ_0 HAS IMPROPER DIMENSION'
! IF(SIZE(OZ_U ).NE.NDIMR) STOP 'VPROD_MK: UZ_0 HAS IMPROPER DIMENSION'

! !$OMP PARALLEL DO DEFAULT(SHARED) &
! !$OMP& PRIVATE(A1,A2,A3,C1,C2,C3,B1,B2,B3,D1,D2,D3)
! DO NN=1,NR
!     A1=RUR_0(NN)
!     A2=RUP_0(NN)
!     A3=UZ_0 (NN)
!     C1=RUR_0(NN)
!     C2=RUP_0(NN)
!     C3=UZ_0 (NN)
!     B1=REAL (ROR_U(NN))
!     B2=REAL (ROP_U(NN))
!     B3=REAL (OZ_U(NN))
!     D1=AIMAG(ROR_U(NN))
!     D2=AIMAG(ROP_U(NN))
!     D3=AIMAG(OZ_U(NN))
!     ROR_U(NN)=CMPLX(A2*B3-A3*B2, C2*D3-C3*D2,P8)
!     ROP_U(NN)=CMPLX(A3*B1-A1*B3, C3*D1-C1*D3,P8)
!     OZ_U(NN)=CMPLX(A1*B2-A2*B1, C1*D2-C2*D1,P8)&
!                 /TFM%R(NN)**2.D0
! ENDDO
! !$OMP END PARALLEL DO

! RETURN
! END SUBROUTINE VPROD_MK
! ! ======================================================================

! ! PROJECT
! SUBROUTINE PROJECT_MK(RUR,RUP,UZ,PSI,CHI)
! !=======================================================================
! COMPLEX(P8),DIMENSION(:),INTENT(IN):: RUR,RUP,UZ
! COMPLEX(P8),DIMENSION(:),INTENT(INOUT):: PSI,CHI

! REAL(P8),DIMENSION(:,:),ALLOCATABLE:: V,D,T,PFD
! COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: W1,W2

! INTEGER:: KK,MM,NN,MV,N,I
! REAL(P8)::KV
! REAL(P8),DIMENSION(NRH):: FFF

! ALLOCATE( V(NRCHOP,NRH) )
! ALLOCATE( D(NRCHOP,NRH) )
! ALLOCATE( T(NRCHOP,NRH) )
! ALLOCATE( PFD(NRH,NRCHOP+1) )
! ALLOCATE( W1(NRCHOP), W2(NRCHOP) )

! W1=0
! FFF=(1-TFM%X(:NRH)**2)/TFM%W(:NRH)

! MV=M(2)
! NN=NRCHOPS(2)
! PFD(:NRH,:NN+1)=BANMUL( TFM%PF(:NRH,:NN+1,2), &
!         BAND_LEG_XXDX(NN+1,MV,TFM%NORM(:,2)),2 )
! DO I=1,NN
!     N=MAX(1,MV+I-1)
!     V(I,:NRH)= TFM%PF(:NRH,I,2)/(N*(N+1)*FFF)
!     D(I,:NRH)= PFD(:NRH,I)/(N*(N+1)*FFF)
!     T(I,:NRH)= TFM%PF(:NRH,I,2)*TFM%W(:NRH)   
! ENDDO

! CALL EOMUL_MK(V(:NN,:),RUR(:NR),W1(:NN))
! CALL OEMUL_MK(D(:NN,:),RUP(:NR),W2(:NN))

! PSI(:NN) = -(IU*MV)*W1(:NN)-W2(:NN)
! !if (MPI_GLB_RANK.eq.0) call mcat(W2)
! CALL OEMUL_MK(D(:NN,:),RUR(:NR),W1(:NN))

! CALL EOMUL_MK(V(:NN,:),RUP(:NR),W2(:NN))
! CALL EOMUL_MK(T(:NN,:),UZ(:NR),CHI(:NN))

! KV=AK(2,2)
! CHI(:NN)=(IU*KV)*W1(:NN)+MV*KV*W2(:NN)-CHI(:NN)

! DEALLOCATE( V,D,T,PFD )
! DEALLOCATE( W1,W2 )

! PSI(NRCHOPS(2)+1:) = CMPLX(0.D0,0.D0)
! CHI(NRCHOPS(2)+1:) = CMPLX(0.D0,0.D0)

! RETURN
! END SUBROUTINE PROJECT_MK
! !=======================================================================
! SUBROUTINE EOMUL_MK(A,B,C)
! !=======================================================================
! REAL(P8),DIMENSION(:,:):: A
! COMPLEX(P8),DIMENSION(:):: B,C

! COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: BE,BO
! INTEGER:: NI,NJ,NJH

! NI  = SIZE(A,1)
! NJH = SIZE(A,2)
! NJ  = SIZE(B,1)

! IF(NJH*2 .NE. NJ) THEN
! IF (MPI_RANK.EQ.0) THEN
!     WRITE(*,*) 'EOMUL: SIZE MISMATCH.'
!     WRITE(*,*) 'NJ=',NJ,' NJH=',NJH
! ENDIF
! STOP
! ENDIF

! ALLOCATE( BE(NJH) )
! ALLOCATE( BO(NJH) )

! BE = ( B(1:NJH) + B(NJ:NJH+1:-1) )
! BO = ( B(1:NJH) - B(NJ:NJH+1:-1) )

! C(1::2)= A(1::2,:) .MUL. BE

! IF (NI .GT. 1) THEN
! C(2::2)= A(2::2,:) .MUL. BO
! ENDIF

! DEALLOCATE( BE,BO )

! RETURN
! END SUBROUTINE EOMUL_MK
! !=======================================================================
! SUBROUTINE OEMUL_MK(A,B,C)
! !=======================================================================
! REAL(P8),DIMENSION(:,:):: A
! COMPLEX(P8),DIMENSION(:):: B,C

! COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: BE,BO
! INTEGER:: NI,NJ,NJH

! NI  = SIZE(A,1)
! NJH = SIZE(A,2)
! NJ  = SIZE(B,1)

! IF(NJH*2 .NE. NJ) THEN
! IF (MPI_RANK.EQ.0) THEN
!     WRITE(*,*) 'OEMUL: SIZE MISMATCH.'
!     WRITE(*,*) 'NJ=',NJ,' NJH=',NJH
! ENDIF
! STOP
! ENDIF

! ALLOCATE( BE(NJH), BO(NJH) )

! BE = ( B(1:NJH) + B(NJ:NJH+1:-1) )
! BO = ( B(1:NJH) - B(NJ:NJH+1:-1) )

! C(1::2)= A(1::2,:) .MUL. BO

! C(2::2)= A(2::2,:) .MUL. BE

! DEALLOCATE( BE,BO )

! RETURN
! END SUBROUTINE OEMUL_MK
! ! ======================================================================

! ! RTRAN
! SUBROUTINE RTRAN_MK(A,IS)
! !=======================================================================
! IMPLICIT NONE
! COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: A
! INTEGER,INTENT(IN):: IS

! COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: B,BE,BO
! INTEGER:: NN,I,MM

! ! PHYSICAL TO MAPPED LEGENDRE SPACE:
! IF(IS.EQ.-1) THEN                                                  

! ALLOCATE(B(NRCHOPDIM))
! B = 0.D0

! ALLOCATE( BE(NRH) )
! ALLOCATE( BO(NRH) )

! DO I=1,NRH
!     BE(I) = (A(I)+A(NR-I+1))*TFM%W(I)
!     BO(I) = (A(I)-A(NR-I+1))*TFM%W(I)
! ENDDO

! NN = NRCHOPS(2)

! ! TFM%PF(RadialCollocPts, RadialModes, AzimuthalModes)
! IF (NN.GE.1) THEN
! B(1:NN:2)=TRANSPOSE(TFM%PF(:NRH,1:NN:2,2)) .MUL. BE
! ENDIF

! IF (NN.GE.2) THEN
! B(2:NN:2)=TRANSPOSE(TFM%PF(:NRH,2:NN:2,2)) .MUL. BO
! ENDIF

! DEALLOCATE(A)
! ALLOCATE(A(NRCHOPDIM))
! A = B
! DEALLOCATE( BE,BO,B )

! ! MAPPED LEGENDRE SPACE TO PHYSICAL:
! ELSE                                                               

! ALLOCATE(B(NDIMR))
! B = 0.D0

! ALLOCATE( BE(NRH) ) !NXCHOPDIM) )
! ALLOCATE( BO(NRH) ) !NXCHOPDIM) )

! NN = NRCHOPS(2)

! IF(NN.GE.1) THEN
! BE= TFM%PF(:NRH,1:NN:2,2) .MUL. A(1:NN:2)!NXCHOPDIM)
! ELSE
! BE=0
! ENDIF

! IF(NN.GE.2) THEN
! BO= TFM%PF(:NRH,2:NN:2,2) .MUL. A(2:NN:2)!NXCHOPDIM)
! ELSE
! BO=0
! ENDIF

! B(1:NRH)= BE(:)+BO(:)
! B(NR:NRH+1:-1)= BE(:)-BO(:)

! B(NRCHOPS(2)+1:) = CMPLX(0.D0,0.D0)

! DEALLOCATE(A)
! ALLOCATE(A(NDIMR))
! A = B
! DEALLOCATE( BE,BO,B )
! ENDIF

! RETURN
! END SUBROUTINE RTRAN_MK
! ! ======================================================================

! ! INIT
! SUBROUTINE INIT_LOOP(RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0)
! ! ======================================================================
! IMPLICIT NONE
! REAL(P8),DIMENSION(:):: RUR_0,RUP_0,UZ_0,ROR_0,ROP_0,OZ_0
! TYPE(SCALAR):: A, B, C   ! EXTRA SCALAR-TYPE VARIABLES
! TYPE(SCALAR):: D, E, F   ! EXTRA SCALAR-TYPE VARIABLES
! INTEGER     :: I, J, K   ! EXTRA INTEGER VARIABLES (FOR ITER.)
! REAL(P8)    :: X, Y, Z   ! EXTRA REAL(P8) VARIABLES

! ! DEALLOCATE(NRCHOPS)
! ! DEALLOCATE(NTCHOPS)
! ! DEALLOCATE(AK)
! ! DEALLOCATE(M)
! ! DEALLOCATE(TFM%R,TFM%W,TFM%X, &
! ! TFM%LN,TFM%NORM,TFM%PF,TFM%AT0,TFM%AT1, &
! ! TFM%Z,TFM%THR,TFM%THI,TFM%TH)

! ! CALL LEGINIT()

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!! STAGE 1. GENERATE PSI0 !!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALL ALLOCATE(A,PPP_SPACE)

! ! INITIALIZE
! A%E = 0.D0

! DO I=1,QPAIR%N
!     ! INPUT OMEGA_Z (=2Q*EXP(-R**2)) INTO GAUSDATA
!     ! IN ORDER TO OBTAIN PSI0, WE USE THE FORMULA
!     ! OMEGA_Z = - DELSQH PSI0      <-->
!     ! PSI0    = DELSQH^(-1)(-OMEGA_Z).
!     GAUSDATA%XC = QPAIR%X(I)
!     GAUSDATA%YC = QPAIR%Y(I)
!     GAUSDATA%Q  = QPAIR%Q(I)*2.
!     GAUSDATA%B  = 1.D0
!     GAUSDATA%DX = QPAIR%DX(I)
!     GAUSDATA%DY = QPAIR%DY(I)
!     GAUSDATA%K  = QPAIR%K(I)
!     GAUSDATA%DP = QPAIR%DP(I)
!     GAUSDATA%R  = QPAIR%R(I)
!     GAUSDATA%DR = QPAIR%DR(I)
!     GAUSDATA%RK = QPAIR%RK(I)
!     GAUSDATA%RP = QPAIR%RP(I)
    
!     ! CREATE A GAUSSIAN VORTEX MULTIPLIED BY (1-X)^(-2)
!     ! NOW A%E HAS THE FIELD F(R,T,Z) = 2Q*EXP(-B*R**2)*(1-X)^(-2)
!     ! WHERE X = (R^2-ELL^2)/(R^2+ELL^2)
!     CALL LAY(A,GAUS,2)
    
!     ! CHECK WHETHER THE GAUSSIAN IS ACTUALLY GAUSSIAN
!     ! FOR DEBUGGIN PURPOSES. TURN IT OFF IN PRACTICE
!     ! CALL CHECKGAUSSIAN(A,GAUSDATA%XC,GAUSDATA%YC,2)
! ENDDO

! ! (1-X)^(-2)*OMEGA_Z (PPP) --> (1-X)^(-2)*OMEGA_Z (FFF)
! ! THE FUNCTION IS TRUNCATED INTO THE FINITE BASIS FUNCTION SET,
! ! CAUSING THE INFORMATION LOSS DUE TO APPROXIMATION, WHICH IS
! ! INEVITABLE UNLESS WE USE THE INFINITE-SIZE FULL BASIS SET
! ! STILL, THE TRANSFORMED COEFFICIENT MATRIX IN FFF SPACE
! ! IS THE BEST APPROXIMATION OF THE INPUT FIELD IN PPP SPACE.
! CALL TOFF(A)

! ! (1-X)^(-2)*OMEGA_Z (FFF) --> -(1-X)^(-2)*OMEGA_Z (FFF)
! A%E=-A%E

! ! MUTE A LOGTERM
! A%LN=0

! ! ADD RANDOM NOISE IF THE NOISE OPTION IS ON
! IF(QPAIR%NOISE .EQ. 1) THEN
!     CALL ALLOCATE(B)
!     B=A
!     CALL NOISE(A,0.02_P8)
!     CALL TOFP(A)
!     CALL TOFP(B)
!     CALL FLATTEN(A,B)
!     CALL TOFF(A)
!     CALL DEALLOCATE(B)
! ENDIF

! ! -(1-X)^(-2)*OMEGA_Z (FFF) -> DELSQH^(-1)(-OMEGA_Z) = PSI0 (FFF)
! ! NOTE: IDELSQH = ((1-X)^(-2) DELSQH)^(-1) = DELSQH^(-1)*((1-X)^(-2))^(-1)
! ! CREATES INVERSE HORIZONTAL DEL-SQUARE OPERATOR WITH EXTRA (1-X)^2
! ! B NOW HAS A LOG TERM ASSOCIATED WITH THE R-DERIVATIVE 
! ! IN THE INVERSE DEL SQ.
! CALL ALLOCATE(B)
! CALL IDELSQH(A,B)

! CALL DEALLOCATE(A)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!! STAGE 2. GENERATE CHI0 !!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALL ALLOCATE(A,PPP_SPACE)

! ! INITIALIZE
! A%E=0

! DO I=1,QPAIR%N
!     ! INPUT V_Z (=H*EXP(-B*R**2)) INTO GAUSDATA
!     ! IN ORDER TO OBTAIN CHI0, WE USE THE FORMULA
!     ! V_Z = - DELSQH CHI0      <-->
!     ! CHI0    = DELSQH^(-1)(-V_Z).
!     GAUSDATA%XC = QPAIR%X(I)
!     GAUSDATA%YC = QPAIR%Y(I)
!     GAUSDATA%Q  = QPAIR%H(I)
!     GAUSDATA%B  = QPAIR%B(I)
!     GAUSDATA%DX = QPAIR%DX(I)
!     GAUSDATA%DY = QPAIR%DY(I)
!     GAUSDATA%K  = QPAIR%K(I)
!     GAUSDATA%DP = QPAIR%DP(I)
!     GAUSDATA%R  = QPAIR%R(I)
!     GAUSDATA%DR = QPAIR%DR(I)
!     GAUSDATA%RK = QPAIR%RK(I)
!     GAUSDATA%RP = QPAIR%RP(I)

!     ! CREATE A GAUSSIAN VORTEX MULTIPLIED BY (1-X)^(-2)
!     ! NOW A%E HAS THE FIELD F(R,T,Z) = B*EXP(-R**2)*(1-X)^(-2)
!     ! WHERE X = (R^2-ELL^2)/(R^2+ELL^2)
!     CALL LAY(A,GAUS,2)

!     ! CHECK WHETHER THE GAUSSIAN IS ACTUALLY GAUSSIAN
!     ! FOR DEBUGGIN PURPOSES. TURN IT OFF IN PRACTICE
!     ! CALL CHECKGAUSSIAN(A,GAUSDATA%XC,GAUSDATA%YC,2)
! ENDDO

! ! (1-X)^(-2)*V_Z (PPP) --> (1-X)^(-2)*V_Z (FFF)
! ! THE FUNCTION IS TRUNCATED INTO THE FINITE BASIS FUNCTION SET,
! ! CAUSING THE INFORMATION LOSS DUE TO APPROXIMATION, WHICH IS
! ! INEVITABLE UNLESS WE USE THE INFINITE-SIZE FULL BASIS SET
! ! STILL, THE TRANSFORMED COEFFICIENT MATRIX IN FFF SPACE
! ! IS THE BEST APPROXIMATION OF THE INPUT FIELD IN PPP SPACE.
! CALL TOFF(A)

! ! (1-X)^(-2)*V_Z (FFF) --> -(1-X)^(-2)*V_Z (FFF)
! A%E=-A%E

! ! MUTE A LOGTERM
! A%LN=0

! ! ADD RANDOM NOISE IF THE NOISE OPTION IS ON
! IF(QPAIR%NOISE.EQ.1) THEN
!     CALL ALLOCATE(C)
!     C=A
!     CALL NOISE(A,0.02_P8)
!     CALL TOFP(A)
!     CALL TOFP(C)
!     CALL FLATTEN(A,C)
!     CALL TOFF(A)
!     CALL DEALLOCATE(C)
! ENDIF

! ! -(1-X)^(-2)*V_Z (FFF) -> DELSQH^(-1)(-V_Z) = CHI0 (FFF)
! ! NOTE: IDELSQH = ((1-X)^(-2) DELSQH)^(-1) = DELSQH^(-1)*((1-X)^(-2))^(-1)
! ! CREATES INVERSE HORIZONTAL DEL-SQUARE OPERATOR WITH EXTRA (1-X)^2
! ! B NOW HAS A LOG TERM ASSOCIATED WITH THE R-DERIVATIVE 
! ! IN THE INVERSE DEL SQ.
! CALL ALLOCATE(C)
! CALL IDELSQH(A,C)

! CALL DEALLOCATE(A)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!! STAGE 3.  TO VEL & VOR !!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! B: PSI0
! ! C: CHI0

! CALL ALLOCATE(D)
! CALL ALLOCATE(E)
! CALL ALLOCATE(F)

! CALL PC2VEL(B,C,D,E,F)
! CALL RTRAN(D,1)
! CALL RTRAN(E,1)
! CALL RTRAN(F,1)
! RUR_0(1:NR) = D%E(:NR,1,1)
! RUP_0(1:NR) = E%E(:NR,1,1)
!  UZ_0(1:NR) = F%E(:NR,1,1)

! CALL DEALLOCATE(D)
! CALL DEALLOCATE(E)
! CALL DEALLOCATE(F)

! CALL ALLOCATE(D)
! CALL ALLOCATE(E)
! CALL ALLOCATE(F)

! CALL PC2VOR(B,C,D,E,F)
! CALL RTRAN(D,1)
! CALL RTRAN(E,1)
! CALL RTRAN(F,1)
! ROR_0(1:NR) = D%E(:NR,1,1)
! ROP_0(1:NR) = E%E(:NR,1,1)
!  OZ_0(1:NR) = F%E(:NR,1,1)

! ! !DEBUG:
! ! IF (MPI_GLB_RANK.EQ.2) THEN
! !  CALL MCAT(B%E(:,1,1))
! !  CALL MCAT(C%E(:,1,1))
! ! ENDIF

! CALL DEALLOCATE(D)
! CALL DEALLOCATE(E)
! CALL DEALLOCATE(F)

! CALL DEALLOCATE(B)
! CALL DEALLOCATE(C)

! RETURN
! END SUBROUTINE INIT_LOOP
! ! ======================================================================

! SUBROUTINE EIGENDECOMPOSE(A,EV,ER,EL)
! !=======================================================================
! ! [USAGE]:                                   
! ! CALCULATE THE EIGENDECOMPOSITION OF A GIVEN SQUARE MATRIX A (A=RDR^-1)
! ! IN ADDITION TO THE CALCULATION OF THE RIGHT EIGENVECTOR MATRIX R,
! ! THE LEFT EIGENVECTOR MATRIX L S.T. L*A = DL* CAN BE COMPUTED OPTINALLY
! ! [PARAMETERS]:                
! ! A >> N X N COMPLEX NON-HERMITIAN (GENERAL) MATRIX
! ! EV >> N COMPLEX EIGENVALUES IN N-VECTOR
! ! ER >> RIGHT EIGENVECTOR MATRIX SATISFYING AR = R X DIAG(EV)
! ! EL >> (OPTIONAL) LEFT EIGENVECTOR MATRIX SATISFYING L*A = DIAG(EV)X L* 
! !                  NORMALIZATION IS PERFORMED TO SATISFY EL* X ER = I_N 
! ! [DEPENDENCIES]:
! ! 1. ZGGEV(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! ! [UPDATES]:
! ! CODED BY SANGJOON LEE @ APR 24 2021
! !=======================================================================
! IMPLICIT NONE
! COMPLEX(P8), DIMENSION(:,:), INTENT(IN) :: A
! COMPLEX(P8), DIMENSION(:), INTENT(INOUT) :: EV
! COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT) :: ER
! COMPLEX(P8), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: EL

! COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(A,2)) :: W1,W2

! COMPLEX(P8), DIMENSION(SIZE(A,1)):: AEIG,BEIG
! COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(A,2)):: VV
! COMPLEX(P8), DIMENSION(2*SIZE(A,1)):: WK
! REAL(P8)   , DIMENSION(8*SIZE(A,1)):: RWK

! INTEGER(P4) :: I,NI,LW,INFO

! IF(SIZE(A,1).NE.SIZE(A,2)) THEN
! WRITE (*,*) 'EIGENDECOMPOSE: MATRIX NOT SQUARE'
! WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
! STOP
! ENDIF

! IF(SIZE(A,1).GT.SIZE(EV)) THEN
! WRITE (*,*) 'EIGENDECOMPOSE: EV VECTOR SIZE IS TOO SMALL'
! WRITE (*,*) SIZE(A,1),', BUT SIZE(EV) = ',SIZE(EV)
! STOP
! ENDIF

! NI=SIZE(A,1)

! W1=A
! W2=EYE(NI)

! LW = SIZE(WK)
! ER = CMPLX(0.D0,0.D0,P8)
! IF ( PRESENT(EL) ) EL = CMPLX(0.D0,0.D0,P8)

! IF ( PRESENT(EL) ) THEN
! CALL ZGGEV('V','V',NI,W1,NI,W2,NI,AEIG,BEIG, &
!                         EL,NI,ER,NI,WK,LW,RWK,INFO)
! DO I = 1, NI
!     EL(:,I) = EL(:,I)/DOT_PRODUCT(EL(:,I),ER(:,I))
! ENDDO
! ELSE
! CALL ZGGEV('N','V',NI,W1,NI,W2,NI,AEIG,BEIG, &
!                         VV,NI,ER,NI,WK,LW,RWK,INFO)
! ENDIF

! IF (INFO .NE. 0) THEN
! WRITE (*,*) 'EIGVEC2: EIGENVECTOR GENERATION FAILURE.'
! STOP
! ENDIF

! DO I = 1, NI
! IF(ABS(BEIG(I)).EQ.0.D0) THEN
!     EV(I) = HUGE(0.D0)                                             ! ASSIGN INFINITY TO AVOID THE DIVISION-BY-ZERO ERROR
! ELSE
!     EV(I) = AEIG(I)/BEIG(I)
! ENDIF
! ENDDO

! RETURN
! END SUBROUTINE EIGENDECOMPOSE
! ! ======================================================================

END PROGRAM EVP_MPI