MODULE MOD_FFT ! LEVEL 2.5 MODULE
! FFT ROUTINES & MPI OPERATIONS
    USE OMP_LIB
    USE MPI
    USE MOD_MISC, ONLY : P4,P8,PI,IU,ITOA3,CP8_SIZE,MCAT                   ! LEVEL 0
    USE MOD_EIG                                                        ! LEVEL 1
    USE MOD_LIN_LEGENDRE                                               ! LEVEL 1
    USE MOD_SCALAR3                                                    ! LEVEL 2
    IMPLICIT NONE
    PRIVATE
! ======================================================================
! =========================== FFT STRATEGY =============================
! ======================================================================
! N1: SUBCOMM_1
! N2: SUBCOMM_2
!
! HORFFT:
! PPP: NDIMR/N1, NDIMTH, NDIMX/N2 << PPP
! >>FFT
! PFP: NDIMR/N1, NDIMTH, NDIMX/N2 << PFP
!
! VERFFT:
! PFP: NDIMR/N1, NDIMTH, NDIMX/N2
! >>CHOP
! PFP: NDIMR/N1, NTCHOPDIM, NDIMX/N2 -> TYPE_PFP0
! >>EXCHANGE
! PFP: NDIMR/N1, NTCHOPDIM/N2, NDIMX -> TYPE_PFP1
! >>FFT
! PFF: NDIMR/N1, NTCHOPDIM/N2, NDIMX
! >>CHOP
! PFF: NDIMR/N1, NTCHOPDIM/N2, NXCHOPDIM -> TYPE_PFF0
! >>EXCHANGE
! PFF: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1 -> TYPE_PFF1

! RTRAN:
! PFF: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1 << PFF
! >>RTRAN
! FFF: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1
! >>CHOP
! FFF: NRCHOPDIM, NTCHOPDIM/N2, NXCHOPDIM/N1 << FFF
!
!=======================================================================
!======================== PUBLIC DECLARATION ===========================
!=======================================================================

! =========================== FFT UTILITIES ============================
! INITIALIZE ALL SPECTRAL METHOD PARAMETERS AND SET THE TFM KIT:
PUBLIC:: LEGINIT
! ALLOCATE/DEALLOCATE MEMORY TO VARIABLES OF "TYPE(SCALAR)":
PUBLIC:: ALLOCATE, DEALLOCATE                                      ! ALLOCATE(A,SPACE) ALLOCATES MEMORY FOR THE ???_SPACE (??? = PPP, PFP, PFF, FFF)
! WRAPPER FOR "EQUAL" OPERATOR (COPY0) AMONG TYPE(SCALAR) VARS.:
PUBLIC:: ASSIGNMENT(=)
! PERFORM THE FFT IN THETA (=HORFFT) AND Z (=VERFFT) DIRECTIONS:
PUBLIC:: HORFFT
PUBLIC:: VERFFT
! TRANSFORMS AMONG PPP, PFP, PFF, AND FFF SPACES:
PUBLIC:: RTRAN                                                     ! FOR EXAMPLE, RTRAN(A,1) TRANSFORMS FROM FFF_SPACE(FUNCTION) TO PFF_SPACE(PHYSICAL).
PUBLIC:: TOFF,TOFP
    
! =========================== MPI UTILITIES ============================
! CREATE CARTESIAN SUB-GRID COMMUNICATORS:
PUBLIC:: SUBCOMM_CART, DECOMPOSE
! ASSEMBLE LOCAL ARRAY INTO GLOBAL ARRAY:
PUBLIC:: MASSEMBLE
! LOAD/SAVE A MATRIX WITH TYPE(SCALAR) FROM/INTO A SPECIFIED FILE
PUBLIC:: MLOAD,MSAVE                                               ! SIMLIAR TO MLOAD & MSAVE IN MOD_MISC FOR GENERAL MATRICES, BUT FOR TYPE(SCALAR)

! ========================== MOD_FFT PRIVATE ===========================
! ! DECOMPOSE A 1D DOMAIN:
! PUBLIC:: DECOMPOSE
! ! EXCHANGE ARRAY DIMENSIONS:
! PUBLIC:: EXCHANGE
! ! CREATE SUBARRAY DATATYPE FOR EACH PROC THAT SPLITS THE ORIGINAL AR
! ! RAY INTO NPROC PRATS ALONG SPECIFIED DIM:
! PUBLIC:: SUBARRAY

! =========================== UTILITY FUNCS ============================
PUBLIC:: local_size, local_index, local_proc, count_proc

!=======================================================================
!============================ INTERFACES ===============================
!=======================================================================
    INTERFACE ALLOCATE
        MODULE PROCEDURE SALLOC
    END INTERFACE

    INTERFACE DEALLOCATE
        MODULE PROCEDURE SFREE
    END INTERFACE

    INTERFACE ASSIGNMENT(=)
        MODULE PROCEDURE COPY0
    END INTERFACE

    INTERFACE EXCHANGE
        MODULE PROCEDURE EXCHANGE_3DCOMPLEX
        MODULE PROCEDURE EXCHANGE_3DCOMPLEX_FAST
    END INTERFACE

    INTERFACE MLOAD
    MODULE PROCEDURE MLOAD0
    END INTERFACE

    INTERFACE MSAVE
        MODULE PROCEDURE MSAVE0
    END INTERFACE

CONTAINS
!=======================================================================
!============================ SUBROUTINES ==============================
!=======================================================================
! SUBROUTINE LEGINIT(NRIN,NTHIN,NXIN,NRCHOPIN,NTCHOPIN,NXCHOPIN,&
!     ZLENIN,ELLIN,MKLINKIN,MINCIN)
SUBROUTINE LEGINIT(MPI_COMM_INPUT, M_INPUT)
!=======================================================================
! [USAGE]: 
! INITIALIZE ALL SPECTRAL METHOD PARAMETERS AND 
! SET THE TRANSFORM KIT VIA TFM VARIABLE
! [PARAMETERS]:
! NRIN >> INTEGER VALUE USED FOR NR (radial)
! NTHIN >> INTEGER VALUE USED FOR NTH
! NXIN >> INTEGER VALUE USED FOR NX (axial)
! NRCHOPIN >> INTEGER VALUE USED FOR NRCHOP
! NTCHOPIN >> INTEGER VALUE USED FOR NTCHOP
! NXCHOPIN >> INTEGER VALUE USED FOR NXCHOP
! ZLENIN >> REAL VALUE USED FOR ZLEN
! ELLIN >> REAL VALUE USED FOR ELL
! MKLINKIN >> (OPTIONAL) INTEGER VALUE USED FOR MKLINK (DEFAULT: 0)
! MINCIN >> (OPTIONAL) INTEGER VALUE USED FOR MKLINK (DEFAULT: 1)
! [NOTE]:
! SEE ANNOTATIONS IN THE 'PARAMETERS' SECTION OF THIS MODULE TO
! COMPREHEND WHICH INPUT ASSIGNED TO WHICH PUBLIC PARAMETER AND WHAT
! THE MEANING OF EACH PARAMETER IS
! [DEPENDENCIES]:
! 1. LEG_ZERO(~)  @ MOD_LIN_LEGENDRE
! 2. LEG_NORM(~)  @ MOD_LIN_LEGENDRE
! 3. LEG_TBL(~)   @ MOD_LIN_LEGENDRE
! 4. EXSET(~)     @ MOD_SCALAR3
! [NOTES]:
! SUBCOMM_CART WILL REODER THE PROC INDEX
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    ! INTEGER,INTENT(IN):: NRIN,NTHIN,NXIN,NRCHOPIN,NTCHOPIN,NXCHOPIN
    ! INTEGER,INTENT(IN),OPTIONAL:: MKLINKIN,MINCIN

    INTEGER :: MM,KK,KV,I
    ! REAL(P8):: ZLENIN,ELLIN

    INTEGER, DIMENSION(:), ALLOCATABLE:: SUB_GROUPS
    ! INTEGER, DIMENSION(3):: SIZE_PFP0,SIZE_PFP1,SIZE_PFF0,SIZE_PFF1

    INTEGER, OPTIONAL:: M_INPUT,MPI_COMM_INPUT

    IF (PRESENT(MPI_COMM_INPUT)) THEN
        MPI_COMM_IVP = MPI_COMM_INPUT
    ELSE
        MPI_COMM_IVP = MPI_COMM_WORLD
    ENDIF
    CALL MPI_COMM_RANK(MPI_COMM_IVP,MPI_RANK,IERR)

    NRH = NR/2
    ELL2 = ELL**2.D0

    NDIMR  = NR+3
    NDIMTH = NTH+1
    NDIMX  = NX
    IF(NDIMTH.EQ.2) NDIMTH=1

    ! NRCHOPDIM  = MIN(NRCHOP+4,NR)       ! For Init_ene
    NRCHOPDIM  = MIN(NRCHOP+6,NR)         ! Original: +3
    NTCHOPDIM  = MIN(NTCHOP,NTH)
    NXCHOPDIM  = NXCHOP*2-1
    NXCHOPH    = NX-NXCHOP+2                                           ! CHOP LOCATION FOR CONJUGATE SIDE
    IF (NXCHOP==1) NXCHOPH=1

    IF ((2*NXCHOP>NX).AND.((NX>1).OR.(NXCHOP>1))) THEN
        IF (MPI_RANK.EQ.0) THEN
            WRITE(*,*) 'ERROR: SCALAR3() -- NXCHOP MUST BE < NX/2,'
            WRITE(*,*) '                    UNLESS NX = NXCHOP = 1'
        ENDIF
      STOP
    ENDIF

    IF (.NOT.(ALLOCATED(M))) ALLOCATE(M(NTCHOPDIM))
    M = (/ (MINC*(MM-1),MM=1,NTCHOPDIM) /)
    IF (PRESENT(M_INPUT).AND.(NTCHOPDIM.EQ.2)) M(2) = M_INPUT

    IF (.NOT.(ALLOCATED(NRCHOPS))) THEN
        ALLOCATE(NRCHOPS(NTH))
        ALLOCATE(NTCHOPS(NR))
    ENDIF
    CALL CHOPSET(0)

    IF (.NOT.(ALLOCATED(AK))) ALLOCATE( AK(NTCHOPDIM,NXCHOPDIM) )
    DO MM=1,NTCHOPDIM
      DO KK=1,NXCHOPDIM
        KV=KK-1 + MKLINK*(MM-1)
        IF(KK.GT.NXCHOP) THEN
          KV=-(NXCHOPDIM-KK+1)
        ENDIF
        AK(MM,KK) = 2*PI/ZLEN*KV
      ENDDO
    ENDDO

    IF (.NOT.(ASSOCIATED(TFM%R))) THEN
        ALLOCATE( TFM%R(NR) )                                              ! START SETTING UP TRANSFORM PARAMETERS INTO THE TFM VARIABLE
        ALLOCATE( TFM%dR(NR))
        ALLOCATE( TFM%W(NR) )
        ALLOCATE( TFM%X(NR) )
        ALLOCATE( TFM%LN(NR) )
        ALLOCATE( TFM%NORM(NRCHOPDIM+14,NTCHOPDIM) )
        ALLOCATE( TFM%LOGNORM(NRCHOPDIM+14,NTCHOPDIM) )
        ALLOCATE( TFM%PF(NRH,NRCHOPDIM+1,NTCHOPDIM) )
        ALLOCATE( TFM%AT0(NRCHOPDIM) )
        ALLOCATE( TFM%AT1(NRCHOPDIM) )
        ALLOCATE( TFM%THR(MINC*NTH+1) )
        ALLOCATE( TFM%THI(MINC*NTH+1) )
        ALLOCATE( TFM%TH(2*MINC*NTH+1) )
    ENDIF

    CALL LEG_ZERO(TFM%X,TFM%R,ELL,TFM%W)
    TFM%dR(1:NR-1) = TFM%R(2:NR)-TFM%R(1:NR-1)
    TFM%dR(NR) = TFM%dR(NR-1)

    !ALLOCATE( TFM%LN(NR) )
    TFM%LN = -LOG(1-TFM%X)

    !ALLOCATE( TFM%NORM(NRCHOPDIM+14,NTCHOPDIM) )
    TFM%NORM = LEG_NORM(NRCHOPDIM+14, M )   ! BECOME INACCURATE AT HIGHER SPECTRAL MODES AND SHOULD BE REPLACED BY TFM%LOGNORM INSTEAD

    !ALLOCATE( TFM%LOGNORM(NRCHOPDIM+14,NTCHOPDIM) )
    TFM%LOGNORM = LEG_LOG_NORM(NRCHOPDIM+14, M )

    !ALLOCATE( TFM%PF(NRH,NRCHOPDIM+1,NTCHOPDIM) )
    ! TFM%PF = LEG_TBL(TFM%X(:NRH),NRCHOPDIM+1, M ,TFM%NORM)
    TFM%PF = LEG_TBL(TFM%X(:NRH),NRCHOPDIM+1, M ,TFM%LOGNORM)

    !ALLOCATE( TFM%AT0(NRCHOPDIM) )
    !ALLOCATE( TFM%AT1(NRCHOPDIM) )
    ! TFM%AT0 = LEG_TBL( -1.D0,NRCHOPDIM, 0 , TFM%NORM(:,1))
    ! TFM%AT1 = LEG_TBL(  1.D0,NRCHOPDIM, 0 , TFM%NORM(:,1))
    TFM%AT0 = LEG_TBL( -1.D0,NRCHOPDIM, 0 , TFM%LOGNORM(:,1))
    TFM%AT1 = LEG_TBL(  1.D0,NRCHOPDIM, 0 , TFM%LOGNORM(:,1))

    ! TFM%NMAX = MAX(NTH,NX)*2
    ! ALLOCATE( TFM%EX(TFM%NMAX) )
    ! CALL EXSET( TFM%EX, TFM%NMAX )

    IF(MKLINK.EQ.0) THEN
        IF (.NOT.(ASSOCIATED(TFM%Z))) ALLOCATE( TFM%Z(NX+1) )
        TFM%Z = (/ (ZLEN*(KK-1)/NX, KK=1,NX+1) /)
    ELSE IF(MKLINK.EQ.1 .OR. MKLINK.EQ.-1) THEN
        IF (.NOT.(ASSOCIATED(TFM%Z))) ALLOCATE( TFM%Z(2*NTH+1) )
        TFM%Z = (/ (ZLEN*(KK-1)/(2*NTH), KK=1,2*NTH+1) /)
    ELSE
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'LEGINIT: INVALID MKLINK. MKLINK=',MKLINK
      STOP
    ENDIF

    ! ALLOCATE( TFM%THR(MINC*NTH+1) )
    ! ALLOCATE( TFM%THI(MINC*NTH+1) )
    TFM%THR = (/ (2*PI/MINC/NTH*(I-1), I=1,MINC*NTH+1) /)
    TFM%THI = (/ (2*PI/MINC/NTH*(I-0.5D0), I=1,MINC*NTH+1) /)

    ! ALLOCATE( TFM%TH(2*MINC*NTH+1) )
    TFM%TH = (/ (PI/MINC/NTH*(I-1), I=1,2*MINC*NTH+1) /)

    ! FORM MPI COMMUNICATOR GROUPS:
    IF (PRESENT(M_INPUT)) THEN ! FREE PREVIOUSLY USED COMMUNICATORS
        IF (SUBCOMM_1.NE.0) CALL MPI_COMM_FREE(SUBCOMM_1,IERR)
        IF (SUBCOMM_2.NE.0) CALL MPI_COMM_FREE(SUBCOMM_2,IERR)
    ENDIF
    CALL SUBCOMM_CART(MPI_COMM_IVP, 2, SUB_GROUPS)
! ============= FROM THIS POINT, PROC INDEX IS REORDERED ===============
    IF (NXCHOPDIM.GE.NTCHOPDIM) THEN
        SUBCOMM_1 = SUB_GROUPS(1) ! FFF: NX/N1
        SUBCOMM_2 = SUB_GROUPS(2) ! FFF: NT/N2
    ELSE
        SUBCOMM_1 = SUB_GROUPS(2)
        SUBCOMM_2 = SUB_GROUPS(1)       
    ENDIF
    DEALLOCATE(SUB_GROUPS)
    ! WRITE(*,*) MPI_RANK,':',NXCHOPDIM,count_proc(SUBCOMM_1),NTCHOPDIM,count_proc(SUBCOMM_2)

    ! FORM MPI SUBARRAY DATA TYPES:
    SIZE_PFP0 = (/local_size(NDIMR,SUBCOMM_1), &
                                    NTCHOPDIM, &
                  local_size(NDIMX,SUBCOMM_2) /)
    TYPE_PFP0 = create_new_type3D(SUBCOMM_2, SIZE_PFP0, 2, MPI_DOUBLE_COMPLEX)

    SIZE_PFP1 = (/local_size(NDIMR,SUBCOMM_1), &
              local_size(NTCHOPDIM,SUBCOMM_2), &
                                        NDIMX /)
    TYPE_PFP1 = create_new_type3D(SUBCOMM_2, SIZE_PFP1, 3, MPI_DOUBLE_COMPLEX)

    SIZE_PFF0 = (/local_size(NDIMR,SUBCOMM_1), &
              local_size(NTCHOPDIM,SUBCOMM_2), &
                                    NXCHOPDIM /)
    TYPE_PFF0 = create_new_type3D(SUBCOMM_1, SIZE_PFF0, 3, MPI_DOUBLE_COMPLEX)
    !WRITE(*,*) LBOUND(TYPE_PFF0,1)

    SIZE_PFF1 = (/                      NDIMR, &
              local_size(NTCHOPDIM,SUBCOMM_2), &
              local_size(NXCHOPDIM,SUBCOMM_1) /)
    TYPE_PFF1 = create_new_type3D(SUBCOMM_1, SIZE_PFF1, 1, MPI_DOUBLE_COMPLEX)

    ! WRITE(*,311) MPI_RANK,local_proc(SUBCOMM_1),local_proc(SUBCOMM_2),SIZE_PFP0,SIZE_PFP1,SIZE_PFF0,SIZE_PFF1
    ! 311  FORMAT(I3,'(',I3,'-',I3,'):',4(I3,' ',I3,' ',I3,' - '))

    RETURN
    END SUBROUTINE LEGINIT
!=======================================================================
SUBROUTINE SALLOC(A,SP)
!=======================================================================
! [USAGE]: 
! ALLOCATE AN ARRAY AT A%E WHERE A IS SCALAR-TYPE VARIABLE WITH SP TAG
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES
! SP >> (OPTIONAL) SPACE TAG TO BE ATTACHED. IF NONE, FFF IS ASSIGNED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! [NOTE]:
! PPP: NDIMR/N1, NDIMTH, NDIMX/N2
! PFP: NDIMR/N1, NDIMTH, NDIMX/N2
! PFF: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1
! FFF: NRCHOPDIM, NTCHOPDIM/N2, NXCHOPDIM/N1
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
        IMPLICIT NONE
        TYPE(SCALAR):: A
        INTEGER,OPTIONAL:: SP

        INTEGER:: MM

        IF(PRESENT(SP)) THEN
        A%SPACE = SP
        ELSE
        A%SPACE = FFF_SPACE
        ENDIF

        IF(ASSOCIATED(A%E)) THEN 
        NULLIFY ( A%E )
        ENDIF

        ! RESET INDEX
        A%INR = 0
        A%INTH = 0
        A%INX = 0

        ! IF(A%SPACE .EQ. FFF_SPACE) THEN
        ! ALLOCATE( A%E(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM),STAT=IERR )
        ! ELSE
        ! ALLOCATE( A%E(NDIMR,NDIMTH,NDIMX),STAT=IERR )
        ! ENDIF

        SELECT CASE(A%SPACE)
        CASE (PFF_SPACE)
            ALLOCATE( A%E(NDIMR    , local_size(NTCHOPDIM,SUBCOMM_2), local_size(NXCHOPDIM,SUBCOMM_1)),STAT=IERR )
            A%INTH = local_index(NTCHOPDIM,SUBCOMM_2)
            A%INX  = local_index(NXCHOPDIM,SUBCOMM_1)

        CASE (FFF_SPACE)
            ALLOCATE( A%E(NRCHOPDIM, local_size(NTCHOPDIM,SUBCOMM_2), local_size(NXCHOPDIM,SUBCOMM_1)),STAT=IERR )
            A%INTH = local_index(NTCHOPDIM,SUBCOMM_2)
            A%INX  = local_index(NXCHOPDIM,SUBCOMM_1)

        ! PPP_SPACE OR PFP_SPACE:
        CASE DEFAULT 
            ALLOCATE( A%E(local_size(NDIMR,SUBCOMM_1), NDIMTH, local_size(NDIMX,SUBCOMM_2)),STAT=IERR )
            A%INR  = local_index(NDIMR,SUBCOMM_1)
            A%INX  = local_index(NDIMX,SUBCOMM_2)

        END SELECT

        IF (IERR.NE.0) THEN
            WRITE(*,*) "SALLOC: ERROR OCCURED - COULD NOT ALLOCATE ARRAY."
        STOP
        ENDIF

        RETURN
        END SUBROUTINE SALLOC
!=======================================================================
    SUBROUTINE SFREE(A)
!=======================================================================
! [USAGE]: 
! FREE AN ARRAY WHERE A SCLAR-TYPE VARIABLE A INDICATES VIA A%E
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES
! [NOTES]:
! DEALLOCATE TARGET OF A%E
! NULLIFY POINTER A%E
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
        IMPLICIT NONE
        TYPE(SCALAR):: A

        IF(.NOT. ASSOCIATED( A%E )) THEN
        WRITE(*,*) 'SFREE:TRIED TO DEALLOCATE THOUGH NOT ASSOCIATED.'
        STOP
        ENDIF

        DEALLOCATE( A%E )
        NULLIFY ( A%E )

        RETURN
        END SUBROUTINE SFREE
!=======================================================================
SUBROUTINE COPY0(B,A)
!=======================================================================
! [USAGE]: 
! COPY THE SCALAR-TYPE VARIABLE A TO B
! WRAPPED BY 'EQUAL(=)' SYMBOL, THUS ABLE TO BE CALLED VIA 'B = A'
! [PARAMETERS]:
! A >> ORIGINAL VARIABLE 
! B >> SCALAR-TYPE VARIABLE TO BE COPIED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
        IMPLICIT NONE
        TYPE(SCALAR),INTENT(INOUT):: B
        TYPE(SCALAR),INTENT(IN):: A
        INTEGER:: MM

        IF(.NOT.ASSOCIATED(B%E)) THEN
            WRITE(*,*) 'COPY0: NO MEMORY FOR LHS.'
            STOP
        ENDIF

        ! IF(SIZE(B%E)<SIZE(A%E)) THEN
        ! CALL DEALLOCATE(B)
        ! CALL ALLOCATE(B,A%SPACE)
        ! ENDIF
        IF(B%SPACE.NE.A%SPACE) THEN
            CALL DEALLOCATE(B)
            CALL ALLOCATE(B,A%SPACE)
        ENDIF

        ! B%SPACE = A%SPACE
        B%LN = A%LN
        B%INR = A%INR
        B%INTH = A%INTH
        B%INX = A%INX

        ! IF(A%SPACE .EQ. FFF_SPACE) THEN
        !     DO MM=1,NTCHOP
        !         B%E(:NRCHOPS(MM),MM,:) = A%E(:NRCHOPS(MM),MM,:)
        !     ENDDO
        ! CALL CHOPDO(B)
        ! ELSE
        !     B%E(:SIZE(A%E,1),:SIZE(A%E,2),:SIZE(A%E,3)) = A%E
        !     B%E(:,:,SIZE(A%E,3)+1:)=0
        !     B%E(:,SIZE(A%E,2)+1:,:SIZE(A%E,3))=0
        !     B%E(SIZE(A%E,1)+1:,:SIZE(A%E,2),:SIZE(A%E,3))=0
        ! ENDIF
        IF(A%SPACE .EQ. FFF_SPACE) THEN
            B%E = A%E
        CALL CHOPDO(B)
        ELSE
            B%E(:SIZE(A%E,1),:SIZE(A%E,2),:SIZE(A%E,3)) = A%E
            B%E(:,:,SIZE(A%E,3)+1:)=0
            B%E(:,SIZE(A%E,2)+1:,:SIZE(A%E,3))=0
            B%E(SIZE(A%E,1)+1:,:SIZE(A%E,2),:SIZE(A%E,3))=0
        ENDIF

        RETURN
        END SUBROUTINE COPY0
!=======================================================================
SUBROUTINE HORFFT(A,IS)
!=======================================================================
! [USAGE]: 
! PERFORM HORIZONTAL FFT WITH RESPECT TO AZIMUTHAL (THETA) DIRECTION
! EXPECTED TO BE CALLED WHEN DATA IS IN PPP SPACE
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES
! IS >> FORWARD OR BACKWARD FFT (1: BACKWARD, -1: FORWARD)
! [NOTES]:
! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
!    CALL RTRAN(A,1) (FFF -> PFF)
!    CALL VERFFT(A,1) (PFF -> PFP)
!    CALL HORFFT(A,1) (PFP -> PPP)
! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
!    CALL HORFFT(A,-1) (PPP -> PFP)
!    CALL VERFFT(A,-1) (PFP -> PFF)
!    CALL RTRAN(A,-1) (PFF -> FFF)
! [DEPENDENCIES]
! 1. DFFTW_PLAN_DFT                           @ FFTW3 
! 2. DFFTW_EXECUTE_DFT                        @ FFTW3
! 3. DFFTW_DESTROY_PLAN                       @ FFTW3
! 4. FFTW_FORWARD,FFTW_BACKWARD,FFTW_ESTIMATE @ FFTW3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! FFT SUBROUTINES ARE REPLACED WITH FFTW3 LIBRARY @ NOV NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    TYPE(SCALAR):: A
    INTEGER,INTENT(IN):: IS

    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: B
    REAL(P8),DIMENSION(:),ALLOCATABLE:: C
    INTEGER:: II,KK,MM,RSIZE,XSIZE
    REAL(P8):: FAC
    INTEGER(P8):: PLAN

    ! PPP: NDIMR/N1, NDIMTH, NDIMX/N2
    ! PFP: NDIMR/N1, NDIMTH, NDIMX/N2

    RSIZE = SIZE(A%E,1)
    XSIZE = SIZE(A%E,3)
    CALL CHOPDO(A) 

    ! PHYSICAL TO FOURIER SPACE:
    IF(IS .EQ.-1) THEN                                                 
    IF(A%SPACE.NE.PPP_SPACE)THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'HORFFT: NOT IN PHYSICAL SPACE'
        STOP
    ENDIF
    IF(NTH.EQ.1) THEN
        A%SPACE=PFP_SPACE
        RETURN
    ENDIF

    ALLOCATE( B(NDIMTH) )! for 2*NTH, R2CFFT gives NTH+1
    ALLOCATE( C(2*NTH ) )
    B = A%E(1,:NDIMTH,1)
    DO II=1,NTH
        C(2*II-1) = REAL(B(II))
        C(2*II  ) = AIMAG(B(II))
    ENDDO
    CALL DFFTW_PLAN_DFT_R2C_1D(PLAN,2*NTH,C,B,FFTW_ESTIMATE)

    DEALLOCATE( B )
    DEALLOCATE( C )

    ALLOCATE(B(NDIMTH))
    ALLOCATE(C(2*NTH))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B,C,KK,II,MM) COLLAPSE(2)
    DO KK=1,XSIZE !NX                                                ! THIS DO LOOP IS PARALLELIZABLE 
        DO II=1,RSIZE !NR 

        B = A%E(II,:NDIMTH,KK)
        DO MM = 1,NTH
            C(2*MM-1) = REAL(B(MM))
            C(2*MM  ) = AIMAG(B(MM))
        ENDDO
        CALL DFFTW_EXECUTE_DFT_R2C(PLAN,C,B)
        A%E(II,:NDIMTH,KK) = B/(2*NTH)

        ENDDO
    ENDDO
!$OMP END PARALLEL DO
    DEALLOCATE( B )
    DEALLOCATE( C )

    CALL DFFTW_DESTROY_PLAN(PLAN)

    A%E(:,NTCHOP+1:,:) = CMPLX(0.D0,0.D0)                            ! CHOPPING IN THETA DIRECTION

    A%E(:,1,:) = CMPLX(REAL(A%E(:,1,:)),0.D0)                        ! MAKE SURE TO ELIMINATE COMPLEX RESIDUAL DUE TO MACHINE ERROR

    A%SPACE = PFP_SPACE                                              ! CHANGE SPACE TAG FROM PPP -> PFP


    ! FOURIER TO PHYSICAL SPACE:
    ELSE                                                               
    IF(A%SPACE.NE.PFP_SPACE)THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'HORFFT: NOT IN PFP SPACE'
        STOP
    ENDIF
    IF(NTH.EQ.1) THEN
        A%SPACE=PPP_SPACE
        RETURN
    ENDIF

    A%E(:,1,:) = CMPLX(REAL(A%E(:,1,:)),0.D0)                        ! MAKE SURE TO ELIMINATE COMPLEX RESIDUAL DUE TO MACHINE ERROR

    A%E(:,NTCHOP+1:,:) = CMPLX(0.D0,0.D0)                            ! CHOPPING IN THETA DIRECTION

    ALLOCATE( B(NDIMTH) )
    ALLOCATE( C(2*NTH ) )
    B = A%E(1,:NDIMTH,1)
    C = 0.D0
    CALL DFFTW_PLAN_DFT_C2R_1D(PLAN,2*NTH,B,C,FFTW_ESTIMATE)
                                        
    DEALLOCATE( B )
    DEALLOCATE( C )

    ALLOCATE( B(NDIMTH) )
    ALLOCATE( C(2*NTH ) )
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B,C,KK,II,MM) COLLAPSE(2)
    DO KK=1,XSIZE !NX                                                ! THIS DO LOOP IS PARALLELIZABLE 
        DO II=1,RSIZE !NR

        B = A%E(II,:NDIMTH,KK)
        C = 0.D0
        CALL DFFTW_EXECUTE_DFT_C2R(PLAN,B,C)
        DO MM = 1,NTH
            A%E(II,MM,KK) = CMPLX(C(2*MM-1),C(2*MM),P8)
        ENDDO
        A%E(II,NDIMTH,KK) = CMPLX(0.D0, 0.D0)

        ENDDO
    ENDDO
!$OMP END PARALLEL DO
    DEALLOCATE( B )
    DEALLOCATE( C )

    CALL DFFTW_DESTROY_PLAN(PLAN)

    A%SPACE = PPP_SPACE                                              ! CHANGE SPACE TAG FROM PFP -> PPP

    ENDIF

    CALL CHOPDO(A)
    RETURN
    END SUBROUTINE HORFFT
!=======================================================================
    SUBROUTINE VERFFT(A,IS)
!=======================================================================
! [USAGE]: 
! PERFORM VERTICAL FFT WITH RESPECT TO AXIAL (Z) DIRECTION
! EXPECTED TO BE CALLED WHEN DATA IS IN PFP SPACE
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES
! IS >> FORWARD OR BACKWARD FFT (1: BACKWARD, -1: FORWARD)
! [NOTES]:
! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
!    CALL RTRAN(A,1) (FFF -> PFF)
!    CALL VERFFT(A,1) (PFF -> PFP)
!    CALL HORFFT(A,1) (PFP -> PPP)
! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
!    CALL HORFFT(A,-1) (PPP -> PFP)
!    CALL VERFFT(A,-1) (PFP -> PFF)
!    CALL RTRAN(A,-1) (PFF -> FFF)
! [DEPENDENCIES]
! 1. DFFTW_PLAN_DFT                           @ FFTW3 
! 2. DFFTW_EXECUTE_DFT                        @ FFTW3
! 3. DFFTW_DESTROY_PLAN                       @ FFTW3
! 4. FFTW_FORWARD,FFTW_BACKWARD,FFTW_ESTIMATE @ FFTW3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! FFT SUBROUTINES ARE REPLACED WITH FFTW3 LIBRARY @ NOV NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    TYPE(SCALAR),INTENT(INOUT):: A
    INTEGER,INTENT(IN):: IS

    COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: B
    INTEGER:: JJ,KK
    REAL(P8):: FAC
    INTEGER(P8):: PLAN
    ! TYPE(C_PTR):: PLAN2

    COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: A_PFP0,A_PFP1,A_PFF0,A_PFF1
    REAL(P8):: LN

    ! PFP: NDIMR/N1, NDIMTH, NDIMX/N2
    ! PFF: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1

    ! MPI PREP =============================================================
    LN = A%LN
    CALL CHOPDO(A)
    ! ALLOCATE: TWO UTILITY ARRAYS IN PFP0 AND PFP1
    ALLOCATE(A_PFP0(SIZE_PFP0(1),SIZE_PFP0(2),SIZE_PFP0(3)))
    ALLOCATE(A_PFP1(SIZE_PFP1(1),SIZE_PFP1(2),SIZE_PFP1(3)))
    ! ALLOCATE: TWO UTILITY ARRAYS IN PFF0 AND PFF1
    ALLOCATE(A_PFF0(SIZE_PFF0(1),SIZE_PFF0(2),SIZE_PFF0(3)))
    ALLOCATE(A_PFF1(SIZE_PFF1(1),SIZE_PFF1(2),SIZE_PFF1(3)))  


    ! PHYSICAL TO FORUIER SPACE: ===========================================
    IF(IS .EQ.-1) THEN                                                 
    IF(A%SPACE.NE.PFP_SPACE) THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'VERFFT: NOT IN PFP SPACE'
        STOP
    ENDIF
    IF(NX.EQ.1) THEN
        A%SPACE=PFF_SPACE
        RETURN
    ENDIF
 
    ! CHOP: NDIMTH -> NTCHOPDIM (NDIMR/N1, NTCHOPDOM, NDIMX/N2)
    A_PFP0 = A%E(:,:NTCHOPDIM,:)
    ! EXCHANGE: A_PFP0 -> A_PFP1(NDIMR/N1, NTCHOPDIM/N2, NDIMX)
    CALL EXCHANGE_3DCOMPLEX_FAST(SUBCOMM_2,A_PFP0,TYPE_PFP0,A_PFP1,TYPE_PFP1)
    DEALLOCATE(A_PFP0)

    ! CREATE FFT PLAN
    ALLOCATE( B(NX) )
    ! B = A%E(1,1,:NX)
    B = A_PFP1(1,1,:NX)
    CALL DFFTW_PLAN_DFT_1D(PLAN,NX,B,B,FFTW_FORWARD,&
                                        FFTW_ESTIMATE)
    DEALLOCATE( B )

    ALLOCATE( B(NX) )
    ! FFT DO LOOP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B,JJ,KK) COLLAPSE(2)
    DO JJ=1,SIZE(A_PFP1,2) !NTCHOP                                  !!!!! THIS DO LOOP IS PARALLELIZABLE 
        DO KK=1,SIZE(A_PFP1,1) !NR
        !B = A%E(KK,JJ,:NX)
        B = A_PFP1(KK,JJ,:NX)
        CALL DFFTW_EXECUTE_DFT(PLAN,B,B)
        !A%E(KK,JJ,:NX) = B/NX
        A_PFP1(KK,JJ,:NX) = B/NX
        ENDDO
    ENDDO
!$OMP END PARALLEL DO
    DEALLOCATE( B )
    CALL DFFTW_DESTROY_PLAN(PLAN)

    ! A%E(:NR,:NTCHOP,NXCHOP+1:NXCHOPH-1) = CMPLX(0.D0,0.D0)           ! CHOPPING IN X DIRECTION

    ! CHOP: NDIMX -> NXCHOPDIM (NDIMR/N1, NTCHOPDIM/N2, NXCHOPDIM)
    A_PFF0(:,:,:NXCHOP) = A_PFP1(:,:,:NXCHOP)
    A_PFF0(:,:,NXCHOP+1:) = A_PFP1(:,:,NXCHOPH:)
    ! EXCHANGE: A_PFF0 -> A_PFF1(NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1)
    CALL EXCHANGE_3DCOMPLEX_FAST(SUBCOMM_1,A_PFF0,TYPE_PFF0,A_PFF1,TYPE_PFF1)
    DEALLOCATE(A_PFP1,A_PFF0)

    ! OUTPUT A BACK (NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1)
    CALL DEALLOCATE(A)
    CALL ALLOCATE(A, PFF_SPACE) ! RESET A%INR, INTH, INX
    A%E = A_PFF1
    A%LN = LN
    ! CALL CHOPDO(A)
    DEALLOCATE(A_PFF1)
    ! A%SPACE = PFF_SPACE                                              ! CHANGE SPACE TAG FROM PFP -> FFP


    ! FOURIER TO PHYSICAL ==================================================
    ELSE                                                               
    IF(A%SPACE.NE.PFF_SPACE) THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'VERFFT: NOT IN PFF SPACE'
        STOP
    ENDIF
    IF(NX.EQ.1) THEN
        A%SPACE=PFP_SPACE
        RETURN
    ENDIF
    
    ! A%E -> A_PFF1(NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1)
    A_PFF1 = A%E
    ! EXCHANGE: A_PFF1 -> A_PFF0(NDIMR/N1, NTCHOPDIM/N2, NXCHOPDIM)
    CALL EXCHANGE_3DCOMPLEX_FAST(SUBCOMM_1,A_PFF1,TYPE_PFF1,A_PFF0,TYPE_PFF0)
    ! UNCHOP: NXCHOPDIM -> NDIMX(NDIMR/N1, NTCHOPDIM/N2, NDIMX)
    A_PFP1 = CMPLX(0.D0,0.D0)
    A_PFP1(:,:,:NXCHOP) = A_PFF0(:,:,:NXCHOP)
    A_PFP1(:,:,NXCHOPH:) = A_PFF0(:,:,NXCHOP+1:)
    !A%E(:NR,:NTCHOP,NXCHOP+1:NXCHOPH-1) = CMPLX(0.D0,0.D0)           ! CHOPPING IN X DIRECTION
    A_PFP1(:,:,NXCHOP+1:NXCHOPH-1) = CMPLX(0.D0,0.D0)
    DEALLOCATE(A_PFF0,A_PFF1)

    ! CREATE FFT PLAN
    ALLOCATE( B(NX) )
    !B = A%E(1,1,:NX)
    B = A_PFP1(1,1,:NX)
    CALL DFFTW_PLAN_DFT_1D(PLAN,NX,B,B,FFTW_BACKWARD,&
                                        FFTW_ESTIMATE)
    !PLAN2 = FFTW_PLAN_DFT_1D(NX,B,B,FFTW_FORWARD,FFTW_ESTIMATE)
    DEALLOCATE( B )

    ALLOCATE( B(NX) )
    ! FFT DO LOOP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B) COLLAPSE(2)
    DO JJ=1,SIZE(A_PFP1,2) !NTCHOP                                   !!!!! THIS DO LOOP IS PARALLELIZABLE 
        DO KK=1,SIZE(A_PFP1,1) !NR
        !B = A%E(KK,JJ,:NX)
        B = A_PFP1(KK,JJ,:NX)
        !if (MPI_RANK.eq.0) write(*,*) JJ,KK,B(1)
        CALL DFFTW_EXECUTE_DFT(PLAN,B,B)
        !if (MPI_RANK.eq.0) write(*,*) 'afterDFFTW',B(1)
        !A%E(KK,JJ,:NX) = B
        A_PFP1(KK,JJ,:NX) = B
        ENDDO
    ENDDO
!$OMP END PARALLEL DO
    DEALLOCATE( B )
    CALL DFFTW_DESTROY_PLAN(PLAN)

    ! EXCHANGE: A_PFP1 -> A_PFP0(NDIMR/N1, NTCHOPDIM, NDIMX/N2)
    CALL EXCHANGE_3DCOMPLEX_FAST(SUBCOMM_2,A_PFP1,TYPE_PFP1,A_PFP0,TYPE_PFP0)
    DEALLOCATE(A_PFP1)
    
    ! OUTPUT A BACK (NDIMR/N1, NDIMTH, NDIMX/N2)
    CALL DEALLOCATE(A)
    CALL ALLOCATE(A, PFP_SPACE) ! RESET A%INR, INTH, INX
    ! UNCHOP: NTCHOPDIM -> NDIMTH (NDIMR/N1, NDIMTH, NDIMX/N2)
    A%E(:,:NTCHOPDIM,:) = A_PFP0
    A%E(:,NTCHOPDIM+1:,:) = CMPLX(0.D0,0.D0)
    A%LN = LN
    DEALLOCATE(A_PFP0)
    A%SPACE = PFP_SPACE                                              ! CHANGE SPACE TAG FROM FFP -> PFP
    
    ENDIF

    CALL CHOPDO(A)
    RETURN
    END SUBROUTINE VERFFT
!=======================================================================
    SUBROUTINE RTRAN(A,IS)
!=======================================================================
! [USAGE]: 
! PERFORM MAPPED LEGENDRE TRANSFORM WITH RESPECT TO X (OR R) DIRECTION
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES
! IS >> FORWARD OR BACKWARD FFT (1: BACKWARD, -1: FORWARD)
! [NOTES]:
! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
!    CALL RTRAN(A,1) (FFF -> PFF)
!    CALL VERFFT(A,1) (PFF -> PFP)
!    CALL HORFFT(A,1) (PFP -> PPP)
! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
!    CALL HORFFT(A,-1) (PPP -> PFP)
!    CALL VERFFT(A,-1) (PFP -> PFF)
!    CALL RTRAN(A,-1) (PFF -> FFF)
! [DEPENDENCIES]
! 1. ALLOCATE(~) @ MOD_SCALAR3
! 2. CHOPDO(~) @ MOD_SCALAR3
! 3. DEALLOCATE(~)  @ MOD_SCALAR3
! 4. OPERATOR(.MUL.) @ MOD_EIG
! [UPDATES]:
! FFT SUBROUTINES ARE REPLACED WITH FFTW3 LIBRARY @ NOV NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR),INTENT(INOUT):: A
    INTEGER,INTENT(IN):: IS

    TYPE(SCALAR):: B
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: BE,BO
    INTEGER:: NN,I,MM
    INTEGER:: XSIZE, THSIZE

    ! MPI PREP
    XSIZE = SIZE(A%E,3)
    THSIZE = SIZE(A%E,2)

    ! PHYSICAL TO MAPPED LEGENDRE SPACE:
    IF(IS.EQ.-1) THEN                                                  
    IF(A%SPACE .NE. PFF_SPACE) THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'RTRAN: NOT IN PFF_SPACE.'
        STOP
    ENDIF

    CALL ALLOCATE(B,FFF_SPACE)
    CALL CHOPDO(B)

    ! ALLOCATE( BE(NRH,NXCHOPDIM) )
    ! ALLOCATE( BO(NRH,NXCHOPDIM) )
    ! ALLOCATE( BE(NRH,XSIZE) )
    ! ALLOCATE( BO(NRH,XSIZE) )
    ALLOCATE( BE(XSIZE,NRH) )
    ALLOCATE( BO(XSIZE,NRH) )

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NN,BE,BO,MM,I)
    DO MM=1, THSIZE !NTCHOP

        ! NN = NRCHOPS(MM+A%INTH)
        NN = NRCHOPS(MM+B%INTH)
        IF (NN.LT.1) CYCLE

        DO I=1,NRH
            ! BE(I,:NXCHOP) = (A%E(I,MM,:NXCHOP)+A%E(NR-I+1,MM,:NXCHOP))&
            !                 *TFM%W(I)
            ! BO(I,:NXCHOP) = (A%E(I,MM,:NXCHOP)-A%E(NR-I+1,MM,:NXCHOP))&
            !                 *TFM%W(I)

            ! IF(NXCHOP .NE. 1) THEN
            !     BE(I,NXCHOP+1:) = (A%E(I,MM,NXCHOPH:NX)+A%E(NR-I+1,MM,&
            !                         NXCHOPH:NX))*TFM%W(I)  
            !     BO(I,NXCHOP+1:) = (A%E(I,MM,NXCHOPH:NX)-A%E(NR-I+1,MM,&
            !                         NXCHOPH:NX))*TFM%W(I)
            ! ENDIF

            ! A%E: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1
            ! BE(I,:) = (A%E(I,MM,:)+A%E(NR-I+1,MM,:))*TFM%W(I)
            ! BO(I,:) = (A%E(I,MM,:)-A%E(NR-I+1,MM,:))*TFM%W(I)
            BE(:,I) = (A%E(I,MM,:)+A%E(NR-I+1,MM,:))*TFM%W(I)
            BO(:,I) = (A%E(I,MM,:)-A%E(NR-I+1,MM,:))*TFM%W(I)

        ENDDO

        

        ! TFM%PF(RadialCollocPts, RadialModes, AzimuthalModes)
        ! IF(NN.GE.1) THEN
        ! B%E(1:NN:2,MM,:)=TRANSPOSE(TFM%PF(:NRH,1:NN:2,MM+A%INTH)) .MUL. BE
        ! ENDIF
        B%E(1:NN:2,MM,:)=TRANSPOSE(BE .MUL. TFM%PF(:NRH,1:NN:2,MM+B%INTH))
        
        IF (NN.GE.2) THEN
        ! B%E(2:NN:2,MM,:)=TRANSPOSE(TFM%PF(:NRH,2:NN:2,MM+A%INTH)) .MUL. BO
        B%E(2:NN:2,MM,:)=TRANSPOSE(BO .MUL. TFM%PF(:NRH,2:NN:2,MM+B%INTH))
        ENDIF

    ENDDO
!$OMP END PARALLEL DO

    CALL DEALLOCATE(A)

    A%E => B%E
    A%INR = B%INR
    A%INTH = B%INTH
    A%INX = B%INX
    NULLIFY(B%E)
    A%SPACE=FFF_SPACE

    DEALLOCATE( BE,BO )


    ! MAPPED LEGENDRE SPACE TO PHYSICAL:
    ELSE                                                               
    IF(A%SPACE .NE. FFF_SPACE) THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'RTRAN: NOT IN FFF_SPACE'
        STOP
    ENDIF

    CALL ALLOCATE(B,PFF_SPACE)
    CALL CHOPDO(B)
    B%LN=0

    ALLOCATE( BE(NRH,XSIZE) ) !NXCHOPDIM) )
    ALLOCATE( BO(NRH,XSIZE) ) !NXCHOPDIM) )

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NN,BE,BO,MM)
    DO MM=1,THSIZE !NTCHOP
        NN = NRCHOPS(MM+A%INTH)

        IF(NN.GE.1) THEN
        BE= TFM%PF(:NRH,1:NN:2,MM+A%INTH) .MUL. A%E(1:NN:2,MM,:)!NXCHOPDIM)
        ELSE
        BE=0
        ENDIF

        IF(NN.GE.2) THEN
        BO= TFM%PF(:NRH,2:NN:2,MM+A%INTH) .MUL. A%E(2:NN:2,MM,:)!NXCHOPDIM)
        ELSE
        BO=0
        ENDIF

        ! B%E(1:NRH,MM,1:NXCHOP)= BE(:,:NXCHOP)+BO(:,:NXCHOP)

        ! IF(NXCHOP .NE. 1) THEN
        ! B%E(1:NRH,MM,NXCHOPH:NX)= BE(:,NXCHOP+1:)+BO(:,NXCHOP+1:)
        ! B%E(NR:NRH+1:-1,MM,NXCHOPH:NX)= BE(:,NXCHOP+1:)-BO(:,NXCHOP+1:)
        ! ENDIF

        ! B%E(NR:NRH+1:-1,MM,:NXCHOP)= BE(:,:NXCHOP)-BO(:,:NXCHOP)

        B%E(1:NRH,MM,:)= BE(:,:)+BO(:,:)
        B%E(NR:NRH+1:-1,MM,:)= BE(:,:)-BO(:,:)

    ENDDO
!$OMP END PARALLEL DO
    
    ! B%E(:,NTCHOP+1:,:)=0
    ! B%E(:,:NTCHOP,NXCHOP+1:NXCHOPH-1)=0
    B%E(:,MAX(1,NTCHOP-A%INTH+1):,:) = CMPLX(0.D0,0.D0)

    IF(A%INTH.EQ.0 .AND. A%INX.EQ.0) THEN
        IF(A%LN.NE.0.0) THEN
            WRITE(*,*) 'RTRAN:TO PHYSICAL SPACE THOUGH LOGTERM IS NONZERO.'
            WRITE(*,*) 'LOGTERM=',A%LN
            B%E(1:NR,1,1)=B%E(1:NR,1,1)+A%LN*TFM%LN(1:NR)
        ENDIF
    ENDIF

    CALL DEALLOCATE(A)

    A%E => B%E
    NULLIFY(B%E)
    A%SPACE = PFF_SPACE
    A%INR = B%INR
    A%INTH = B%INTH
    A%INX = B%INX

    DEALLOCATE( BE,BO )
    ENDIF
    
    CALL CHOPDO(A)
    RETURN
    END SUBROUTINE RTRAN
!=======================================================================
    SUBROUTINE RTRAN2(A,IS) ! RTRAN FOR TEST PURPOSE
!=======================================================================
! [USAGE]: 
! INSTEAD OF USING MATRIX-MATRIX MULTIPLICATION, PERFORM MATRIX-VECTOR
! MULTIPLICATION IN FOR LOOP, WHICH FORCES LAPACK'S ZGEMM TO PRODUCE IDE
! -NTICAL RESULTS (NO MACHINE ROUND-OFF DIFFERENCES) WITH DIFFERENT NUMB
! -ER OF PROCESSORS.
! PERFORM MAPPED LEGENDRE TRANSFORM WITH RESPECT TO X (OR R) DIRECTION
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES
! IS >> FORWARD OR BACKWARD FFT (1: BACKWARD, -1: FORWARD)
! [NOTES]:
! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
!    CALL RTRAN(A,1) (FFF -> PFF)
!    CALL VERFFT(A,1) (PFF -> PFP)
!    CALL HORFFT(A,1) (PFP -> PPP)
! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
!    CALL HORFFT(A,-1) (PPP -> PFP)
!    CALL VERFFT(A,-1) (PFP -> PFF)
!    CALL RTRAN(A,-1) (PFF -> FFF)
! [DEPENDENCIES]
! 1. ALLOCATE(~) @ MOD_SCALAR3
! 2. CHOPDO(~) @ MOD_SCALAR3
! 3. DEALLOCATE(~)  @ MOD_SCALAR3
! 4. OPERATOR(.MUL.) @ MOD_EIG
! [UPDATES]:
! FFT SUBROUTINES ARE REPLACED WITH FFTW3 LIBRARY @ NOV NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR),INTENT(INOUT):: A
    INTEGER,INTENT(IN):: IS

    TYPE(SCALAR):: B
    COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: BE,BO
    INTEGER:: NN,I,MM
    INTEGER:: XSIZE, THSIZE

    ! MPI PREP
    XSIZE = SIZE(A%E,3)
    THSIZE = SIZE(A%E,2)

    ! PHYSICAL TO MAPPED LEGENDRE SPACE:
    IF(IS.EQ.-1) THEN                                                  
    IF(A%SPACE .NE. PFF_SPACE) THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'RTRAN: NOT IN PFF_SPACE.'
        STOP
    ENDIF

    CALL ALLOCATE(B,FFF_SPACE)
    CALL CHOPDO(B)

    ! ALLOCATE( BE(NRH,NXCHOPDIM) )
    ! ALLOCATE( BO(NRH,NXCHOPDIM) )
    ALLOCATE( BE(NRH,XSIZE) )
    ALLOCATE( BO(NRH,XSIZE) )

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NN,BE,BO,MM,I)
    DO MM=1, THSIZE !NTCHOP

        NN = NRCHOPS(MM+A%INTH)
        IF(NN.LT.1) CYCLE

        DO I=1,NRH

            ! A%E: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1
            BE(I,:) = (A%E(I,MM,:)+A%E(NR-I+1,MM,:))*TFM%W(I)
            BO(I,:) = (A%E(I,MM,:)-A%E(NR-I+1,MM,:))*TFM%W(I)

        ENDDO

        ! TFM%PF(RadialCollocPts, RadialModes, AzimuthalModes)
        DO I=1,XSIZE
            B%E(1:NN:2,MM,I:I)=TRANSPOSE(TFM%PF(:NRH,1:NN:2,MM+A%INTH)) .MUL. BE(:,I:I)
            
            IF (NN.GE.2) THEN
            B%E(2:NN:2,MM,I:I)=TRANSPOSE(TFM%PF(:NRH,2:NN:2,MM+A%INTH)) .MUL. BO(:,I:I)
            ENDIF
        ENDDO

    ENDDO
!$OMP END PARALLEL DO

    CALL DEALLOCATE(A)

    A%E => B%E
    A%INR = B%INR
    A%INTH = B%INTH
    A%INX = B%INX
    NULLIFY(B%E)
    A%SPACE=FFF_SPACE

    DEALLOCATE( BE,BO )


    ! MAPPED LEGENDRE SPACE TO PHYSICAL:
    ELSE                                                               
    IF(A%SPACE .NE. FFF_SPACE) THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'RTRAN: NOT IN FFF_SPACE'
        STOP
    ENDIF

    CALL ALLOCATE(B,PFF_SPACE)
    CALL CHOPDO(B)
    B%LN=0

    ALLOCATE( BE(NRH,XSIZE) ) !NXCHOPDIM) )
    ALLOCATE( BO(NRH,XSIZE) ) !NXCHOPDIM) )

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NN,BE,BO,MM,I)
    DO MM=1,THSIZE !NTCHOP
        NN = NRCHOPS(MM+A%INTH)

        IF(NN.GE.1) THEN
            DO I=1,XSIZE
                BE(:,I:I)= TFM%PF(:NRH,1:NN:2,MM+A%INTH) .MUL. A%E(1:NN:2,MM,I:I)!NXCHOPDIM)
            ENDDO
        ELSE
        BE=0
        ENDIF

        IF(NN.GE.2) THEN
            DO I=1,XSIZE
                BO(:,I:I)= TFM%PF(:NRH,2:NN:2,MM+A%INTH) .MUL. A%E(2:NN:2,MM,I:I)!NXCHOPDIM)
            ENDDO
        ELSE
        BO=0
        ENDIF

        ! B%E(1:NRH,MM,1:NXCHOP)= BE(:,:NXCHOP)+BO(:,:NXCHOP)

        ! IF(NXCHOP .NE. 1) THEN
        ! B%E(1:NRH,MM,NXCHOPH:NX)= BE(:,NXCHOP+1:)+BO(:,NXCHOP+1:)
        ! B%E(NR:NRH+1:-1,MM,NXCHOPH:NX)= BE(:,NXCHOP+1:)-BO(:,NXCHOP+1:)
        ! ENDIF

        ! B%E(NR:NRH+1:-1,MM,:NXCHOP)= BE(:,:NXCHOP)-BO(:,:NXCHOP)

        B%E(1:NRH,MM,:)= BE(:,:)+BO(:,:)
        B%E(NR:NRH+1:-1,MM,:)= BE(:,:)-BO(:,:)

    ENDDO
!$OMP END PARALLEL DO
    
    ! B%E(:,NTCHOP+1:,:)=0
    ! B%E(:,:NTCHOP,NXCHOP+1:NXCHOPH-1)=0
    B%E(:,MAX(1,NTCHOP-A%INTH+1):,:) = CMPLX(0.D0,0.D0)

    IF(A%INTH.EQ.0 .AND. A%INX.EQ.0) THEN
        IF(A%LN.NE.0.0) THEN
            WRITE(*,*) 'RTRAN:TO PHYSICAL SPACE THOUGH LOGTERM IS NONZERO.'
            WRITE(*,*) 'LOGTERM=',A%LN
            B%E(1:NR,1,1)=B%E(1:NR,1,1)+A%LN*TFM%LN(1:NR)
        ENDIF
    ENDIF

    CALL DEALLOCATE(A)

    A%E => B%E
    NULLIFY(B%E)
    A%SPACE = PFF_SPACE
    A%INR = B%INR
    A%INTH = B%INTH
    A%INX = B%INX

    DEALLOCATE( BE,BO )
    ENDIF
    
    CALL CHOPDO(A)
    RETURN
    END SUBROUTINE RTRAN2
!=======================================================================
    SUBROUTINE TOFF(A)
!=======================================================================
! [USAGE]: 
! TRANSFORM A FROM PPP SPACE TO FFF SPACE
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES IN PPP SPACE
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR):: A
    INTEGER:: N

    N = A%SPACE
    IF (N.EQ.PPP_SPACE) THEN
    CALL HORFFT(A,-1)
    CALL VERFFT(A,-1)
    CALL RTRAN(A,-1)
    ELSEIF (N.EQ.PFP_SPACE) THEN
    CALL VERFFT(A,-1)
    CALL RTRAN(A,-1)
    ELSEIF (N.EQ.PFF_SPACE) THEN
    CALL RTRAN(A,-1)
    ENDIF

    RETURN
    END SUBROUTINE TOFF
!=======================================================================
    SUBROUTINE TOFP(A)
!=======================================================================
! [USAGE]: 
! TRANSFORM A FROM FFF SPACE TO PPP SPACE
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES IN FFF SPACE
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR):: A
    INTEGER:: N

    N = A%SPACE
    IF (N.EQ.FFF_SPACE) THEN
    CALL RTRAN(A,1)
    CALL VERFFT(A,1)
    CALL HORFFT(A,1)
    ELSEIF (N.EQ.PFF_SPACE) THEN
    CALL VERFFT(A,1)
    CALL HORFFT(A,1)
    ELSEIF (N.EQ.PFP_SPACE) THEN
    CALL HORFFT(A,1)
    ENDIF
    
    RETURN
    END SUBROUTINE TOFP
! ======================================================================
    subroutine DECOMPOSE(NSIZE,NPROCS,PROC_NUM,NSIZE_PROC,INDEX_PROC)
! ======================================================================
! [USAGE]:
! DECOMPOSE 1D DOMAIN INTO ALL PROCESSORS
! [PARAMETERS]:
! NSIZE >> TOT NUM OF ELEMENTS ALONG THAT DIMENSION 
! NPROCS >> TOT NUM OF PROCS ALONG THAT DIMENSION
! PROC_NUM >> RANK OF THAT PROCESSOR
! NSIZE_PROC >> LOC NUM OF ELEMENTS ALONG THAT DIMENSION
! INDEX_PROC >> START INDEX OF THAT PROCESSOR ALONG THAT DIMENSION
! [NOTE]:
! INDEX STARTS FROM 0
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: NSIZE, NPROCS, PROC_NUM, NSIZE_PROC, INDEX_PROC
    INTEGER:: Q,R

    Q = NSIZE/NPROCS
    R = MOD(NSIZE,NPROCS)
    IF (R.GT.PROC_NUM) THEN
        NSIZE_PROC = Q+1
        INDEX_PROC = NSIZE_PROC*PROC_NUM
    ELSE
        NSIZE_PROC = Q
        INDEX_PROC = NSIZE_PROC*PROC_NUM + R
    ENDIF

    ! NSIZE_PROC = Q + (R.GT.PROC_NUM) ! DOESN'T WORK
    ! INDEX_PROC = Q * PROC_NUM + MIN(R,PROC_NUM)

    RETURN
    end subroutine DECOMPOSE
! ======================================================================
    subroutine SUBARRAY(ELEMENT_DATA_TYPE,NDIM,DATASIZE_PROC,DIM,NPROCS,NEW_DATA_TYPE)
! ======================================================================
! [USAGE]: CREATE SUBARRAY DATATYPES FOR EACH PROC. ORIGINAL ARRAY IS DE
! COMPOSED ALONG DIM INTO SUBARRAYS.
! [PARAMETERS]: 
! OLD_DATA_TYPE >> ORIGINAL MPI ELEMENTARY DATATYPE OF THE ARRAY
! NDIM >> NUMBER OF DIMENSIONS OF THE LOCAL ARRAY
! DATASIZE_PROC >> ORIGINAL DATASIZE OF THE LOCAL ARRAY
! DIM >> DIMENSION ALONG WHICH THE LOCAL ARRAY IS FURTHER DECOMPOSED
! NPROCS >> NUMBER OF PROCS ALONG THAT DIM
! NEW_DATA_TYPE >> MPI DERIVED DATATYPE FOR ALL PROCS (DIM = NPROCS)
! [NOTE]:
! 1. PROC INDEX STARTS FROM 0
! 2. DIMN INDEX STARTS FROM 1
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: ELEMENT_DATA_TYPE, NDIM, DIM, NPROCS
    INTEGER,DIMENSION(NDIM):: DATASIZE_PROC, SUBDSIZE_PROC, SUBSTARTS_PROC
    INTEGER,DIMENSION(0:NPROCS-1):: NEW_DATA_TYPE

    INTEGER:: I 

    DO I = 1,NDIM ! MEMO >>> MAY NEED TO CHANGE THIS TO 0~NDIM-1
        SUBDSIZE_PROC(I) = DATASIZE_PROC(I)
        SUBSTARTS_PROC(I) = 0
    ENDDO

    ! ALONG THAT DIM
    DO I = 0,NPROCS-1
        ! DECOMPOSE A DOMAIN OF DATASIZE_PROC(DIM) INTO NPROCS
        ! -> DECOMPOSED DOMAIN SIZE ALONG THAT DIM: NSIZE_PROC(DIM)
        ! -> OTHER DIMS KEEP THE ORIGINAL DOMAIN SIZE
        ! -> DECOMPOSED DOMAIN START INDEX ALONG THAT DIM: SUBSTARTS_PROC(DIM)
        ! -> OTHER DIMS KEEP THE ORIGINAL START INDEX
        CALL DECOMPOSE(DATASIZE_PROC(DIM),NPROCS,I,SUBDSIZE_PROC(DIM),SUBSTARTS_PROC(DIM))
        ! CREATE A NDIM-DIMENSION SUBARRAY DATATYPE FOR EACH PROC
        CALL MPI_TYPE_CREATE_SUBARRAY(NDIM,DATASIZE_PROC,SUBDSIZE_PROC,SUBSTARTS_PROC, & 
            MPI_ORDER_FORTRAN,ELEMENT_DATA_TYPE,NEW_DATA_TYPE(I),IERR)
        CALL MPI_TYPE_COMMIT(NEW_DATA_TYPE(I),IERR)
    ENDDO

    RETURN
    END subroutine SUBARRAY
! ======================================================================
    subroutine EXCHANGE_3DCOMPLEX(COMM,DATA_SIZE_PROC_OLD,ARRAY_PROC_OLD,DIM_OLD &
                                      ,DATA_SIZE_PROC_NEW,ARRAY_PROC_NEW,DIM_NEW)
! ======================================================================
! [USAGE]:
! FOR PROCESSOR GROUP 'COMM', SWAP DIM_OLD AND DIM_NEW.
! [PARAMETERS]:
! COMM >> PROCESSOR GROUP
! DATA_SIZE_PROC >> SIZE OF THE LOCAL ARRAY
! ARRAY_PROC >> LOCAL ARRAY TO BE SWAPPED
! DIM >> DIMENSIONS TO BE SWAPPED
! [NOTE]
! 1. USE THIS ONLY IF THE DATATYPE ONLY NEEDS TO BE USED ONCE
! 2. PROC INDEX STARTS FROM 0
! 3. DIMN INDEX STARTS FROM 1
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: NDIM = 3
    INTEGER:: COMM, DIM_OLD, DIM_NEW
    INTEGER,DIMENSION(3),INTENT(IN):: DATA_SIZE_PROC_OLD, DATA_SIZE_PROC_NEW
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(IN):: ARRAY_PROC_OLD
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(INOUT):: ARRAY_PROC_NEW

    INTEGER:: I , NPROC_COMM
    INTEGER,DIMENSION(:),ALLOCATABLE:: DATA_TYPE_OLD, DATA_TYPE_NEW, counts, displs

    ! DETERMINE THE NUMBER OF PROCS IN PROCESSOR GROUP: COMM
    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,IERR)
    ALLOCATE(DATA_TYPE_OLD(0:NPROC_COMM-1))
    ALLOCATE(DATA_TYPE_NEW(0:NPROC_COMM-1))
    ALLOCATE(counts(0:NPROC_COMM-1))
    ALLOCATE(displs(0:NPROC_COMM-1))

    ! CREATE SUBARRAY DATATYPES
    ! OLD LOCAL ARRAY: ARRAY_OLD OF SIZE_OLD IS DECOMPOSED INTO SUBARRAY ALONG DIM_OLD
    CALL SUBARRAY(MPI_DOUBLE_COMPLEX,NDIM,DATA_SIZE_PROC_OLD,DIM_OLD,NPROC_COMM,DATA_TYPE_OLD)
    ! NEW LOCAL ARRAY: ARRAY_NEW OF SIZE_NEW IS DECOMPOSED INTO SUBARRAY ALONG DIM_NEW
    CALL SUBARRAY(MPI_DOUBLE_COMPLEX,NDIM,DATA_SIZE_PROC_NEW,DIM_NEW,NPROC_COMM,DATA_TYPE_NEW)

    ! SWAP SUBARRAYS
    ! EACH SUBARRAY IS COUNTED AS ONE UNIT
    DO I = 0,NPROC_COMM-1
        counts(I) = 1; ! SWAP ONE SUBARRAY A TIME
        displs(I) = 0; ! DIRECTLY START FROM THE FIRST MEMORY LOC OF THE ARRAY
    enddo

    CALL MPI_ALLTOALLW(ARRAY_PROC_OLD, counts, displs, DATA_TYPE_OLD, &
                    ARRAY_PROC_NEW, counts, displs, DATA_TYPE_NEW, COMM, IERR)

    ! FREE THE DATATYPE
    DO I = 0,NPROC_COMM-1
        CALL MPI_TYPE_FREE(DATA_TYPE_OLD(I),IERR)
        CALL MPI_TYPE_FREE(DATA_TYPE_NEW(I),IERR)
    ENDDO

    DEALLOCATE(DATA_TYPE_OLD, DATA_TYPE_NEW, counts, displs)

    end subroutine EXCHANGE_3DCOMPLEX
! ======================================================================
    subroutine EXCHANGE_3DCOMPLEX_FAST(COMM,ARRAY_PROC_OLD,DATA_TYPE_OLD &
                                           ,ARRAY_PROC_NEW,DATA_TYPE_NEW)
! ======================================================================
! [USAGE]:
! FOR PROCESSOR GROUP 'COMM', SWAP DIM_OLD AND DIM_NEW.
! [PARAMETERS]:
! COMM >> PROCESSOR GROUP
! DATA_SIZE_PROC >> SIZE OF THE LOCAL ARRAY
! ARRAY_PROC >> LOCAL ARRAY TO BE SWAPPED
! DIM >> DIMENSIONS TO BE SWAPPED
! [NOTE]
! 1. USE THIS ONLY IF THE DATATYPE ONLY NEEDS TO BE USED ONCE
! 2. PROC INDEX STARTS FROM 0
! 3. DIMN INDEX STARTS FROM 1
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: NDIM = 3
    INTEGER:: COMM
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(IN):: ARRAY_PROC_OLD
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(INOUT):: ARRAY_PROC_NEW
    INTEGER,DIMENSION(:),INTENT(IN):: DATA_TYPE_OLD, DATA_TYPE_NEW

    INTEGER:: I , NPROC_COMM
    INTEGER,DIMENSION(:),ALLOCATABLE:: counts, displs, DATA_TYPE_OLD_COPY

    ! CHECK
    IF (RANK(ARRAY_PROC_OLD).NE.NDIM) THEN
    ! ARRAY_PROC_OLD must have NDIM ranks
    ! i.e NDIM = 3, ARRAY_PROC_OLD(:,:,:)
        CALL MPI_COMM_RANK(MPI_COMM_IVP, MPI_RANK, IERR)
        WRITE(*,*) 'ERROR in GLOBAL PROC#'//ITOA3(MPI_RANK)//': EXCHANGE - INCOMPATIBLE INPUT ARRAY DIMS'
    ENDIF
    ! IF (LBOUND(DATA_TYPE_OLD,1).NE.0) THEN
    ! ! DATA_TYPE_OLD(NPROC_COMM) should count from 0
    ! ! i.e DATA_TYPE_OLD(0), ... , DATA_TYPE_OLD(NPROC_COMM-1)
    !     CALL MPI_COMM_RANK(MPI_COMM_IVP,PROC_NUM,IERR)
    !     IF (PROC_NUM.EQ.0) THEN
    !         WRITE(*,*) 'ERROR in GLOBAL PROC#'//ITOA3(PROC_NUM)//': EXCHANGE - SUBARRAY DATATYPE MUST START FROM 0'
    !     ENDIF
    ! ENDIF

    ! SWAP SUBARRAYS
    ! EACH SUBARRAY IS COUNTED AS ONE UNIT
    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,IERR)
    ALLOCATE(counts(0:NPROC_COMM-1))
    ALLOCATE(displs(0:NPROC_COMM-1))
    ALLOCATE(DATA_TYPE_OLD_COPY(0:NPROC_COMM-1))
    DATA_TYPE_OLD_COPY(0:NPROC_COMM-1) = DATA_TYPE_OLD(:)

    DO I = 0,NPROC_COMM-1
    counts(I) = 1; ! SWAP ONE SUBARRAY A TIME
    displs(I) = 0; ! DIRECTLY START FROM THE FIRST MEMORY LOC OF THE ARRAY
    enddo

    CALL MPI_ALLTOALLW(ARRAY_PROC_OLD, counts, displs, DATA_TYPE_OLD, &
    ARRAY_PROC_NEW, counts, displs, DATA_TYPE_NEW, COMM, IERR)

    DEALLOCATE(counts, displs)

    end subroutine EXCHANGE_3DCOMPLEX_FAST
! ======================================================================
    subroutine SUBCOMM_CART(COMM,NDIM,SUBCOMMS)
! ======================================================================
! [USAGE]:
! CREATE COMMUNICATORS (SUBCOMMS(I)) FOR EACH DIMENSION(I) THAT HAS CART
! ESIAN TOPOLOGY. 
! [EXAMPLE]: 
! FOR N1xN2xN3 = NPROC_COMM (NDIM = 3) PROCESSORS,
! SUBCOMMS(1) >> N2XN3 GROUPS EACH HAS N1 PROCS
! SUBCOMMS(2) >> N1XN3 GROUPS EACH HAS N2 PROCS
! [NOTES]:
! EA PROC CALCULATES ITS CORRESPONDING SUBCOMMS. PROCS ON THE SAME LINE
! IN I-TH DIM WILL SHARE THE SAME SUBCOMMS(I)
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM
    INTEGER:: NDIM
    INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT)::SUBCOMMS

    INTEGER:: COMM_CART, NPROC_COMM, I 
    INTEGER,DIMENSION(1:NDIM):: DIMS
    LOGICAL,DIMENSION(1:NDIM):: PERIODS, REMDIMS

    DIMS = 0
    PERIODS = .FALSE.
    REMDIMS = .FALSE.
    ALLOCATE(SUBCOMMS(NDIM))

    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,IERR)

    ! CREATES A DIVISION OF PROCESSORS IN A CARTESIAN GRID
    ! DIMS: NUMBER OF PROCESSORS IN EACH DIM
    CALL MPI_DIMS_CREATE(NPROC_COMM,NDIM,DIMS,IERR)
    ! write(*,*) MPI_RANK,'-',dims(1),'x',dims(2)

    ! MAKES PROCESSOR GROUP THAT ATTACHES THE CARTESIAN TOPOLOGY
    ! COMM: ORIGINAL PROCESSOR GROUP
    ! NDIM: NUMBER OF DIMS IN CARTESIAN GRID
    ! DIMS: NUMBER OF PROCESSORS IN EACH DIM
    ! PERI: WHETHER EA DIM IS PERIODIC
    ! REOR: REORDER OR NOT
    ! COMM_CART: NEW PROCESSOR GROUP
    CALL MPI_CART_CREATE(COMM,NDIM,DIMS,PERIODS,.TRUE.,COMM_CART,IERR)

    ! CREATES SUB-PROCESSOR GROUP
    ! REMDIMS: TELLS WHICH DIM(S) IS KEPT IN THE SUBGRID
    ! E.X: FOR 3D, REMDIMS = .FALSE.,.TRUE.,.FALSE.
    !      THEN CREATES DIM1xDIM3 SUBGROUPS, EA WITH DIM2 PROCESSORS
    DO I = 1,NDIM
        REMDIMS(I) = .TRUE.
        CALL MPI_CART_SUB(COMM_CART,REMDIMS,SUBCOMMS(I),IERR)
        REMDIMS(I) = .FALSE.
    ENDDO

    ! FOR 2D GRID ONLY:
    ! MAKE SURE SUBCOMMS(1) IS ALWAYS BIGGER THAN SUBCOMMS(2)
    IF (NDIM.EQ.2) THEN
        IF (count_proc(SUBCOMMS(1)).LT.count_proc(SUBCOMMS(2))) THEN
            IERR = SUBCOMMS(1)
            SUBCOMMS(1) = SUBCOMMS(2)
            SUBCOMMS(2) = IERR
        ENDIF
    ENDIF
    ! write(*,*) MPI_RANK,'-',count_proc(SUBCOMMS(1)),'x',count_proc(SUBCOMMS(2)), &
    !     '-',local_proc(SUBCOMMS(1)), &
    !     '-',local_proc(SUBCOMMS(2))

    ! FUTURE EDIT NOTE:
    ! FOR NTH = 1 OR NX = 1, MAKE SUBCOMMS(I) HAVE ONE PROC IN EA GROUP
    ! FOR SLAB DECOMP, MAKE SUBCOMMS(I) HAVE ONE PROC IN EA GROUP
    ! (USE MPI_COMM_SPLIT)

    CALL MPI_COMM_FREE(COMM_CART,IERR)      

    end subroutine SUBCOMM_CART
! ======================================================================
    SUBROUTINE MASSEMBLE(LOCAL_ARRAY, GLOBAL_ARRAY, axis)
! ======================================================================
! [USAGE]: 
! ASSEMBLE COMPLEX LOCAL ARRAYS INTO GLOBAL ARRAY IN PROC#0
! [PARAMETERS]:
! axis >> axis along which the domain is NOT chopped    
! MPI_Datatype >> element type of local array (e.g. MPI_INTEGER)
! [UPDATES]:
! CODED BY JINGE WANG @ SEP 19 2021
! [NOTE1]:
! SUBCOMMS_L,axis,SUBCOMMS_R
! or axis,SUBCOMMS_L,SUBCOMMS_R
! DOES NOT work with SUBCOMMS_L,SUBCOMMS_R,axis (not needed in IVP)
! [NOTE2]:
! USE public SUBCOMM_1 and SUBCOMM_2. Need to change their definition to
! match the actual usage.
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(IN):: LOCAL_ARRAY
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(INOUT):: GLOBAL_ARRAY
    integer,INTENT(IN):: axis 

    INTEGER:: II, JJ, KK, N1_glb, N2_glb, N3_glb, N1_loc, N2_loc, N3_loc, ELEMENT_SIZE 
    INTEGER(KIND=MPI_ADDRESS_KIND):: EXTEND_SIZE
    INTEGER:: SUBCOMM_L, SUBCOMM_R, MPI_Datatype
    INTEGER:: SUBARRAY_TYPE, SUBARRAY_TYPE_resized, DISPLACEMENT_loc, RECVCOUNT_loc
    integer,DIMENSION(:),ALLOCATABLE:: DISPLACEMENT, RECVCOUNT
    COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: LOCAL_ARRAY2, SUBGLOBAL_ARRAY, GLOBAL_ARRAY_COPY

    N1_glb = SIZE(GLOBAL_ARRAY,1)
    N2_glb = SIZE(GLOBAL_ARRAY,2)
    N3_glb = SIZE(GLOBAL_ARRAY,3)

    N1_loc = SIZE(LOCAL_ARRAY,1)
    N2_loc = SIZE(LOCAL_ARRAY,2)
    N3_loc = SIZE(LOCAL_ARRAY,3)

    MPI_Datatype = MPI_DOUBLE_COMPLEX
    IF (CP8_SIZE.EQ.0) THEN
    CALL MPI_TYPE_SIZE(MPI_Double_Complex,CP8_SIZE,IERR)
    ENDIF
    ELEMENT_SIZE = CP8_SIZE

    
    IF (axis.EQ.1) THEN ! axis,SUBCOMMS_L,SUBCOMMS_R => FFF,PFF
        SUBCOMM_L = SUBCOMM_2
        SUBCOMM_R = SUBCOMM_1

        ALLOCATE(SUBGLOBAL_ARRAY(N1_glb,N2_glb,N3_loc))
        ALLOCATE(LOCAL_ARRAY2(N1_glb,N3_loc,N2_loc))

        LOCAL_ARRAY2 = RESHAPE(LOCAL_ARRAY,SHAPE(LOCAL_ARRAY2),ORDER = [1,3,2])

        CALL MPI_TYPE_VECTOR(N3_loc, N1_glb, N2_glb*N1_glb, &
                            MPI_Datatype, SUBARRAY_TYPE, IERR)
        EXTEND_SIZE = ELEMENT_SIZE*N1_glb
        RECVCOUNT_loc = N2_loc
        DISPLACEMENT_loc = local_index(N2_glb,SUBCOMM_L)

    ELSEIF (axis.EQ.2) THEN ! SUBCOMMS_L,axis,SUBCOMMS_R => PPP,PFP
        SUBCOMM_L = SUBCOMM_1 ! >>> SHOULD BE SUBCOMM_2
        SUBCOMM_R = SUBCOMM_2 ! >>> SHOULD BE SUBCOMM_1

        ALLOCATE(SUBGLOBAL_ARRAY(N2_glb,N1_glb,N3_loc))
        ALLOCATE(LOCAL_ARRAY2(N2_glb,N3_loc,N1_loc))

        DO KK = 1,N1_loc
            DO JJ = 1,N3_loc
                DO II = 1,N2_glb
                    LOCAL_ARRAY2(II,JJ,KK) = LOCAL_ARRAY(KK,II,JJ)
                ENDDO
            ENDDO
        ENDDO

        CALL MPI_TYPE_VECTOR(N3_loc, N2_glb, N2_glb*N1_glb, &
                            MPI_Datatype, SUBARRAY_TYPE, IERR)
        EXTEND_SIZE = ELEMENT_SIZE*N2_glb
        RECVCOUNT_loc = N1_loc
        DISPLACEMENT_loc = local_index(N1_glb,SUBCOMM_L)

    ELSE
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'MASSEMBLE: undefined for the current axis'
        RETURN
    ENDIF

    ! CREATE ELEMENTAL SUBARRAY OF LOCAL_ARRAY
    CALL MPI_TYPE_CREATE_RESIZED(SUBARRAY_TYPE, 0, EXTEND_SIZE, SUBARRAY_TYPE_resized, IERR)
    CALL MPI_TYPE_COMMIT(SUBARRAY_TYPE_resized,IERR)

    ! IN EA SUBCOMM_L, GATHER LOCAL_ARRAY TO FORM SUBGLOBAL_ARRAY in SUBCOMM_L's #0 PROC
    ALLOCATE(RECVCOUNT(0:count_proc(SUBCOMM_L)-1))
    ALLOCATE(DISPLACEMENT(0:count_proc(SUBCOMM_L)-1))
    CALL MPI_ALLGATHER(RECVCOUNT_loc,1,MPI_INTEGER,RECVCOUNT,1,MPI_INTEGER,SUBCOMM_L,IERR)
    CALL MPI_ALLGATHER(DISPLACEMENT_loc,1,MPI_INTEGER,DISPLACEMENT,1,MPI_INTEGER,SUBCOMM_L,IERR)
    CALL MPI_GATHERV(LOCAL_ARRAY2, N1_loc*N2_loc*N3_loc, MPI_Datatype, &
                    SUBGLOBAL_ARRAY, RECVCOUNT, DISPLACEMENT, SUBARRAY_TYPE_resized, &
                    0, SUBCOMM_L, IERR)

    ! IN THE SUBCOMM_R THAT CONTAINS SUBCOMM_Ls' #0 PROCs: 
    ! GATHER SUBGLOBAL_ARRAY TO FORM GLOBAL_ARRAY in SUBCOMM_R's #0 proc (=> GLOBAL #0 PROC)
    IF (local_proc(SUBCOMM_L).EQ.0) THEN

        DEALLOCATE(RECVCOUNT,DISPLACEMENT)
        ALLOCATE(RECVCOUNT(0:count_proc(SUBCOMM_R)-1))
        ALLOCATE(DISPLACEMENT(0:count_proc(SUBCOMM_R)-1))
        RECVCOUNT_loc = N1_glb*N2_glb*N3_loc
        DISPLACEMENT_loc = N1_glb*N2_glb*local_index(N3_glb,SUBCOMM_R)
        CALL MPI_ALLGATHER(RECVCOUNT_loc,1,MPI_INTEGER,RECVCOUNT,1,MPI_INTEGER,SUBCOMM_R,IERR)
        CALL MPI_ALLGATHER(DISPLACEMENT_loc,1,MPI_INTEGER,DISPLACEMENT,1,MPI_INTEGER,SUBCOMM_R,IERR)

        IF (axis.EQ.1) THEN
        
            CALL MPI_GATHERV(SUBGLOBAL_ARRAY, N1_glb*N2_glb*N3_loc, MPI_Datatype, &
                            GLOBAL_ARRAY, RECVCOUNT, DISPLACEMENT, MPI_Datatype, &
                            0, SUBCOMM_R, IERR)

            ! DEBUG:
            ! CALL MCAT(SUBGLOBAL_ARRAY(1,:,:))

        ELSEIF (axis.EQ.2) THEN

            ALLOCATE(GLOBAL_ARRAY_COPY(N2_glb,N1_glb,N3_glb))
            CALL MPI_GATHERV(SUBGLOBAL_ARRAY, N1_glb*N2_glb*N3_loc, MPI_Datatype, &
                            GLOBAL_ARRAY_COPY, RECVCOUNT, DISPLACEMENT, MPI_Datatype, &
                            0, SUBCOMM_R, IERR)
            GLOBAL_ARRAY = RESHAPE(GLOBAL_ARRAY_COPY,SHAPE(GLOBAL_ARRAY),ORDER = [2,1,3]) ! NEED TO REORDER
            DEALLOCATE(GLOBAL_ARRAY_COPY)

            ! DEBUG:
            ! CALL MCAT(SUBGLOBAL_ARRAY(:,3,:))        

        ENDIF    
    ENDIF                 

    CALL MPI_TYPE_FREE(SUBARRAY_TYPE_resized,IERR)
    DEALLOCATE(LOCAL_ARRAY2,SUBGLOBAL_ARRAY,RECVCOUNT,DISPLACEMENT)

    END SUBROUTINE MASSEMBLE
! ======================================================================
    SUBROUTINE MSAVE0(A,FN,GLB)
!=======================================================================
! [USAGE]: 
! WRAPPER OF MSAVEX. PASSING ARGUMENTS A AND FN INTO MSAVEX
! [PARAMETERS]:
! A >> SCALAR-TYPE VARIABLE TO SAVE
! FN >> FILENAME
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR):: A
    CHARACTER(LEN=*),INTENT(IN):: FN
    LOGICAL, OPTIONAL:: GLB

    IF (PRESENT(GLB)) THEN
        CALL MSAVEX(A,FN)
        WRITE(*,*) 'WRITTEN TO ',TRIM(FN)
    ELSE
        CALL MPISAVEX(A,FN)
    ENDIF

    RETURN
    END SUBROUTINE MSAVE0
!=======================================================================
    SUBROUTINE MSAVEX(A,FN)
!=======================================================================
! [USAGE]: 
! SAVE THE GLOBAL SCALAR-TYPE VARIABLE A INTO FN
! [PARAMETERS]:
! A >> SCALAR-TYPE VARIABLE TO SAVE
! FN >> FILENAME
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
! MUST BE CALLED IN ROOT PROC
!=======================================================================
    IMPLICIT NONE
    TYPE(SCALAR):: A
    CHARACTER(LEN=*):: FN

    INTEGER:: STATUS,MM,KK,I,PROC_NUM

    CALL MPI_COMM_RANK(MPI_COMM_IVP, PROC_NUM, IERR)
    IF (PROC_NUM.NE.0) THEN
        WRITE(*,*) 'MSAVEX_IN_SCALAR3: NOT IN ROOT PROC'
        STOP
    ENDIF

    OPEN(UNIT=7,FILE=FN,STATUS='UNKNOWN',&
        FORM='UNFORMATTED',IOSTAT=STATUS)

    IF(STATUS.NE.0) THEN
        WRITE(*,*) 'MSAVEX_IN_SCALAR3: FAILED TO OPEN ',TRIM(FN)
        STOP
    ENDIF

    WRITE(7) NR,NTH,NX
    WRITE(7) NRCHOP,NTCHOP,NXCHOP
    WRITE(7) A%SPACE,A%LN
    WRITE(7) ZLEN,ELL
    WRITE(7) MINC,MKLINK
    WRITE(7) SIZE(A%E,1),SIZE(A%E,2),SIZE(A%E,3)

    ! IF(A%SPACE.EQ.FFF_SPACE) THEN
    ! DO MM=1,NTCHOP
    !     WRITE(7) A%E(:NRCHOPS(MM),MM,:NXCHOPDIM)
    ! ENDDO
    ! ELSE
    ! WRITE(7) A%E(:NR,:NTH,:NX)
    ! ENDIF

    WRITE(7) A%E(:,:,:)

    CLOSE(7)

    RETURN
    END SUBROUTINE MSAVEX
!=======================================================================
    SUBROUTINE MPISAVEX(local_scalar,FN)
! ======================================================================
! [USAGE]: 
! SAVE THE LOCAL SCALAR-TYPE VARIABLE local_scalar INTO FN
! [PARAMETERS]:
! local_scalar >> SCALAR-TYPE VARIABLE TO SAVE
! FN >> FILENAME
! [UPDATES]:
! WRITTEN BY JINGE WANG @ SEP 22 2021
! SAVE TWO FILES:
! >>> FN.info for scalar info
! >>> FN for local_scalar%E
!=======================================================================
    TYPE(SCALAR):: local_scalar
    CHARACTER(LEN=*):: FN

    INTEGER:: PROC_NUM 
    INTEGER:: ARRAYSIZE, FILEBLK, LOCAL_ARRAYSIZE
    INTEGER:: STATUS, MPIFILE, INDEX
    INTEGER(KIND=MPI_ADDRESS_KIND):: DISPLACEMENT

    ! 1. save basic info:
    CALL MPI_COMM_RANK(MPI_COMM_IVP, PROC_NUM, IERR)
    IF (PROC_NUM.EQ.0) THEN
    
        OPEN(UNIT=7,FILE=TRIM(FN)//'.info',STATUS='UNKNOWN',&
        FORM='UNFORMATTED',IOSTAT=STATUS)
    
        IF(STATUS.NE.0) THEN
            WRITE(*,*) 'MPISAVEX_IN_SCALAR3: ROOT PROC FAILED TO OPEN ',TRIM(FN),'.info'
            STOP
        ENDIF
    
        WRITE(7) NDIMR,NDIMTH,NDIMX
        WRITE(7) NRCHOPDIM,NTCHOPDIM,NXCHOPDIM
        WRITE(7) local_scalar%SPACE,local_scalar%LN
        WRITE(7) ZLEN,ELL
        WRITE(7) MINC,MKLINK
    
        CLOSE(7)
        WRITE(*,*) 'WRITTEN TO ',TRIM(FN)
    ENDIF

    ! 2. determine element size:
    ! PPP: NDIMR/N1, NDIMTH, NDIMX/N2
    ! PFP: NDIMR/N1, NDIMTH, NDIMX/N2
    ! PFF: NDIMR, NTCHOPDIM/N2, NXCHOPDIM/N1
    ! FFF: NRCHOPDIM, NTCHOPDIM/N2, NXCHOPDIM/N1
    SELECT CASE(local_scalar%SPACE)
    CASE (PFF_SPACE)
        ARRAYSIZE = NDIMR * CEILING(REAL(NTCHOPDIM)/count_proc(SUBCOMM_2)) * CEILING(REAL(NXCHOPDIM)/count_proc(SUBCOMM_1))
    CASE (FFF_SPACE)
        ARRAYSIZE = NRCHOPDIM * CEILING(REAL(NTCHOPDIM)/count_proc(SUBCOMM_2)) * CEILING(REAL(NXCHOPDIM)/count_proc(SUBCOMM_1))
    ! PPP_SPACE OR PFP_SPACE:
    CASE DEFAULT
        ARRAYSIZE = CEILING(REAL(NDIMR)/count_proc(SUBCOMM_1)) * NDIMTH * CEILING(REAL(NDIMX)/count_proc(SUBCOMM_2))
    END SELECT
    LOCAL_ARRAYSIZE = SIZE(local_scalar%E,1)*SIZE(local_scalar%E,2)*SIZE(local_scalar%E,3)
    
    ! 3. determine Processor index
    INDEX = local_proc(SUBCOMM_2) + local_proc(SUBCOMM_1) * count_proc(SUBCOMM_2)
    
    ! 4. determine element size of MPI_DOUBLE_COMPLEX
    IF (CP8_SIZE.EQ.0) THEN
        CALL MPI_TYPE_SIZE(MPI_Double_Complex,CP8_SIZE,IERR)
    ENDIF
    
    ! 5. create FILEBLK (derived datatype)
    call MPI_TYPE_CONTIGUOUS(ARRAYSIZE, MPI_DOUBLE_COMPLEX, FILEBLK, IERR)
    call MPI_TYPE_COMMIT(FILEBLK, IERR)
    
    ! 4. save A%E in each proc
    DISPLACEMENT = INDEX*ARRAYSIZE*CP8_SIZE
    CALL MPI_FILE_OPEN(MPI_COMM_IVP, FN, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, MPIFILE, IERR)
    CALL MPI_FILE_SET_VIEW(MPIFILE, DISPLACEMENT, MPI_DOUBLE_COMPLEX, FILEBLK, &
                        'NATIVE', MPI_INFO_NULL, IERR)
    CALL MPI_FILE_WRITE_ALL(MPIFILE, local_scalar%E, LOCAL_ARRAYSIZE, MPI_DOUBLE_COMPLEX, STATUS, IERR)
    CALL MPI_FILE_CLOSE(MPIFILE, IERR)
    CALL MPI_TYPE_FREE(FILEBLK, IERR)
    
    ! DEBUG:
    ! IF (PROC_NUM.EQ.0) THEN
    !     CALL MCAT(REAL(LOCAL_ARRAYS(:,1,:)))
    ! ENDIF

    END subroutine MPISAVEX
! ======================================================================
    SUBROUTINE MLOAD0(FN,A,GLB)
!=======================================================================
! [USAGE]: 
! WRAPPER OF MLOADX. PASSING ARGUMENTS FN AND A INTO MLOADX
! [PARAMETERS]:
! FN >> FILENAME
! A >> SCALAR-TYPE VARIABLE WHERE THE LOADED VALUES FROM FN ARE STROED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    CHARACTER(LEN=*):: FN
    TYPE(SCALAR):: A
    LOGICAL, OPTIONAL:: GLB

    IF (PRESENT(GLB)) THEN
        CALL MLOADX(FN,A)
        WRITE(*,*) 'READ FROM ',TRIM(FN)
    ELSE
        CALL MPILOADX(FN,A)
    ENDIF

    RETURN
    END SUBROUTINE MLOAD0
!=======================================================================
    SUBROUTINE MLOADX(FN,A)
!=======================================================================
! [USAGE]: 
! READ VALUES FROM FN AND THEN STORE THEM INTO A SCALR-TYPE VARIABLE A
! [PARAMETERS]:
! FN >> FILENAME
! A >> SCALAR-TYPE VARIABLE WHERE THE LOADED VALUES FROM FN ARE STROED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
    IMPLICIT NONE
    CHARACTER(LEN=*):: FN
    TYPE(SCALAR):: A

    INTEGER:: STATUS,MM
    INTEGER:: IZ,ITH,IX
    INTEGER:: IZCHOP,ITCHOP,IXCHOP
    REAL(P8):: IZLEN,IELL
    INTEGER:: IINC,IKLINK,IOF
    INTEGER:: SIZE1,SIZE2,SIZE3

    OPEN(UNIT=7,FILE=FN,STATUS='OLD',&
        FORM='UNFORMATTED',IOSTAT=STATUS)

    IF(STATUS.NE.0) THEN
    WRITE(*,*) 'MLOADX_IN_SCALAR3: FAILED TO OPEN ',TRIM(FN)
    STOP
    ENDIF

    READ(7) IZ,ITH,IX
    READ(7) IZCHOP,ITCHOP,IXCHOP
    READ(7) A%SPACE,A%LN
    READ(7) IZLEN,IELL
    READ(7) IINC,IKLINK
    READ(7) SIZE1,SIZE2,SIZE3

    ! IF(A%SPACE.EQ.FFF_SPACE) THEN
    ! IOF = IZCHOP-NRCHOP
    ! CALL CHOPSET(IOF)
    ! A%E=0
    ! DO MM=1,ITCHOP
    ! READ(7) A%E(:NRCHOPS(MM),MM,:NXCHOPDIM)
    ! ENDDO
    ! CALL CHOPSET(-IOF)
    ! CALL CHOPDO(A)
    ! ELSE
    ! READ(7) A%E(:IZ,:ITH,:IX)
    ! ENDIF

    NULLIFY(A%E)
    ALLOCATE(A%E(SIZE1,SIZE2,SIZE3))
    A%E = 0.D0
    READ(7) A%E(:,:,:)
    A%INR = 0
    A%INTH = 0
    A%INX = 0

    CLOSE(7)

    IF(IZLEN.NE.ZLEN) THEN
    WRITE(*,*) 'MLOADX: ZLEN INCONSISTENT.'
    WRITE(*,*) 'PROG =',ZLEN
    WRITE(*,*) 'FILE =',IZLEN
    ENDIF

    IF(IELL.NE.ELL) THEN
    WRITE(*,*) 'MLOADX: ELL INCONSISTENT.'
    WRITE(*,*) 'PROG =',ELL
    WRITE(*,*) 'FILE =',IELL
    ENDIF

    IF(MINC.NE.IINC) THEN
    WRITE(*,*) 'MLOADX: MINC INCONSISTENT.'
    WRITE(*,*) 'PROG =',MINC
    WRITE(*,*) 'FILE =',IINC
    ENDIF

    IF(MKLINK.NE.IKLINK) THEN
    WRITE(*,*) 'MLOADX: MKLINK INCONSISTENT.'
    WRITE(*,*) 'PROG =',MKLINK
    WRITE(*,*) 'FILE =',IKLINK
    ENDIF

    IF (NRCHOP.NE.IZCHOP .OR. NTCHOP.NE.ITCHOP .OR. NXCHOP.NE.IXCHOP)&
    THEN
    WRITE(*,*) 'MLOADX: CHOPPING DATA INCONSISTENT.'
    WRITE(*,*) 'PROG NR,NT,NXCHOP=',NRCHOP,NTCHOP,NXCHOP
    WRITE(*,*) 'FILE NR,NT,NXCHOP=',IZCHOP,ITCHOP,IXCHOP
    ENDIF

    RETURN
    END SUBROUTINE MLOADX
! ======================================================================
    subroutine MPILOADX(FN,local_scalar)
! ======================================================================
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
        TYPE(SCALAR):: local_scalar
        CHARACTER(LEN=*):: FN

        INTEGER:: PROC_NUM 
        INTEGER:: ARRAYSIZE, FILEBLK, LOCAL_ARRAYSIZE
        INTEGER:: STATUS, MPIFILE, INDEX
        INTEGER(KIND=MPI_ADDRESS_KIND):: DISPLACEMENT

        INTEGER:: INDIMR, INDIMTH, INDIMX, INRCHOPDIM, INTCHOPDIM, INXCHOPDIM
        INTEGER:: IMINC, IMKLINK
        REAL(P8):: IZLEN, IELL
        
        ! 1. root proc collects all info and bcast to all procs
        CALL MPI_COMM_RANK(MPI_COMM_IVP, PROC_NUM, IERR)
        IF (PROC_NUM.EQ.0) THEN
        
            OPEN(UNIT=7,FILE=TRIM(FN)//'.info',STATUS='OLD',&
            FORM='UNFORMATTED',IOSTAT=STATUS)
        
            IF(STATUS.NE.0) THEN
                WRITE(*,*) 'MSAVEX_IN_SCALAR3: ROOT PROC FAILED TO READ ',TRIM(FN),'.info'
                STOP
            ENDIF
        
            READ(7) INDIMR,INDIMTH,INDIMX
            READ(7) INRCHOPDIM,INTCHOPDIM,INXCHOPDIM
            READ(7) local_scalar%SPACE,local_scalar%LN
            READ(7) IZLEN,IELL
            READ(7) IMINC,IMKLINK
            
            CLOSE(7)
            WRITE(*,*) 'MPIREAD FROM ',TRIM(FN)
    
            ! Check consistency
            IF(IZLEN.NE.ZLEN) THEN
                WRITE(*,*) 'MLOADX: ZLEN INCONSISTENT.'
                WRITE(*,*) 'PROG =',ZLEN
                WRITE(*,*) 'FILE =',IZLEN
            ENDIF
    
            IF(IELL.NE.ELL) THEN
                WRITE(*,*) 'MLOADX: ELL INCONSISTENT.'
                WRITE(*,*) 'PROG =',ELL
                WRITE(*,*) 'FILE =',IELL
            ENDIF
    
            IF(MINC.NE.IMINC) THEN
                WRITE(*,*) 'MLOADX: MINC INCONSISTENT.'
                WRITE(*,*) 'PROG =',MINC
                WRITE(*,*) 'FILE =',IMINC
            ENDIF
    
            IF(MKLINK.NE.IMKLINK) THEN
                WRITE(*,*) 'MLOADX: MKLINK INCONSISTENT.'
                WRITE(*,*) 'PROG =',MKLINK
                WRITE(*,*) 'FILE =',IMKLINK
            ENDIF
            ! ! DEBUG:
            ! WRITE(*,*) IZLEN
            ! WRITE(*,*) IELL
            ! WRITE(*,*) INDIMR,INDIMTH,INDIMX
            ! WRITE(*,*) local_scalar%SPACE,local_scalar%LN
        ENDIF
        CALL MPI_BCAST(INDIMR,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
        CALL MPI_BCAST(INDIMTH,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
        CALL MPI_BCAST(INDIMX,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
        CALL MPI_BCAST(INRCHOPDIM,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
        CALL MPI_BCAST(INTCHOPDIM,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
        CALL MPI_BCAST(INXCHOPDIM,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
        CALL MPI_BCAST(local_scalar%SPACE,1,MPI_INTEGER,0,MPI_COMM_IVP,IERR)
    
        ! 2. determine element size
        IF (ASSOCIATED(local_scalar%E)) NULLIFY(local_scalar%E)
        local_scalar%INR = 0
        local_scalar%INTH = 0
        local_scalar%INX = 0

        SELECT CASE(local_scalar%SPACE)
        CASE (PFF_SPACE)
            ARRAYSIZE = INDIMR * CEILING(REAL(INTCHOPDIM)/count_proc(SUBCOMM_2)) * CEILING(REAL(INXCHOPDIM)/count_proc(SUBCOMM_1))
            ALLOCATE(local_scalar%E(INDIMR,local_size(INTCHOPDIM,SUBCOMM_2),local_size(INXCHOPDIM,SUBCOMM_1)))
            local_scalar%INTH = local_index(INTCHOPDIM,SUBCOMM_2)
            local_scalar%INX  = local_index(INXCHOPDIM,SUBCOMM_1)

        CASE (FFF_SPACE)
            ARRAYSIZE = INRCHOPDIM * CEILING(REAL(INTCHOPDIM)/count_proc(SUBCOMM_2)) * CEILING(REAL(INXCHOPDIM)/count_proc(SUBCOMM_1))
            ALLOCATE(local_scalar%E(INRCHOPDIM,local_size(INTCHOPDIM,SUBCOMM_2),local_size(INXCHOPDIM,SUBCOMM_1)))
            local_scalar%INTH = local_index(INTCHOPDIM,SUBCOMM_2)
            local_scalar%INX  = local_index(INXCHOPDIM,SUBCOMM_1)

        ! PPP_SPACE OR PFP_SPACE:
        CASE DEFAULT
            ARRAYSIZE = CEILING(REAL(INDIMR)/count_proc(SUBCOMM_1)) * INDIMTH * CEILING(REAL(INDIMX)/count_proc(SUBCOMM_2))
            ALLOCATE(local_scalar%E(local_size(INDIMR,SUBCOMM_1),INDIMTH,local_size(INDIMX,SUBCOMM_2)))
            local_scalar%INR  = local_index(INDIMR,SUBCOMM_1)
            local_scalar%INX  = local_index(INDIMX,SUBCOMM_2)

        END SELECT
        LOCAL_ARRAYSIZE = SIZE(local_scalar%E,1)*SIZE(local_scalar%E,2)*SIZE(local_scalar%E,3)
        
        ! 3. determine Processor index
        INDEX = local_proc(SUBCOMM_2) + local_proc(SUBCOMM_1) * count_proc(SUBCOMM_2)
        
        ! 4. determine element size of MPI_DOUBLE_COMPLEX
        IF (CP8_SIZE.EQ.0) THEN
            CALL MPI_TYPE_SIZE(MPI_Double_Complex,CP8_SIZE,IERR)
        ENDIF
        
        ! 5. create FILEBLK (derived datatype)
        call MPI_TYPE_CONTIGUOUS(ARRAYSIZE, MPI_DOUBLE_COMPLEX, FILEBLK, IERR)
        call MPI_TYPE_COMMIT(FILEBLK, IERR)
        
        ! 6. save A%E in each proc
        CALL MPI_FILE_OPEN(MPI_COMM_IVP, FN, MPI_MODE_RDONLY, MPI_INFO_NULL, MPIFILE, IERR)
        DISPLACEMENT = INDEX*ARRAYSIZE*CP8_SIZE
        CALL MPI_FILE_SET_VIEW(MPIFILE, DISPLACEMENT, MPI_DOUBLE_COMPLEX, FILEBLK, &
                            'NATIVE', MPI_INFO_NULL, IERR)
        CALL MPI_FILE_READ_ALL(MPIFILE, local_scalar%E, LOCAL_ARRAYSIZE, MPI_DOUBLE_COMPLEX, STATUS, IERR)
        CALL MPI_FILE_CLOSE(MPIFILE, IERR)
        CALL MPI_TYPE_FREE(FILEBLK, IERR)

    END subroutine MPILOADX
! ======================================================================

! ======================================================================
!                           UTILITY FUNCTIONS                           
! ======================================================================
    function local_size(NSIZE,COMM)
! ======================================================================
! [USAGE]:
! CALCULATE THE LOCAL SIZE
! [NOTE]:
! 1. ASSUMES THE PROCESSOR GROUP COMM IS ALIGNED 1D.
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM, NSIZE
    INTEGER:: local_size

    INTEGER:: NPROC_COMM, RANK_COMM, INDEX_COMM 

    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,IERR)
    CALL MPI_COMM_RANK(COMM,RANK_COMM,IERR)
    CALL DECOMPOSE(NSIZE,NPROC_COMM,RANK_COMM,local_size,INDEX_COMM)

    end function local_size
! ======================================================================
function local_index(NSIZE,COMM)
! ======================================================================
! [USAGE]:
! CALCULATE THE LOCAL INDEX
! [NOTE]:
! 1. ASSUMES THE PROCESSOR GROUP COMM IS ALIGNED 1D.
! 2. INDEX COUNTS FROM 0
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM, NSIZE
    INTEGER:: local_index

    INTEGER:: NPROC_COMM, RANK_COMM, SIZE_COMM 

    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,IERR)
    CALL MPI_COMM_RANK(COMM,RANK_COMM,IERR)
    CALL DECOMPOSE(NSIZE,NPROC_COMM,RANK_COMM,SIZE_COMM,local_index)

    end function local_index
! ======================================================================
function local_proc(COMM)
! ======================================================================
! [USAGE]:
! CALCULATE THE LOCAL INDEX
! [NOTE]:
! 1. ASSUMES THE PROCESSOR GROUP COMM IS ALIGNED 1D.
! 2. INDEX COUNTS FROM 0
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM
    INTEGER:: local_proc

    CALL MPI_COMM_RANK(COMM,local_proc,IERR)

    end function local_proc
! ======================================================================
    function count_proc(COMM)
! ======================================================================
! [USAGE]:
! COUNT THE NUMBER OF PROCS IN THAT PROC_GROUP
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM, count_proc 

    CALL MPI_COMM_SIZE(COMM, count_proc, IERR)

    end function count_proc
! ======================================================================
    function create_new_type3D(COMM,DATA_SIZE_PROC_OLD,DIM_OLD,DATA_TYPE)
! ======================================================================
! [USAGE]:
! CREATE SUBARRAY DATATYPE FOR 3D ARRAY OF DATA_SIZE_PROC_OLD STORED IN
! COMMUNICATOR GROUP: COMM
! [NOTE]:
! COMPANION FUNC FOR EXCHANGE_3DCOMPLEX_FAST
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: NDIM = 3
    INTEGER:: COMM, DIM_OLD, DATA_TYPE
    INTEGER,DIMENSION(3),INTENT(IN):: DATA_SIZE_PROC_OLD
    INTEGER,DIMENSION(:),ALLOCATABLE:: create_new_type3D

    INTEGER:: I , NPROC_COMM

    ! DETERMINE THE NUMBER OF PROCS IN PROCESSOR GROUP: COMM
    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,IERR)
    ALLOCATE(create_new_type3D(0:NPROC_COMM-1))

    ! CREATE SUBARRAY DATATYPES
    ! OLD LOCAL ARRAY: ARRAY_OLD OF SIZE_OLD IS DECOMPOSED INTO SUBARRAY ALONG DIM_OLD
    CALL SUBARRAY(DATA_TYPE,NDIM,DATA_SIZE_PROC_OLD,DIM_OLD,NPROC_COMM,create_new_type3D)

    end function create_new_type3D
! ======================================================================
END MODULE MOD_FFT