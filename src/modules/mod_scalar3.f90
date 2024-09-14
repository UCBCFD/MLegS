MODULE MOD_SCALAR3 ! LEVEL 2 MODULE
! TYPE SETUP: SCALAR & TRANSFORM
      USE omp_lib
      USE MPI
      USE MOD_MISC, ONLY : P4,P8,PI,IU, CP8_SIZE                         ! LEVEL 0
!XUSE USE MOD_FD                                                         ! LEVEL 1
      USE MOD_EIG                                                        ! LEVEL 1
      USE MOD_LIN_LEGENDRE                                               ! LEVEL 1
!XUSE USE MOD_BANDMAT                                                    ! LEVEL 1
      IMPLICIT NONE
      PRIVATE
!=======================================================================
!============================ PARAMETERS ===============================
!=======================================================================
      ! 1.1> DIMENSIONS
      INTEGER, PUBLIC:: NR,NRH                                           ! # OF COLLOCATION POINTS (= SPECTRAL COEFFICIENTS IF IN FUNC. SPACE) IN THE RADIAL DIRECTION R (OR X WHEN UNMAPPED); NRH = NR/2
      INTEGER, PUBLIC:: NX                                               ! # OF COLLOCATION POINTS (= SPECTRAL COEFFICIENTS IF IN FUNC. SPACE) IN THE AXIAL DIRECTION Z
      INTEGER, PUBLIC:: NTH                                              ! # OF COLLOCATION POINTS (= SPECTRAL COEFFICIENTS IF IN FUNC. SPACE) IN THE AZIMUTHAL DIRECTION THETA
      INTEGER, PUBLIC:: NRCHOP,NTCHOP                                    ! CHOP INDICES IN RADIAL AND TH FOR TRUCATION IN ORDER TO LET THE MAIMUM DEGREE OF THE BASIS BE LIMITED.
      INTEGER, PUBLIC,DIMENSION(:),ALLOCATABLE:: NRCHOPS,NTCHOPS         ! CHOP INDICES IN RADIAL AND TH FOR EACH N OF P^M_N
      INTEGER, PUBLIC:: NXCHOP,NXCHOPH                                   ! CHOP INDEX IN Z DETERMINING THE MAXIMUM DEGREE USED IN THIS SPECTRAL ANALYSIS) ; NXCHOPH = CHOPPING LOCATION FOR THE CONJUGATE SIDE
      REAL(P8),PUBLIC:: ZLEN,ZLEN0                                       ! PERIOD IN THE AXIAL DIRECTION (NOTE: Z-PERIODICITY IS ASSUMED SO THAT WE USE FOURIER SPECTRAL METHOD IN Z. THE EVP CODE CHANGES ZLEN FOR FAST CALCULATION, SO THE ORIGINAL ZLEN IS SAVED AS ZLEN0 FOR ENERGY CALCULATION.)
      REAL(P8),PUBLIC:: ELL,ELL2                                         ! THE MAP PARAMTER ELL (MATSUSHIMA AND MARCUS, 1997); ELL2 = ELL**2
      INTEGER, PUBLIC:: MKLINK                                           ! MKLINK = 1 OR -1 IF M AND K ARE LINKED TOGETHER. OTHERWISE, SET 0 AS DEFAULT.
      INTEGER, PUBLIC:: MINC                                             ! THE NUMBER OF SYMMETRY OCCURRING IN THE AZIMUTHAL DIRECTION (E.G., SEEING 4 IDENTICAL PATTERNS IN THETA, THEN MINC = 4). DEFAULT IS 1 (NO INNER SYMMETRY).
      INTEGER, PUBLIC:: NDIMR,NDIMTH,NDIMX                               ! # OF DIMENSIONS IN EACH DIRECTION
      INTEGER, PUBLIC:: NRCHOPDIM,NTCHOPDIM,NXCHOPDIM                    ! # OF DIMENSIONS IN EACH DIRECTION AFTER CHOPPED
      INTEGER, PUBLIC,DIMENSION(:),ALLOCATABLE:: M                       ! ACCESSIBLE AZIMUTHAL WAVENUMBERS
      REAL(P8),PUBLIC,DIMENSION(:,:),ALLOCATABLE:: AK                    ! ACCESSIBLE AXIAL WAVENUMBERS GIVEN LEGENDRE ORDER N AND AZIMUTHAL WAVENUMBER M
      ! 1.2> MPI-SPECIFIC
      INTEGER, PUBLIC:: MPI_RANK, MPI_COMM_IVP, MPI_THREAD_MODE
      ! INTEGER, PUBLIC:: SUBCOMM_1, SUBCOMM_2                             ! SUBCOMMUNICATOR GROUPS FOR MPI EXCHANGE
      INTEGER, PUBLIC:: SUBCOMM_1 = 0, SUBCOMM_2 = 0                     ! SUBCOMMUNICATOR GROUPS FOR MPI EXCHANGE
      INTEGER, PUBLIC,DIMENSION(:),ALLOCATABLE:: &
                        TYPE_PFP0,TYPE_PFP1,TYPE_PFF0,TYPE_PFF1          ! SUBARRAY DATA TYPES FOR MPI EXCHANGE
      INTEGER, PUBLIC,DIMENSION(3):: &
                        SIZE_PFP0,SIZE_PFP1,SIZE_PFF0,SIZE_PFF1          ! SUBARRAY DATA TYPES FOR MPI EXCHANGE
      INTEGER, PUBLIC:: IERR                                             ! MPI ERROR STATUS

      ! 2> TYPES
      PUBLIC:: SCALAR                                                    ! ANY PHYSICAL QUANTITY IN INTEREST WILL BE DECLARED THORUGH THIS TYPE (U,V,W,P,ETC.)
      TYPE SCALAR
      COMPLEX(P8),DIMENSION(:,:,:),POINTER:: E                           ! PHYSICAL QUANTITY (OR FUNCTION COEFFICIENT IF IN FUNCTION SPACE) IN EACH NODE; ORDER CONVENTION = (Z, THETA, X)
      REAL(P8)                            :: LN                          ! LOCARITHMIC TERM UTILIZED FOR THE TOROIDAL-POLOIDAL DECOMPOSITION (EQ. 54 OF MATSUSHIMA AND MARCUS (1997))
      INTEGER                             :: SPACE                       ! TAG SHOWING IN WHICH SPACE THIS SCALAR-TYPE VALUE IS. SEE 4> SPACE TAG 
      INTEGER                             :: INR,INTH,INX                ! LOCAL ARRAY STARTING INDEX
      END TYPE

      ! 3> TRANSFORM KIT
      PUBLIC:: TRANSFORM                                                 ! PHY. -> FUNC. OR FUNC -> PHY. TRANSFORMATION
      TYPE TRANSFORM
      REAL(P8),   DIMENSION(:),    POINTER:: R,dR                        ! COLLOCATION POINTS IN R (EVALUATED VIA GAUSS-LEGENDRE QUADRATURE)
      REAL(P8),   DIMENSION(:),    POINTER:: X                           ! UNMAPPED COLLOCATION POINTS IN R (X=(R^2-ELL^2)/(R^2+ELL^2))
      REAL(P8),   DIMENSION(:),    POINTER:: W                           ! GAUSS-LEGENDRE WEIGHT. W(I) IS FOR X(I)
      REAL(P8),   DIMENSION(:,:),  POINTER:: NORM,LOGNORM                ! NORMALIZATION FACTORS OF P^M_N
      REAL(P8),   DIMENSION(:,:,:),POINTER:: PF                          ! LEGENDRE POLYNOMIAL VALUE TABLE FOR P^M_N(X_I) GIVEN M, N AND I
      REAL(P8),   DIMENSION(:),    POINTER:: AT0,AT1                     ! LEGENDRE POLYNOMIAL EVALUATION AT X=-1(R=0 -> AT0) AND X=1(R=INF -> AT1)
      REAL(P8),   DIMENSION(:),    POINTER:: LN                          ! MINUS LOG VALUE OF X 
      REAL(P8),   DIMENSION(:),    POINTER:: Z                           ! EQUISPACED COLLOCATION POINTS IN Z (FOR FOURIER SPECTRAL ANALYSIS)
      REAL(P8),   DIMENSION(:),    POINTER:: THR,THI,TH                  ! EQUISPACED COLLOCATION POINTS IN THETA (FOR FOURIER SPECTRAL ANALYSIS)
      ! COMPLEX(P8),DIMENSION(:),    POINTER:: EX                          ! STORES EXP((I*2*PI/N)*J) FOR FFT
      ! INTEGER(P8):: NMAX                                                 ! SET TO BE MAX(NTH,NX)*2
      END TYPE
      TYPE(TRANSFORM),PUBLIC:: TFM                                       ! UNIVERSAL TRANSFORM KIT DECLARED IN THIS CLASS ACROSS THE CODE

      ! 4> SPACE TAG
      INTEGER,PUBLIC,PARAMETER:: PPP_SPACE=0                             ! (X (OR R), THETA, Z) = (PHYS.          , PHYS.          , PHYS.          )
      INTEGER,PUBLIC,PARAMETER:: FFF_SPACE=1                             ! (X (OR R), THETA, Z) = (FUNC. (MAP_LEG), FUNC. (FOURIER), FUNC. (FOURIER))
      INTEGER,PUBLIC,PARAMETER:: PFP_SPACE=2                             ! (X (OR R), THETA, Z) = (PHYS.          , FUNC. (FOURIER), PHYS.          )
      INTEGER,PUBLIC,PARAMETER:: PFF_SPACE=3                             ! (X (OR R), THETA, Z) = (PHYS.          , FUNC. (FOURIER), FUNC. (FOURIER))
!=======================================================================
!======================== PUBLIC DECLARATION ===========================
!=======================================================================
      ! CHOP BASED ON THE INITIALZED CHOP VALUES IN EACH DIRECTION
      PUBLIC:: CHOPSET
      PUBLIC:: CHOPDO
      ! RETURNS THE FUNCTION VALUE AT INFINITY(CALCAT1) OR ORIGIN(0) (ONLY BE CALLED BY ROOT PROC)
      PUBLIC:: CALCAT1,CALCAT0                                           ! A VALUE IS RETURNED FOR EACH K.
      ! INTEGRATION OF A FUNCTION OVER A FULL-(~) OR HALF-(~H) DOMAIN (ONLY BE CALLED BY ROOT PROC)
      PUBLIC:: INTEG                                                     ! AVAILABLE ONLY WHEN SPACE = FFF
      PUBLIC:: INTEGH                                                    ! AVAILABLE ONLY WHEN SPACE = PFF
      ! PRODUCT & INTEGRATE F=A*B*(1-X)**2. OVER THE DOMAIN
      PUBLIC:: PRODCT                                                    ! AVAILABLE ONLY WHEN SPACE = PFF
      
      ! ! PERFORM THE FFT IN THETA (=HORFFT) AND Z (=VERFFT) DIRECTIONS
      ! PUBLIC:: HORFFT
      ! PUBLIC:: VERFFT
      ! ! TRANSFORMS AMONG PPP, PFP, PFF, AND FFF SPACES
      ! PUBLIC:: RTRAN                                                     ! FOR EXAMPLE, RTRAN(A,1) TRANSFORMS FROM FFF_SPACE(FUNCTION) TO PFF_SPACE(PHYSICAL).
      ! ! AUTOMATICALLY DO PROPER TRANSFORM BASED ON THE CURRENT SPACE
      ! PUBLIC:: TOFF,TOFP

CONTAINS
!=======================================================================
!============================ SUBROUTINES ==============================
!=======================================================================
      SUBROUTINE EXSET(EX, NPTS)
!=======================================================================
! [USAGE]: 
! CALCULATE EXP((i*2*PI/NPTS)*J) INTENSIVELY USED FOR FOURIER SPACE
! WHERE J RANGES FROM -NPTS/2 TO NPTS/2
! [PARAMETERS]:
! EX >> AN ARRAY WHERE THE EXPONENTS ARE STORED
! NPTS >> THE SIZE OF THE ARRAY EX
! [UPDATES]:
! RE-CODED IN MODERN FORTRAN FORMAT BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      INTEGER(P8):: NPTS
      COMPLEX(P8),DIMENSION(1:NPTS):: EX

      INTEGER(P8):: NH,J
      REAL(P8)   :: W

      W = 2*PI/FLOAT(NPTS)
      NH = NPTS/2

      DO J = 1,NH
        EX(J) = EXP(CMPLX(0.D0,DFLOAT(J-1)*W))
        EX(J+NH) = EXP(CMPLX(0.D0,DFLOAT(1-J)*W))
      ENDDO

      RETURN
      END SUBROUTINE EXSET
!=======================================================================
      SUBROUTINE CHOPSET(IOF)
!=======================================================================
! [USAGE]: 
! SETTING UP THE CHOP LOCATIONS FOR AZIMUTHAL AND AXIAL FUNCTION COEFFS.
! [PARAMETERS]:
! IOF >> OFFSET OF NRCHOP
! [NOTES]:
! SEE THE BELOW DIAGRAM TO UNDERSTAND HOW CHOPPING IS DONE.
!     NR-> +----------------------------------------+---------+
! NRCHOP-> +--+                                     +         +
! (NRDEG)  +  +-----+              INACTIVE         +    I    +
!          +        +-----+                         +    N    +
!          +              +-----+                   +    A    +
!          +                    +-----+ GRAD=-1     +    C    +
!          +                          +-----+       +    T    +
!          +                                +-----+ +    I    +
!          +         ACTIVE COEFFICIENTS          +-+    V    +
!          +                                        +    E    +
!          +                                        +         +
!          +----------------------------------------+---------+
!                                                   ↑         ↑ 
! [UPDATES]:                                   NTCHOP       NTH
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
      INTEGER:: IOF
      INTEGER:: MM,NN,NRDEG
      NRCHOP = NRCHOP + IOF
      NRDEG  = NRCHOP 
      IF(NRCHOP.GT.NRCHOPDIM) THEN
        IF (MPI_RANK.EQ.0) THEN
        WRITE(*,*) 'CHOPSET: NRCHOP TOO LARGE.'
        WRITE(*,*) '  NRCHOP=',NRCHOP,'  NRCHOPDIM=',NRCHOPDIM
        ENDIF
      STOP
      ENDIF

      DO MM=1,NTCHOPDIM
        NRCHOPS(MM) = MAX(MIN(NRCHOP,NRDEG-M(MM)),0)
      ENDDO
      NRCHOPS(NTCHOPDIM+1:)=0
      DO NN=1,NR
        NTCHOPS(NN) = MAX(MIN(NTCHOP,(NRDEG-NN+MINC)/MINC),0)
      ENDDO

      RETURN
      END SUBROUTINE CHOPSET
!=======================================================================
      SUBROUTINE CHOPDO(A)
!=======================================================================
! [USAGE]: 
! PERFORM CHOPPING FUNCTION SPACE COEFFFICIENTS.
! [PARAMETERS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN FUNCTION SPACE INFORMATION
!      HENCE, EXPECTED TO BE A%SPACE = FFF_SPACE
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
      TYPE(SCALAR):: A
      INTEGER::MM,NTCHOP_LOC,NR_LOC
   
      ! FFF_SPACE:
      IF (A%SPACE.EQ.FFF_SPACE) THEN
        NTCHOP_LOC = NTCHOP-A%INTH
        ! IF (NTCHOP_LOC.GT.0) THEN
        DO MM = 1,MIN(SIZE(A%E,2),NTCHOP_LOC)
          A%E(NRCHOPS(A%INTH+MM)+1:,MM,:) = CMPLX(0.D0,0.D0)
        ENDDO
        ! ENDIF
        A%E(:,MAX(1,NTCHOP_LOC+1):,:) = CMPLX(0.D0,0.D0)
        ! DO MM=1,NTCHOP
        !   A%E(NRCHOPS(MM)+1:,MM,:)= CMPLX(0.D0,0.D0)
        ! ENDDO
        ! A%E(:,NTCHOP+1:,:) = CMPLX(0.D0,0.D0)

      ! PFF_SPACE:
      ELSEIF (A%SPACE.EQ.PFF_SPACE) THEN
        A%E(NR+1,:,:)  = 0.D0

      ! PPP_SPACE OR PFP_SPACE:
      ELSE
        A%E(:,NTH+1,:) = 0.D0
        NR_LOC = NR - A%INR
        IF (NR_LOC .LT. SIZE(A%E,1)) THEN
          A%E(NR_LOC+1:,:,:) = 0.D0
        ENDIF
        ! A%E(:,NTH+1,:) = 0.D0
        ! A%E(NR+1,:,:)  = 0.D0
      ENDIF
   
      RETURN
      END SUBROUTINE CHOPDO
!=======================================================================
!       SUBROUTINE HORFFT(A,IS)
! !=======================================================================
! ! [USAGE]: 
! ! PERFORM HORIZONTAL FFT WITH RESPECT TO AZIMUTHAL (THETA) DIRECTION
! ! EXPECTED TO BE CALLED WHEN DATA IS IN PPP SPACE
! ! [PARAMETERS]:
! ! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES
! ! IS >> FORWARD OR BACKWARD FFT (1: BACKWARD, -1: FORWARD)
! ! [NOTES]:
! ! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
! !    CALL RTRAN(A,1) (FFF -> PFF)
! !    CALL VERFFT(A,1) (PFF -> PFP)
! !    CALL HORFFT(A,1) (PFP -> PPP)
! ! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
! !    CALL HORFFT(A,-1) (PPP -> PFP)
! !    CALL VERFFT(A,-1) (PFP -> PFF)
! !    CALL RTRAN(A,-1) (PFF -> FFF)
! ! [DEPENDENCIES]
! ! 1. DFFTW_PLAN_DFT                           @ FFTW3 
! ! 2. DFFTW_EXECUTE_DFT                        @ FFTW3
! ! 3. DFFTW_DESTROY_PLAN                       @ FFTW3
! ! 4. FFTW_FORWARD,FFTW_BACKWARD,FFTW_ESTIMATE @ FFTW3
! ! [UPDATES]:
! ! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! ! FFT SUBROUTINES ARE REPLACED WITH FFTW3 LIBRARY @ NOV NOV 11 2020
! !=======================================================================
!       IMPLICIT NONE
!       INCLUDE 'fftw3.f'
!       TYPE(SCALAR):: A
!       INTEGER,INTENT(IN):: IS

!       COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: B
!       REAL(P8),DIMENSION(:),ALLOCATABLE:: C
!       INTEGER:: II,KK,MM
!       REAL(P8):: FAC
!       INTEGER(P8):: PLAN

!       IF(IS .EQ.-1) THEN                                                 ! PHYSICAL TO FOURIER SPACE
!         IF(A%SPACE.NE.PPP_SPACE)THEN
!           WRITE(*,*) 'HORFFT: NOT IN PHYSICAL SPACE'
!         ENDIF
!         IF(NTH.EQ.1) THEN
!           A%SPACE=PFP_SPACE
!           RETURN
!         ENDIF

!         ALLOCATE( B(NDIMTH) )
!         ALLOCATE( C(2*NTH ) )
!         B = A%E(1,:NDIMTH,1)
!         DO II=1,NTH
!           C(2*II-1) = REAL(B(II))
!           C(2*II  ) = AIMAG(B(II))
!         ENDDO
!         CALL DFFTW_PLAN_DFT_R2C_1D(PLAN,2*NTH,C,B,FFTW_ESTIMATE)

!         DEALLOCATE( B )
!         DEALLOCATE( C )

!         ALLOCATE(B(NDIMTH))
!         ALLOCATE(C(2*NTH))
! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B,C)

!         DO II=1,NR                                                       ! THIS DO LOOP IS PARALLELIZABLE 
!           DO KK=1,NX
!             B = A%E(II,:NDIMTH,KK)
!             DO MM = 1,NTH
!               C(2*MM-1) = REAL(B(MM))
!               C(2*MM  ) = AIMAG(B(MM))
!             ENDDO
!             CALL DFFTW_EXECUTE_DFT_R2C(PLAN,C,B)
!             A%E(II,:NDIMTH,KK) = B/(2*NTH)
!           ENDDO
!         ENDDO
! !$OMP END PARALLEL DO
!         DEALLOCATE( B )
!         DEALLOCATE( C )

!         CALL DFFTW_DESTROY_PLAN(PLAN)

!         A%E(:,NTCHOP+1:,:) = CMPLX(0.D0,0.D0)                            ! CHOPPING IN THETA DIRECTION

!         A%E(:,1,:) = CMPLX(REAL(A%E(:,1,:)),0.D0)                        ! MAKE SURE TO ELIMINATE COMPLEX RESIDUAL DUE TO MACHINE ERROR

!         A%SPACE = PFP_SPACE                                              ! CHANGE SPACE TAG FROM PPP -> PFP

!       ELSE                                                               ! FOURIER TO PHYSICAL SPACE
!         IF(A%SPACE.NE.PFP_SPACE)THEN
!           WRITE(*,*) 'HORFFT: NOT IN PFP SPACE'
!         ENDIF
!         IF(NTH.EQ.1) THEN
!           A%SPACE=PPP_SPACE
!           RETURN
!         ENDIF

!         A%E(:,1,:) = CMPLX(REAL(A%E(:,1,:)),0.D0)                        ! MAKE SURE TO ELIMINATE COMPLEX RESIDUAL DUE TO MACHINE ERROR

!         A%E(:,NTCHOP+1:,:) = CMPLX(0.D0,0.D0)                            ! CHOPPING IN THETA DIRECTION

!         ALLOCATE( B(NDIMTH) )
!         ALLOCATE( C(2*NTH ) )
!         B = A%E(1,:NDIMTH,1)
!         C = 0.D0
!         CALL DFFTW_PLAN_DFT_C2R_1D(PLAN,2*NTH,B,C,FFTW_ESTIMATE)
                                            
!         DEALLOCATE( B )
!         DEALLOCATE( C )

!         ALLOCATE( B(NDIMTH) )
!         ALLOCATE( C(2*NTH ) )
! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B,C)
!         DO II=1,NR                                                       ! THIS DO LOOP IS PARALLELIZABLE 
!           DO KK=1,NX
!             B = A%E(II,:NDIMTH,KK)
!             C = 0.D0
!             CALL DFFTW_EXECUTE_DFT_C2R(PLAN,B,C)
!             DO MM = 1,NTH
!               A%E(II,MM,KK) = CMPLX(C(2*MM-1),C(2*MM),P8)
!             ENDDO
!             A%E(II,NDIMTH,KK) = CMPLX(0.D0, 0.D0)
!           ENDDO
!         ENDDO
! !$OMP END PARALLEL DO
!         DEALLOCATE( B )
!         DEALLOCATE( C )

!         CALL DFFTW_DESTROY_PLAN(PLAN)

!         A%SPACE = PPP_SPACE                                              ! CHANGE SPACE TAG FROM PFP -> PPP

!       ENDIF

!       RETURN
!       END SUBROUTINE HORFFT
! !=======================================================================
!       SUBROUTINE VERFFT(A,IS)
! !=======================================================================
! ! [USAGE]: 
! ! PERFORM VERTICAL FFT WITH RESPECT TO AXIAL (Z) DIRECTION
! ! EXPECTED TO BE CALLED WHEN DATA IS IN PFP SPACE
! ! [PARAMETERS]:
! ! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES
! ! IS >> FORWARD OR BACKWARD FFT (1: BACKWARD, -1: FORWARD)
! ! [NOTES]:
! ! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
! !    CALL RTRAN(A,1) (FFF -> PFF)
! !    CALL VERFFT(A,1) (PFF -> PFP)
! !    CALL HORFFT(A,1) (PFP -> PPP)
! ! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
! !    CALL HORFFT(A,-1) (PPP -> PFP)
! !    CALL VERFFT(A,-1) (PFP -> PFF)
! !    CALL RTRAN(A,-1) (PFF -> FFF)
! ! [DEPENDENCIES]
! ! 1. DFFTW_PLAN_DFT                           @ FFTW3 
! ! 2. DFFTW_EXECUTE_DFT                        @ FFTW3
! ! 3. DFFTW_DESTROY_PLAN                       @ FFTW3
! ! 4. FFTW_FORWARD,FFTW_BACKWARD,FFTW_ESTIMATE @ FFTW3
! ! [UPDATES]:
! ! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! ! FFT SUBROUTINES ARE REPLACED WITH FFTW3 LIBRARY @ NOV NOV 11 2020
! !=======================================================================
!       IMPLICIT NONE
!       INCLUDE 'fftw3.f'
!       TYPE(SCALAR):: A
!       INTEGER,INTENT(IN):: IS

!       COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: B
!       INTEGER:: JJ,KK
!       REAL(P8):: FAC
!       INTEGER(P8):: PLAN

!       IF(IS .EQ.-1) THEN                                                 ! PHYSICAL TO FORUIER SPACE
!         IF(A%SPACE.NE.PFP_SPACE) THEN
!           WRITE(*,*) 'VERFFT: NOT IN PFP SPACE'
!         ENDIF
!         IF(NX.EQ.1) THEN
!           A%SPACE=PFF_SPACE
!           RETURN
!         ENDIF

!         ALLOCATE( B(NX) )
!         B = A%E(1,1,:NX)
!         CALL DFFTW_PLAN_DFT_1D(PLAN,NX,B,B,FFTW_FORWARD,&
!                                           FFTW_ESTIMATE)
!         DEALLOCATE( B )

!         ALLOCATE( B(NX) )
! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B)

!         DO JJ=1,NTCHOP                                                   !!!!! THIS DO LOOP IS PARALLELIZABLE 
!           DO KK=1,NR
!             B = A%E(KK,JJ,:NX)
!             CALL DFFTW_EXECUTE_DFT(PLAN,B,B)
!             A%E(KK,JJ,:NX) = B/NX
!           ENDDO
!         ENDDO
! !$OMP END PARALLEL DO
!         DEALLOCATE( B )

!         CALL DFFTW_DESTROY_PLAN(PLAN)
  
!         A%E(:NR,:NTCHOP,NXCHOP+1:NXCHOPH-1) = CMPLX(0.D0,0.D0)           ! CHOPPING IN X DIRECTION
  
!         A%SPACE = PFF_SPACE                                              ! CHANGE SPACE TAG FROM PFP -> FFP

!       ELSE                                                               ! FOURIER TO PHYSICAL
!         IF(A%SPACE.NE.PFF_SPACE) THEN
!           WRITE(*,*) 'VERFFT: NOT IN PFF SPACE'
!         ENDIF
!         IF(NX.EQ.1) THEN
!           A%SPACE=PFP_SPACE
!           RETURN
!         ENDIF

!         A%E(:NR,:NTCHOP,NXCHOP+1:NXCHOPH-1) = CMPLX(0.D0,0.D0)           ! CHOPPING IN X DIRECTION

!         ALLOCATE( B(NTH) )
!         B = A%E(1,1,:NX)
!         CALL DFFTW_PLAN_DFT_1D(PLAN,NX,B,B,FFTW_BACKWARD,&
!                                            FFTW_ESTIMATE)
!         DEALLOCATE( B )

!         ALLOCATE( B(NX) )
! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(B)
!         DO JJ=1,NTCHOP                                                   !!!!! THIS DO LOOP IS PARALLELIZABLE 
!           DO KK=1,NR
!             B = A%E(KK,JJ,:NX)
!             CALL DFFTW_EXECUTE_DFT(PLAN,B,B)
!             A%E(KK,JJ,:NX) = B
!           ENDDO
!         ENDDO
! !$OMP END PARALLEL DO
!         DEALLOCATE( B )

!         CALL DFFTW_DESTROY_PLAN(PLAN)
  
!         A%SPACE = PFP_SPACE                                              ! CHANGE SPACE TAG FROM FFP -> PFP
      
!       ENDIF

!       RETURN
!       END SUBROUTINE VERFFT
! !=======================================================================
!       SUBROUTINE RTRAN(A,IS)
! !=======================================================================
! ! [USAGE]: 
! ! PERFORM MAPPED LEGENDRE TRANSFORM WITH RESPECT TO X (OR R) DIRECTION
! ! [PARAMETERS]:
! ! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES
! ! IS >> FORWARD OR BACKWARD FFT (1: BACKWARD, -1: FORWARD)
! ! [NOTES]:
! ! 1. WHEN GOING FROM FFF SPACE TO PPP SPACE, EXPECTED CALLING SEQUENCE =
! !    CALL RTRAN(A,1) (FFF -> PFF)
! !    CALL VERFFT(A,1) (PFF -> PFP)
! !    CALL HORFFT(A,1) (PFP -> PPP)
! ! 2. WHEN GOING FROM PPP SPACE TO FFF SPACE, EXPECTED CALLING SEQUENCE =
! !    CALL HORFFT(A,-1) (PPP -> PFP)
! !    CALL VERFFT(A,-1) (PFP -> PFF)
! !    CALL RTRAN(A,-1) (PFF -> FFF)
! ! [DEPENDENCIES]
! ! 1. SALLOC(~) @ MOD_SCALAR3
! ! 2. CHOPDO(~) @ MOD_SCALAR3
! ! 3. SFREE(~)  @ MOD_SCALAR3
! ! 4. OPERATOR(.MUL.) @ MOD_EIG
! ! [UPDATES]:
! ! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! ! FFT SUBROUTINES ARE REPLACED WITH FFTW3 LIBRARY @ NOV NOV 11 2020
! !=======================================================================
!       IMPLICIT NONE
!       TYPE(SCALAR),INTENT(INOUT):: A
!       INTEGER,INTENT(IN):: IS

!       TYPE(SCALAR):: B
!       COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: BE,BO
!       INTEGER:: NN,I,MM

!       IF(IS.EQ.-1) THEN                                                  ! PHYSICAL TO MAPPED LEGENDRE SPACE
!         IF(A%SPACE .NE. PFF_SPACE) THEN
!           WRITE(*,*) 'RTRAN: NOT IN PFF_SPACE.'
!         ENDIF

!         CALL SALLOC(B,FFF_SPACE)
!         CALL CHOPDO(B)

!         ALLOCATE( BE(NRH,NXCHOPDIM) )
!         ALLOCATE( BO(NRH,NXCHOPDIM) )

! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NN,BE,BO)
!         DO MM=1,NTCHOP

!           DO I=1,NRH
!             BE(I,:NXCHOP) = (A%E(I,MM,:NXCHOP)+A%E(NR-I+1,MM,:NXCHOP))&
!                             *TFM%W(I)
!             BO(I,:NXCHOP) = (A%E(I,MM,:NXCHOP)-A%E(NR-I+1,MM,:NXCHOP))&
!                             *TFM%W(I)

!             IF(NXCHOP .NE. 1) THEN
!               BE(I,NXCHOP+1:) = (A%E(I,MM,NXCHOPH:NX)+A%E(NR-I+1,MM,&
!                                  NXCHOPH:NX))*TFM%W(I)  
!               BO(I,NXCHOP+1:) = (A%E(I,MM,NXCHOPH:NX)-A%E(NR-I+1,MM,&
!                                  NXCHOPH:NX))*TFM%W(I)
!             ENDIF
!           ENDDO

!           NN = NRCHOPS(MM)

!           IF(NN.GE.1) THEN
!             B%E(1:NN:2,MM,:)=TRANSPOSE(TFM%PF(:NRH,1:NN:2,MM)) .MUL. BE
!           ENDIF
          
!           IF (NN.GE.2) THEN
!             B%E(2:NN:2,MM,:)=TRANSPOSE(TFM%PF(:NRH,2:NN:2,MM)) .MUL. BO
!           ENDIF

!         ENDDO
! !$OMP END PARALLEL DO

!         CALL DEALLOCATE(A)

!         A%E => B%E
!         NULLIFY(B%E)
!         A%SPACE=FFF_SPACE

!         DEALLOCATE( BE,BO )
!       ELSE                                                               ! MAPPED LEGENDRE SPACE TO PHYSICAL
!         IF(A%SPACE .NE. FFF_SPACE) THEN
!           WRITE(*,*) 'RTRAN: NOT IN FFF_SPACE'
!         ENDIF

!         CALL SALLOC(B,PPP_SPACE)
!         CALL CHOPDO(B)
!         B%LN=0

!         ALLOCATE( BE(NRH,NXCHOPDIM) )
!         ALLOCATE( BO(NRH,NXCHOPDIM) )

! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NN,BE,BO)
!         DO MM=1,NTCHOP
!           NN = NRCHOPS(MM)

!           IF(NN.GE.1) THEN
!             BE= TFM%PF(:NRH,1:NN:2,MM) .MUL. A%E(1:NN:2,MM,:NXCHOPDIM)
!           ELSE
!             BE=0
!           ENDIF

!           IF(NN.GE.2) THEN
!             BO= TFM%PF(:NRH,2:NN:2,MM) .MUL. A%E(2:NN:2,MM,:NXCHOPDIM)
!           ELSE
!             BO=0
!           ENDIF

!           B%E(1:NRH,MM,1:NXCHOP)= BE(:,:NXCHOP)+BO(:,:NXCHOP)

!           IF(NXCHOP .NE. 1) THEN
!             B%E(1:NRH,MM,NXCHOPH:NX)= BE(:,NXCHOP+1:)+BO(:,NXCHOP+1:)
!             B%E(NR:NRH+1:-1,MM,NXCHOPH:NX)= BE(:,NXCHOP+1:)-BO(:,NXCHOP+1:)
!           ENDIF
  
!           B%E(NR:NRH+1:-1,MM,:NXCHOP)= BE(:,:NXCHOP)-BO(:,:NXCHOP)
!         ENDDO
! !$OMP END PARALLEL DO
        
!         B%E(:,NTCHOP+1:,:)=0
!         B%E(:,:NTCHOP,NXCHOP+1:NXCHOPH-1)=0
  
!         IF(A%LN.NE.0.0) THEN
!           WRITE(*,*) 'RTRAN:TO PHYSICAL SPACE THOUGH LOGTERM IS NONRERO.'
!           WRITE(*,*) 'LOGTERM=',A%LN
!           B%E(1:NR,1,1)=B%E(1:NR,1,1)+A%LN*TFM%LN(1:NR)
!         ENDIF

!         CALL SFREE(A)

!         A%E => B%E
!         NULLIFY(B%E)
!         A%SPACE = PFF_SPACE
  
!         DEALLOCATE( BE,BO )
!       ENDIF
      
!       RETURN
!       END SUBROUTINE RTRAN
!=======================================================================
!       SUBROUTINE TOFF(A)
! !=======================================================================
! ! [USAGE]: 
! ! TRANSFORM A FROM PPP SPACE TO FFF SPACE
! ! [PARAMETERS]:
! ! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES IN PPP SPACE
! ! [UPDATES]:
! ! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! !=======================================================================
!       IMPLICIT NONE
!       TYPE(SCALAR):: A
!       INTEGER:: N

!       N = A%SPACE
!       IF (N.EQ.PPP_SPACE) THEN
!         CALL HORFFT(A,-1)
!         CALL VERFFT(A,-1)
!         CALL RTRAN(A,-1)
!       ELSEIF (N.EQ.PFP_SPACE) THEN
!         CALL VERFFT(A,-1)
!         CALL RTRAN(A,-1)
!       ELSEIF (N.EQ.PFF_SPACE) THEN
!         CALL RTRAN(A,-1)
!       ENDIF

!       RETURN
!       END SUBROUTINE TOFF
! !=======================================================================
!       SUBROUTINE TOFP(A)
! !=======================================================================
! ! [USAGE]: 
! ! TRANSFORM A FROM FFF SPACE TO PPP SPACE
! ! [PARAMETERS]:
! ! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES IN FFF SPACE
! ! [UPDATES]:
! ! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! !=======================================================================
!       IMPLICIT NONE
!       TYPE(SCALAR):: A
!       INTEGER:: N

!       N = A%SPACE
!       IF (N.EQ.FFF_SPACE) THEN
!         CALL RTRAN(A,1)
!         CALL VERFFT(A,1)
!         CALL HORFFT(A,1)
!       ELSEIF (N.EQ.PFF_SPACE) THEN
!         CALL VERFFT(A,1)
!         CALL HORFFT(A,1)
!       ELSEIF (N.EQ.PFP_SPACE) THEN
!         CALL HORFFT(A,1)
!       ENDIF
      
!       RETURN
!       END SUBROUTINE TOFP
!=======================================================================


      
!=======================================================================
!============================ FUNCTIONS ================================
!=======================================================================
      FUNCTION CALCAT0(A)
!=======================================================================
! [USAGE]:
! RETURNS THE FUNCTION VALUE AT ORIGIN (R = 0)
! [INPUTS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES IN FFF SPACE
! [OUTPUTS]:
! CALCAT0 >> FUNCTION VALUE AT ORIGIN (R = 0)
! [NOTE]:
! SHOULD ONLY CALLED BY PROCS THAT HAVE M=0 MODE
!=======================================================================
      TYPE(SCALAR),INTENT(IN):: A

      COMPLEX(P8),DIMENSION(SIZE(A%E,3)):: CALCAT0

      IF (A%INTH.NE.0) THEN
        WRITE(*,*) 'CALCAT0: INPUT MUST CONTAIN M=0 MODE'
      ENDIF

      CALCAT0 = TRANSPOSE(A%E(:NRCHOP,1,:)).MUL. TFM%AT0(:NRCHOP)

      RETURN
      END FUNCTION CALCAT0
!=======================================================================
      FUNCTION CALCAT1(A)
!=======================================================================
! [USAGE]:
! RETURNS THE FUNCTION VALUE AT INFINITY (R = INF)
! [INPUTS]:
! A >> TYPE(SCALAR) VARIABLE CONTAIN (X,THETA,Z) VALUES IN FFF SPACE
! [OUTPUTS]:
! CALCAT0 >> FUNCTION VALUE AT ORIGIN (R = INF)
! [NOTE]:
! SHOULD ONLY CALLED BY PROCS THAT HAVE M=0 MODE
! SAFE TO CALL IN EA LOCALARRAY WITH INTH = 0
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
      TYPE(SCALAR),INTENT(IN):: A

      COMPLEX(P8),DIMENSION(SIZE(A%E,3)):: CALCAT1

      IF (A%INTH.NE.0) THEN
        WRITE(*,*) 'CALCAT1: INPUT MUST CONTAIN M=0 MODE'
      ENDIF

      CALCAT1 = TRANSPOSE(A%E(:NRCHOP,1,:)) .MUL. TFM%AT1(:NRCHOP)

      RETURN
      END FUNCTION CALCAT1
!=======================================================================
      FUNCTION INTEG(F)
!=======================================================================
! [USAGE]:
! INTEGRATION OF A FUNCTION OVER A DOMAIN 
!     0  <  R   < INFTY
!     0  <  PHI < 2*PI
!     0  <  Z   < ZLEN
! CALL IN F/F SPACE
! INTEGRATED FUNCTION G IS IN THE FORM G = F*(1-MU)^2
! [INPUTS]:
! F >> A FUNTION INPUT FOR INTEGRATION (ACTUAL INTEGRATION IS DONE ON G)
! [OUTPUTS]:
! INTEG >> INTEGRATION EQUAL TO INT_[FULL_RANGE](G(R,PHI,Z) DR DPHI DZ)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 18 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR),INTENT(IN):: F
      REAL(P8)::INTEG

      IF(F%SPACE.NE.FFF_SPACE) THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'INTEG: NOT IN FFF_SPACE'
        STOP
      ENDIF      

      IF ((F%INTH.EQ.0).AND.(F%INX.EQ.0)) THEN
        IF(F%LN.NE.0.0) THEN
          WRITE(*,*) 'INTEG: LOGTERM NOT ZERO'
          WRITE(*,*) 'LOGTERM=',F%LN
        ENDIF
        INTEG = 4*PI*ZLEN*ELL2*REAL(F%E(1,1,1))*TFM%NORM(1,1)
      ELSE
        INTEG = 0
      ENDIF
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, INTEG, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_IVP, IERR)

      RETURN
      END FUNCTION INTEG
!=======================================================================
      FUNCTION PRODCT(A,B)
!=======================================================================
! [USAGE]:
! CALCULATE THE PRODUCT AND INTEGRATE OVER THE DOMAIN
! CALL IN R-PHYSICAL / PHI,Z-FOURIER SPACE
! WHAT IS ACTUALLY INTEGRATED IS G = A*B*(1-MU)^2
! [INPUTS]:
! A >> FIRST SCALAR-TYPE VARIABLE FOR INTEGRATION
! B >> SECOND SCALAR-TYPE VARIABLE FOR INTEGRATION
! [OUTPUTS]:
! PRODCT >> INTEGRATION EQUAL TO INT_[FULL_RANGE](G(R,PHI,Z) DR DPHI DZ)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 18 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR),INTENT(IN):: A,B
      REAL(P8):: PRODCT

      COMPLEX(P8):: PROD(NR)
      INTEGER:: MM,KK

      IF(A%SPACE.NE.PFF_SPACE .OR. B%SPACE.NE.PFF_SPACE) THEN
        IF (MPI_RANK.EQ.0) THEN
          WRITE(*,*) 'PRODCT:NOT IN PFF_SPACE'
          WRITE(*,*) 'A%SPACE,B%SPACE=',A%SPACE,B%SPACE
        ENDIF
        STOP
      ENDIF

      IF(A%LN.NE.0.0 .OR. B%LN.NE.0.0) THEN
        WRITE(*,*) 'PRODCT:LOGTERM NOT ZERO'
      ENDIF

      ! PROD=0
      ! DO KK = 1,NX
      !   IF(KK.GT.NXCHOP .AND. KK.LT.NXCHOPH) CYCLE
      !   PROD = PROD+A%E(:NR,1,KK)*CONJG(B%E(:NR,1,KK))
      ! ENDDO

      ! DO KK = 1,NX
      !   IF(KK.GT.NXCHOP .AND. KK.LT.NXCHOPH) CYCLE
      !   DO MM = 2,NTCHOP
      !     PROD = PROD + 2*(REAL(A%E(:NR,MM,KK))*REAL(B%E(:NR,MM,KK)) &
      !                 + AIMAG(A%E(:NR,MM,KK))*AIMAG(B%E(:NR,MM,KK)))
      !   ENDDO
      ! ENDDO

      PROD=0
      DO KK = 1,SIZE(A%E,3)
        DO MM = 1,SIZE(A%E,2) !NTCHOP
          IF (MM+A%INTH.EQ.1) THEN ! M(1) = 0
          PROD = PROD+A%E(:NR,1,KK)*CONJG(B%E(:NR,1,KK))
          ELSE
          PROD = PROD + 2*(REAL(A%E(:NR,MM,KK))*REAL(B%E(:NR,MM,KK)) &
                      + AIMAG(A%E(:NR,MM,KK))*AIMAG(B%E(:NR,MM,KK)))
          ENDIF
        ENDDO
      ENDDO
      PRODCT = SUM((PROD*TFM%W)*TFM%PF(1,1,1))

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, PRODCT, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_IVP, IERR)
      PRODCT = 4*PI*ZLEN*ELL2*PRODCT*TFM%NORM(1,1)

      RETURN
      END FUNCTION PRODCT
!=======================================================================
      FUNCTION INTEGH(F,AIN,BIN)
!=======================================================================
! [USAGE]:
! INTEGRATION OF A FUNCTION OVER A PART OF THE FULL DOMAIN 
!     0   <  R   < INFTY
!     AIN <  PHI < BIN
!     0   <  Z   < ZLEN
! CALL IN P/F SPACE, MEANING THAT F IS IN HAT FORM
! INTEGRATED FUNCTION G IS IN THE FORM G = F*(1-MU)^2
! [INPUTS]:
! F >> A FUNTION INPUT FOR INTEGRATION (ACTUAL INTEGRATION IS DONE ON G)
! AIN >> (OPTIONAL) INFIMUM OF PHI. DEFAULT IS -PI/2.
! BIN >> (OPTIONAL) SUPREMUM OF PHI. DEFAULT IS PI/2.
! [OUTPUTS]:
! INTEGH >> INTEGRATION EQUAL TO INT_[PART_RANGE](G(R,PHI,Z) DR DPHI DZ)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 18 2020
! MPI-ED BY JINGE WANG @ SEP 29 2021
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR),INTENT(IN):: F
      REAL(P8),INTENT(IN),OPTIONAL:: AIN, BIN
      REAL(P8):: INTEGH

      COMPLEX(P8):: FAB
      ! COMPLEX(P8),PARAMETER:: IU = (0.0D0,1.0D0) !Declared in mod_misc.f90
      REAL(P8):: A,B
      INTEGER:: MM,M0

      IF(F%SPACE.NE.PFF_SPACE) THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'INTEGH: NOT IN PFF_SPACE'
        STOP
      ENDIF

      IF(PRESENT(AIN)) THEN
        A=AIN
      ELSE
        A=-PI/2                                                          ! DEFAULT: RIGHT HAND SIDE HALF PLANE
      ENDIF

      IF(PRESENT(BIN)) THEN
        B=BIN
      ELSE
        B=PI/2
      ENDIF

      ! INTEGH = (B-A)*SUM(TFM%W*F%E(:NR,1,1))

      ! DO MM=2,NTCHOP
      !   M0=M(MM)
      !   FAB = (EXP(IU*M0*B)-EXP(IU*M0*A))/(IU*M0)
      !   INTEGH = INTEGH+2*REAL(FAB*SUM(TFM%W*F%E(:NR,MM,1)))
      ! ENDDO

      ! INTEGH = INTEGH*ZLEN*ELL2

      INTEGH = 0.D0
      IF (F%INX.EQ.0) THEN
        !INTEGH = (B-A)*SUM(TFM%W*F%E(:NR,1,1))
        DO MM=1,SIZE(F%E,2) !NTCHOP
          M0=M(MM+F%INTH) ! M(1) = 0
          IF (M0.EQ.0) THEN
            INTEGH = INTEGH + (B-A)*SUM(TFM%W*F%E(:NR,1,1))
          ELSE
            FAB = (EXP(IU*M0*B)-EXP(IU*M0*A))/(IU*M0)
            INTEGH = INTEGH + 2*REAL(FAB*SUM(TFM%W*F%E(:NR,MM,1)))
          ENDIF
        ENDDO
      ENDIF

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, INTEGH, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_IVP, IERR)
      INTEGH = INTEGH*ZLEN*ELL2

      RETURN
      END FUNCTION INTEGH
!=======================================================================

END MODULE MOD_SCALAR3