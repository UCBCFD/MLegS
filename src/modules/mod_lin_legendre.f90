MODULE MOD_LIN_LEGENDRE ! LEVEL 1 MODULE
! LEGENDRE OPERATORS (MATRIX OR BANDMATRIX) & LEGENDRE UTILITIES
! ALL MATRIX OPERATORS ARE FOR A SINGLE M. INDEPDENT OF K.
  USE omp_lib
  USE MPI
  USE MOD_MISC, ONLY: P4, P8, PI, IU                                   ! LEVEL 0
  IMPLICIT NONE
  PRIVATE
!=======================================================================
!============================ PARAMETERS ===============================
!=======================================================================
  REAL(P8), DIMENSION(1:129) :: FACVAL = 0.D0                        ! STORE FACT1 VALUES ONCE EXECUTED
!=======================================================================
!======================== PUBLIC DECLARATION ===========================
!=======================================================================
  ! FIND ZEROS AND WEIGHTS
  PUBLIC :: LEG_ZERO
  ! FIND NORMALIZATION FACTORS OF LEGENDRE POLYNOMIALS P^M_N(MU)
  PUBLIC :: LEG_NORM, LEG_LOG_NORM
  ! FIND VALUE TABLE OF ASSOCIATED LEGENDRE POLYNOMIALS P^M_N(MU)
  PUBLIC :: LEG_TBL
  ! BASIC OPERATORS: (1-X)* <==> 2R^2/(R^2+L^2)*
  PUBLIC :: LEG_XM
  PUBLIC :: BAND_LEG_XM
  PUBLIC :: LOGLEG_XM
  PUBLIC :: BAND_LOGLEG_XM
  ! BASIC OPERATORS: (1+X)* <==> 2L^2/(R^2+L^2)*
  PUBLIC :: LEG_XP
  PUBLIC :: BAND_LEG_XP
  PUBLIC :: LOGLEG_XP
  PUBLIC :: BAND_LOGLEG_XP
  ! BASIC OPERATORS: X*     <==> (R^2-L^2)/(R^2+L^2)*
  PUBLIC :: LEG_X
  PUBLIC :: BAND_LEG_X
  PUBLIC :: LOGLEG_X
  PUBLIC :: BAND_LOGLEG_X
  ! BASIC OPERATORS: (1-X**2)*D/DX <==> R*D/DR
  PUBLIC :: LEG_XXDX
  PUBLIC :: BAND_LEG_XXDX
  PUBLIC :: LOGLEG_XXDX
  PUBLIC :: BAND_LOGLEG_XXDX
  ! BASIC OPERATORS: DEL^2_PERP (1/R*D/DR(R*D/DR) - M^2/R^2)
  PUBLIC :: LEG_RAT_DEL2H
  PUBLIC :: BAND_LEG_RAT_DEL2H
  PUBLIC :: LOGLEG_RAT_DEL2H
  PUBLIC :: BAND_LOGLEG_RAT_DEL2H
  ! BASIC OPERATORS: DEL^2 (1/R*D/DR(R*D/DR) - M^2/R^2 - AK^2)
  PUBLIC :: LEG_RAT_DEL2
  PUBLIC :: BAND_LEG_RAT_DEL2
  PUBLIC :: LOGLEG_RAT_DEL2
  PUBLIC :: BAND_LOGLEG_RAT_DEL2
!=======================================================================
!============================ INTERFACES ===============================
!=======================================================================
  INTERFACE LEG_NORM
    MODULE PROCEDURE LEG_NORM1                                          ! NORM OF P^M_N FOR ONE M
    MODULE PROCEDURE LEG_NORMT                                          ! NORM OF P^M_N FOR MULTIPLE M'S
  END INTERFACE

  INTERFACE LEG_LOG_NORM
    MODULE PROCEDURE LEG_LOG_NORM1                                      ! LOG(NORM) OF P^M_N FOR ONE M
    MODULE PROCEDURE LEG_LOG_NORMT                                      ! LOG(NORM) OF P^M_N FOR MULTIPLE M'S
  END INTERFACE

  INTERFACE LEG_TBL
    MODULE PROCEDURE LEG_TBL1                                           ! FOR FIXED M, COMPUTE P^M_N VAL @ POINT X
    MODULE PROCEDURE LEG_TBL2                                           ! FOR FIXED M, COMPUTE P^M_N VAL AT X(:)
    MODULE PROCEDURE LEG_TBL3                                           ! MAKE FULL TABLE FOR MULTIPLE M'S
  END INTERFACE

  INTERFACE FACT
    MODULE PROCEDURE FACT1                                              ! GET (2M-1)!! FOR A FIXED M
    MODULE PROCEDURE FACTA                                              ! CALCULATE ALL (2M-1)!! VALUES
  END INTERFACE
CONTAINS
!=======================================================================
!============================ SUBROUTINES ==============================
!=======================================================================
  SUBROUTINE GAUSS_LEGENDRE(X1, X2, X, W, N)
!=======================================================================
! [USAGE]:
! GET THE ABSCISSAS AND WEIGHTS OF GAUSS-LEGENDRE N-POINT QUADRATURES
! FROM NUMERICAL RECIPES IN FORTRAN90 (2ND ED.) (Press et al., 1996)
! [PARAMETERS]:
! X1,X2 >> LOWER AND UPPER BOUNDS
! X >> ABSCISSAS (OR COLLOCATION POINTS) FOR P^M_N
! W >> WEIGHTS FOR P^M_N FOR EACH ZERO (* P^M_N: GENERAL LEGENDRE POLY.)
! N >> # OF QUADRATURES TO BE CALCULATED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    USE MOD_MISC, ONLY: PI
    IMPLICIT NONE
    REAL(P8), INTENT(IN)   :: X1, X2                                    ! LOWER & UPPER BOUNDS
    REAL(P8), DIMENSION(:) :: X, W                                      ! VALUE & WEIGHT (1:N)
    INTEGER(P4)             :: N                                        ! TOTAL NUMBER OF QUATRATURES

    REAL(P8), PARAMETER    :: EPS = 1.D-15                              ! ERROR LIMIT
    REAL(P8)                :: XM, XL, P1, P2, P3, Z, Z1, PP
    INTEGER(P4)             :: M, I, J

    M = CEILING((N+1)/2.D0)
    XM = (X2+X1)/2.D0
    XL = (X2-X1)/2.D0

    Z1 = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Z,P1,P2,P3,PP,Z1)
    DO I = 1, M                                                         ! REFINEMENT USING NEWTON'S METHOD
      Z = COS(PI*(I-.25D0)/(N+.5D0))
      DO WHILE (ABS(Z-Z1) .GT. EPS)
        P1 = 1.D0
        P2 = 0.D0
        DO J = 1, N
          P3 = P2
          P2 = P1
          P1 = ((2*J-1)*Z*P2-(J-1)*P3)/J
        END DO
        PP = N*(Z*P1-P2)/(Z*Z-1)
        Z1 = Z
        Z = Z1-P1/PP
      END DO
      X(I) = XM-XL*Z
      X(N+1-I) = XM+XL*Z
      W(I) = 2*XL/((1-Z*Z)*PP*PP)
      W(N+1-I) = W(I)
    END DO
!$OMP END PARALLEL DO

    RETURN
  END SUBROUTINE GAUSS_LEGENDRE
!=======================================================================
  SUBROUTINE LEG_ZERO(X, R, ELL, W)
!=======================================================================
! [USAGE]:
! CALCULATE THE COLLOCATION POINTS AND WEIGHTS AFTER RATIONAL MAPPING
! OF THE LEGENDRE POLYNOMIALS (P^M_L_N(R) WHERE 0 <= R < +INFTY)
! (* P^M_L_N: RAT. LEGENDRE FUNCTION = P^M_N([R^2-ELL^2]/[R^2+ELL^2]))
! [PARAMETERS]:
! X >> ABSCISSAS (OR COLLOCATION POINTS) FOR P^M_N
! R >> COLLOCATION POINTS FOR RATIONAL MAPPING
! ELL >> THE MAP PARAMETER OPTIMIZING THE CONVERGENCE OF THE EXPANSION
! W >> WEIGHTS FOR P^M_N(X) = P^M_L_N(R) FOR EACH ZERO
!      (* P^M_N: GENERAL LEGENDRE POLYNOMIAL)
! [DEPENDENCIES]:
! 1. GAUSS_LEGENDRE(X1,X2,X,W,N) SUBROUTINE @ MOD_LIN_LEGENDRE
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    REAL(P8), DIMENSION(:) :: X                                         ! ABSCISSAS OF LEGENDRE POLYNOMIALS
    REAL(P8), DIMENSION(:), OPTIONAL :: R                               ! COLLOCATION POINTS FOR RATIONAL MAPPING
    REAL(P8), OPTIONAL :: ELL                                           ! MAPPING PARAMETER
    REAL(P8), DIMENSION(:), OPTIONAL :: W                               ! WEIGHTS FOR QUADRATURE

    REAL(P8), DIMENSION(SIZE(X)):: WK

    IF (PRESENT(R) .NEQV. PRESENT(ELL)) THEN                            ! BOTH PARAM. VARS. MUST BE PRESENT
      WRITE(*,*) 'ERR_LEG_ZEROS: BOTH RATIONAL MAPPING COLLOC. PTS.'&
                  'ARRAY & PARAMETER VAR. MUST BE DECLARED IN ADVANCE'&
                  'TO OBTAIN THE MAPPED COLLOCATION POINTS'
      STOP 'INSUFFICIENT ASSIGNMENT'
    END IF

    IF (PRESENT(W)) THEN
      CALL GAUSS_LEGENDRE(-1.D0, 1.D0, X, W, SIZE(X))
    ELSE
      CALL GAUSS_LEGENDRE(-1.D0, 1.D0, X, WK, SIZE(X))
    END IF

    R = ELL*SQRT((1+X)/(1-X))                                           ! EQUATION (2) (MATSUSHIMA & MARCUS, 1997)

    RETURN
  END SUBROUTINE LEG_ZERO
!=======================================================================
!============================ FUNCTIONS ================================
!=======================================================================
  FUNCTION LEG_NORM1(NDIM, M) ! FAMILY OF LEG_NORM
!=======================================================================
! [USAGE]:
! FIND THE NORMALIZATION FACTOR N^M_N FOR ONE M
! EX) F = LEG_NORM(64, M) <-- SCALAR INPUT M
! [INPUTS]:
! NDIM >> (MAXIMUM) DEGREE OF THE ASSOCIATED LEGENDRE POLYNOMIAL
! M >> ORDER OF THE ASSOCIATED LEGENDRE POLYNOMIAL
! [OUTPUTS]:
! LEG_NORM1 >> NORMALIZATION FACTOR TABLE, SIZE: NDIM X 1
! [DEPENDENCIES]:
! 1. LEG_NORMT(NDIM,MS) FUNCTION @ MOD_LIN_LEGENDRE
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NDIM                                     ! NUMBER OF COEFFICIENTS
    INTEGER(P4), INTENT(IN) :: M                                        ! SINGLE DEGREE
    REAL(P8), DIMENSION(1:NDIM) :: LEG_NORM1                            ! RETURNED ARRAY

    LEG_NORM1 = (/LEG_NORMT(NDIM, (/ABS(M)/))/)

    RETURN
  END FUNCTION LEG_NORM1
!=======================================================================
  FUNCTION LEG_NORMT(NDIM, MS) ! FAMILY OF LEG_NORM
!=======================================================================
! [USAGE]:
! FIND THE NORMALIZATION FACTORS N^M_N FOR MULTIPLE M'S IN MS
! EX) F = LEG_NORM(64, (/0, MS/)) <-- ARRAY INPUT (/0, M/) = MS
! [NOTES]:
! 1. CHECK THE FOLLOWING RECURSIVE FORMULAS
! P^L_L  (X) = (-1)**L*(2*L-1)!!*(1-X**2)**(L/2)
! P^L_L+1(X) = X*(2*L+1)*P^L_L(X)
! P^K_L(X)   = 1/(L-K)*[(2*L-1)*X*P^K_L-1(X) - (L+K-1)*P^K_L-2(X)]
! TO OBTAIN THE NORMALIZATION FACTORS AS FOLLOWS:
! N^M_N = SQRT(K(2N+1)*(N-|M|!)/(N+|M|)!) WHERE K = 1 FOR M=0; 2 OTHER.
! 2. IN THIS FUNCTION RECIPROCAL IS OUTPUT OF N^M_N.
! GIVEN N = NN-1 & M = MS(MM), LEG_NORMT(NN,MM) * P^M_N(R)
! IS AN NORMALIZED FUNCTION FOR THE ORDER M AND DEGREE N!
! [INPUTS]:
! NDIM >> THE (MAXIMUM) DEGREE OF THE ASSOCIATED LEGENDRE POLYNOMIAL
! MS >> "VECTOR" OF THE ORDER OF THE ASSOCIATED LEGENDRE POLYNOMIAL
! [OUTPUTS]:
! LEG_NORMT >> NORMALIZATION FACTOR TABLE, SIZE: NDIM X SIZE(MS)
!              FOR EACH M, CALCULATES
!              (/N^M_|M| , N^M_|M|+1 , ... , N^M_|M|+NDIM-1/)
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NDIM                                     ! NUMBER OF COEFFICIENTS
    INTEGER(P4), DIMENSION(:), INTENT(IN) :: MS                         ! MULTIPLE DEGREES
    REAL(P8), DIMENSION(1:NDIM, 1:SIZE(MS)) :: LEG_NORMT                ! RETURNED ARRAY

    INTEGER(P4) :: NN, MM, N, M, ME
    REAL(P8), DIMENSION(:), ALLOCATABLE :: WK

    ME = MAXVAL(ABS(MS))
    ALLOCATE (WK(0:ME))

    WK(0) = 0.5D0

    DO M = 1, ME
      WK(M) = WK(M-1)*(2*M+1.D0)/(2*M*(2*M-1.D0)**2)                    ! EQUATION (13) (MATSUSHIMA & MARCUS, 1997) WITHOUT SQRT (VALID UNTIL M ~ 85)
    END DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(M)
    DO MM = 1, SIZE(MS)
      M = ABS(MS(MM))
      LEG_NORMT(1, MM) = WK(M)*((10.D0**2)**M)                          ! 10**2M TO AVOID UNDERFLOW
    END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(M,N)
    DO MM = 1, SIZE(MS)
      M = ABS(MS(MM))
      DO NN = 2, NDIM
        N = M+(NN-1)
        LEG_NORMT(NN, MM) = LEG_NORMT(NN-1, MM)*(2*N+1.D0)/(2*N-1.D0) &
                            *(N-M)/(N+M)
      END DO
    END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(M)
    DO MM = 1, SIZE(MS)
      M = ABS(MS(MM))
      LEG_NORMT(:, MM) = SQRT(LEG_NORMT(:, MM))*1.D1**(-M)              ! TAKE SQRT. 10**-M TO RECOVER UNDERFLOW PROTECTION
    END DO
!$OMP END PARALLEL DO

    DEALLOCATE (WK)

    RETURN
  END FUNCTION LEG_NORMT
!=======================================================================
  FUNCTION LEG_LOG_NORM1(NDIM, M) ! FAMILY OF LEG_LOG_NORM
!=======================================================================
! [USAGE]:
! FIND THE LOG OF THE NORMALIZATION FACTOR N^M_N FOR ONE M
! EX) F = LEG_LOG_NORM(64, M) <-- SCALAR INPUT M
! [INPUTS]:
! NDIM >> (MAXIMUM) DEGREE OF THE ASSOCIATED LEGENDRE POLYNOMIAL
! M >> ORDER OF THE ASSOCIATED LEGENDRE POLYNOMIAL
! [OUTPUTS]:
! LEG_LOG_NORM1 >> LOG NORMALIZATION FACTOR TABLE, SIZE: NDIM X 1
! [DEPENDENCIES]:
! 1. LEG_LOG_NORMT(NDIM,MS) FUNCTION @ MOD_LIN_LEGENDRE
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ MAR 14 2021
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NDIM                                     ! NUMBER OF COEFFICIENTS
    INTEGER(P4), INTENT(IN) :: M                                        ! SINGLE DEGREE
    REAL(P8), DIMENSION(1:NDIM) :: LEG_LOG_NORM1                        ! RETURNED ARRAY

    LEG_LOG_NORM1 = (/LEG_LOG_NORMT(NDIM, (/ABS(M)/))/)

    RETURN
  END FUNCTION LEG_LOG_NORM1
!=======================================================================
  FUNCTION LEG_LOG_NORMT(NDIM, MS) ! FAMILY OF LEG_LOG_NORM
!=======================================================================
! [USAGE]:
! FIND THE LOG OF THE NORMALIZATION FACTORS N^M_N FOR MULTIPLE M'S IN MS
! EX) F = LEG_LOG_NORM(64, (/0, MS/)) <-- ARRAY INPUT (/0, M/) = MS
! [NOTES]:
! 1. CHECK THE FOLLOWING RECURSIVE FORMULAS
! P^L_L  (X) = (-1)**L*(2*L-1)!!*(1-X**2)**(L/2)
! P^L_L+1(X) = X*(2*L+1)*P^L_L(X)
! P^K_L(X)   = 1/(L-K)*[(2*L-1)*X*P^K_L-1(X) - (L+K-1)*P^K_L-2(X)]
! TO OBTAIN THE NORMALIZATION FACTORS AS FOLLOWS:
! N^M_N = SQRT(K(2N+1)*(N-|M|!)/(N+|M|)!) WHERE K = 1 FOR M=0; 2 OTHER.
! 2. IN THIS FUNCTION RECIPROCAL IS OUTPUT OF N^M_N.
! GIVEN N = NN-1 & M = MS(MM), LEG_NORMT(NN,MM) * P^M_N(R)
! IS AN NORMALIZED FUNCTION FOR THE ORDER M AND DEGREE N!
! [INPUTS]:
! NDIM >> THE (MAXIMUM) DEGREE OF THE ASSOCIATED LEGENDRE POLYNOMIAL
! MS >> "VECTOR" OF THE ORDER OF THE ASSOCIATED LEGENDRE POLYNOMIAL
! [OUTPUTS]:
! LEG_LOG_NORMT >> LOG NORMALIZATION FACTOR TABLE, SIZE: NDIM X SIZE(MS)
!              FOR EACH M, CALCULATES
!              (/LOG(N^M_|M|) , ... , LOG(N^M_|M|+NDIM-1)/)
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ MAR 14, 2021
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NDIM                                     ! NUMBER OF COEFFICIENTS
    INTEGER(P4), DIMENSION(:), INTENT(IN) :: MS                         ! MULTIPLE DEGREES
    REAL(P8), DIMENSION(1:NDIM, 1:SIZE(MS)) :: LEG_LOG_NORMT            ! RETURNED ARRAY

    INTEGER(P4) :: NN, MM, N, M, ME
    REAL(P8), DIMENSION(:), ALLOCATABLE :: WK_LOG

    ME = MAXVAL(ABS(MS))
    ALLOCATE (WK_LOG(0:ME))

    WK_LOG(0) = LOG(0.5D0)

    DO M = 1, ME
      WK_LOG(M) = WK_LOG(M-1)+LOG(2*M+1.D0)-LOG(2*M*(2*M-1.D0)**2)      ! EQUATION (13) (MATSUSHIMA & MARCUS, 1997) WITHOUT SQRT
    END DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(M)
    DO MM = 1, SIZE(MS)
      M = ABS(MS(MM))
      LEG_LOG_NORMT(1, MM) = WK_LOG(M)
    END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(M,N)
    DO MM = 1, SIZE(MS)
      M = ABS(MS(MM))
      DO NN = 2, NDIM
        N = M+(NN-1)
        LEG_LOG_NORMT(NN, MM) = LEG_LOG_NORMT(NN-1, MM) &
                                +LOG(2*N+1.D0)-LOG(2*N-1.D0) &
                                +LOG(1.D0*(N-M))-LOG(1.D0*(N+M))
      END DO
    END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(M)
    DO MM = 1, SIZE(MS)
      LEG_LOG_NORMT(:, MM) = .5D0*LEG_LOG_NORMT(:, MM)
    END DO
!$OMP END PARALLEL DO

    DEALLOCATE (WK_LOG)

    RETURN
  END FUNCTION LEG_LOG_NORMT
! ======================================================================
  FUNCTION LEG_TBL1(X, NE, M, LOGNORM) ! FAMILY OF LEG_TBL
!=======================================================================
! [USAGE]:
! TABULATE THE LEGENDRE POLYNOMIALS AT ONE POINT X FOR A FIXED M
! [INPUTS]:
! X >> ONE RADIAL POINT TO EVALUATE
! NE >> MAX. DEGREE OF THE ASSOCIATED LEGENDRE POLYNOMIAL
! M >> ORDER OF THE ASSOCIATED LEGENDRE POLYNOMIAL TO EVALUATE
! LOGNORM >> (OPTIONAL) TABLE OF LN(NORMALIZATION FACTORS) FROM LEG_NORM
! [OUTPUTS]:
! LEG_TBL1 >> THE LEGENDRE POLY. COEFFCIENT AT X & M FROM N = 1 TO NE
!             (/P^M_L_|M|(X), P^M_L_|M|+1(X) , ..., P^M_L_|M|+NE-1(X)/)
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE WHEN LOGNORM IS INPUT
! 2. LEG_TBL3(X,NE,MS,LOGNORM) FUNCTION @ MOD_LIN_LEGENDRE
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    REAL(P8), INTENT(IN)    :: X                                        ! POINT TO EVALUATE
    INTEGER(P4), INTENT(IN) :: NE                                       ! MAXIMUM DEGREE
    INTEGER(P4), INTENT(IN) :: M                                        ! SINGLE ORDER
    REAL(P8), DIMENSION(:), INTENT(IN), OPTIONAL :: LOGNORM             ! TABLE OF NORMALIZATION FACTORS FROM LEG_NORM
    REAL(P8), DIMENSION(1:NE) :: LEG_TBL1                               ! RETURNED LEGENDRE POLYNOMIAL TABLE

    REAL(P8), DIMENSION(1:NE, 1) :: LOGNORM1                            ! ARRAY TO STORE NORMALIZATION FACTORS FOR M

    IF (PRESENT(LOGNORM)) THEN
      LOGNORM1(:, 1) = LOGNORM(:NE)
      LEG_TBL1 = (/LEG_TBL3((/X/), NE, (/ABS(M)/), LOGNORM1)/)
    ELSE
      LEG_TBL1 = (/LEG_TBL3((/X/), NE, (/ABS(M)/))/)
    END IF

    RETURN
  END FUNCTION LEG_TBL1
!=======================================================================
  FUNCTION LEG_TBL2(X, NE, M, LOGNORM) ! FAMILY OF LEG_TBL
!=======================================================================
! [USAGE]:
! TABULATE THE LEGENDRE POLYNOMIALS FOR A FIXED M AT MULTIPLE X'S
! [INPUTS]:
! X >> MULTIPLE POINTS (EX. COLLOC. PTS.) TO EVALUATE
! NE >> MAX. DEGREE OF THE ASSOCIATED LEGENDRE POLYNOMIAL
! M >> ORDER OF THE ASSOCIATED LEGENDRE POLYNOMIAL TO EVALUATE
! LOGNORM >> (OPTIONAL) TABLE OF LN(NORMALIZATION FACTORS) FROM LEG_NORM
! [OUTPUTS]:
! LEG_TBL2 >> THE LEGENDRE POLY. COEFFCIENT AT X & M FROM N = 1 TO NE
!             FOR EACH X,
!             (/P^M_L_|M|(X), P^M_L_|M|+1(X) , ..., P^M_L_|M|+NE-1(X)/)
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE WHEN LOGNORM IS INPUT
! 2. LEG_TBL3(X,NE,MS,LOGNORM) FUNCTION @ MOD_LIN_LEGENDRE
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    REAL(P8), DIMENSION(:), INTENT(IN) :: X                             ! POINT TO EVALUATE
    INTEGER(P4), INTENT(IN) :: NE                                       ! MAXIMUM DEGREE
    INTEGER(P4), INTENT(IN) :: M                                        ! SINGLE ORDER
    REAL(P8), DIMENSION(:), INTENT(IN), OPTIONAL :: LOGNORM             ! TABLE OF NORMALIZATION FACTORS FROM LEG_NORM
    REAL(P8), DIMENSION(SIZE(X), NE) :: LEG_TBL2                        ! RETURNED LEGENDRE POLYNOMIAL TABLE

    IF (PRESENT(LOGNORM)) THEN
      LEG_TBL2 = RESHAPE(LEG_TBL3(X, NE, (/ABS(M)/), &
                                  RESHAPE(LOGNORM, (/SIZE(LOGNORM), 1/))), (/SIZE(X), NE/))
    ELSE
      LEG_TBL2 = RESHAPE(LEG_TBL3(X, NE, (/ABS(M)/)), (/SIZE(X), NE/))
    END IF

    RETURN
  END FUNCTION LEG_TBL2
!=======================================================================
!       FUNCTION LEG_TBL3(X,NE,MS,NORM) ! FAMILY OF LEG_TBL
! !=======================================================================
! ! [USAGE]:
! ! MAKE A FULL TABLE OF THE ASSOCIATED LEGENDRE POLYNOMIALS
! ! LEG_TBL3(RadialCollocPts,RadialModes,AzimuthalModes)
! ! [INPUTS]:
! ! X >> MULTIPLE POINTS (EX. COLLOC. PTS.) TO EVALUATE
! ! NE >> MAX. DEGREE OF THE ASSOCIATED LEGENDRE POLYNOMIAL
! ! MS >> SET OF ORDERS OF THE ASSOCIATED LEGENDRE POLYNOMIAL TO EVALUATE
! ! NORM >> (OPTIONAL) TABLE OF THE NORMALIZATION FACTORS FROM LEG_NORM
! ! [OUTPUTS]:
! ! LEG_TBL3 >> THE LEGENDRE POLY. COEFFCIENT AT X & M FROM N = 1 TO NE
! !             FOR EACH X AND M,
! !             (/P^M_L_|M|(X), P^M_L_|M|+1(X) , ..., P^M_L_|M|+NE-1(X)/)
! ! [DEPENDENCIES]:
! ! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! ! 2. FACT1(M) FUNCTION @ MOD_LIN_LEGENDRE
! ! [UPDATES]:
! ! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! ! RE-CODED BY SANGJOON LEE @ OCT 26 2020
! !=======================================================================
!       IMPLICIT NONE
!       REAL(P8) , DIMENSION(:), INTENT(IN):: X                            ! POINT TO EVALUATE
!       INTEGER(P4) , INTENT(IN) :: NE                                     ! MAXIMUM DEGREE
!       INTEGER(P4) , INTENT(IN), DIMENSION(:) :: MS                       ! SET OF ORDERS
!       REAL(P8) , DIMENSION(:,:), INTENT(IN), OPTIONAL :: NORM            ! TABLE OF NORMALIZATION FACTORS FROM LEG_NORM
!       REAL(P8) , DIMENSION(SIZE(X),NE,SIZE(MS)) :: LEG_TBL3              ! RETURNED LEGENDRE POLYNOMIAL TABLE

!       INTEGER(P4) :: M, ME, MM, NN, N

!       ME = SIZE(MS)

!       DO MM=1,ME                                                         ! ORTHOGONAL-UNNORMALIZED POLYNOMIAL TABLE
!         M = ABS(MS(MM))                                                  ! STARTING VALUES OF LEG_TBL3
!         IF(M.EQ.0) THEN
!           LEG_TBL3(:,1,MM) = 1.D0
!         ELSE
!           LEG_TBL3(:,1,MM) = ((-1.)**M)*FACT1(M)*(1-X**2)**(0.5D0*M)     ! P^L_L(X) = (-1)**L*(2*L-1)!!*(1-X**2)**(L/2)
!         ENDIF
!       ENDDO

!       DO MM=1,ME
!         M = ABS(MS(MM))
!         LEG_TBL3(:,2,MM) = ((2*M+1)*X)*LEG_TBL3(:,1,MM)                  ! P^L_L+1(X) = X*(2*L+1)*P^L_L(X)
!       ENDDO

!       DO MM=1,ME                                                         ! USE OF THE RECURRENCE RELATIONSHIP
!         M = ABS(MS(MM))
!         DO NN=3,NE
!           N = M+(NN-1)
!           LEG_TBL3(:,NN,MM) = 1.D0/(N-M)*(   &                           ! P^K_L(X) = 1/(L-K)*[(2*L-1)*X*P^K_L-1(X) - (L+K-1)*P^K_L-2(X)]
!           LEG_TBL3(:,NN-1,MM)*((2*N-1)*X)  &
!           -LEG_TBL3(:,NN-2,MM)*(N+M-1))
!         ENDDO
!       ENDDO

!       IF(PRESENT(NORM)) THEN                                             ! ORTHONORMAL-NORMALIZED POLYNOMIAL TABLE
!         IF(SIZE(NORM,1).LT.NE .OR. SIZE(NORM,2).LT.ME) THEN
!           WRITE(*,*) 'NORM SIZE IS NOT COMPATIBLE. IT IS TOO SMALL.'
!           WRITE(*,*) 'NORM=',SIZE(NORM,1),'X',SIZE(NORM,2)
!           STOP
!         ENDIF
!         DO MM=1,ME
!           DO NN=1,NE
!             LEG_TBL3(:,NN,MM)=LEG_TBL3(:,NN,MM)*NORM(NN,MM)
!           ENDDO
!         ENDDO
!       ENDIF

!       RETURN
!       END FUNCTION LEG_TBL3
!=======================================================================
  FUNCTION LEG_TBL3(X, NE, MS, LOGNORM) ! FAMILY OF LEG_TBL
!=======================================================================
! [USAGE]:
! MAKE A FULL TABLE OF THE ASSOCIATED LEGENDRE POLYNOMIALS
! [INPUTS]:
! X >> MULTIPLE POINTS (EX. COLLOC. PTS.) TO EVALUATE
! NE >> MAX. DEGREE OF THE ASSOCIATED LEGENDRE POLYNOMIAL
! MS >> SET OF ORDERS OF THE ASSOCIATED LEGENDRE POLYNOMIAL TO EVALUATE
! LOGNORM >> (OPTIONAL) TABLE OF LN(NORMALIZATION FACTORS) FROM LEG_NORM
! [OUTPUTS]:
! LEG_TBL3 >> THE LEGENDRE POLY. COEFFCIENT AT X & M FROM N = 1 TO NE
!             FOR EACH X AND M,
!             (/P^M_L_|M|(X), P^M_L_|M|+1(X) , ..., P^M_L_|M|+NE-1(X)/)
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE WHEN LOGNORM IS INPUT
! 2. FACT1(M) FUNCTION @ MOD_LIN_LEGENDRE
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    USE FMVALS
    USE FMZM                                                              ! MULTIPLE PRECISION (FOR ENABLING HIGHER M'S)
    IMPLICIT NONE
    REAL(P8), DIMENSION(:), INTENT(IN):: X                                ! POINT TO EVALUATE
    INTEGER(P4), INTENT(IN) :: NE                                         ! MAXIMUM DEGREE
    INTEGER(P4), INTENT(IN), DIMENSION(:) :: MS                           ! SET OF ORDERS
    REAL(P8), DIMENSION(:, :), INTENT(IN), OPTIONAL :: LOGNORM            ! TABLE OF NORMALIZATION FACTORS FROM LEG_NORM
    REAL(P8), DIMENSION(SIZE(X), NE, SIZE(MS)) :: LEG_TBL3                ! RETURNED LEGENDRE POLYNOMIAL TABLE

    INTEGER(P4) :: M, ME, MM, NN, N, XX
    TYPE(FM), ALLOCATABLE, DIMENSION(:, :) :: FM_LEG_TBL3                 ! MULTIPLE PRECISION ARRAY

    CALL FM_SET(50)
    ME = SIZE(MS)

    DO MM = 1, ME
      ALLOCATE (FM_LEG_TBL3(SIZE(X), NE))

      M = ABS(MS(MM))                                                     ! STARTING VALUES OF LEG_TBL3
      IF (M .EQ. 0) THEN
        FM_LEG_TBL3(:, 1) = TO_FM(1.D0)                                   ! HERE LEG_TBL3 STORES P^L_L(X)/(2*L-1)!!
      ELSE
        DO XX = 1, SIZE(X)
          FM_LEG_TBL3(XX, 1) = (TO_FM(-1.D0)**M)*SQRT(TO_FM(1.D0) &
                                                      -TO_FM(X(XX))**2)**M  ! P^L_L(X)/(2*L-1)!! = (-1)**L*(1-X**2)**(L/2)
        END DO
      END IF

      DO XX = 1, SIZE(X)                                                  ! P^L_L+1(X) = X*(2*L+1)*P^L_L(X)  <===>
        FM_LEG_TBL3(XX, 2) = X(XX)*FM_LEG_TBL3(XX, 1)                     ! P^L_L+1(X) / (2*L+1)!! = X * P^L_L(X) / (2L-1)!!
      END DO

      DO NN = 3, NE
        N = M+(NN-1)
        DO XX = 1, SIZE(X)                                                ! USE OF THE RECURRENCE RELATIONSHIP
          FM_LEG_TBL3(XX, NN) = TO_FM(1.D0)/(N-M)*( &                     ! P^K_L(X) = 1/(L-K)*[(2*L-1)*X*P^K_L-1(X) - (L+K-1)*P^K_L-2(X)] <===>
                                FM_LEG_TBL3(XX, NN-1)*X(XX) &             ! P^K_L(X) / (2*L-1)!! = 1/(L-K)*[X*P^K_L-1(X)/(2*L-3)!! - (L+K-1)/(2*L-1)/(2*L-3)*P^K_L-2(X)/(2*L-5)!!]
                                -FM_LEG_TBL3(XX, NN-2)*(N+M-1)/(2*N-1)/(2*N-3))
        END DO
      END DO

      IF (PRESENT(LOGNORM)) THEN                                          ! ORTHONORMAL-NORMALIZED POLYNOMIAL TABLE
        IF (SIZE(LOGNORM, 1) .LT. NE .OR. SIZE(LOGNORM, 2) .LT. ME) THEN
          WRITE (*, *) 'LOGNORM SIZE IS NOT COMPATIBLE. IT IS TOO SMALL.'
          WRITE (*, *) 'LOGNORM=', SIZE(LOGNORM, 1), 'X', SIZE(LOGNORM, 2)
          STOP
        END IF

        DO NN = 1, NE
          FM_LEG_TBL3(:, NN) = FM_LEG_TBL3(:, NN) &
                               *EXP(TO_FM(LOGNORM(NN, MM)) &
                                    +TO_FM(LOGFACT(M+NN-1)))
        END DO
      ELSE
        DO NN = 1, NE
          DO XX = 1, SIZE(X)
            FM_LEG_TBL3(XX, NN) = FM_LEG_TBL3(XX, NN) &
                                  *EXP(TO_FM(LOGFACT(M+NN-1)))            ! MULTIPLYING (2M-1)!!
          END DO
        END DO
      END IF

      DO NN = 1, NE
        DO XX = 1, SIZE(X)
          LEG_TBL3(XX, NN, MM) = TO_DP(FM_LEG_TBL3(XX, NN))
        END DO
      END DO

      CALL FM_DEALLOCATE(FM_LEG_TBL3)
      DEALLOCATE (FM_LEG_TBL3)

    END DO

    RETURN
  END FUNCTION LEG_TBL3
! ======================================================================
  FUNCTION FACT1(M) ! FAMILY OF FACT
!=======================================================================
! [USAGE]:
! COMPUTE (2M)! / (2^M * M!) OR EQUIVALENTLY (2M-1)!! (SEMI-FACTORIAL)
! [INPUTS]:
! M >> FACTORIAL ARGUMENT
! [OUTPUTS]:
! FACT1 >> (2M)! / (2^M * M!) OR EQUIVALENTLY (2M-1)!!
!          TO AVOID OVERFLOW, M SET TO BE LESS THAN 129
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
! REMOVAL OF "SAVE" ARGUMENT FOR STABLE USE, S.LEE, OCT 27 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN)  :: M
    REAL(P8) :: FACT1

    INTEGER(P4) :: I

    IF (M .GT. SIZE(FACVAL)) THEN
      WRITE (*, *) 'FACT1: ARGUMENT OUT OF RANGE.'
      STOP
    END IF

    IF (FACVAL(1) .NE. 1.D0) THEN                                       ! COMPUTE ALL FACTORIALS AT FIRST TIME (INITIALLY ALL ARE 0.D0)
      FACVAL = FACT()
      FACT1 = FACVAL(ABS(M))
    ELSE
      FACT1 = FACVAL(ABS(M))
    END IF

    RETURN
  END FUNCTION FACT1
!=======================================================================
  FUNCTION FACTA() ! FAMILY OF FACT
!=======================================================================
! [USAGE]:
! COMPUTE (2M)! / (2^M * M!) OR EQUIVALENTLY (2M-1)!! FROM M = 1 TO 150
! [INPUTS]:
! N/A
! [OUTPUTS]:
! FACT1 >> (2M)! / (2^M * M!) OR EQUIVALENTLY (2M-1)!! FROM M =1 TO 150
!          TO AVOID OVERFLOW, M SET TO BE LESS THAN 129
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    REAL(P8), DIMENSION(SIZE(FACVAL)) :: FACTA

    INTEGER(P4) :: I

    FACTA(1) = 1.D0
    DO I = 2, SIZE(FACVAL)
      FACTA(I) = FACTA(I-1)*(2*I-1)
    END DO

    RETURN
  END FUNCTION FACTA
!=======================================================================
  FUNCTION LOGFACT(M) ! FAMILY OF FACT
!=======================================================================
! [USAGE]:
! COMPUTE LOG((2M)! / (2^M * M!)) OR EQUIVALENTLY LOG((2M-1)!!)
! [INPUTS]:
! M >> FACTORIAL ARGUMENT
! [OUTPUTS]:
! LOGFACT >> LOG((2M)! / (2^M * M!)) OR EQUIVALENTLY LOG((2M-1)!!)
! [UPDATES]:
! CODED BY SANGJOON LEE @ MAR 14 2021
!=======================================================================
    IMPLICIT NONE
    REAL(P8) :: LOGFACT
    INTEGER :: I, M

    LOGFACT = LOG_GAMMA(2*M+1.D0)-M*LOG(2.D0)-LOG_GAMMA(M+1.D0)

    RETURN
  END FUNCTION LOGFACT
!=======================================================================
  FUNCTION LEG_XM(NI, NJ, M, NORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY (1-X) IN FUNCTION (LEGENDRE) SPACE IN A MATRIX
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X) = SUM_N=|M|^INFTY A^M_N P^M_N(X)
! COEFFICIENTS OF THE FUNCTION (1-X)F(X) B^M_N AS
! (1-X)F(X) = SUM_N=|M|^INFTY B^M_N P^M_N(X)
! CAN BE EVAULATED BY "B^M_N = LEG_XM * A^M_N"
! 2. (MATSUSHIMA & MARCUS, 1997)
! GIVEN THE ORDER M,
! B^M_N = -(N-ABS(M))/(2*N-1)A^M_N-1+A^M_N-(N+ABS(M)+1)/(2*N+3)A^M_N+1
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! NJ >> # OF COLUMNS OF THE OPERATOR
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! NORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! LEG_XM >> OPERATOR TO MULTIPLY (1-X) IN FUNCTION SPACE, SIZE: NI X NJ
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, NJ, M
    REAL(P8), DIMENSION(1:NI, 1:NJ) :: LEG_XM
    REAL(P8), DIMENSION(:), OPTIONAL :: NORM

    INTEGER(P4) :: NN, N, I, J, AM

    AM = ABS(M)

    LEG_XM = 0.D0                                                       ! ASSIGN ZERO MATRIX TO LEG_XM

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) LEG_XM(NN, NN-1) = -(N-AM)/(2*N-1.D0)            ! SUBDIAGONAL COMPONENTS
      IF (NN .LE. NJ) LEG_XM(NN, NN) = 1.D0                             ! MAIN DIAGONAL COMPONENTS
      IF (NN+1 .LE. NJ) LEG_XM(NN, NN+1) = -(N+AM+1.D0)/(2*N+3.D0)      ! SUPERDIAGONAL COMPONENTS
    END DO

    IF (PRESENT(NORM)) THEN
      IF (SIZE(NORM) .LT. MAX(NI, NJ)) THEN
        WRITE (*, *) 'LEG_XM: SIZE(NORM) TOO SMALL'
        WRITE (*, *) 'SIZE(NORM)=', SIZE(NORM)
        WRITE (*, *) 'NI=', NI, '  NJ=', NJ
        STOP
      END IF

      DO I = 1, NI
        DO J = 1, NJ
          LEG_XM(I, J) = LEG_XM(I, J)*(NORM(J)/NORM(I))
        END DO
      END DO
    END IF

    RETURN
  END FUNCTION LEG_XM
!=======================================================================
  FUNCTION LOGLEG_XM(NI, NJ, M, LOGNORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY (1-X) IN FUNCTION (LEGENDRE) SPACE IN A MATRIX
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X) = SUM_N=|M|^INFTY A^M_N P^M_N(X)
! COEFFICIENTS OF THE FUNCTION (1-X)F(X) B^M_N AS
! (1-X)F(X) = SUM_N=|M|^INFTY B^M_N P^M_N(X)
! CAN BE EVAULATED BY "B^M_N = LOGLEG_XM * A^M_N"
! 2. (MATSUSHIMA & MARCUS, 1997)
! GIVEN THE ORDER M,
! B^M_N = -(N-ABS(M))/(2*N-1)A^M_N-1+A^M_N-(N+ABS(M)+1)/(2*N+3)A^M_N+1
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! NJ >> # OF COLUMNS OF THE OPERATOR
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! LOGNORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! LOGLEG_XM >> OPERATOR TO MULTIPLY (1-X) IN FUNCTION SPACE, SIZE: NI X NJ
! [DEPENDENCIES]:
! 1. LOGLEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'LOGNORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, NJ, M
    REAL(P8), DIMENSION(1:NI, 1:NJ) :: LOGLEG_XM
    REAL(P8), DIMENSION(:), OPTIONAL :: LOGNORM

    INTEGER(P4) :: NN, N, I, J, AM

    AM = ABS(M)

    LOGLEG_XM = 0.D0                                                    ! ASSIGN ZERO MATRIX TO LOGLEG_XM

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) LOGLEG_XM(NN, NN-1) = -(N-AM)/(2*N-1.D0)         ! SUBDIAGONAL COMPONENTS
      IF (NN .LE. NJ) LOGLEG_XM(NN, NN) = 1.D0                          ! MAIN DIAGONAL COMPONENTS
      IF (NN+1 .LE. NJ) LOGLEG_XM(NN, NN+1) = -(N+AM+1.D0)/(2*N+3.D0)   ! SUPERDIAGONAL COMPONENTS
    END DO

    IF (PRESENT(LOGNORM)) THEN
      IF (SIZE(LOGNORM) .LT. MAX(NI, NJ)) THEN
        WRITE (*, *) 'LOGLEG_XM: SIZE(LOGNORM) TOO SMALL'
        WRITE (*, *) 'SIZE(LOGNORM)=', SIZE(LOGNORM)
        WRITE (*, *) 'NI=', NI, '  NJ=', NJ
        STOP
      END IF

      DO I = 1, NI
        DO J = 1, NJ
          ! LOGLEG_XM(I,J)=LOGLEG_XM(I,J)*(NORM(J)/NORM(I))
          LOGLEG_XM(I, J) = LOGLEG_XM(I, J)*EXP(LOGNORM(J)-LOGNORM(I))
        END DO
      END DO
    END IF

    RETURN
  END FUNCTION LOGLEG_XM
!=======================================================================
  FUNCTION BAND_LEG_XM(NI, M, NORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY (1-X) IN FUNCTION (LEGENDRE) SPACE
! BASICALLY SAME AS LEG_XM
! HOWEVER, ONLY NON-TRIVIAL COMPONENTS (TRI-DIAGONALS) ARE COMPUTED
! AND SHOWN AS A TRIPLE BAND MATRIX (NI X 3) FORM FOR EFFICIENCY
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! NORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! BAND_LEG_XM >> BAND MATRIX STYLE OPER. TO MULTIPLY (1-X) IN FUNCTION
!               SPACE, SIZE: NIX3 (1: SUBDIAG., 2: DIAG., 3: SUPERDIAG.)
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, M
    REAL(P8), DIMENSION(NI, 3) :: BAND_LEG_XM
    REAL(P8), DIMENSION(:), OPTIONAL:: NORM

    INTEGER :: NN, N, I, J, AM

    AM = ABS(M)
    BAND_LEG_XM = 0.D0

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) BAND_LEG_XM(NN, 1) = -(N-AM)/(2*N-1.D0)          ! SUBDIAGONAL COMPONENTS
      IF (NN .LE. NI) BAND_LEG_XM(NN, 2) = 1.D0                         ! MAIN DIAGONAL COMPONENTS
      IF (NN+1 .LE. NI) BAND_LEG_XM(NN, 3) = -(N+AM+1.D0)/(2*N+3.D0)    ! SUPERDIAGONAL COMPONENTS
    END DO

    IF (PRESENT(NORM)) THEN
      IF (SIZE(NORM) .LT. NI) THEN
        WRITE (*, *) 'BAND_LEG_XM: SIZE(NORM) TOO SMALL'
        WRITE (*, *) 'SIZE(NORM)=', SIZE(NORM)
        WRITE (*, *) 'NI=', NI
        STOP
      END IF

      DO I = 1, NI
        IF (I .GT. 1) THEN
          BAND_LEG_XM(I, 1) = BAND_LEG_XM(I, 1)*(NORM(I-1)/NORM(I))
        END IF
        IF (I .LT. NI) THEN
          BAND_LEG_XM(I, 3) = BAND_LEG_XM(I, 3)*(NORM(I+1)/NORM(I))
        END IF
      END DO
    END IF

    RETURN
  END FUNCTION BAND_LEG_XM
!=======================================================================
  FUNCTION BAND_LOGLEG_XM(NI, M, LOGNORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY (1-X) IN FUNCTION (LEGENDRE) SPACE
! BASICALLY SAME AS LOGLEG_XM
! HOWEVER, ONLY NON-TRIVIAL COMPONENTS (TRI-DIAGONALS) ARE COMPUTED
! AND SHOWN AS A TRIPLE BAND MATRIX (NI X 3) FORM FOR EFFICIENCY
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! LOGNORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! BAND_LOGLEG_XM >> BAND MATRIX STYLE OPER. TO MULTIPLY (1-X) IN FUNCTION
!               SPACE, SIZE: NIX3 (1: SUBDIAG., 2: DIAG., 3: SUPERDIAG.)
! [DEPENDENCIES]:
! 1. LOGLEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'LOGNORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, M
    REAL(P8), DIMENSION(NI, 3) :: BAND_LOGLEG_XM
    REAL(P8), DIMENSION(:), OPTIONAL:: LOGNORM

    INTEGER :: NN, N, I, J, AM

    AM = ABS(M)
    BAND_LOGLEG_XM = 0.D0

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) BAND_LOGLEG_XM(NN, 1) = -(N-AM)/(2*N-1.D0)       ! SUBDIAGONAL COMPONENTS
      IF (NN .LE. NI) BAND_LOGLEG_XM(NN, 2) = 1.D0                      ! MAIN DIAGONAL COMPONENTS
      IF (NN+1 .LE. NI) BAND_LOGLEG_XM(NN, 3) = -(N+AM+1.D0)/(2*N+3.D0) ! SUPERDIAGONAL COMPONENTS
    END DO

    !   IF(PRESENT(NORM)) THEN
    !     IF(SIZE(NORM) .LT. NI) THEN
    !       WRITE(*,*) 'BAND_LOGLEG_XM: SIZE(NORM) TOO SMALL'
    !       WRITE(*,*) 'SIZE(NORM)=',SIZE(NORM)
    !       WRITE(*,*) 'NI=',NI
    !       STOP
    !     ENDIF

    !     DO I = 1, NI
    !       IF(I .GT. 1) THEN
    !         BAND_LOGLEG_XM(I,1)=BAND_LOGLEG_XM(I,1)*(NORM(I-1)/NORM(I))
    !       ENDIF
    !       IF(I .LT. NI) THEN
    !         BAND_LOGLEG_XM(I,3)=BAND_LOGLEG_XM(I,3)*(NORM(I+1)/NORM(I))
    !       ENDIF
    !     ENDDO
    !   ENDIF

    IF (NI .LE. 1) RETURN

    IF (PRESENT(LOGNORM)) THEN
      IF (SIZE(LOGNORM) .LT. NI) THEN
        WRITE (*, *) 'BAND_LOGLEG_XM: SIZE(LOGNORM) TOO SMALL'
        WRITE (*, *) 'SIZE(LOGNORM)=', SIZE(LOGNORM)
        WRITE (*, *) 'NI=', NI
        STOP
      END IF

      BAND_LOGLEG_XM(2:NI, 1) = BAND_LOGLEG_XM(2:NI, 1)*EXP(LOGNORM(1:NI-1)-LOGNORM(2:NI))
      BAND_LOGLEG_XM(1:NI-1, 3) = BAND_LOGLEG_XM(1:NI-1, 3)*EXP(LOGNORM(2:NI)-LOGNORM(1:NI-1))

    END IF

    RETURN
  END FUNCTION BAND_LOGLEG_XM
!=======================================================================
  FUNCTION LEG_XP(NI, NJ, M, NORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY (1+X) IN FUNCTION (LEGENDRE) SPACE IN A MATRIX
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X) = SUM_N=|M|^INFTY A^M_N P^M_N(X)
! COEFFICIENTS OF THE FUNCTION (1+X)F(X) B^M_N AS
! (1+X)F(X) = SUM_N=|M|^INFTY B^M_N P^M_N(X)
! CAN BE EVAULATED BY "B^M_N = LEG_XP * A^M_N"
! 2. (MATSUSHIMA & MARCUS, 1997)
! GIVEN THE ORDER M,
! B^M_N = (N-ABS(M))/(2*N-1)A^M_N-1+A^M_N+(N+ABS(M)+1)/(2*N+3)A^M_N+1
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! NJ >> # OF COLUMNS OF THE OPERATOR
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! NORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! LEG_XP >> OPERATOR TO MULTIPLY (1+X) IN FUNCTION SPACE, SIZE: NI X NJ
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, NJ, M
    REAL(P8), DIMENSION(1:NI, 1:NJ) :: LEG_XP
    REAL(P8), DIMENSION(:), OPTIONAL :: NORM

    INTEGER(P4) :: NN, N, I, J, AM

    AM = ABS(M)

    LEG_XP = 0.D0                                                       ! ASSIGN ZERO MATRIX TO LEG_XP

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) LEG_XP(NN, NN-1) = (N-AM)/(2*N-1.D0)             ! SUBDIAGONAL COMPONENTS
      IF (NN .LE. NJ) LEG_XP(NN, NN) = 1.D0                             ! MAIN DIAGONAL COMPONENTS
      IF (NN+1 .LE. NJ) LEG_XP(NN, NN+1) = (N+AM+1.D0)/(2*N+3.D0)       ! SUPERDIAGONAL COMPONENTS
    END DO

    IF (PRESENT(NORM)) THEN
      IF (SIZE(NORM) .LT. MAX(NI, NJ)) THEN
        WRITE (*, *) 'LEG_XP: SIZE(NORM) TOO SMALL'
        WRITE (*, *) 'SIZE(NORM)=', SIZE(NORM)
        WRITE (*, *) 'NI=', NI, '  NJ=', NJ
        STOP
      END IF

      DO I = 1, NI
        DO J = 1, NJ
          LEG_XP(I, J) = LEG_XP(I, J)*(NORM(J)/NORM(I))
        END DO
      END DO
    END IF

    RETURN
  END FUNCTION LEG_XP
!=======================================================================
  FUNCTION LOGLEG_XP(NI, NJ, M, LOGNORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY (1+X) IN FUNCTION (LEGENDRE) SPACE IN A MATRIX
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X) = SUM_N=|M|^INFTY A^M_N P^M_N(X)
! COEFFICIENTS OF THE FUNCTION (1+X)F(X) B^M_N AS
! (1+X)F(X) = SUM_N=|M|^INFTY B^M_N P^M_N(X)
! CAN BE EVAULATED BY "B^M_N = LOGLEG_XP * A^M_N"
! 2. (MATSUSHIMA & MARCUS, 1997)
! GIVEN THE ORDER M,
! B^M_N = (N-ABS(M))/(2*N-1)A^M_N-1+A^M_N+(N+ABS(M)+1)/(2*N+3)A^M_N+1
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! NJ >> # OF COLUMNS OF THE OPERATOR
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! LOGNORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! LOGLEG_XP >> OPERATOR TO MULTIPLY (1+X) IN FUNCTION SPACE, SIZE: NI X NJ
! [DEPENDENCIES]:
! 1. LOGLEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'LOGNORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, NJ, M
    REAL(P8), DIMENSION(1:NI, 1:NJ) :: LOGLEG_XP
    REAL(P8), DIMENSION(:), OPTIONAL :: LOGNORM

    INTEGER(P4) :: NN, N, I, J, AM

    AM = ABS(M)

    LOGLEG_XP = 0.D0                                                    ! ASSIGN ZERO MATRIX TO LOGLEG_XP

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) LOGLEG_XP(NN, NN-1) = (N-AM)/(2*N-1.D0)          ! SUBDIAGONAL COMPONENTS
      IF (NN .LE. NJ) LOGLEG_XP(NN, NN) = 1.D0                          ! MAIN DIAGONAL COMPONENTS
      IF (NN+1 .LE. NJ) LOGLEG_XP(NN, NN+1) = (N+AM+1.D0)/(2*N+3.D0)    ! SUPERDIAGONAL COMPONENTS
    END DO

    IF (PRESENT(LOGNORM)) THEN
      IF (SIZE(LOGNORM) .LT. MAX(NI, NJ)) THEN
        WRITE (*, *) 'LOGLEG_XP: SIZE(LOGNORM) TOO SMALL'
        WRITE (*, *) 'SIZE(LOGNORM)=', SIZE(LOGNORM)
        WRITE (*, *) 'NI=', NI, '  NJ=', NJ
        STOP
      END IF

      DO I = 1, NI
        DO J = 1, NJ
          ! LOGLEG_XP(I,J)=LOGLEG_XP(I,J)*(NORM(J)/NORM(I))
          LOGLEG_XP(I, J) = LOGLEG_XP(I, J)*EXP(LOGNORM(J)-LOGNORM(I))
        END DO
      END DO
    END IF

    RETURN
  END FUNCTION LOGLEG_XP
!=======================================================================
  FUNCTION BAND_LEG_XP(NI, M, NORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY (1+X) IN FUNCTION (LEGENDRE) SPACE
! BASICALLY SAME AS LEG_XP
! HOWEVER, ONLY NON-TRIVIAL COMPONENTS (TRI-DIAGONALS) ARE COMPUTED
! AND SHOWN AS A TRIPLE BAND MATRIX (NI X 3) FORM FOR EFFICIENCY
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! NORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! BAND_LEG_XP >> BAND MATRIX STYLE OPER. TO MULTIPLY (1+X) IN FUNCTION
!               SPACE, SIZE: NIX3 (1: SUBDIAG., 2: DIAG., 3: SUPERDIAG.)
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, M
    REAL(P8), DIMENSION(NI, 3) :: BAND_LEG_XP
    REAL(P8), DIMENSION(:), OPTIONAL:: NORM

    INTEGER :: NN, N, I, J, AM

    AM = ABS(M)
    BAND_LEG_XP = 0.D0

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) BAND_LEG_XP(NN, 1) = (N-AM)/(2*N-1.D0)           ! SUBDIAGONAL COMPONENTS
      IF (NN .LE. NI) BAND_LEG_XP(NN, 2) = 1.D0                         ! MAIN DIAGONAL COMPONENTS
      IF (NN+1 .LE. NI) BAND_LEG_XP(NN, 3) = (N+AM+1.D0)/(2*N+3.D0)     ! SUPERDIAGONAL COMPONENTS
    END DO

    IF (PRESENT(NORM)) THEN
      IF (SIZE(NORM) .LT. NI) THEN
        WRITE (*, *) 'BAND_LEG_XP: SIZE(NORM) TOO SMALL'
        WRITE (*, *) 'SIZE(NORM)=', SIZE(NORM)
        WRITE (*, *) 'NI=', NI
        STOP
      END IF

      DO I = 1, NI
        IF (I .GT. 1) THEN
          BAND_LEG_XP(I, 1) = BAND_LEG_XP(I, 1)*(NORM(I-1)/NORM(I))
        END IF
        IF (I .LT. NI) THEN
          BAND_LEG_XP(I, 3) = BAND_LEG_XP(I, 3)*(NORM(I+1)/NORM(I))
        END IF
      END DO
    END IF

    RETURN
  END FUNCTION BAND_LEG_XP
!=======================================================================
  FUNCTION BAND_LOGLEG_XP(NI, M, LOGNORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY (1+X) IN FUNCTION (LEGENDRE) SPACE
! BASICALLY SAME AS LOGLEG_XP
! HOWEVER, ONLY NON-TRIVIAL COMPONENTS (TRI-DIAGONALS) ARE COMPUTED
! AND SHOWN AS A TRIPLE BAND MATRIX (NI X 3) FORM FOR EFFICIENCY
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! LOGNORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! BAND_LOGLEG_XP >> BAND MATRIX STYLE OPER. TO MULTIPLY (1+X) IN FUNCTION
!               SPACE, SIZE: NIX3 (1: SUBDIAG., 2: DIAG., 3: SUPERDIAG.)
! [DEPENDENCIES]:
! 1. LOGLEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'LOGNORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, M
    REAL(P8), DIMENSION(NI, 3) :: BAND_LOGLEG_XP
    REAL(P8), DIMENSION(:), OPTIONAL:: LOGNORM

    INTEGER :: NN, N, I, J, AM

    AM = ABS(M)
    BAND_LOGLEG_XP = 0.D0

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) BAND_LOGLEG_XP(NN, 1) = (N-AM)/(2*N-1.D0)        ! SUBDIAGONAL COMPONENTS
      IF (NN .LE. NI) BAND_LOGLEG_XP(NN, 2) = 1.D0                      ! MAIN DIAGONAL COMPONENTS
      IF (NN+1 .LE. NI) BAND_LOGLEG_XP(NN, 3) = (N+AM+1.D0)/(2*N+3.D0)  ! SUPERDIAGONAL COMPONENTS
    END DO

    !   IF(PRESENT(NORM)) THEN
    !     IF(SIZE(NORM) .LT. NI) THEN
    !       WRITE(*,*) 'BAND_LOGLEG_XP: SIZE(NORM) TOO SMALL'
    !       WRITE(*,*) 'SIZE(NORM)=',SIZE(NORM)
    !       WRITE(*,*) 'NI=',NI
    !       STOP
    !     ENDIF

    !     DO I = 1, NI
    !       IF(I .GT. 1) THEN
    !         BAND_LOGLEG_XP(I,1)=BAND_LOGLEG_XP(I,1)*(NORM(I-1)/NORM(I))
    !       ENDIF
    !       IF(I .LT. NI) THEN
    !         BAND_LOGLEG_XP(I,3)=BAND_LOGLEG_XP(I,3)*(NORM(I+1)/NORM(I))
    !       ENDIF
    !     ENDDO
    !   ENDIF

    IF (NI .LE. 1) RETURN

    IF (PRESENT(LOGNORM)) THEN
      IF (SIZE(LOGNORM) .LT. NI) THEN
        WRITE (*, *) 'BAND_LOGLEG_XP: SIZE(LOGNORM) TOO SMALL'
        WRITE (*, *) 'SIZE(LOGNORM)=', SIZE(LOGNORM)
        WRITE (*, *) 'NI=', NI
        STOP
      END IF

      BAND_LOGLEG_XP(2:NI, 1) = BAND_LOGLEG_XP(2:NI, 1)*EXP(LOGNORM(1:NI-1)-LOGNORM(2:NI))
      BAND_LOGLEG_XP(1:NI-1, 3) = BAND_LOGLEG_XP(1:NI-1, 3)*EXP(LOGNORM(2:NI)-LOGNORM(1:NI-1))

    END IF

    RETURN
  END FUNCTION BAND_LOGLEG_XP
!=======================================================================
  FUNCTION LEG_X(NI, NJ, M, NORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY X IN FUNCTION (LEGENDRE) SPACE IN A MATRIX
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X) = SUM_N=|M|^INFTY A^M_N P^M_N(X)
! COEFFICIENTS OF THE FUNCTION X*F(X) B^M_N AS
! X*F(X) = SUM_N=|M|^INFTY B^M_N P^M_N(X)
! CAN BE EVAULATED BY "B^M_N = LEG_X * A^M_N"
! 2. (MATSUSHIMA & MARCUS, 1997)
! GIVEN THE ORDER M,
! B^M_N = (N-ABS(M))/(2*N-1)A^M_N-1+(N+ABS(M)+1)/(2*N+3)A^M_N+1
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! NJ >> # OF COLUMNS OF THE OPERATOR
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! NORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! LEG_X >> OPERATOR TO MULTIPLY X IN FUNCTION SPACE, SIZE: NI X NJ
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, NJ, M
    REAL(P8), DIMENSION(1:NI, 1:NJ) :: LEG_X
    REAL(P8), DIMENSION(:), OPTIONAL :: NORM

    INTEGER(P4) :: NN, N, I, J, AM

    AM = ABS(M)

    LEG_X = 0.D0                                                        ! ASSIGN ZERO MATRIX TO LEG_X

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) LEG_X(NN, NN-1) = (N-AM)/(2*N-1.D0)              ! SUBDIAGONAL COMPONENTS
      IF (NN+1 .LE. NJ) LEG_X(NN, NN+1) = (N+AM+1.D0)/(2*N+3.D0)        ! SUPERDIAGONAL COMPONENTS
    END DO

    IF (PRESENT(NORM)) THEN
      IF (SIZE(NORM) .LT. MAX(NI, NJ)) THEN
        WRITE (*, *) 'LEG_X: SIZE(NORM) TOO SMALL'
        WRITE (*, *) 'SIZE(NORM)=', SIZE(NORM)
        WRITE (*, *) 'NI=', NI, '  NJ=', NJ
        STOP
      END IF

      DO I = 1, NI
        DO J = 1, NJ
          LEG_X(I, J) = LEG_X(I, J)*(NORM(J)/NORM(I))
        END DO
      END DO
    END IF

    RETURN
  END FUNCTION LEG_X
!=======================================================================
  FUNCTION LOGLEG_X(NI, NJ, M, LOGNORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY X IN FUNCTION (LEGENDRE) SPACE IN A MATRIX
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X) = SUM_N=|M|^INFTY A^M_N P^M_N(X)
! COEFFICIENTS OF THE FUNCTION X*F(X) B^M_N AS
! X*F(X) = SUM_N=|M|^INFTY B^M_N P^M_N(X)
! CAN BE EVAULATED BY "B^M_N = LOGLEG_X * A^M_N"
! 2. (MATSUSHIMA & MARCUS, 1997)
! GIVEN THE ORDER M,
! B^M_N = (N-ABS(M))/(2*N-1)A^M_N-1+(N+ABS(M)+1)/(2*N+3)A^M_N+1
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! NJ >> # OF COLUMNS OF THE OPERATOR
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! LOGNORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! LOGLEG_X >> OPERATOR TO MULTIPLY X IN FUNCTION SPACE, SIZE: NI X NJ
! [DEPENDENCIES]:
! 1. LOGLEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'LOGNORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, NJ, M
    REAL(P8), DIMENSION(1:NI, 1:NJ) :: LOGLEG_X
    REAL(P8), DIMENSION(:), OPTIONAL :: LOGNORM

    INTEGER(P4) :: NN, N, I, J, AM

    AM = ABS(M)

    LOGLEG_X = 0.D0                                                     ! ASSIGN ZERO MATRIX TO LOGLEG_X

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) LOGLEG_X(NN, NN-1) = (N-AM)/(2*N-1.D0)           ! SUBDIAGONAL COMPONENTS
      IF (NN+1 .LE. NJ) LOGLEG_X(NN, NN+1) = (N+AM+1.D0)/(2*N+3.D0)     ! SUPERDIAGONAL COMPONENTS
    END DO

    IF (PRESENT(LOGNORM)) THEN
      IF (SIZE(LOGNORM) .LT. MAX(NI, NJ)) THEN
        WRITE (*, *) 'LOGLEG_X: SIZE(LOGNORM) TOO SMALL'
        WRITE (*, *) 'SIZE(LOGNORM)=', SIZE(LOGNORM)
        WRITE (*, *) 'NI=', NI, '  NJ=', NJ
        STOP
      END IF

      DO I = 1, NI
        DO J = 1, NJ
          ! LOGLEG_X(I,J)=LOGLEG_X(I,J)*(NORM(J)/NORM(I))
          LOGLEG_X(I, J) = LOGLEG_X(I, J)*EXP(LOGNORM(J)-LOGNORM(I))
        END DO
      END DO
    END IF

    RETURN
  END FUNCTION LOGLEG_X
!=======================================================================
  FUNCTION BAND_LEG_X(NI, M, NORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY X IN FUNCTION (LEGENDRE) SPACE
! BASICALLY SAME AS LEG_X
! HOWEVER, ONLY NON-TRIVIAL COMPONENTS (TRI-DIAGONALS) ARE COMPUTED
! AND SHOWN AS A TRIPLE BAND MATRIX (NI X 3) FORM FOR EFFICIENCY
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! NORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! BAND_LEG_X >> BAND MATRIX STYLE OPERATOR TO MULTIPLY X IN FUNCTION
!               SPACE, SIZE: NIX3 (1: SUBDIAG., 2: DIAG., 3: SUPERDIAG.)
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, M
    REAL(P8), DIMENSION(NI, 3) :: BAND_LEG_X
    REAL(P8), DIMENSION(:), OPTIONAL:: NORM

    INTEGER :: NN, N, I, J, AM

    AM = ABS(M)
    BAND_LEG_X = 0.D0

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) BAND_LEG_X(NN, 1) = (N-AM)/(2*N-1.D0)            ! SUBDIAGONAL COMPONENTS
      IF (NN+1 .LE. NI) BAND_LEG_X(NN, 3) = (N+AM+1.D0)/(2*N+3.D0)      ! SUPERDIAGONAL COMPONENTS
    END DO

    IF (PRESENT(NORM)) THEN
      IF (SIZE(NORM) .LT. NI) THEN
        WRITE (*, *) 'BAND_LEG_X: SIZE(NORM) TOO SMALL'
        WRITE (*, *) 'SIZE(NORM)=', SIZE(NORM)
        WRITE (*, *) 'NI=', NI
        STOP
      END IF

      DO I = 1, NI
        IF (I .GT. 1) THEN
          BAND_LEG_X(I, 1) = BAND_LEG_X(I, 1)*(NORM(I-1)/NORM(I))
        END IF
        IF (I .LT. NI) THEN
          BAND_LEG_X(I, 3) = BAND_LEG_X(I, 3)*(NORM(I+1)/NORM(I))
        END IF
      END DO
    END IF

    RETURN
  END FUNCTION BAND_LEG_X
!=======================================================================
  FUNCTION BAND_LOGLEG_X(NI, M, LOGNORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY X IN FUNCTION (LEGENDRE) SPACE
! BASICALLY SAME AS LOGLEG_X
! HOWEVER, ONLY NON-TRIVIAL COMPONENTS (TRI-DIAGONALS) ARE COMPUTED
! AND SHOWN AS A TRIPLE BAND MATRIX (NI X 3) FORM FOR EFFICIENCY
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! LOGNORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! BAND_LOGLEG_X >> BAND MATRIX STYLE OPERATOR TO MULTIPLY X IN FUNCTION
!               SPACE, SIZE: NIX3 (1: SUBDIAG., 2: DIAG., 3: SUPERDIAG.)
! [DEPENDENCIES]:
! 1. LOGLEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'LOGNORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, M
    REAL(P8), DIMENSION(NI, 3) :: BAND_LOGLEG_X
    REAL(P8), DIMENSION(:), OPTIONAL:: LOGNORM

    INTEGER :: NN, N, I, J, AM

    AM = ABS(M)
    BAND_LOGLEG_X = 0.D0

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) BAND_LOGLEG_X(NN, 1) = (N-AM)/(2*N-1.D0)         ! SUBDIAGONAL COMPONENTS
      IF (NN+1 .LE. NI) BAND_LOGLEG_X(NN, 3) = (N+AM+1.D0)/(2*N+3.D0)   ! SUPERDIAGONAL COMPONENTS
    END DO

    !   IF(PRESENT(NORM)) THEN
    !     IF(SIZE(NORM) .LT. NI) THEN
    !       WRITE(*,*) 'BAND_LOGLEG_X: SIZE(NORM) TOO SMALL'
    !       WRITE(*,*) 'SIZE(NORM)=',SIZE(NORM)
    !       WRITE(*,*) 'NI=',NI
    !       STOP
    !     ENDIF

    !     DO I = 1, NI
    !       IF(I .GT. 1) THEN
    !         BAND_LOGLEG_X(I,1)=BAND_LOGLEG_X(I,1)*(NORM(I-1)/NORM(I))
    !       ENDIF
    !       IF(I .LT. NI) THEN
    !         BAND_LOGLEG_X(I,3)=BAND_LOGLEG_X(I,3)*(NORM(I+1)/NORM(I))
    !       ENDIF
    !     ENDDO
    !   ENDIF

    IF (NI .LE. 1) RETURN

    IF (PRESENT(LOGNORM)) THEN
      IF (SIZE(LOGNORM) .LT. NI) THEN
        WRITE (*, *) 'BAND_LOGLEG_X: SIZE(LOGNORM) TOO SMALL'
        WRITE (*, *) 'SIZE(LOGNORM)=', SIZE(LOGNORM)
        WRITE (*, *) 'NI=', NI
        STOP
      END IF

      BAND_LOGLEG_X(2:NI, 1) = BAND_LOGLEG_X(2:NI, 1)*EXP(LOGNORM(1:NI-1)-LOGNORM(2:NI))
      BAND_LOGLEG_X(1:NI-1, 3) = BAND_LOGLEG_X(1:NI-1, 3)*EXP(LOGNORM(2:NI)-LOGNORM(1:NI-1))

    END IF

    RETURN
  END FUNCTION BAND_LOGLEG_X
!=======================================================================
  FUNCTION LEG_XXDX(NI, NJ, M, NORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY (1-X**2)D/DX IN FUNCTION SPACE IN A MATRIX
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X) = SUM_N=|M|^INFTY A^M_N P^M_N(X)
! COEFFICIENTS OF THE FUNCTION (1-X**2)DF(X)/DX B^M_N AS
! (1-X**2)DF(X)/DX = SUM_N=|M|^INFTY B^M_N P^M_N(X)
! CAN BE EVAULATED BY "B^M_N = LEG_XXDX * A^M_N"
! 2. (MATSUSHIMA & MARCUS, 1997)
! GIVEN THE ORDER M,
! B^M_N =(N-1)(N-ABS(M))/(2*N-1)A^M_N-1+(N+2)(N+ABS(M)+1)/(2*N+3)A^M_N+1
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! NJ >> # OF COLUMNS OF THE OPERATOR
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! NORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! LEG_XXDX >> OPERATOR (1-X**2)D/DX IN FUNCTION SPACE, SIZE:NI X NJ
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, NJ, M
    REAL(P8), DIMENSION(1:NI, 1:NJ) :: LEG_XXDX
    REAL(P8), DIMENSION(:), OPTIONAL :: NORM

    INTEGER(P4) :: NN, N, I, J, AM

    AM = ABS(M)

    LEG_XXDX = 0.D0                                                     ! ASSIGN ZERO MATRIX TO LEG_XXDX

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) THEN
        LEG_XXDX(NN, NN-1) = -(N-1.D0)*(N-AM)/(2*N-1.D0)                ! SUBDIAGONAL COMPONENTS
      END IF
      IF (NN+1 .LE. NJ) THEN
        LEG_XXDX(NN, NN+1) = (N+2.D0)*(N+AM+1.D0)/(2*N+3.D0)            ! SUPERDIAGONAL COMPONENTS
      END IF
    END DO

    IF (PRESENT(NORM)) THEN
      IF (SIZE(NORM) .LT. MAX(NI, NJ)) THEN
        WRITE (*, *) 'LEG_XXDX: SIZE(NORM) TOO SMALL'
        WRITE (*, *) 'SIZE(NORM)=', SIZE(NORM)
        WRITE (*, *) 'NI=', NI, '  NJ=', NJ
        STOP
      END IF

      DO I = 1, NI
        DO J = 1, NJ
          LEG_XXDX(I, J) = LEG_XXDX(I, J)*(NORM(J)/NORM(I))
        END DO
      END DO
    END IF

    RETURN
  END FUNCTION LEG_XXDX
!=======================================================================
  FUNCTION LOGLEG_XXDX(NI, NJ, M, LOGNORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY (1-X**2)D/DX IN FUNCTION SPACE IN A MATRIX
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X) = SUM_N=|M|^INFTY A^M_N P^M_N(X)
! COEFFICIENTS OF THE FUNCTION (1-X**2)DF(X)/DX B^M_N AS
! (1-X**2)DF(X)/DX = SUM_N=|M|^INFTY B^M_N P^M_N(X)
! CAN BE EVAULATED BY "B^M_N = LOGLEG_XXDX * A^M_N"
! 2. (MATSUSHIMA & MARCUS, 1997)
! GIVEN THE ORDER M,
! B^M_N =(N-1)(N-ABS(M))/(2*N-1)A^M_N-1+(N+2)(N+ABS(M)+1)/(2*N+3)A^M_N+1
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! NJ >> # OF COLUMNS OF THE OPERATOR
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! LOGNORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! LOGLEG_XXDX >> OPERATOR (1-X**2)D/DX IN FUNCTION SPACE, SIZE:NI X NJ
! [DEPENDENCIES]:
! 1. LOGLEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'LOGNORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, NJ, M
    REAL(P8), DIMENSION(1:NI, 1:NJ) :: LOGLEG_XXDX
    REAL(P8), DIMENSION(:), OPTIONAL :: LOGNORM

    INTEGER(P4) :: NN, N, I, J, AM

    AM = ABS(M)

    LOGLEG_XXDX = 0.D0                                                  ! ASSIGN ZERO MATRIX TO LOGLEG_XXDX

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) THEN
        LOGLEG_XXDX(NN, NN-1) = -(N-1.D0)*(N-AM)/(2*N-1.D0)             ! SUBDIAGONAL COMPONENTS
      END IF
      IF (NN+1 .LE. NJ) THEN
        LOGLEG_XXDX(NN, NN+1) = (N+2.D0)*(N+AM+1.D0)/(2*N+3.D0)         ! SUPERDIAGONAL COMPONENTS
      END IF
    END DO

    IF (PRESENT(LOGNORM)) THEN
      IF (SIZE(LOGNORM) .LT. MAX(NI, NJ)) THEN
        WRITE (*, *) 'LOGLEG_XXDX: SIZE(LOGNORM) TOO SMALL'
        WRITE (*, *) 'SIZE(LOGNORM)=', SIZE(LOGNORM)
        WRITE (*, *) 'NI=', NI, '  NJ=', NJ
        STOP
      END IF

      DO I = 1, NI
        DO J = 1, NJ
          ! LOGLEG_XXDX(I,J)=LOGLEG_XXDX(I,J)*(NORM(J)/NORM(I))
          LOGLEG_XXDX(I, J) = LOGLEG_XXDX(I, J)*EXP(LOGNORM(J)-LOGNORM(I))
        END DO
      END DO
    END IF

    RETURN
  END FUNCTION LOGLEG_XXDX
!=======================================================================
  FUNCTION BAND_LEG_XXDX(NI, M, NORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY (1-X**2)D/DX IN FUNCTION (LEGENDRE) SPACE
! BASICALLY SAME AS LEG_XXDX
! HOWEVER, ONLY NON-TRIVIAL COMPONENTS (TRI-DIAGONALS) ARE COMPUTED
! AND SHOWN AS A TRIPLE BAND MATRIX (NI X 3) FORM FOR EFFICIENCY
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! NORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! BAND_LEG_XXDX >> BAND MATRIX STYLE OPER. (1-X**2)D/DX IN FUNCTION
!               SPACE, SIZE: NIX3 (1: SUBDIAG., 2: DIAG., 3: SUPERDIAG.)
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, M
    REAL(P8), DIMENSION(NI, 3) :: BAND_LEG_XXDX
    REAL(P8), DIMENSION(:), OPTIONAL:: NORM

    INTEGER :: NN, N, I, J, AM

    AM = ABS(M)
    BAND_LEG_XXDX = 0.D0

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) THEN
        BAND_LEG_XXDX(NN, 1) = -(N-1.D0)*(N-AM)/(2*N-1.D0)              ! SUBDIAGONAL COMPONENTS
      END IF
      IF (NN+1 .LE. NI) THEN
        BAND_LEG_XXDX(NN, 3) = (N+2.D0)*(N+AM+1.D0)/(2*N+3.D0)          ! SUPERDIAGONAL COMPONENTS
      END IF
    END DO

    IF (PRESENT(NORM)) THEN
      IF (SIZE(NORM) .LT. NI) THEN
        WRITE (*, *) 'BAND_LEG_XXDX: SIZE(NORM) TOO SMALL'
        WRITE (*, *) 'SIZE(NORM)=', SIZE(NORM)
        WRITE (*, *) 'NI=', NI
        STOP
      END IF

      DO I = 1, NI
        IF (I .GT. 1) THEN
          BAND_LEG_XXDX(I, 1) = BAND_LEG_XXDX(I, 1)*(NORM(I-1)/NORM(I))
        END IF
        IF (I .LT. NI) THEN
          BAND_LEG_XXDX(I, 3) = BAND_LEG_XXDX(I, 3)*(NORM(I+1)/NORM(I))
        END IF
      END DO
    END IF

    RETURN
  END FUNCTION BAND_LEG_XXDX
!=======================================================================
  FUNCTION BAND_LOGLEG_XXDX(NI, M, LOGNORM)
!=======================================================================
! [USAGE]:
! OPERATOR TO MULTIPLY (1-X**2)D/DX IN FUNCTION (LEGENDRE) SPACE
! BASICALLY SAME AS LOGLEG_XXDX
! HOWEVER, ONLY NON-TRIVIAL COMPONENTS (TRI-DIAGONALS) ARE COMPUTED
! AND SHOWN AS A TRIPLE BAND MATRIX (NI X 3) FORM FOR EFFICIENCY
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! LOGNORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! BAND_LOGLEG_XXDX >> BAND MATRIX STYLE OPER. (1-X**2)D/DX IN FUNCTION
!               SPACE, SIZE: NIX3 (1: SUBDIAG., 2: DIAG., 3: SUPERDIAG.)
! [DEPENDENCIES]:
! 1. LOGLEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'LOGNORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, M
    REAL(P8), DIMENSION(NI, 3) :: BAND_LOGLEG_XXDX
    REAL(P8), DIMENSION(:), OPTIONAL:: LOGNORM

    INTEGER :: NN, N, I, J, AM

    AM = ABS(M)
    BAND_LOGLEG_XXDX = 0.D0

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-1 .GE. 1) THEN
        BAND_LOGLEG_XXDX(NN, 1) = -(N-1.D0)*(N-AM)/(2*N-1.D0)           ! SUBDIAGONAL COMPONENTS
      END IF
      IF (NN+1 .LE. NI) THEN
        BAND_LOGLEG_XXDX(NN, 3) = (N+2.D0)*(N+AM+1.D0)/(2*N+3.D0)       ! SUPERDIAGONAL COMPONENTS
      END IF
    END DO

    !   IF(PRESENT(NORM)) THEN
    !     IF(SIZE(NORM) .LT. NI) THEN
    !       WRITE(*,*) 'BAND_LOGLEG_XXDX: SIZE(NORM) TOO SMALL'
    !       WRITE(*,*) 'SIZE(NORM)=',SIZE(NORM)
    !       WRITE(*,*) 'NI=',NI
    !       STOP
    !     ENDIF

    !     DO I = 1, NI
    !       IF(I .GT. 1) THEN
    !         BAND_LOGLEG_XXDX(I,1)=BAND_LOGLEG_XXDX(I,1)*(NORM(I-1)/NORM(I))
    !       ENDIF
    !       IF(I .LT. NI) THEN
    !         BAND_LOGLEG_XXDX(I,3)=BAND_LOGLEG_XXDX(I,3)*(NORM(I+1)/NORM(I))
    !       ENDIF
    !     ENDDO
    !   ENDIF

    IF (NI .LE. 1) RETURN

    IF (PRESENT(LOGNORM)) THEN
      IF (SIZE(LOGNORM) .LT. NI) THEN
        WRITE (*, *) 'BAND_LOGLEG_XXDX: SIZE(LOGNORM) TOO SMALL'
        WRITE (*, *) 'SIZE(LOGNORM)=', SIZE(LOGNORM)
        WRITE (*, *) 'NI=', NI
        STOP
      END IF

      BAND_LOGLEG_XXDX(2:NI, 1) = BAND_LOGLEG_XXDX(2:NI, 1)*EXP(LOGNORM(1:NI-1)-LOGNORM(2:NI))
      BAND_LOGLEG_XXDX(1:NI-1, 3) = BAND_LOGLEG_XXDX(1:NI-1, 3)*EXP(LOGNORM(2:NI)-LOGNORM(1:NI-1))

    END IF

    RETURN
  END FUNCTION BAND_LOGLEG_XXDX
!=======================================================================
  FUNCTION LEG_RAT_DEL2H(NI, NJ, M, ELL, NORM)
!=======================================================================
! [USAGE]:
! HORIZONTAL LAPLACIAN OPERATOR (DEL^2_PERP) IN FUNCTION SPACE
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X(R),PHI,Z) = SUM_N=|M|^INFTY A^M_N P^M_N(X(R))*EXP(I*M*PHI+I*AK*Z)
! COEFFICIENTS OF THE FUNCTION DEL2_PERP(F(X)) B^M_N AS
! DEL2_PERP(F(R)) = SUM_N=|M|^INFTY B^M_N P^M_N(X(R), PHI, Z)
!                                                  *EXP(I*M*PHI+I*AK*Z)
! CAN BE EVAULATED BY "B^M_N = LEG_RAT_DEL2H * A^M_N"
! 2. TO FIND THE RECURRENCE RELATIONSHIP,
! SEE EQUATION (107) IN THE APPENDIX OF MATSUSHIMA & MARCUS (1997)
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! NJ >> # OF COLUMNS OF THE OPERATOR
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! ELL >> THE MAP PARAMETER OPTIMIZING THE CONVERGENCE OF THE EXPANSION
! NORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! LEG_RAT_DEL2H >> HORIZONTAL LAPLACIAN OPERATOR, SIZE:NI X NJ
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, NJ, M
    REAL(P8), INTENT(IN) :: ELL
    REAL(P8), DIMENSION(NI, NJ) :: LEG_RAT_DEL2H
    REAL(P8), DIMENSION(:), OPTIONAL :: NORM

    INTEGER(P4) :: NN, N, I, J, AM
    REAL(P8) :: ELL2
    REAL(P8), PARAMETER :: S = 0.D0                                     ! IF S .NE. 0, THEN THE OPERATOR BECOMES
    ! [(1+X)**(-S)][DEL^2_PERP(M)][(1-X)**S]
    AM = ABS(M)
    LEG_RAT_DEL2H = 0
    ELL2 = ELL**2.

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-2 .GE. 1) LEG_RAT_DEL2H(NN, NN-2) = -(N-AM-1.D0)*(N-AM) &
                                                 *(N-2+S)*(N-1+S)/(2*N-3.)/(2*N-1.D0)/ELL2
      IF (NN-1 .GE. 1) LEG_RAT_DEL2H(NN, NN-1) = 2.D0*N*(N-AM)*(N-1+S) &
                                                 /(2*N-1.D0)/ELL2
      LEG_RAT_DEL2H(NN, NN) = (-2.D0*N*(N+1.D0)*(3.D0*N*N+3*N-AM**2-2) &
                               +2*S*(S-2)*(N*N+N+M*M-1.D0))/(2*N-1.D0)/(2*N+3.D0)/ELL2
      IF (NN+1 .LE. NJ) LEG_RAT_DEL2H(NN, NN+1) = 2.D0*(N+1) &
                                                  *(N+AM+1.D0)*(N+2-S)/(2*N+3.D0)/ELL2
      IF (NN+2 .LE. NJ) LEG_RAT_DEL2H(NN, NN+2) = -(N+AM+1.D0) &
                                                  *(N+AM+2.D0)*(N+3-S)*(N+2-S)/(2*N+3.D0)/(2*N+5.D0)/ELL2
    END DO

    IF (PRESENT(NORM)) THEN
      IF (SIZE(NORM) < MAX(NI, NJ)) THEN
        WRITE (*, *) 'LEG_RAT_DEL2H: SIZE(NORM) TOO SMALL'
        WRITE (*, *) 'SIZE(NORM)=', SIZE(NORM)
        WRITE (*, *) 'NI=', NI, '  NJ=', NJ
        STOP
      END IF

      DO I = 1, NI
        DO J = 1, NJ
          LEG_RAT_DEL2H(I, J) = LEG_RAT_DEL2H(I, J)*(NORM(J)/NORM(I))
        END DO
      END DO
    END IF

    RETURN
  END FUNCTION LEG_RAT_DEL2H
!=======================================================================
  FUNCTION LOGLEG_RAT_DEL2H(NI, NJ, M, ELL, LOGNORM)
!=======================================================================
! [USAGE]:
! HORIZONTAL LAPLACIAN OPERATOR (DEL^2_PERP) IN FUNCTION SPACE
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X(R),PHI,Z) = SUM_N=|M|^INFTY A^M_N P^M_N(X(R))*EXP(I*M*PHI+I*AK*Z)
! COEFFICIENTS OF THE FUNCTION DEL2_PERP(F(X)) B^M_N AS
! DEL2_PERP(F(X(R))) = SUM_N=|M|^INFTY B^M_N P^M_N(X(R), PHI, Z)
!                                                  *EXP(I*M*PHI+I*AK*Z)
! CAN BE EVAULATED BY "B^M_N = LOGLEG_RAT_DEL2H * A^M_N"
! 2. TO FIND THE RECURRENCE RELATIONSHIP,
! SEE EQUATION (107) IN THE APPENDIX OF MATSUSHIMA & MARCUS (1997)
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! NJ >> # OF COLUMNS OF THE OPERATOR
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! ELL >> THE MAP PARAMETER OPTIMIZING THE CONVERGENCE OF THE EXPANSION
! LOGNORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! LOGLEG_RAT_DEL2H >> HORIZONTAL LAPLACIAN OPERATOR, SIZE:NI X NJ
! [DEPENDENCIES]:
! 1. LOGLEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'LOGNORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, NJ, M
    REAL(P8), INTENT(IN) :: ELL
    REAL(P8), DIMENSION(NI, NJ) :: LOGLEG_RAT_DEL2H
    REAL(P8), DIMENSION(:), OPTIONAL :: LOGNORM

    INTEGER(P4) :: NN, N, I, J, AM
    REAL(P8) :: ELL2
    REAL(P8), PARAMETER :: S = 0.D0                                     ! IF S .NE. 0, THEN THE OPERATOR BECOMES
    ! [(1+X)**(-S)][DEL^2_PERP(M)][(1-X)**S]
    AM = ABS(M)
    LOGLEG_RAT_DEL2H = 0
    ELL2 = ELL**2.

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-2 .GE. 1) LOGLEG_RAT_DEL2H(NN, NN-2) = -(N-AM-1.D0)*(N-AM) &
                                                    *(N-2+S)*(N-1+S)/(2*N-3.)/(2*N-1.D0)/ELL2
      IF (NN-1 .GE. 1) LOGLEG_RAT_DEL2H(NN, NN-1) = 2.D0*N*(N-AM)*(N-1+S) &
                                                    /(2*N-1.D0)/ELL2
      LOGLEG_RAT_DEL2H(NN, NN) = (-2.D0*N*(N+1.D0)*(3.D0*N*N+3*N-AM**2-2) &
                                  +2*S*(S-2)*(N*N+N+M*M-1.D0))/(2*N-1.D0)/(2*N+3.D0)/ELL2
      IF (NN+1 .LE. NJ) LOGLEG_RAT_DEL2H(NN, NN+1) = 2.D0*(N+1) &
                                                     *(N+AM+1.D0)*(N+2-S)/(2*N+3.D0)/ELL2
      IF (NN+2 .LE. NJ) LOGLEG_RAT_DEL2H(NN, NN+2) = -(N+AM+1.D0) &
                                                     *(N+AM+2.D0)*(N+3-S)*(N+2-S)/(2*N+3.D0)/(2*N+5.D0)/ELL2
    END DO

    IF (PRESENT(LOGNORM)) THEN
      IF (SIZE(LOGNORM) < MAX(NI, NJ)) THEN
        WRITE (*, *) 'LOGLEG_RAT_DEL2H: SIZE(LOGNORM) TOO SMALL'
        WRITE (*, *) 'SIZE(LOGNORM)=', SIZE(LOGNORM)
        WRITE (*, *) 'NI=', NI, '  NJ=', NJ
        STOP
      END IF

      DO I = 1, NI
        DO J = 1, NJ
          ! LOGLEG_RAT_DEL2H(I,J)=LOGLEG_RAT_DEL2H(I,J)*(NORM(J)/NORM(I))
          LOGLEG_RAT_DEL2H(I, J) = LOGLEG_RAT_DEL2H(I, J)*EXP(LOGNORM(J)-LOGNORM(I))
        END DO
      END DO
    END IF

    RETURN
  END FUNCTION LOGLEG_RAT_DEL2H
!=======================================================================
  FUNCTION BAND_LEG_RAT_DEL2H(NI, M, ELL, NORM)
!=======================================================================
! [USAGE]:
! HORIZONTAL LAPLACIAN OPERATOR (DEL^2_PERP) IN FUNCTION SPACE
! BASICALLY SAME AS LEG_RAT_DEL2H
! HOWEVER, ONLY NON-TRIVIAL COMPONENTS (PENTA-DIAGONALS) ARE COMPUTED
! AND SHOWN AS A TRIPLE BAND MATRIX (NI X 5) FORM FOR EFFICIENCY
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! ELL >> THE MAP PARAMETER OPTIMIZING THE CONVERGENCE OF THE EXPANSION
! NORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! BAND_LEG_RAT_DEL2H >> BAND MATRIX STYLE OPER. OF HORIZONTAL LAPLACIAN
!                       SIZE: NI X 5
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, M
    REAL(P8), INTENT(IN) :: ELL
    REAL(P8), DIMENSION(NI, 5) :: BAND_LEG_RAT_DEL2H
    REAL(P8), DIMENSION(:), OPTIONAL :: NORM
    INTEGER(P4) :: NN, N, I, J, AM
    REAL(P8) :: ELL2
    REAL(P8), PARAMETER :: S = 0.D0                                     ! IF S .NE. 0, THEN THE OPERATOR BECOMES
    ! [(1+X)**(-S)][DEL^2_PERP(M)][(1-X)**S]
    AM = ABS(M)
    BAND_LEG_RAT_DEL2H = 0
    ELL2 = ELL**2.

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-2 .GE. 1) THEN
        BAND_LEG_RAT_DEL2H(NN, 1) = -(N-AM-1.D0)*(N-AM)*(N-2+S) &
                                    *(N-1+S)/(2*N-3.)/(2*N-1.D0)/ELL2
      END IF
      IF (NN-1 .GE. 1) THEN
        BAND_LEG_RAT_DEL2H(NN, 2) = 2.D0*N*(N-AM)*(N-1+S) &
                                    /(2*N-1.D0)/ELL2
      END IF
      BAND_LEG_RAT_DEL2H(NN, 3) = (-2.D0*N*(N+1.D0) &
                                   *(3.D0*N*N+3*N-AM**2-2) &
                                   +2*S*(S-2)*(N*N+N+M*M-1.D0))/(2*N-1.D0)/(2*N+3.D0)/ELL2
      IF (NN+1 .LE. NI) THEN
        BAND_LEG_RAT_DEL2H(NN, 4) = 2.D0*(N+1)*(N+AM+1.D0) &
                                    *(N+2-S)/(2*N+3.D0)/ELL2
      END IF
      IF (NN+2 .LE. NI) THEN
        BAND_LEG_RAT_DEL2H(NN, 5) = -(N+AM+1.D0) &
                                    *(N+AM+2.D0)*(N+3-S)*(N+2-S)/(2*N+3.D0)/(2*N+5.D0)/ELL2
      END IF
    END DO

    IF (PRESENT(NORM)) THEN
      IF (SIZE(NORM) .LT. NI) THEN
        WRITE (*, *) 'BAND_LEG_RAT_DEL2H: SIZE(NORM) TOO SMALL'
        WRITE (*, *) 'SIZE(NORM)=', SIZE(NORM)
        WRITE (*, *) 'NI=', NI
        STOP
      END IF

      DO I = 1, NI
        IF (I .GT. 2) BAND_LEG_RAT_DEL2H(I, 1) = &
          BAND_LEG_RAT_DEL2H(I, 1)*(NORM(I-2)/NORM(I))
        IF (I .GT. 1) BAND_LEG_RAT_DEL2H(I, 2) = &
          BAND_LEG_RAT_DEL2H(I, 2)*(NORM(I-1)/NORM(I))
        IF (I .LT. NI) BAND_LEG_RAT_DEL2H(I, 4) = &
          BAND_LEG_RAT_DEL2H(I, 4)*(NORM(I+1)/NORM(I))
        IF (I .LT. NI-1) BAND_LEG_RAT_DEL2H(I, 5) = &
          BAND_LEG_RAT_DEL2H(I, 5)*(NORM(I+2)/NORM(I))
      END DO
    END IF

    RETURN
  END FUNCTION BAND_LEG_RAT_DEL2H
!=======================================================================
  FUNCTION BAND_LOGLEG_RAT_DEL2H(NI, M, ELL, LOGNORM)
!=======================================================================
! [USAGE]:
! HORIZONTAL LAPLACIAN OPERATOR (DEL^2_PERP) IN FUNCTION SPACE
! BASICALLY SAME AS LOGLEG_RAT_DEL2H
! HOWEVER, ONLY NON-TRIVIAL COMPONENTS (PENTA-DIAGONALS) ARE COMPUTED
! AND SHOWN AS A TRIPLE BAND MATRIX (NI X 5) FORM FOR EFFICIENCY
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! ELL >> THE MAP PARAMETER OPTIMIZING THE CONVERGENCE OF THE EXPANSION
! LOGNORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! BAND_LOGLEG_RAT_DEL2H >> BAND MATRIX STYLE OPER. OF HORIZONTAL LAPLACIAN
!                       SIZE: NI X 5
! [DEPENDENCIES]:
! 1. LOGLEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'LOGNORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, M
    REAL(P8), INTENT(IN) :: ELL
    REAL(P8), DIMENSION(NI, 5) :: BAND_LOGLEG_RAT_DEL2H
    REAL(P8), DIMENSION(:), OPTIONAL :: LOGNORM
    INTEGER(P4) :: NN, N, I, J, AM
    REAL(P8) :: ELL2
    REAL(P8), PARAMETER :: S = 0.D0                                     ! IF S .NE. 0, THEN THE OPERATOR BECOMES
    ! [(1+X)**(-S)][DEL^2_PERP(M)][(1-X)**S]
    AM = ABS(M)
    BAND_LOGLEG_RAT_DEL2H = 0
    ELL2 = ELL**2.

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-2 .GE. 1) THEN
        BAND_LOGLEG_RAT_DEL2H(NN, 1) = -(N-AM-1.D0)*(N-AM)*(N-2+S) &
                                       *(N-1+S)/(2*N-3.)/(2*N-1.D0)/ELL2
      END IF
      IF (NN-1 .GE. 1) THEN
        BAND_LOGLEG_RAT_DEL2H(NN, 2) = 2.D0*N*(N-AM)*(N-1+S) &
                                       /(2*N-1.D0)/ELL2
      END IF
      BAND_LOGLEG_RAT_DEL2H(NN, 3) = (-2.D0*N*(N+1.D0) &
                                      *(3.D0*N*N+3*N-AM**2-2) &
                                      +2*S*(S-2)*(N*N+N+M*M-1.D0))/(2*N-1.D0)/(2*N+3.D0)/ELL2
      IF (NN+1 .LE. NI) THEN
        BAND_LOGLEG_RAT_DEL2H(NN, 4) = 2.D0*(N+1)*(N+AM+1.D0) &
                                       *(N+2-S)/(2*N+3.D0)/ELL2
      END IF
      IF (NN+2 .LE. NI) THEN
        BAND_LOGLEG_RAT_DEL2H(NN, 5) = -(N+AM+1.D0) &
                                       *(N+AM+2.D0)*(N+3-S)*(N+2-S)/(2*N+3.D0)/(2*N+5.D0)/ELL2
      END IF
    END DO

    !   IF(PRESENT(NORM)) THEN
    !     IF(SIZE(NORM) .LT. NI) THEN
    !       WRITE(*,*) 'BAND_LOGLEG_RAT_DEL2H: SIZE(NORM) TOO SMALL'
    !       WRITE(*,*) 'SIZE(NORM)=',SIZE(NORM)
    !       WRITE(*,*) 'NI=',NI
    !       STOP
    !     ENDIF

    !     DO I = 1, NI
    !       IF(I.GT.2) BAND_LOGLEG_RAT_DEL2H(I,1)= &
    !         BAND_LOGLEG_RAT_DEL2H(I,1)*(NORM(I-2)/NORM(I))
    !       IF(I.GT.1) BAND_LOGLEG_RAT_DEL2H(I,2)= &
    !         BAND_LOGLEG_RAT_DEL2H(I,2)*(NORM(I-1)/NORM(I))
    !       IF(I.LT.NI) BAND_LOGLEG_RAT_DEL2H(I,4)= &
    !         BAND_LOGLEG_RAT_DEL2H(I,4)*(NORM(I+1)/NORM(I))
    !       IF(I.LT.NI-1) BAND_LOGLEG_RAT_DEL2H(I,5)= &
    !         BAND_LOGLEG_RAT_DEL2H(I,5)*(NORM(I+2)/NORM(I))
    !     ENDDO
    !   ENDIF

    IF (NI .LE. 1) RETURN

    IF (PRESENT(LOGNORM)) THEN

      IF (SIZE(LOGNORM) .LT. NI) THEN
        WRITE (*, *) 'BAND_LOGLEG_RAT_DEL2H: SIZE(LOGNORM) TOO SMALL'
        WRITE (*, *) 'SIZE(LOGNORM)=', SIZE(LOGNORM)
        WRITE (*, *) 'NI=', NI
        STOP
      END IF

      BAND_LOGLEG_RAT_DEL2H(2:NI, 2) = &
        BAND_LOGLEG_RAT_DEL2H(2:NI, 2)*EXP(LOGNORM(1:NI-1)-LOGNORM(2:NI))

      BAND_LOGLEG_RAT_DEL2H(1:NI-1, 4) = &
        BAND_LOGLEG_RAT_DEL2H(1:NI-1, 4)*EXP(LOGNORM(2:NI)-LOGNORM(1:NI-1))

      IF (NI .LE. 2) RETURN

      BAND_LOGLEG_RAT_DEL2H(3:NI, 1) = &
        BAND_LOGLEG_RAT_DEL2H(3:NI, 1)*EXP(LOGNORM(1:NI-2)-LOGNORM(3:NI))

      BAND_LOGLEG_RAT_DEL2H(1:NI-2, 5) = &
        BAND_LOGLEG_RAT_DEL2H(1:NI-2, 5)*EXP(LOGNORM(3:NI)-LOGNORM(1:NI-2))

    END IF

    RETURN
  END FUNCTION BAND_LOGLEG_RAT_DEL2H
!=======================================================================
  FUNCTION LEG_RAT_DEL2(NI, NJ, M, AK, ELL, NORM)
!=======================================================================
! [USAGE]:
! LAPLACIAN OPERATOR (DEL^2) IN FUNCTION SPACE
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X(R),PHI,Z) = SUM_N=|M|^INFTY A^M_N P^M_N(X(R))*EXP(I*M*PHI+I*AK*Z)
! COEFFICIENTS OF THE FUNCTION DEL2(F(X)) B^M_N AS
! DEL2(F(X(R))) = SUM_N=|M|^INFTY B^M_N P^M_N(X(R), PHI, Z)
!                                                  *EXP(I*M*PHI+I*AK*Z)
! CAN BE EVAULATED BY "B^M_N = LEG_RAT_DEL2 * A^M_N"
! 2. TO FIND THE RECURRENCE RELATIONSHIP,
! SEE EQUATION (107) IN THE APPENDIX OF MATSUSHIMA & MARCUS (1997)
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! NJ >> # OF COLUMNS OF THE OPERATOR
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! ELL >> THE MAP PARAMETER OPTIMIZING THE CONVERGENCE OF THE EXPANSION
! NORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! LEG_RAT_DEL2 >> HORIZONTAL LAPLACIAN OPERATOR, SIZE:NI X NJ
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, NJ, M
    REAL(P8), INTENT(IN) :: AK
    REAL(P8), INTENT(IN) :: ELL
    REAL(P8), DIMENSION(NI, NJ) :: LEG_RAT_DEL2
    REAL(P8), DIMENSION(:), OPTIONAL :: NORM

    INTEGER(P4) :: NN, N, I, J, AM
    REAL(P8) :: ELL2
    REAL(P8), PARAMETER :: S = 0.D0                                    ! IF S .NE. 0, THEN THE OPERATOR BECOMES
    ! [(1+X)**(-S)][DEL^2(M)][(1-X)**S]
    AM = ABS(M)
    LEG_RAT_DEL2 = 0
    ELL2 = ELL**2.

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-2 .GE. 1) LEG_RAT_DEL2(NN, NN-2) = -(N-AM-1.D0)*(N-AM) &
                                                *(N-2+S)*(N-1+S)/(2*N-3.)/(2*N-1.D0)/ELL2
      IF (NN-1 .GE. 1) LEG_RAT_DEL2(NN, NN-1) = 2.D0*N*(N-AM)*(N-1+S) &
                                                /(2*N-1.D0)/ELL2
      LEG_RAT_DEL2(NN, NN) = (-2.D0*N*(N+1.D0)*(3.D0*N*N+3*N-AM**2-2) &
                              +2*S*(S-2)*(N*N+N+M*M-1.D0))/(2*N-1.D0)/(2*N+3.D0)/ELL2
      IF (NN+1 .LE. NJ) LEG_RAT_DEL2(NN, NN+1) = 2.D0*(N+1) &
                                                 *(N+AM+1.D0)*(N+2-S)/(2*N+3.D0)/ELL2
      IF (NN+2 .LE. NJ) LEG_RAT_DEL2(NN, NN+2) = -(N+AM+1.D0) &
                                                 *(N+AM+2.D0)*(N+3-S)*(N+2-S)/(2*N+3.D0)/(2*N+5.D0)/ELL2
    END DO

    DO NN = 1, NI
      LEG_RAT_DEL2(NN, NN) = LEG_RAT_DEL2(NN, NN)-AK**2.
    END DO

    IF (PRESENT(NORM)) THEN
      IF (SIZE(NORM) < MAX(NI, NJ)) THEN
        WRITE (*, *) 'LEG_RAT_DEL2: SIZE(NORM) TOO SMALL'
        WRITE (*, *) 'SIZE(NORM)=', SIZE(NORM)
        WRITE (*, *) 'NI=', NI, '  NJ=', NJ
        STOP
      END IF

      DO I = 1, NI
        DO J = 1, NJ
          LEG_RAT_DEL2(I, J) = LEG_RAT_DEL2(I, J)*(NORM(J)/NORM(I))
        END DO
      END DO
    END IF

    RETURN
  END FUNCTION LEG_RAT_DEL2
!=======================================================================
  FUNCTION LOGLEG_RAT_DEL2(NI, NJ, M, AK, ELL, LOGNORM)
!=======================================================================
! [USAGE]:
! LAPLACIAN OPERATOR (DEL^2) IN FUNCTION SPACE
! [NOTES]:
! 1. FOR A FUNCTION F(X) EXPENDED WITH COEFFCIEINT A^M_N AS
! F(X(R),PHI,Z) = SUM_N=|M|^INFTY A^M_N P^M_N(X(R))*EXP(I*M*PHI+I*AK*Z)
! COEFFICIENTS OF THE FUNCTION DEL2(F(X)) B^M_N AS
! DEL2(F(X(R))) = SUM_N=|M|^INFTY B^M_N P^M_N(X(R), PHI, Z)
!                                                  *EXP(I*M*PHI+I*AK*Z)
! CAN BE EVAULATED BY "B^M_N = LOGLEG_RAT_DEL2 * A^M_N"
! 2. TO FIND THE RECURRENCE RELATIONSHIP,
! SEE EQUATION (107) IN THE APPENDIX OF MATSUSHIMA & MARCUS (1997)
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! NJ >> # OF COLUMNS OF THE OPERATOR
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! ELL >> THE MAP PARAMETER OPTIMIZING THE CONVERGENCE OF THE EXPANSION
! LOGNORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! LOGLEG_RAT_DEL2 >> HORIZONTAL LAPLACIAN OPERATOR, SIZE:NI X NJ
! [DEPENDENCIES]:
! 1. LOGLEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'LOGNORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, NJ, M
    REAL(P8), INTENT(IN) :: AK
    REAL(P8), INTENT(IN) :: ELL
    REAL(P8), DIMENSION(NI, NJ) :: LOGLEG_RAT_DEL2
    REAL(P8), DIMENSION(:), OPTIONAL :: LOGNORM

    INTEGER(P4) :: NN, N, I, J, AM
    REAL(P8) :: ELL2
    REAL(P8), PARAMETER :: S = 0.D0                                     ! IF S .NE. 0, THEN THE OPERATOR BECOMES
    ! [(1+X)**(-S)][DEL^2(M)][(1-X)**S]
    AM = ABS(M)
    LOGLEG_RAT_DEL2 = 0
    ELL2 = ELL**2.

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-2 .GE. 1) LOGLEG_RAT_DEL2(NN, NN-2) = -(N-AM-1.D0)*(N-AM) &
                                                   *(N-2+S)*(N-1+S)/(2*N-3.)/(2*N-1.D0)/ELL2
      IF (NN-1 .GE. 1) LOGLEG_RAT_DEL2(NN, NN-1) = 2.D0*N*(N-AM)*(N-1+S) &
                                                   /(2*N-1.D0)/ELL2
      LOGLEG_RAT_DEL2(NN, NN) = (-2.D0*N*(N+1.D0)*(3.D0*N*N+3*N-AM**2-2) &
                                 +2*S*(S-2)*(N*N+N+M*M-1.D0))/(2*N-1.D0)/(2*N+3.D0)/ELL2
      IF (NN+1 .LE. NJ) LOGLEG_RAT_DEL2(NN, NN+1) = 2.D0*(N+1) &
                                                    *(N+AM+1.D0)*(N+2-S)/(2*N+3.D0)/ELL2
      IF (NN+2 .LE. NJ) LOGLEG_RAT_DEL2(NN, NN+2) = -(N+AM+1.D0) &
                                                    *(N+AM+2.D0)*(N+3-S)*(N+2-S)/(2*N+3.D0)/(2*N+5.D0)/ELL2
    END DO

    DO NN = 1, NI
      LOGLEG_RAT_DEL2(NN, NN) = LOGLEG_RAT_DEL2(NN, NN)-AK**2.
    END DO

    IF (PRESENT(LOGNORM)) THEN
      IF (SIZE(LOGNORM) < MAX(NI, NJ)) THEN
        WRITE (*, *) 'LOGLEG_RAT_DEL2: SIZE(LOGNORM) TOO SMALL'
        WRITE (*, *) 'SIZE(LOGNORM)=', SIZE(LOGNORM)
        WRITE (*, *) 'NI=', NI, '  NJ=', NJ
        STOP
      END IF

      DO I = 1, NI
        DO J = 1, NJ
          ! LOGLEG_RAT_DEL2(I,J)=LOGLEG_RAT_DEL2(I,J)*(NORM(J)/NORM(I))
          LOGLEG_RAT_DEL2(I, J) = LOGLEG_RAT_DEL2(I, J)*EXP(LOGNORM(J)-LOGNORM(I))
        END DO
      END DO
    END IF

    RETURN
  END FUNCTION LOGLEG_RAT_DEL2
!=======================================================================
  FUNCTION BAND_LEG_RAT_DEL2(NI, M, AK, ELL, NORM)
!=======================================================================
! [USAGE]:
! LAPLACIAN OPERATOR (DEL^2) IN FUNCTION SPACE
! BASICALLY SAME AS LEG_RAT_DEL2
! HOWEVER, ONLY NON-TRIVIAL COMPONENTS (PENTA-DIAGONALS) ARE COMPUTED
! AND SHOWN AS A TRIPLE BAND MATRIX (NI X 5) FORM FOR EFFICIENCY
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! AK >> THE AXIAL WAVE NUMBER
! ELL >> THE MAP PARAMETER OPTIMIZING THE CONVERGENCE OF THE EXPANSION
! NORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! BAND_LEG_RAT_DEL2 >> BAND MATRIX STYLE OPER. OF HORIZONTAL LAPLACIAN
!                       SIZE: NI X 5
! [DEPENDENCIES]:
! 1. LEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'NORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, M
    REAL(P8), INTENT(IN) :: AK
    REAL(P8), INTENT(IN) :: ELL
    REAL(P8), DIMENSION(NI, 5) :: BAND_LEG_RAT_DEL2
    REAL(P8), DIMENSION(:), OPTIONAL :: NORM
    INTEGER(P4) :: NN, N, I, J, AM
    REAL(P8) :: ELL2
    REAL(P8), PARAMETER :: S = 0.D0                                     ! IF S .NE. 0, THEN THE OPERATOR BECOMES
    ! [(1+X)**(-S)][DEL^2(M)][(1-X)**S]
    AM = ABS(M)
    BAND_LEG_RAT_DEL2 = 0
    ELL2 = ELL**2.

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-2 .GE. 1) THEN
        BAND_LEG_RAT_DEL2(NN, 1) = -(N-AM-1.D0)*(N-AM)*(N-2+S) &
                                   *(N-1+S)/(2*N-3.)/(2*N-1.D0)/ELL2
      END IF
      IF (NN-1 .GE. 1) THEN
        BAND_LEG_RAT_DEL2(NN, 2) = 2.D0*N*(N-AM)*(N-1+S) &
                                   /(2*N-1.D0)/ELL2
      END IF
      BAND_LEG_RAT_DEL2(NN, 3) = (-2.D0*N*(N+1.D0) &
                                  *(3.D0*N*N+3*N-AM**2-2) &
                                  +2*S*(S-2)*(N*N+N+M*M-1.D0))/(2*N-1.D0)/(2*N+3.D0)/ELL2
      IF (NN+1 .LE. NI) THEN
        BAND_LEG_RAT_DEL2(NN, 4) = 2.D0*(N+1)*(N+AM+1.D0) &
                                   *(N+2-S)/(2*N+3.D0)/ELL2
      END IF
      IF (NN+2 .LE. NI) THEN
        BAND_LEG_RAT_DEL2(NN, 5) = -(N+AM+1.D0) &
                                   *(N+AM+2.D0)*(N+3-S)*(N+2-S)/(2*N+3.D0)/(2*N+5.D0)/ELL2
      END IF
    END DO
    BAND_LEG_RAT_DEL2(:, 3) = BAND_LEG_RAT_DEL2(:, 3)-AK**2.

    IF (PRESENT(NORM)) THEN
      IF (SIZE(NORM) .LT. NI) THEN
        WRITE (*, *) 'BAND_LEG_RAT_DEL2: SIZE(NORM) TOO SMALL'
        WRITE (*, *) 'SIZE(NORM)=', SIZE(NORM)
        WRITE (*, *) 'NI=', NI
        STOP
      END IF

      DO I = 1, NI
        IF (I .GT. 2) BAND_LEG_RAT_DEL2(I, 1) = &
          BAND_LEG_RAT_DEL2(I, 1)*(NORM(I-2)/NORM(I))
        IF (I .GT. 1) BAND_LEG_RAT_DEL2(I, 2) = &
          BAND_LEG_RAT_DEL2(I, 2)*(NORM(I-1)/NORM(I))
        IF (I .LT. NI) BAND_LEG_RAT_DEL2(I, 4) = &
          BAND_LEG_RAT_DEL2(I, 4)*(NORM(I+1)/NORM(I))
        IF (I .LT. NI-1) BAND_LEG_RAT_DEL2(I, 5) = &
          BAND_LEG_RAT_DEL2(I, 5)*(NORM(I+2)/NORM(I))
      END DO
    END IF

    RETURN
  END FUNCTION BAND_LEG_RAT_DEL2
!=======================================================================
  FUNCTION BAND_LOGLEG_RAT_DEL2(NI, M, AK, ELL, LOGNORM)
!=======================================================================
! [USAGE]:
! LAPLACIAN OPERATOR (DEL^2) IN FUNCTION SPACE
! BASICALLY SAME AS LOGLEG_RAT_DEL2
! HOWEVER, ONLY NON-TRIVIAL COMPONENTS (PENTA-DIAGONALS) ARE COMPUTED
! AND SHOWN AS A TRIPLE BAND MATRIX (NI X 5) FORM FOR EFFICIENCY
! [INPUTS]:
! NI >> # OF ROWS OF THE OPERATOR (HIGHEST DEGREE TO BE CONSIDERED)
! M >> THE ORDER OF LEGENDRE POLYNOMIALS TO OPERATE ON
! AK >> THE AXIAL WAVE NUMBER
! ELL >> THE MAP PARAMETER OPTIMIZING THE CONVERGENCE OF THE EXPANSION
! LOGNORM >> OPTIONAL. IF PRESENT, NORMALIZED OPERATOR CAN BE RETURNED
! [OUTPUTS]:
! BAND_LOGLEG_RAT_DEL2 >> BAND MATRIX STYLE OPER. OF HORIZONTAL LAPLACIAN
!                       SIZE: NI X 5
! [DEPENDENCIES]:
! 1. LOGLEG_NORM(NDIM,M) FUNCTION @ MOD_LIN_LEGENDRE IF 'LOGNORM' IS INPUT
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA
! RE-CODED BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
    IMPLICIT NONE
    INTEGER(P4), INTENT(IN) :: NI, M
    REAL(P8), INTENT(IN) :: AK
    REAL(P8), INTENT(IN) :: ELL
    REAL(P8), DIMENSION(NI, 5) :: BAND_LOGLEG_RAT_DEL2
    REAL(P8), DIMENSION(:), OPTIONAL :: LOGNORM
    INTEGER(P4) :: NN, N, I, J, AM
    REAL(P8) :: ELL2
    REAL(P8), PARAMETER :: S = 0.D0                                     ! IF S .NE. 0, THEN THE OPERATOR BECOMES
    ! [(1+X)**(-S)][DEL^2(M)][(1-X)**S]
    AM = ABS(M)
    BAND_LOGLEG_RAT_DEL2 = 0
    ELL2 = ELL**2.

    DO NN = 1, NI
      N = AM+(NN-1)
      IF (NN-2 .GE. 1) THEN
        BAND_LOGLEG_RAT_DEL2(NN, 1) = -(N-AM-1.D0)*(N-AM)*(N-2+S) &
                                      *(N-1+S)/(2*N-3.)/(2*N-1.D0)/ELL2
      END IF
      IF (NN-1 .GE. 1) THEN
        BAND_LOGLEG_RAT_DEL2(NN, 2) = 2.D0*N*(N-AM)*(N-1+S) &
                                      /(2*N-1.D0)/ELL2
      END IF
      BAND_LOGLEG_RAT_DEL2(NN, 3) = (-2.D0*N*(N+1.D0) &
                                     *(3.D0*N*N+3*N-AM**2-2) &
                                     +2*S*(S-2)*(N*N+N+M*M-1.D0))/(2*N-1.D0)/(2*N+3.D0)/ELL2
      IF (NN+1 .LE. NI) THEN
        BAND_LOGLEG_RAT_DEL2(NN, 4) = 2.D0*(N+1)*(N+AM+1.D0) &
                                      *(N+2-S)/(2*N+3.D0)/ELL2
      END IF
      IF (NN+2 .LE. NI) THEN
        BAND_LOGLEG_RAT_DEL2(NN, 5) = -(N+AM+1.D0) &
                                      *(N+AM+2.D0)*(N+3-S)*(N+2-S)/(2*N+3.D0)/(2*N+5.D0)/ELL2
      END IF
    END DO
    BAND_LOGLEG_RAT_DEL2(:, 3) = BAND_LOGLEG_RAT_DEL2(:, 3)-AK**2.

    !   IF(PRESENT(NORM)) THEN
    !     IF(SIZE(NORM) .LT. NI) THEN
    !       WRITE(*,*) 'BAND_LOGLEG_RAT_DEL2: SIZE(NORM) TOO SMALL'
    !       WRITE(*,*) 'SIZE(NORM)=',SIZE(NORM)
    !       WRITE(*,*) 'NI=',NI
    !       STOP
    !     ENDIF

    !     DO I = 1, NI
    !       IF(I.GT.2) BAND_LOGLEG_RAT_DEL2(I,1)= &
    !         BAND_LOGLEG_RAT_DEL2(I,1)*(NORM(I-2)/NORM(I))
    !       IF(I.GT.1) BAND_LOGLEG_RAT_DEL2(I,2)= &
    !         BAND_LOGLEG_RAT_DEL2(I,2)*(NORM(I-1)/NORM(I))
    !       IF(I.LT.NI) BAND_LOGLEG_RAT_DEL2(I,4)= &
    !         BAND_LOGLEG_RAT_DEL2(I,4)*(NORM(I+1)/NORM(I))
    !       IF(I.LT.NI-1) BAND_LOGLEG_RAT_DEL2(I,5)= &
    !         BAND_LOGLEG_RAT_DEL2(I,5)*(NORM(I+2)/NORM(I))
    !     ENDDO
    !   ENDIF

    IF (NI .LE. 1) RETURN

    IF (PRESENT(LOGNORM)) THEN
      IF (SIZE(LOGNORM) .LT. NI) THEN
        WRITE (*, *) 'BAND_LOGLEG_RAT_DEL2: SIZE(LOGNORM) TOO SMALL'
        WRITE (*, *) 'SIZE(LOGNORM)=', SIZE(LOGNORM)
        WRITE (*, *) 'NI=', NI
        STOP
      END IF

      BAND_LOGLEG_RAT_DEL2(2:NI, 2) = &
        BAND_LOGLEG_RAT_DEL2(2:NI, 2)*EXP(LOGNORM(1:NI-1)-LOGNORM(2:NI))

      BAND_LOGLEG_RAT_DEL2(1:NI-1, 4) = &
        BAND_LOGLEG_RAT_DEL2(1:NI-1, 4)*EXP(LOGNORM(2:NI)-LOGNORM(1:NI-1))

      IF (NI .LE. 2) RETURN

      BAND_LOGLEG_RAT_DEL2(3:NI, 1) = &
        BAND_LOGLEG_RAT_DEL2(3:NI, 1)*EXP(LOGNORM(1:NI-2)-LOGNORM(3:NI))

      BAND_LOGLEG_RAT_DEL2(1:NI-2, 5) = &
        BAND_LOGLEG_RAT_DEL2(1:NI-2, 5)*EXP(LOGNORM(3:NI)-LOGNORM(1:NI-2))

    END IF

    RETURN
  END FUNCTION BAND_LOGLEG_RAT_DEL2
!=======================================================================

END MODULE MOD_LIN_LEGENDRE
