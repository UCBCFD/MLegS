MODULE MOD_LAYOUT ! LEVEL 3 MODULE
! FIELD/LAYOUT MANIPULATIONS
      USE OMP_LIB
      USE MPI
      USE MOD_MISC, ONLY : P4,P8,PI,IU,MSAVE,MLOAD                       ! LEVEL 0
!XUSE USE MOD_FD                                                         ! LEVEL 1
      USE MOD_EIG                                                        ! LEVEL 1
!XUSE USE MOD_LIN_LEGENDRE                                               ! LEVEL 1
!XUSE USE MOD_BANDMAT                                                    ! LEVEL 1
      USE MOD_SCALAR3                                                    ! LEVEL 2
      USE MOD_FFT                                                        ! LEVEL 2.5
      IMPLICIT NONE
      PRIVATE
!=======================================================================
!============================ PARAMETERS ===============================
!=======================================================================
      !1> VORTEX PAIR DATA
      PUBLIC:: PAIRDATA
      TYPE PAIRDATA
      INTEGER:: N,NOISE
      REAL(P8),DIMENSION(8):: X,Y,Q,H,B
      REAL(P8),DIMENSION(8):: DX,DY,K,DP
      REAL(P8),DIMENSION(8):: R,DR,RK,RP
      END TYPE
      TYPE(PAIRDATA),PUBLIC:: QPAIR

      !2> GAUSS DATA (FOR UNPERTURBED, EXACT SOLUTION)
      PUBLIC:: GAUS_DATA
      TYPE GAUS_DATA
      REAL(P8):: XC,YC,Q,B
      REAL(P8):: DX,DY,K,DP
      REAL(P8):: R,DR,RK,RP
      END TYPE
      TYPE(GAUS_DATA),PUBLIC:: GAUSDATA                                  ! MADE BY PAIRMAKE
!=======================================================================
!======================== PUBLIC DECLARATION ===========================
!=======================================================================
      ! ADD THE FIELD FUNCTION INTO AN ARRAY
      PUBLIC:: LAY
      ! MULTIPLY THE FIELD FUNCTION INTO AN ARRAY
      PUBLIC:: LAYP
      ! PUT RANDOM NOISE ON WITH A CERTAIN NOISE LEVEL
      PUBLIC:: NOISE
      ! FLATTEN NOISE
      PUBLIC:: FLATTEN
      ! SINUSOIDALLY PERTURBED GAUSSIAN
      PUBLIC:: GAUS
      ! CHECK A GAUSSIAN FIELD
      PUBLIC:: CHECKGAUSSIAN
      ! ! RANDOM NUMBER GENERATOR
      ! PUBLIC:: RAN0 ! THIS LEADS TO SEGMENTATION ERROR
!=======================================================================
!======================== CUSTOM DECLARATION ===========================
!=======================================================================
      ! COS(M*TH) FNC
      PUBLIC:: COSM
CONTAINS
!=======================================================================
!============================ SUBROUTINES ==============================
!=======================================================================
      SUBROUTINE LAY(A,FUNC,IPIN)
!=======================================================================
! [USAGE]: 
! ADD A REAL FIELD FUNCTION INTO AN ARRAY A IN A PPP SPACE
! [PARAMETERS]:
! A >> SCALAR-TYPE INPUT
! FUNC >> FIELD FUNCTION TO ADD (E.G., GAUS @ MOD_LAYOUT. SEE BELOW.)
!         FUNC SHOULD BE A FIELD IN A CYLINDRICAL COORDINATE (PPP SPACE)
! IPIN >> (OPTIONAL) CORRECTION FACTOR OF (1-X)**(-IPIN) IS ADDED.
!         DEFAULT IS 0.
! [DEPENDENCIES]:
! 1. TOFP(~) @ MOD_SCALAR3
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR):: A
      REAL(P8),EXTERNAL:: FUNC
      INTEGER,OPTIONAL:: IPIN

      REAL(P8):: R,Z,F,THR,THI
      INTEGER:: KK,NN,MM,IP

      IF(PRESENT(IPIN)) THEN
        IP = IPIN
      ELSE
        IP = 0
      ENDIF
      
      IF (A%SPACE.NE.PPP_SPACE) THEN
      CALL TOFP(A)                                                       ! TO A PPP SPACE
      ENDIF

      ! PPP: NDIMR/N1, NDIMTH, NDIMX/N2     
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Z,R,F,THR,THI) COLLAPSE(2)
      DO KK=1,SIZE(A%E,3) !NX
        DO MM=1,NTH ! PPP HAS ALL THETA POINTS
          DO NN=1,SIZE(A%E,1) !NR
            Z  = TFM%Z(KK+A%INX)                                         ! COLLOCATION POINTS IN Z FROM 0 TO ZLEN
            R  = TFM%R(NN+A%INR)                                         ! COLLOCATION POINTS IN R FROM 0 TO INFINITY
            F  = (1-TFM%X(NN+A%INR))**IP                                 ! COLLOCATION POINTS IN MAPPED R FROM -1 TO 1
            THR = TFM%THR(MM)                                            ! COLLOCATION POINTS IN THETA, ALTERNATE 
            THI = TFM%THI(MM)                                            ! HALF IN THETAR AND HALF IN THETAI
            A%E(NN,MM,KK)= A%E(NN,MM,KK)+&
                           CMPLX(FUNC(R,THR,Z),FUNC(R,THI,Z),P8)/F
          ENDDO                                                          ! ALL THESE TRANSFORM VECTORS ARE CREATED BY  SUBROUTINE LEGINIT, 
        ENDDO                                                            ! WHICH IS CALLED AT THE END OF READIN (IN READIN.F90)
      ENDDO
!$OMP END PARALLEL DO

      CALL CHOPDO(A)                                                     ! MAKING SURE OUT-OF-RANGE VALUES ARE ZERO
      RETURN
      END SUBROUTINE LAY
!=======================================================================
      SUBROUTINE LAYP(A,FUNC)
!=======================================================================
! [USAGE]: 
! MULTIPLY A REAL FIELD FUNCTION INTO AN ARRAY A IN A PPP SPACE
! [PARAMETERS]:
! A >> SCALAR-TYPE INPUT
! FUNC >> FIELD FUNCTION TO MULTIPLY (E.G., GAUS @ MOD_LAYOUT.SEE BELOW)
!         FUNC SHOULD BE A FIELD IN A CYLINDRICAL COORDINATE (PPP SPACE)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      TYPE(SCALAR):: A
      REAL(P8),EXTERNAL:: FUNC

      REAL(P8):: R,Z,F,THR,THI,AAA,BBB,FR,FI
      INTEGER:: KK,NN,MM

      IF(A%SPACE .NE. PPP_SPACE) THEN
        IF (MPI_RANK.EQ.0) WRITE(*,*) 'LAYP: INPUT IS NOT IN PPP_SPACE.'
        STOP
      ENDIF

      ! PPP: NDIMR/N1, NDIMTH, NDIMX/N2     
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Z,R,THR,THI,AAA,BBB,FR,FI) COLLAPSE(2) 
      DO KK=1,SIZE(A%E,3) !NX
        DO MM=1,NTH ! PPP HAS ALL THETA POINTS
          DO NN=1,SIZE(A%E,1) !NR
            Z  = TFM%Z(KK+A%INX)
            R  = TFM%R(NN+A%INR)
            THR = TFM%THR(MM)
            THI = TFM%THI(MM)
            AAA = REAL(A%E(NN,MM,KK))
            BBB = AIMAG(A%E(NN,MM,KK))
            FR = FUNC(R,THR,Z)
            FI = FUNC(R,THI,Z)
            A%E(NN,MM,KK)= CMPLX(AAA*FR,BBB*FI,P8)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL CHOPDO(A)
      RETURN
      END SUBROUTINE LAYP
!=======================================================================
      SUBROUTINE NOISE(A,ALEVEL)
!=======================================================================
! [USAGE]: 
! PUT RANDOM NOISE ON A SCALAR-TYPE VARIABLE A
! [PARAMETERS]:
! A >> SCALAR-TYPE INPUT WHERE NOISE IS INFUSED
! ALEVEL >> NOISE LEVEL
! [DEPENDENCIES]:
! 1. RAN0(~) @ MOD_LAYOUT
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      TYPE(SCALAR):: A
      REAL(P8):: ALEVEL

      REAL(P8):: FAC,REFVAL,X,Y,AKMAX,RANVAL
      INTEGER:: MM,NN,N,KK,ISEED

      CALL TOFF(A)

      CALL RANDOM_SEED()
      CALL RANDOM_NUMBER(RANVAL)
      ISEED = INT(1.D6*RANVAL)

      REFVAL = ALEVEL*MAXVAL(ABS(A%E))
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,REFVAL,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_IVP,IERR)
      AKMAX = MAXVAL(ABS(AK(:NTCHOP,:NXCHOP))) ! ALL PROCS HAVE FULL AK.

      IF(AKMAX.EQ.0.0) AKMAX=1.D0

      ! FFF: NRCHOPDIM, NTCHOPDIM/N2, NXCHOPDIM/N1
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(N,FAC,X,Y) COLLAPSE(2)
      DO MM=1,MIN(NTCHOP-A%INTH,SIZE(A%E,2)) !NTCHOP
        DO NN=1,NRCHOPS(MM+A%INTH) ! FFF HAS ALL RADIAL POINTS
        N = M(MM+A%INTH)+(NN-1)
          DO KK=1,SIZE(A%E,3) !NXCHOPDIM
            FAC=0.1*EXP(-12*((AK(MM+A%INTH,KK+A%INX)/AKMAX)**2+(FLOAT(N)/NRCHOP)**2))
            X=FAC*RAN0(ISEED)
            Y=FAC*RAN0(ISEED)
            A%E(NN,MM,KK)= A%E(NN,MM,KK) +CMPLX(X,Y)*REFVAL
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL CHOPDO(A)
      RETURN
      END SUBROUTINE NOISE
!=======================================================================
      SUBROUTINE FLATTEN(A,B)
!=======================================================================
! [USAGE]: 
! FLATTEN THE NOISE ACCORDING TO THE MAGNITUDE OF REFERENCE VALUE B
! [PARAMETERS]:
! A >> SCALAR-TYPE INPUT WHERE NOISE IS INFUSED
! B >> SCALAR-TYPE INPUT WHERE THE NOISE LEVEL IS REFERENCED
!      BOTH A AND B MUST BE IN A PPP SPACE
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR):: A,B

      REAL(P8):: REFVAL,THR,THI,REFR,REFI,FACR,FACI,DEVR,DEVI
      INTEGER:: KK,NN,MM

      IF ( A%SPACE.NE.PPP_SPACE .OR. B%SPACE.NE.PPP_SPACE ) THEN
            IF (MPI_RANK.EQ.0) WRITE(*,*) 'FLATTEN: BOTH INPUTS MUST BE IN PPP_SPACE'
            STOP
      ENDIF

      REFVAL = MAXVAL(ABS(REAL(A%E)))
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,REFVAL,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_IVP,IERR)
      
      IF(REFVAL.EQ.0.0) REFVAL=1.D0

      ! PPP: NDIMR/N1, NDIMTH, NDIMX/N2
!$OMP PARALLEL DO DEFAULT(SHARED)&
!$OMP&PRIVATE(THR,THI,REFR,REFI,FACR,FACI,DEVR,DEVI) COLLAPSE(2)
      DO KK=1,SIZE(A%E,3) !NX
        DO MM=1,NTH
          DO NN=1,SIZE(A%E,1) !NR
            THR = TFM%THR(MM)
            THI = TFM%THI(MM)
            REFR = REAL(B%E(NN,MM,KK))
            REFI = AIMAG(B%E(NN,MM,KK))
            FACR = TANH((REFR/REFVAL*20))**2.
            FACI = TANH(REFI/REFVAL*20)**2.
            DEVR = REAL(A%E(NN,MM,KK))-REFR
            DEVI = AIMAG(A%E(NN,MM,KK))-REFI
            A%E(NN,MM,KK)= CMPLX(REFR,REFI,P8)+&
                           CMPLX(DEVR*FACR,DEVI*FACI,P8)
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL CHOPDO(A)
      RETURN
      END SUBROUTINE FLATTEN
!=======================================================================
!       SUBROUTINE CHECKGAUSSIAN(A,XCENTER,YCENTER,IPIN)
! !=======================================================================
! ! [USAGE]: 
! ! CHECK TO SEE IF THERE IS REALLY A GAUSSIAN FIELD
! ! [PARAMETERS]:
! ! A >> SCALAR-TYPE INPUT WHERE NOISE IS INFUSED
! ! XCENTER >> X COORDINATE OF THE VORTEX CENTER
! ! YCENTER >> Y COORDINATE OF THE VORTEX CENTER
! ! IPIN >> SAME AS ONE IN LAY @ MOD_LAYOUT
! ! [DEPENDENCIES]:
! ! 1. (DE)ALLOCATE(~) @ MOD_SCALAR3
! ! 2. MSAVE(~) @ MOD_SCALAR3
! ! [NOTE]:
! ! CHECK TO SEE IF YOU REALLY HAVE A GAUSSIAN FIELD LOCATED 
! ! AT XCENTER, YCENTER
! ! THE FIELD SHOULD CONTAIN AN EXTRA (1-X)^(-IPIN) JUST AS "LAY" DOES.
! ! --- SCALAR "A" SHOULD BE IN PPP SPACE ---
! ! GENERALLY, YOU WOULD WANT TO PASS THE OMEGA_Z FIELD TO CHECK.
! ! [UPDATES]:
! ! RE-CODED BY SANGJOON LEE @ NOV 11 2020
! !=======================================================================
!       IMPLICIT NONE
!       TYPE(SCALAR)::A
!       REAL(P8):: XCENTER, YCENTER

!       TYPE(SCALAR):: B
!       REAL(P8):: RSQ1, RSQ2
!       REAL(P8), DIMENSION(NX):: TMP1, TMP2
!       INTEGER:: MM, NN, KK, IPIN, NRTMP

!       IF(A%SPACE .NE. PPP_SPACE) THEN
!         WRITE(*,*) 'CHECKGAUSSIAN: NOT IN PPP_SPACE.'
!         STOP
!       ENDIF

!       WRITE(*,*) '----------------'
!       WRITE(*,*)  'CHECKING FOR A GAUSSIAN AT',XCENTER,YCENTER

!       CALL ALLOCATE(B,PPP_SPACE)

!       CALL MSAVE(A%E(:,:,1),'omega.dat') ! MEMO>>>NEED TO CHANGE MSAVE FIRST

!       ! PPP: NDIMR/N1, NDIMTH, NDIMX/N2
!       NRTMP = SIZE(A%E,1)    
!       DO NN = 1,MIN(NR-A%INR,SIZE(A%E,1))
!         IF (TFM%R(NN+A%INR).GT.ELL) THEN
!           NRTMP = NN
!           EXIT
!         ENDIF
!       END DO

!       DO NN = 1,NRTMP
!         !IF (TFM%R(NN).LE.ELL) THEN
!         !VALUES OF EXP(R*R) BEYOND R=21 ARE TOO LARGE.
!         DO MM = 1,NTH

!           RSQ1 = (TFM%R(NN+A%INR)*COS(TFM%THR(MM)) - XCENTER)**2 + &
!           & (TFM%R(NN+A%INR)*SIN(TFM%THR(MM)) - YCENTER)**2
    
!           RSQ2 = (TFM%R(NN+A%INR)*COS(TFM%THI(MM)) - XCENTER)**2 + &
!           & (TFM%R(NN+A%INR)*SIN(TFM%THI(MM)) - YCENTER)**2

!           TMP1 = REAL(A%E(NN,MM,:))*EXP(RSQ1)*((1-TFM%X(NN+A%INR))**(IPIN))    ! FINAL R VALUES ARE INFINITE AND SO RESULTS ARE VERY LARGE R'S
!           TMP2 = IMAG(A%E(NN,MM,:))*EXP(RSQ2)*((1-TFM%X(NN+A%INR))**(IPIN))    ! ARE CORRPUTED: USUALLY BEYOND ABOUT R=6 (EXP(36) = 10E+16)

!           B%E(NN,MM,:) = CMPLX(TMP1,TMP2)
!           !IF (NN.GT.(3*NR/4)) EXIT

!         ENDDO
!         !ENDIF
!       ENDDO

!       CALL MSAVE(B%E(:,:,1),'firstPlanes.dat') ! MEMO>>>NEED TO CHANGE MSAVE FIRST

!       CALL DEALLOCATE(B)

!       RETURN
!       END SUBROUTINE CHECKGAUSSIAN
! ======================================================================
      SUBROUTINE CHECKGAUSSIAN(A,XCENTER,YCENTER,IPIN)
!=======================================================================
! [USAGE]: 
! CHECK TO SEE IF THERE IS REALLY A GAUSSIAN FIELD
! [PARAMETERS]:
! A >> SCALAR-TYPE INPUT WHERE NOISE IS INFUSED
! XCENTER >> X COORDINATE OF THE VORTEX CENTER
! YCENTER >> Y COORDINATE OF THE VORTEX CENTER
! IPIN >> SAME AS ONE IN LAY @ MOD_LAYOUT
! [DEPENDENCIES]:
! 1. (DE)ALLOCATE(~) @ MOD_SCALAR3
! 2. MSAVE(~) @ MOD_SCALAR3
! [NOTE]:
! CHECK TO SEE IF YOU REALLY HAVE A GAUSSIAN FIELD LOCATED 
! AT XCENTER, YCENTER
! THE FIELD SHOULD CONTAIN AN EXTRA (1-X)^(-IPIN) JUST AS "LAY" DOES.
! --- SCALAR "A" SHOULD BE IN PPP SPACE ---
! GENERALLY, YOU WOULD WANT TO PASS THE OMEGA_Z FIELD TO CHECK.
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 11 2020
!=======================================================================
      IMPLICIT NONE
      TYPE(SCALAR)::A
      REAL(P8):: XCENTER, YCENTER

      TYPE(SCALAR):: B
      COMPLEX(P8), DIMENSION(NDIMR,NDIMTH,NDIMX):: B_GLB
      REAL(P8):: RSQ1, RSQ2
      REAL(P8), DIMENSION(NX):: TMP1, TMP2
      INTEGER:: MM, NN, KK, IPIN, NRTMP

      IF(A%SPACE .NE. PPP_SPACE) THEN
            IF (MPI_RANK.EQ.0) WRITE(*,*) 'CHECKGAUSSIAN: NOT IN PPP_SPACE.'
            STOP
      ENDIF

      IF (MPI_RANK.EQ.0) THEN
            WRITE(*,*) '----------------'
            WRITE(*,*)  'CHECKING FOR A GAUSSIAN AT',XCENTER,YCENTER
      ENDIF

      CALL ALLOCATE(B,PPP_SPACE)

      !CALL MSAVE(A%E(:,:,1),'omega.dat') ! MEMO>>>NEED TO CHANGE MSAVE FIRST
      CALL MASSEMBLE(A%E, B_GLB, 2)
      IF (MPI_RANK.EQ.0) CALL MSAVE(B_GLB(:,:,1),'omega.dat')

      ! PPP: NDIMR/N1, NDIMTH, NDIMX/N2
      NRTMP = SIZE(A%E,1)    
      DO NN = 1,MIN(NR-A%INR,SIZE(A%E,1))
            IF (TFM%R(NN+A%INR).GT.ELL) THEN
            NRTMP = NN
            EXIT
            ENDIF
      END DO

      B%E = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(RSQ1,RSQ2,TMP1,TMP2) COLLAPSE(2)
      DO NN = 1,NRTMP
            !IF (TFM%R(NN).LE.ELL) THEN
            !VALUES OF EXP(R*R) BEYOND R=21 ARE TOO LARGE.
            DO MM = 1,NTH

            RSQ1 = (TFM%R(NN+A%INR)*COS(TFM%THR(MM)) - XCENTER)**2 + &
            & (TFM%R(NN+A%INR)*SIN(TFM%THR(MM)) - YCENTER)**2
      
            RSQ2 = (TFM%R(NN+A%INR)*COS(TFM%THI(MM)) - XCENTER)**2 + &
            & (TFM%R(NN+A%INR)*SIN(TFM%THI(MM)) - YCENTER)**2

            TMP1 = REAL(A%E(NN,MM,:))*EXP(RSQ1)*((1-TFM%X(NN+A%INR))**(IPIN))    ! FINAL R VALUES ARE INFINITE AND SO RESULTS ARE VERY LARGE R'S
            TMP2 = IMAG(A%E(NN,MM,:))*EXP(RSQ2)*((1-TFM%X(NN+A%INR))**(IPIN))    ! ARE CORRPUTED: USUALLY BEYOND ABOUT R=6 (EXP(36) = 10E+16)

            B%E(NN,MM,:) = CMPLX(TMP1,TMP2)
            !IF (NN.GT.(3*NR/4)) EXIT

            ENDDO
            !ENDIF
      ENDDO
!$OMP END PARALLEL DO

     !CALL MSAVE(B%E(:,:,1),'firstPlanes.dat') ! MEMO>>>NEED TO CHANGE MSAVE FIRST
      CALL MASSEMBLE(B%E, B_GLB, 2)
      IF (MPI_RANK.EQ.0) CALL MSAVE(B%E(:,:,1),'firstPlanes.dat')

      CALL DEALLOCATE(B)
      !DEALLOCATE(B_GLB)

      RETURN
      END SUBROUTINE CHECKGAUSSIAN
!=======================================================================
!============================= FUNCTIONS ===============================
!=======================================================================
      FUNCTION RAN0(IDUM)
!=======================================================================
! [USAGE]:
! CREATE THE UNIFORMLY DISTRIBUTED RANDON NUMBER.
! [INPUTS]:
! IDUM >> RANDOM NUMBER SEED
! [OUTPUTS]:
! RAN0 >> RANDOM NUMBER BETWEEN 0 AND 1 UNIFORMLY DISTRIBUTED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 18 2020
! THIS LEADS TO SEGMENTATION ERROR!!!
!=======================================================================
      INTEGER :: IDUM
      REAL(P8):: RAN0
      INTEGER,PARAMETER:: IA=16807,IM=2147483647,IQ=127773,IR=2836
      INTEGER,PARAMETER:: MASK=123459876
      REAL(P8),PARAMETER:: AM=1.0D0/IM
      INTEGER:: K

      IDUM=IEOR(IDUM,MASK)
      K=IDUM/IQ
      IDUM=IA*(IDUM-K*IQ)-IR*K
      IF(IDUM.LT.0) IDUM=IDUM+IM
      RAN0 = AM*IDUM
      IDUM=IEOR(IDUM,MASK)

      RETURN
      END FUNCTION RAN0
!=======================================================================
      FUNCTION GAUS(R,TH,Z)
!=======================================================================
! [USAGE]:
! PRODUCE THE SINUSOIDALLY PERTURBED GAUSSIAN FIELD
! RETURNS THE VALUE OF VORTICITY OF THE Q-VORTEX IN AXIAL DIRECTION
! [INPUTS]:
! R >> RADIAL COORDINATE
! TH >> AZIMUTHAL COORDINATE
! Z >> AXIAL COORDINATE
! [OUTPUTS]:
! GAUS >> SINUSOIDALLY PERTURBED GAUSSIAN FIELD VALUE
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 18 2020
!=======================================================================
      REAL(P8):: R,TH,Z
      REAL(P8):: GAUS
      
      !REAL(P8)::X,Y,R2,AK,R0,RK,DP,RP,XC,YC !Variable AK was used in parental scope
      REAL(P8)::X,Y,R2,AK1,R0,RK,DP,RP,XC,YC 

      AK1 = 2*PI*GAUSDATA%K /ZLEN                                        !DEFAULT K=1; 2*PI/ZLEN
      RK = 2*PI*GAUSDATA%RK/ZLEN                                         !DEFAULT = 0
      DP = PI/180.*GAUSDATA%DP                                           !DEFAULT = 0
      RP = PI/180.*GAUSDATA%RP                                           !DEFAULT = 0
      X = R*COS(TH)
      Y = R*SIN(TH)
      XC = GAUSDATA%XC + GAUSDATA%DX*SIN(AK1*Z-DP)                       !DEFAULT = 0
      YC = GAUSDATA%YC + GAUSDATA%DY*SIN(AK1*Z-DP)                       !DEFAULT = 0
      R0 = GAUSDATA%R  + GAUSDATA%DR*SIN(RK*Z-RP)                        !DEFAULT = 1, CORESIZE
      R2 = ((X-XC)**2.+(Y-YC)**2.)/R0**2.
      GAUS = GAUSDATA%Q /R0**2.*EXP(-GAUSDATA%B*R2)
      
      RETURN
      END FUNCTION GAUS
!=======================================================================
      FUNCTION COSM(R,TH,Z)
!=======================================================================
! [USAGE]:
! COS(M*TH)*(1-X)
! [INPUTS]:
! R >> RADIAL COORDINATE
! TH >> AZIMUTHAL COORDINATE
! Z >> AXIAL COORDINATE
! [OUTPUTS]:
! COSM >> COS(M*TH)
! [UPDATES]:
!=======================================================================
      REAL(P8):: R,TH,Z
      REAL(P8):: COSM
      
      ! COSM = (COS(2.D0*TH))*R*(2*ELL2)/(ELL2+R**2)
      ! COSM = COS(2.D0*TH)
      COSM = COS(12.D0*TH)+COS(3.D0*TH)+COS(5.D0*Z)+COS(9.D0*Z)
      
      RETURN
      END FUNCTION COSM
!=======================================================================
END MODULE MOD_LAYOUT