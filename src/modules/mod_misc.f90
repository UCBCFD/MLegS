MODULE MOD_MISC ! LEVEL 0 MODULE
! BASIC TOOLS AND UTILITIES
      USE OMP_LIB
      USE MPI
      IMPLICIT NONE
      PRIVATE
!=======================================================================
!============================ PARAMETERS ===============================
!=======================================================================
      INTEGER , PARAMETER, PUBLIC :: P4 = SELECTED_REAL_KIND(P=6,R=37)   ! SINGLE PRECISION
      INTEGER , PARAMETER, PUBLIC :: P8 = SELECTED_REAL_KIND(P=14,R=300) ! DOUBLE PRECISION
      REAL(P8), PARAMETER, PUBLIC :: PI = ACOS(-1.D0)                    ! PI VALUE 3.141..
      COMPLEX(P8),PARAMETER,PUBLIC :: IU = (0.D0,1.D0)                   ! PURE IMAGINARY #
      INTEGER , PUBLIC:: CP8_SIZE = 0
      REAL(P8), PUBLIC:: CFL = 0.D0                                      ! CFL NUMBER

      INTEGER(P4) :: IE, ISSAV                                           ! USED IN READCHAR
      INTEGER(P4) :: II                                                  ! USED IN READCHAR,SPYR,SPYC
      CHARACTER(LEN=256) :: BUF = " "                                    ! USED IN READCHAR
!=======================================================================
!======================== PUBLIC DECLARATION ===========================
!=======================================================================
      ! PRINT THE REAL-TIME CLOCK
      PUBLIC :: PRINT_REAL_TIME
      ! SHOW THE STRUCTURE OF AN INPUT MATRIX
      PUBLIC :: SPY
      ! PRINT THE VALUES OF AN INPUT MATRIX
      PUBLIC :: MCAT
      ! SAVE A MATRIX INTO A SPECIFIED FILE
      PUBLIC :: MSAVE
      ! LOAD A MATRIX FROM SAVED FILES
      PUBLIC :: MLOAD      
      ! CHANGE ASCII STRINGS INTO FLOATING NUMBER
      PUBLIC :: ATOF
      ! CONVERT AN INTEGER/REAL TO CHARACTER(LEN=3,4 or 6)
      PUBLIC :: ITOA3
      PUBLIC :: ITOA4
      PUBLIC :: RTOA6
! ======================================================================

! ====================== ARRAY/MATRIX UTILITIES ========================
      ! LINSPACE AS IN MATLAB
      PUBLIC :: LINSPACE
      ! SORT
      PUBLIC :: SORT
      ! PRODUCE SPECIAL MATRICES
      PUBLIC :: ZEROS                                                     ! ZERO RECTANGULAR MATRIX PRODUCTION
      PUBLIC :: EYE                                                       ! SQUARE IDENTITY MATRIX PRODUCTION
      PUBLIC :: DIAG                                                      ! PRODUCE A DIAGONAL MATRIX (1D INPUT) OR EXTRACT DIAGONAL ELEMENTS (2D INPUT)
! ======================================================================

!=======================================================================
!============================ INTERFACES ===============================
!=======================================================================
      INTERFACE SPY
        MODULE PROCEDURE SPYR
        MODULE PROCEDURE SPYC
      END INTERFACE

      INTERFACE MCAT
        MODULE PROCEDURE MCATI
        MODULE PROCEDURE MCATR
        MODULE PROCEDURE MCATC
        MODULE PROCEDURE MCAT1DR
        MODULE PROCEDURE MCAT1DC
      END INTERFACE

      INTERFACE MSAVE
        MODULE PROCEDURE MSAVEDR
        MODULE PROCEDURE MSAVEDC
        MODULE PROCEDURE MSAVE1R
        MODULE PROCEDURE MSAVE1C
        MODULE PROCEDURE MSAVE3DIMENSIONALDATACOMPLEX
      END INTERFACE

      INTERFACE MLOAD
        MODULE PROCEDURE MLOAD1R
        MODULE PROCEDURE MLOAD1C
        MODULE PROCEDURE MLOADDR
        MODULE PROCEDURE MLOADDC
        MODULE PROCEDURE MLOADR0
        MODULE PROCEDURE MLOADC0
      END INTERFACE

      INTERFACE SORT
        MODULE PROCEDURE SORTR
        MODULE PROCEDURE SORTR_I
        MODULE PROCEDURE SORTR_IO
        MODULE PROCEDURE SORTR_IND
        MODULE PROCEDURE SORTC
        MODULE PROCEDURE SORTC_I
        MODULE PROCEDURE SORTC_IO
        MODULE PROCEDURE SORTC_IND
      END INTERFACE       

      INTERFACE ZEROS
        MODULE PROCEDURE ZEROS1
      END INTERFACE

      INTERFACE EYE
        MODULE PROCEDURE EYE1
      END INTERFACE

      INTERFACE DIAG
        MODULE PROCEDURE DIAG1R
        MODULE PROCEDURE DIAG2R
        MODULE PROCEDURE DIAG1C
        MODULE PROCEDURE DIAG2C
      END INTERFACE
! ======================================================================

CONTAINS
!=======================================================================
!============================ SUBROUTINES ==============================
!=======================================================================
      SUBROUTINE PRINT_REAL_TIME()
!=======================================================================
! [USAGE]: 
! EXPRESS THE DATE AND TIME WHEN THE SUBROUTINE IS CALLED
! [UPDATES]:
! WRITTEN BY SANGJOON LEE @ OCT 26 2020
!=======================================================================
      IMPLICIT NONE
      CHARACTER*8         ::  DATE
      CHARACTER*10        ::  NOW
      CHARACTER*5         ::  ZONE
      INTEGER(P4)         ::  VALS(1:8)

      CALL DATE_AND_TIME(DATE,NOW,ZONE,VALS)                             ! INTRINSIC SUBROUTINE. OUTPUT THE REAL-TIME CLOCK

      WRITE(*,101) VALS(1),VALS(2),VALS(3),VALS(5),VALS(6),VALS(7)
      WRITE(*,*) ''
 101  FORMAT(' @ 'I0.4,'-',I0.2,'-',I0.2,' ',I0.2,':',I0.2,':',I0.2)

      RETURN
      END SUBROUTINE PRINT_REAL_TIME
!=======================================================================
      SUBROUTINE SPYR(A)
!=======================================================================
! [USAGE]: 
! PRINT THE STRUCTURE OF A GIVEN REAL MATRIX A IN A FANCY FORM
! IF ABS. OF AN ELEMENT IS MACHINE-PRECISION 1, PRINT '1'
!               IS EXACTLY 0, PRINT ' '
!               IS IN THE ORDER OF N(.NE.0), PRINT CORSPNDING SYMBOL
!               ACCORDING TO CSR AND TBL.
! [PARAMETERS]:
! A >> THE NI X NJ MATRIX TO BE PRINTED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),DIMENSION(:,:) :: A
      INTEGER(P4) :: NI,NJ
      INTEGER(P4),PARAMETER:: NC=70,NCM=NC-1
      REAL(P8),PARAMETER:: EPS=1.E-11
      CHARACTER(LEN=80) :: B,C
      CHARACTER(LEN=60) :: CSR,TBL
      INTEGER(P4) :: I,J,JB,JE,JC,JJ,L,AT !II
      REAL(P8) :: X

      NI = SIZE(A,1)
      NJ = SIZE(A,2)

      B='----+----1----+----2----+----3----+----4'// &
        '----+----5----+----6----+----7----+----8'

      CSR='LOG10:-2----+----1----+----0----+----1----+----2----'
      TBL='   (-) ________...ihgfedcba@ABCDEFGHI!!!!~~~~~~~~(+)'
      WRITE(*,600) CSR
      WRITE(*,600) TBL
      AT=INDEX(TBL,'@')

      DO J=1,(NJ+NCM)/NC
        JB=NC*J-NCM
        JE=MIN(NC*J,NJ)
        WRITE(*,*) ' '
        WRITE(*,62) B(1:MOD(JE-1,NC)+1)

        DO I=1,NI
          II=1
          DO JJ=JB,JE
            X=ABS(A(I,JJ))
            JC=MOD(JJ-1,NC)+1
            IF(X.GT.1.D0-EPS .AND. X.LT.1.D0+EPS) THEN
              C(JC:JC)='1'
            ELSE IF(X.EQ.0.0) THEN
              C(JC:JC)=' '
            ELSE 
              L=LOG10(X)
              IF(L.LT.-20.) THEN
                C(JC:JC)='_'
              ELSE IF(L.GT.20.) THEN
                C(JC:JC)='~'
              ELSE
                L=L+AT
                C(JC:JC)=TBL(L:L)
              ENDIF
            ENDIF
          ENDDO
          WRITE(*,61) I,C(1:MOD(JE-1,NC)+1)
        ENDDO
      ENDDO

 62   FORMAT(3X,':',A)
 61   FORMAT(I3,':',A)
 600  FORMAT(A60)

      RETURN
      END SUBROUTINE SPYR
!=======================================================================
      SUBROUTINE SPYC(A)
!=======================================================================
! [USAGE]: 
! PRINT THE STRUCTURE OF A GIVEN COMPLEX MATRIX A IN A FANCY FORM
! IF THE ABS. OF AN ELEMENT IS MACHINE-PRECISION 1, PRINT '1'
!               IS EXACTLY 0, PRINT ']'
!               IS IN THE ORDER OF N(.NE.0), PRINT CORSPNDING SYMBOL
!               ACCORDING TO CSR AND TBL.
! [PARAMETERS]:
! A >> THE NI X NJ COMPLEX MATRIX TO BE PRINTED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8),DIMENSION(:,:) :: A
      INTEGER(P4) :: NI,NJ
      INTEGER(P4),PARAMETER:: NC=70,NCM=NC-1
      REAL(P8),PARAMETER:: EPS=1.E-11
      CHARACTER(LEN=80) :: B,C
      CHARACTER(LEN=60) :: CSR,TBL
      INTEGER(P4) :: I,J,JB,JE,JC,JJ,L,AT !II
      REAL(P8) :: X

      NI = SIZE(A,1)
      NJ = SIZE(A,2)

      B='----+----1----+----2----+----3----+----4'// &
        '----+----5----+----6----+----7----+----8'

      CSR='LOG10:-2----+----1----+----0----+----1----+----2----'
      TBL='   (-) ________...ihgfedcba@ABCDEFGHI!!!!~~~~~~~~(+)'
      WRITE(*,600) CSR
      WRITE(*,600) TBL
      AT=INDEX(TBL,'@')

      DO J=1,(NJ+NCM)/NC
        JB=NC*J-NCM
        JE=MIN(NC*J,NJ)
        WRITE(*,*) ' '
        WRITE(*,62) B(1:MOD(JE-1,NC)+1)

        DO I=1,NI
          II=1
          DO JJ=JB,JE
            X=ABS(A(I,JJ))
            JC=MOD(JJ-1,NC)+1
            IF(X.GT.1.D0-EPS .AND. X.LT.1.D0+EPS) THEN
              C(JC:JC)='1'
            ELSE IF(X.EQ.0.0) THEN
              C(JC:JC)=']'
            ELSE 
              L=LOG10(X)
              IF(L.LT.-20.) THEN
                C(JC:JC)='_'
              ELSE IF(L.GT.20.) THEN
                C(JC:JC)='~'
              ELSE
                L=L+AT
                C(JC:JC)=TBL(L:L)
              ENDIF
            ENDIF
          ENDDO
          WRITE(*,61) I,C(1:MOD(JE-1,NC)+1)
        ENDDO
      ENDDO

 62   FORMAT(3X,':',A)
 61   FORMAT(I3,':',A)
 600  FORMAT(A60)

      RETURN
      END SUBROUTINE SPYC
!=======================================================================
SUBROUTINE MCATI(A)
!=======================================================================
! [USAGE]: 
! PRINT A GIVEN REAL MATRIX A IN A ORGANIZED MANNER
! [PARAMETERS]:
! A >> THE NI X NJ REAL MATRIX TO BE PRINTED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      INTEGER,DIMENSION(:,:),INTENT(IN):: A
      INTEGER :: NI,NJ
      INTEGER :: I,J,JB,JE,JJ

      NI = SIZE(A,1)
      NJ = SIZE(A,2)

      DO J=1,(NJ+7)/8
        JB=8*J-7
        JE=MIN(8*J,NJ)
        WRITE(*,*) ' '
        WRITE(*,61) (JJ-1,JJ=JB,JE)
        DO I=1,NI
          WRITE(*,60) I-1, (A(I,JJ),JJ=JB,JE)
        ENDDO
      ENDDO

  61   FORMAT(4X,8I10)
  60   FORMAT(I3,':',1P8I10)

      RETURN
      END SUBROUTINE MCATI
!=======================================================================
      SUBROUTINE MCATR(A)
!=======================================================================
! [USAGE]: 
! PRINT A GIVEN REAL MATRIX A IN A ORGANIZED MANNER
! [PARAMETERS]:
! A >> THE NI X NJ REAL MATRIX TO BE PRINTED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),DIMENSION(:,:),INTENT(IN):: A
      INTEGER :: NI,NJ
      INTEGER :: I,J,JB,JE,JJ

      NI = SIZE(A,1)
      NJ = SIZE(A,2)

      DO J=1,(NJ+7)/8
        JB=8*J-7
        JE=MIN(8*J,NJ)
        WRITE(*,*) ' '
        WRITE(*,61) (JJ,JJ=JB,JE)
        DO I=1,NI
          WRITE(*,60) I, (A(I,JJ),JJ=JB,JE)
        ENDDO
      ENDDO

 61   FORMAT(4X,8I10)
 60   FORMAT(I3,':',1P8E10.2)

      RETURN
      END SUBROUTINE MCATR
!=======================================================================
      SUBROUTINE MCATC(A)
!=======================================================================
! [USAGE]: 
! PRINT A GIVEN COMPLEX MATRIX A IN A ORGANIZED MANNER
! [PARAMETERS]:
! A >> THE NI X NJ COMPLEX MATRIX TO BE PRINTED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8),DIMENSION(:,:),INTENT(IN):: A
      INTEGER :: NI,NJ
      INTEGER :: I,J,JB,JE,JJ

      NI = SIZE(A,1)
      NJ = SIZE(A,2)

      DO J=1,(NJ+3)/4
        JB=4*J-3
        JE=MIN(4*J,NJ)
        WRITE(*,*) ' '
        WRITE(*,61) (JJ,JJ=JB,JE)
        DO I=1,NI
          WRITE(*,60) I,(A(I,JJ),JJ=JB,JE)
        ENDDO
      ENDDO

 61   FORMAT(4X,8I23)
 60   FORMAT(I3,':',1P8(S,E11.3E3,SP,E11.3E3,'i'))

      RETURN
      END SUBROUTINE MCATC
!=======================================================================
      SUBROUTINE MCAT1DR(A)
!=======================================================================
! [USAGE]: 
! PRINT A GIVEN REAL VECTOR, MORE PRECISE THAN MATRIX PRINT SUBR. MCATR
! [PARAMETERS]:
! A >> THE N X 1 REAL VECTOR TO BE PRINTED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),DIMENSION(:),INTENT(IN):: A
      INTEGER:: I

      WRITE(*,60) (I,A(I),I=1,SIZE(A))
 60   FORMAT(I3,':',1PE24.15E3)

      RETURN
      END SUBROUTINE MCAT1DR
!=======================================================================
      SUBROUTINE MCAT1DC(A)
!=======================================================================
! [USAGE]: 
! PRINT A GIVEN COMPLEX VECTOR, MORE PRECISE THAN MCATC SUBROUTINE
! [PARAMETERS]:
! A >> THE N X 1 REAL VECTOR TO BE PRINTED
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8),DIMENSION(:),INTENT(IN):: A
      INTEGER:: I

      WRITE(*,60) (I,A(I),I=1,SIZE(A))
 60   FORMAT(S,I3,':',S,E24.15E3,SP,E24.15E3,'i')

      RETURN
      END SUBROUTINE MCAT1DC
!=======================================================================
      SUBROUTINE MSAVEDR(A,FN)
!=======================================================================
! [USAGE]: 
! STORE A REAL MATRIX VALUE INTO FN
! [PARAMETERS]:
! A >> A REAL MATRIX
! FN >> FILENAME TO STORE A VALUE
! [DEPENDENCIES]:
! 1. MSAVER(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),DIMENSION(:,:),INTENT(IN) :: A
      CHARACTER(LEN=*),INTENT(IN):: FN

      CALL MSAVER(FN,A,SIZE(A,1),SIZE(A,2))

      RETURN
      END SUBROUTINE MSAVEDR
!=======================================================================
      SUBROUTINE MSAVEDC(A,FN)
!=======================================================================
! [USAGE]: 
! STORE A COMPLEX MATRIX VALUE INTO FN
! [PARAMETERS]:
! A >> A COMPLEX MATRIX
! FN >> FILENAME TO STORE A VALUE
! [DEPENDENCIES]:
! 1. MSAVER(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8),DIMENSION(:,:),INTENT(IN) :: A
      CHARACTER(LEN=*),INTENT(IN):: FN

      CALL MSAVEC(FN,A,SIZE(A,1),SIZE(A,2))

      RETURN
      END SUBROUTINE MSAVEDC
!=======================================================================
      SUBROUTINE MSAVE1R(A,FN)
!=======================================================================
! [USAGE]: 
! STORE A REAL VECTOR VALUE INTO FN
! [PARAMETERS]:
! A >> A REAL VECTOR
! FN >> FILENAME TO STORE A VALUE
! [DEPENDENCIES]:
! 1. MSAVER(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),DIMENSION(:),INTENT(IN) :: A
      CHARACTER(LEN=*),INTENT(IN):: FN

      CALL MSAVER(FN,A,SIZE(A),1)

      RETURN
      END SUBROUTINE MSAVE1R
!=======================================================================
      SUBROUTINE MSAVE1C(A,FN)
!=======================================================================
! [USAGE]: 
! STORE A COMPLEX VECTOR VALUE INTO FN
! [PARAMETERS]:
! A >> A COMPLEX VECTOR
! FN >> FILENAME TO STORE A VALUE
! [DEPENDENCIES]:
! 1. MSAVER(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8),DIMENSION(:),INTENT(IN) :: A
      CHARACTER(LEN=*),INTENT(IN):: FN

      CALL MSAVEC(FN,A,SIZE(A),1)

      RETURN
      END SUBROUTINE MSAVE1C
!=======================================================================
      SUBROUTINE MSAVER(FN,A,NI,NJ)
!=======================================================================
! [USAGE]: 
! SAVE A REAL NI X NJ MATRIX VALUE INTO FN
! [PARAMETERS]:
! A >> A REAL MATRIX
! FN >> FILENAME TO STORE A VALUE
! NI >> # OF ROWS OF INPUT A
! NJ >> # OF COLUMNS OF INPUT A
! [DEPENDENCIES]:
! 1. LABOPEN(~) @ MOD_MISC
! 2. LABPR(~)   @ MOD_MISC
! 3. LABCR      @ MOD_MISC
! 4. LABCLOSE   @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      INTEGER :: NI,NJ
      REAL(P8),DIMENSION(NI,NJ):: A
      CHARACTER(LEN=*),INTENT(IN):: FN

      INTEGER :: I,J

      CALL LABOPEN(FN,NJ)

      DO I=1,NI
        DO J=1,NJ
          CALL LABPR(A(I,J))
        ENDDO
        CALL LABCR
      ENDDO

      CALL LABCLOSE

      WRITE(*,60) NI,NJ,TRIM(FN)

 60   FORMAT(I3,'X',I3,' WRITTEN TO ',A)

      RETURN
      END SUBROUTINE MSAVER
!=======================================================================
      SUBROUTINE MSAVEC(FN,A,NI,NJ)
!=======================================================================
! [USAGE]: 
! SAVE A COMPLEX NI X NJ MATRIX VALUE INTO FN
! [PARAMETERS]:
! A >> A COMPLEX MATRIX
! FN >> FILENAME TO STORE A VALUE
! NI >> # OF ROWS OF INPUT A
! NJ >> # OF COLUMNS OF INPUT A
! [DEPENDENCIES]:
! 1. LABOPEN(~) @ MOD_MISC
! 2. LABPR(~)   @ MOD_MISC
! 3. LABCR      @ MOD_MISC
! 4. LABCLOSE   @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      INTEGER :: NI,NJ
      COMPLEX(P8),DIMENSION(NI,NJ):: A
      CHARACTER(LEN=*),INTENT(IN):: FN

      INTEGER :: I,J

      CALL LABOPEN(FN,NJ*2)

      DO I=1,NI
        DO J=1,NJ
          CALL LABPR(REAL(A(I,J)))
          CALL LABPR(AIMAG(A(I,J)))
        ENDDO
        CALL LABCR
      ENDDO

      CALL LABCLOSE

      WRITE(*,60) NI,NJ,TRIM(FN)

 60   FORMAT(I3,'X',I3,'C WRITTEN TO ',A)

      RETURN
      END SUBROUTINE MSAVEC
!=======================================================================
      SUBROUTINE MSAVE3DIMENSIONALDATACOMPLEX(A,FN)
!=======================================================================
! [USAGE]: 
! SAVE A 3-DIMENSIONAL COMPLEX MATRIX VALUE INTO FN
! [PARAMETERS]:
! A >> A COMPLEX MATRIX
! FN >> FILENAME TO STORE A VALUE
! [DEPENDENCIES]:
! 1. MSAVEDO3D(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8),DIMENSION(:,:,:),INTENT(IN) ::A
      CHARACTER(LEN=*),INTENT(IN):: FN
      
      CALL MSAVEDO3D(A,FN,SIZE(A,1),SIZE(A,2),SIZE(A,3))

      RETURN       
      END SUBROUTINE MSAVE3DIMENSIONALDATACOMPLEX
!=======================================================================
      SUBROUTINE MSAVEDO3D(INPUT3D, FILENAME, DIM1, DIM2, DIM3)
!=======================================================================
! [USAGE]: 
! SAVE A 3-DIMENSIONAL COMPLEX MATRIX VALUE INTO FN
! [PARAMETERS]:
! INPUT3D >> A COMPLEX MATRIX
! FILENAME >> FILENAME TO STORE A VALUE
! DIM1 >> 1ST DIMENSION OF INPUT3D
! DIM2 >> 2ND DIMENSION OF INPUT3D
! DIM3 >> 3RD DIMENSION OF INPUT3D
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      ! INTEGER, INTENT(IN):: DIM1, DIM2, DIM3
      ! COMPLEX(P8),DIMENSION(DIM1,DIM2,DIM3),INTENT(IN):: INPUT3D
      ! CHARACTER(LEN=*),INTENT(IN):: FILENAME

      ! INTEGER :: I,J,K

      INTEGER, INTENT(IN):: DIM1, DIM2, DIM3
      COMPLEX(P8),DIMENSION(:,:,:),INTENT(IN):: INPUT3D
      CHARACTER(LEN=*),INTENT(IN):: FILENAME

      INTEGER :: I,J,K

      CALL LABOPEN(FILENAME,DIM2*2)

      DO K=1,DIM3
        DO I=1,DIM1
          DO J=1,DIM2
            CALL LABPR(REAL(INPUT3D(I,J,K)))
            CALL LABPR(AIMAG(INPUT3D(I,J,K)))
          ENDDO
          CALL LABCR
        ENDDO
      ENDDO

      CALL LABCLOSE

      WRITE(*,60) DIM1,DIM2,DIM3,TRIM(FILENAME)

 60   FORMAT(I3,'X',I3,'X',I3,'C WRITTEN TO ',A)

      RETURN
      END SUBROUTINE MSAVEDO3D
!=======================================================================
      SUBROUTINE MLOAD1R(FN,A)
!=======================================================================
! [USAGE]: 
! LOAD A REAL VECTOR VALUE FROM A STORAGE FN
! [PARAMETERS]:
! FN >> FILENAME STORING THE VALUE TO BE LOADED
! A >> ON EXIT, LOADED WITH A VECTOR VALUE STORED IN FN
! [DEPENDENCIES]:
! 1. MLOADR(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN):: FN
      REAL(P8),DIMENSION(:) :: A
      REAL(P8),DIMENSION(SIZE(A),1) :: AA
      INTEGER :: NI,NJ

      NI=SIZE(A,1)
      NJ=1
      CALL MLOADR(FN,AA,NI,NJ)
      A = AA(:,1)

      RETURN
      END SUBROUTINE MLOAD1R
!=======================================================================
      SUBROUTINE MLOAD1C(FN,A)
!=======================================================================
! [USAGE]: 
! LOAD A COPMLEX VECTOR VALUE FROM A STORAGE FN
! [PARAMETERS]:
! FN >> FILENAME STORING THE VALUE TO BE LOADED
! A >> ON EXIT, LOADED WITH A VECTOR VALUE STORED IN FN
! [DEPENDENCIES]:
! 1. MLOADC(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN):: FN
      COMPLEX(P8),DIMENSION(:) :: A
      COMPLEX(P8),DIMENSION(SIZE(A),1) :: AA
      INTEGER :: NI,NJ

      NI=SIZE(A,1)
      NJ=1
      CALL MLOADC(FN,AA,NI,NJ)
      A = AA(:,1)

      RETURN
      END SUBROUTINE MLOAD1C
!=======================================================================
      SUBROUTINE MLOADDR(FN,A)
!=======================================================================
! [USAGE]: 
! LOAD A REAL MATRIX VALUE FROM A STORAGE FN
! [PARAMETERS]:
! FN >> FILENAME STORING THE VALUE TO BE LOADED
! A >> ON EXIT, LOADED WITH A VECTOR VALUE STORED IN FN
! [DEPENDENCIES]:
! 1. MLOADR(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN):: FN
      REAL(P8),DIMENSION(:,:) :: A
      INTEGER :: NI,NJ

      NI=SIZE(A,1)
      NJ=SIZE(A,2)
      CALL MLOADR(FN,A,NI,NJ)

      RETURN
      END SUBROUTINE MLOADDR
!=======================================================================
      SUBROUTINE MLOADDC(FN,A)
!=======================================================================
! [USAGE]: 
! LOAD A COMPLEX MATRIX VALUE FROM A STORAGE FN
! [PARAMETERS]:
! FN >> FILENAME STORING THE VALUE TO BE LOADED
! A >> ON EXIT, LOADED WITH A VECTOR VALUE STORED IN FN
! [DEPENDENCIES]:
! 1. MLOADC(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN):: FN
      COMPLEX(P8),DIMENSION(:,:) :: A
      INTEGER :: NI,NJ

      NI=SIZE(A,1)
      NJ=SIZE(A,2)
      CALL MLOADC(FN,A,NI,NJ)

      RETURN
      END SUBROUTINE MLOADDC
!=======================================================================
      SUBROUTINE MLOADR0(FN,A,NI,NJ)
!=======================================================================
! [USAGE]: 
! LOAD A REAL MATRIX VALUE FROM A STORAGE FN
! [PARAMETERS]:
! FN >> FILENAME STORING THE VALUE TO BE LOADED
! A >> ON EXIT, LOADED WITH A REAL MATRIX VALUE STORED IN FN
! NI >> ON ENTRY, DUMMY VARIABLE
!       ON EXIT, STORES # OF LOADED ROWS 
! NJ >> ON ENTRY, # OF VALUES TO BE READ IN EACH LINE
!       ON EXIT, STORES # OF LOADED COLUMNS 
! [DEPENDENCIES]:
! 1. MLOADR(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN):: FN
      REAL(P8),DIMENSION(:,:) :: A
      INTEGER,INTENT(OUT) :: NI,NJ

      CALL MLOADR(FN,A,NI,NJ)

      RETURN
      END SUBROUTINE MLOADR0
!=======================================================================
      SUBROUTINE MLOADC0(FN,A,NI,NJ)
!=======================================================================
! [USAGE]: 
! LOAD A COPMLEX MATRIX VALUE FROM A STORAGE FN
! [PARAMETERS]:
! FN >> FILENAME STORING THE VALUE TO BE LOADED
! A >> ON EXIT, LOADED WITH A COMPLEX MATRIX VALUE STORED IN FN
! NI >> ON ENTRY, DUMMY VARIABLE
!       ON EXIT, STORES # OF LOADED ROWS 
! NJ >> ON ENTRY, # OF VALUES TO BE READ IN EACH LINE
!       ON EXIT, STORES # OF LOADED COLUMNS 
! [DEPENDENCIES]:
! 1. MLOADC(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN):: FN
      COMPLEX(P8),DIMENSION(:,:) :: A
      INTEGER,INTENT(OUT) :: NI,NJ

      CALL MLOADC(FN,A,NI,NJ)

      RETURN
      END SUBROUTINE MLOADC0
!=======================================================================
      SUBROUTINE MLOADR(FN,A,NI,NJ)
!=======================================================================
! [USAGE]: 
! LOAD A REAL MATRIX VALUE FROM A STORAGE FN
! [PARAMETERS]:
! FN >> FILENAME STORING THE VALUE TO BE LOADED
! A >> ON EXIT, LOADED WITH A REAL MATRIX VALUE STORED IN FN
! NI >> ON ENTRY, DUMMY VARIABLE
!       ON EXIT, STORES # OF LOADED ROWS 
! NJ >> ON ENTRY, # OF VALUES TO BE READ IN EACH LINE
!       ON EXIT, STORES # OF LOADED COLUMNS 
! [DEPENDENCIES]:
! 1. LABOPEN(~) @ MOD_MISC
! 2. LABRD(~)   @ MOD_MISC
! 3. LABCLOSE   @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      INTEGER,INTENT(OUT) :: NI,NJ
      REAL(P8),DIMENSION(:,:):: A
      CHARACTER(LEN=*),INTENT(IN):: FN
      INTEGER :: I,J,IS
      REAL(P8) :: VAL
      
      CALL LABOPEN(FN,NJ)
      I=1
      J=0
      NJ=-1

      DO WHILE(1)
        CALL LABRD(VAL,IS)
        ! VALUE READ NORMALLY
        IF(IS.EQ.0) THEN
          J=J+1
          IF(J > SIZE(A,2).OR. I>SIZE(A,1)) THEN
            WRITE (*,*) 'MLOADR:FILE SIZE TOO LARGE.'
            WRITE (*,*) 'MEM:',SIZE(A,1),'X',SIZE(A,2)
            WRITE (*,*) 'FIL:',I,'X',J
            STOP
          ENDIF
          A(I,J)=VAL
          CYCLE
        ! VALUE READ BUT EOR OCCURED
        ELSE IF(IS.EQ.-1) THEN
          ! >> FIRST LINE ?
          IF(NJ.EQ.-1) THEN
            J=J+1
            IF(J > SIZE(A,2).OR. I>SIZE(A,1)) THEN
              WRITE (*,*) 'MLOADR:FILE SIZE TOO LARGE.'
              WRITE (*,*) 'MEM:',SIZE(A,1),'X',SIZE(A,2)
              WRITE (*,*) 'FIL:',I,'X',J
             STOP
            ENDIF
            A(I,J)=VAL
            I=I+1
            NJ=J
            J=0
            CYCLE
          ELSE IF(J.EQ.NJ-1) THEN
            J=J+1
            IF(J > SIZE(A,2).OR. I>SIZE(A,1)) THEN
              WRITE (*,*) 'MLOADR:FILE SIZE TOO LARGE.'
              WRITE (*,*) 'MEM:',SIZE(A,1),'X',SIZE(A,2)
              WRITE (*,*) 'FIL:',I,'X',J
              STOP
            ENDIF
            A(I,J)=VAL
            I=I+1
            J=0
            CYCLE
          ELSE
            WRITE (*,*) 'MLOADR: EOF AT STRANGE PLACE.'
            NI=I-1
            EXIT
          ENDIF
        ! NO VALUE READ AND EOR OCCURED
        ELSE IF(IS.EQ.-2) THEN
          IF(J.NE.0.AND.NJ.EQ.-1) THEN
            NJ=J
            J=0
            I=I+1
            CYCLE
          ELSE IF(J.EQ.NJ) THEN
            I=I+1
            J=0
            CYCLE
          ELSE IF(J.EQ.0) THEN
            NI=I-1
            EXIT
            ELSE
            WRITE (*,*) 'MLOAD: EOF AT STRANGE PLACE.'
            NI=I-1
            EXIT
          ENDIF
        ! UNEXPECTED CHARACTER READ
        ELSE
          WRITE (*,*) 'MLOAD: STRANGE CHARACTER ?'
          NI=0
          NJ=0
          EXIT
        ENDIF
      ENDDO

      CALL LABCLOSE
      WRITE(*,60) NI,NJ,TRIM(FN)

 60   FORMAT(I3,'X',I3,' READ FROM ',A)

      RETURN
      END SUBROUTINE MLOADR
!=======================================================================
      SUBROUTINE MLOADC(FN,A,NI,NJ)
!=======================================================================
! [USAGE]: 
! LOAD A COPMLEX MATRIX VALUE FROM A STORAGE FN
! [PARAMETERS]:
! FN >> FILENAME STORING THE VALUE TO BE LOADED
! A >> ON EXIT, LOADED WITH A COMPLEX MATRIX VALUE STORED IN FN
! NI >> ON ENTRY, DUMMY VARIABLE
!       ON EXIT, STORES # OF LOADED ROWS 
! NJ >> ON ENTRY, # OF VALUES TO BE READ IN EACH LINE
!       ON EXIT, STORES # OF LOADED COLUMNS 
! [DEPENDENCIES]:
! 1. LABOPEN(~) @ MOD_MISC
! 2. LABRD(~)   @ MOD_MISC
! 3. LABCLOSE   @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      INTEGER,INTENT(OUT) :: NI,NJ
      COMPLEX(P8),DIMENSION(:,:):: A
      CHARACTER(LEN=*),INTENT(IN):: FN
      INTEGER :: I,J,IS
      REAL(P8) :: VAL1,VAL2

      CALL LABOPEN(FN,NJ)
      I=1
      J=0
      NJ=-1

      DO WHILE(1)
        CALL LABRD(VAL1,IS)
        IF(IS.EQ.0) THEN
          CALL LABRD(VAL2,IS)
        ENDIF
        ! VALUE READ NORMALLY
        IF(IS.EQ.0) THEN
          J=J+1
          IF(J>SIZE(A,2).OR. I>SIZE(A,1)) THEN
            WRITE (*,*) 'MLOADC:FILE SIZE TOO LARGE.'
            WRITE (*,*) 'MEM:',SIZE(A,1),'X',SIZE(A,2)
            WRITE (*,*) 'FIL:',I,'X',J
            STOP
          ENDIF
          A(I,J)=CMPLX(VAL1,VAL2,P8)
          CYCLE
        ! VALUE READ BUT EOR OCCURED
        ELSE IF(IS.EQ.-1) THEN
          ! >> FIRST LINE ?
          IF(NJ.EQ.-1) THEN
            J=J+1
            IF(J>SIZE(A,2).OR. I>SIZE(A,1)) THEN
              WRITE (*,*) 'MLOADC:FILE SIZE TOO LARGE.'
              WRITE (*,*) 'MEM:',SIZE(A,1),'X',SIZE(A,2)
              WRITE (*,*) 'FIL:',I,'X',J
              STOP
            ENDIF
            A(I,J)=CMPLX(VAL1,VAL2,P8)
            I=I+1
            NJ=J
            J=0
            CYCLE
          ELSE IF(J.EQ.NJ-1) THEN
            J=J+1
            IF(J>SIZE(A,2).OR. I>SIZE(A,1)) THEN
              WRITE (*,*) 'MLOADC:FILE SIZE TOO LARGE.'
              WRITE (*,*) 'MEM:',SIZE(A,1),'X',SIZE(A,2)
              WRITE (*,*) 'FIL:',I,'X',J
              STOP
            ENDIF
            A(I,J)=CMPLX(VAL1,VAL2,P8)
            I=I+1
            J=0
            CYCLE
          ELSE
            WRITE (*,*) 'MLOAD: EOF AT STRANGE PLACE.'
            NI=I-1
            EXIT
          ENDIF
        ! NO VALUE READ AND EOR OCCURED
        ELSE IF(IS.EQ.-2) THEN
          IF(J.NE.0.AND.NJ.EQ.-1) THEN
            NJ=J
            J=0
            I=I+1
            CYCLE
          ELSE IF(J.EQ.NJ) THEN
            I=I+1
            J=0
            CYCLE
          ELSE IF(J.EQ.0) THEN
            NI=I-1
            EXIT
          ELSE
            WRITE (*,*) 'MLOAD: EOF AT STRANGE PLACE.'
            NI=I-1
            EXIT
          ENDIF
        ! UNEXPECTED CHARACTER
        ELSE
          WRITE (*,*) 'MLOAD: STRANGE CHARACTER ?'
          NI=0
          NJ=0
          EXIT
        ENDIF
      ENDDO

      CALL LABCLOSE
      WRITE(*,60) NI,NJ,TRIM(FN)

      60 FORMAT(I3,'X',I3,'C READ FROM ',A)

      RETURN
      END SUBROUTINE MLOADC
!=======================================================================
      SUBROUTINE LABOPEN(FN,NC)
!=======================================================================
! [USAGE]: 
! OPEN A FILE FN A COPMLEX MATRIX VALUE FROM A STORAGE FN
! [PARAMETERS]:
! FN >> FILENAME STORING THE VALUE TO BE LOADED
! NC >> # OF VALUES TO BE CALLED IN EACH LINE
! [DEPENDENCIES]:
! 1. READCHAR(~) @ MOD_MISC
! 2. ATOF(~)   @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      CHARACTER(LEN=*) :: FN
      INTEGER :: NC

      INTEGER :: IDMY,RECLEN
      CHARACTER(LEN=1):: DMY
      INTEGER:: IS

      RECLEN = 25*(NC+1)
      OPEN(UNIT=11,FILE=FN,STATUS='UNKNOWN',&
           IOSTAT=IS,RECL=RECLEN)

      IF(IS.NE.0) THEN
        WRITE (*,*) 'LABOPEN: CANNOT OPEN',FN
      ENDIF

      ! INITIALIZE
      CALL READCHAR(DMY,IDMY,1)

      RETURN
      END SUBROUTINE LABOPEN
!=======================================================================
      SUBROUTINE LABPR(A)
!=======================================================================
! [USAGE]: 
! WRITE A NUMERIC DATA A INTO THE FILE ASSIGNED TO FILE NO. 11
! [PARAMETERS]:
! A >> A REAL VALUE 
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),INTENT(IN) :: A

      WRITE(UNIT=11,FMT=60,ADVANCE="NO") A
      ! IF YOU MODIFY THE FORMAT(25 COLUMN), MODIFY RECLEN IN LABOPEN

      60 FORMAT(1PE25.16E3)

      RETURN
      END SUBROUTINE LABPR
!=======================================================================
      SUBROUTINE LABRD(A,IS)
!=======================================================================
! [USAGE]: 
! READ A NUMERIC DATA A FROM THE FILE ASSIGNED TO FILE NO. 11
! [PARAMETERS]:
! A >> A REAL VALUE READ FROM THE FILE (11) 
! IS >> IS= 0 : NUMBER READ, NORMAL
!       IS=-1 : NUMBER READ, END OF RECORD/FILE
!       IS=-2 : NUMBER NOT READ, END OF RECORD/FILE
!       IS=-3 : UNEXPECTED CHARACTER ENCOUNTER
! [DEPENDENCIES]:
! 1. READCHAR(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),INTENT(OUT) :: A
      INTEGER :: IS

      CHARACTER(LEN=1) :: C
      CHARACTER(LEN=72) :: STR
      INTEGER :: COUNT

      COUNT=0
      STR=" "

      DO
        CALL READCHAR(C,IS)
        60   FORMAT(A1)
        IF(IS.NE.0) THEN
          A=0.0
          IS=-2
          RETURN
        ELSE IF(INDEX(" , ",C).NE.0) THEN
          CYCLE
        ELSE
          EXIT
        ENDIF
      ENDDO

      COUNT=1
      STR(COUNT:COUNT)=C

      DO
        CALL READCHAR(C,IS)
        IF(IS.NE.0) THEN
          A=ATOF(STR)
          IS=-1
          RETURN
        ELSE IF(INDEX("0123456789.EEDD+-",C).NE.0) THEN
          COUNT=COUNT+1
          STR(COUNT:COUNT)=C
        ELSE IF(INDEX(" , ",C).NE.0) THEN
          A=ATOF(STR)
          IS=0
          RETURN
        ELSE
          WRITE (*,*) 'LABRD: UNEXPECTED CHARACTER : ',C
          IS=-3
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE LABRD
!=======================================================================
      SUBROUTINE READCHAR(C,IS,INIT)
!=======================================================================
! [USAGE]: 
! READ A CHARACTER C
! [PARAMETERS]:
! A >> A REAL VALUE READ FROM THE FILE (11) 
! IS >> NUMBER READ STATUS
! INIT >> IF PRESENT, INITIALIZE THE PARAMETERS (II, IE, ISSAV)
! [DEPENDENCIES]:
! 1. READCHAR(~) @ MOD_MISC
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      CHARACTER(LEN=1) :: C
      INTEGER :: IS
      INTEGER,INTENT(IN),OPTIONAL:: INIT

      IF(PRESENT(INIT)) THEN
        II=0
        IE=0
        ISSAV=0
        BUF=" "
        RETURN
      ENDIF
      
      ! OUTPUT BUFFERED CHARACTER IF NOT EOR
      IF(II.LT.IE) THEN
        II=II+1
        C=BUF(II:II)
        IS=0
        RETURN
      ENDIF
      
      ! END OF RECORD ?
      IF(ISSAV.NE.0)THEN
        IS=ISSAV
        C=' '
        II=0
        IE=0
        ISSAV=0
        RETURN
      ENDIF
      
      ! READ BUFFER
      II=0
      BUF=" "
      READ(UNIT=11,FMT='(A256)',ADVANCE="NO",IOSTAT=ISSAV) BUF
      
      IF (ISSAV.EQ.0) THEN
        IE=256
        II=1
        C=BUF(II:II)
        IS=0
        RETURN
      ELSE  
        IE=LEN(TRIM(BUF))
        IF(IE.EQ.0) THEN
          IS=ISSAV
          C=' '
          II=0
          ISSAV=0
          RETURN
        ENDIF
        II=1
        C=BUF(II:II)
        IS=0
        RETURN
      ENDIF

      RETURN
      END SUBROUTINE READCHAR
!=======================================================================
      SUBROUTINE LABCR()
!=======================================================================
! [USAGE]: 
! WRITE CARRIAGE RETURN INTO FILE NO. 11
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE

      WRITE(11,60) ' '

   60 FORMAT(A1)

      RETURN
      END SUBROUTINE LABCR
!=======================================================================
      SUBROUTINE LABCLOSE()
!=======================================================================
! [USAGE]: 
! CLOSE A FILE ASSIGNED TO FILE NO. 11
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE

      CLOSE(11)

      RETURN
      END SUBROUTINE LABCLOSE
!=======================================================================
SUBROUTINE SORTR(DIR,INPUT)
!=======================================================================
      IMPLICIT NONE
      INTEGER:: DIR,N
      REAL(P8),DIMENSION(:):: INPUT
      REAL(P8),DIMENSION(:),ALLOCATABLE:: INDEX

      N = SIZE(INPUT)
      ALLOCATE(INDEX(N))
      CALL DLASRT(DIR,INPUT,N,INDEX)
      DEALLOCATE(INDEX)

      END SUBROUTINE SORTR
!=======================================================================
SUBROUTINE SORTR_I(DIR,INPUT,INDEX)
!=======================================================================
      IMPLICIT NONE
      INTEGER:: DIR,N
      REAL(P8),DIMENSION(:):: INPUT,INDEX
      REAL(P8),DIMENSION(:),ALLOCATABLE:: OUTPUT
      
      IF(SIZE(INDEX,1).NE.SIZE(INPUT,1)) THEN
        WRITE(*,*) 'SORTR_I: INDEX DIM DOES NOT MATCH INPUT DIM'
        RETURN
      ENDIF

      N = SIZE(INPUT)
      ALLOCATE(OUTPUT(N))
      OUTPUT = INPUT
      CALL DLASRT(DIR,OUTPUT,N,INDEX)
      DEALLOCATE(OUTPUT)

      END SUBROUTINE SORTR_I
!=======================================================================
SUBROUTINE SORTR_IO(DIR,INPUT,INDEX,OUTPUT)
!=======================================================================
      IMPLICIT NONE
      INTEGER:: DIR,N
      REAL(P8),DIMENSION(:):: INPUT,INDEX
      REAL(P8),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: OUTPUT

      IF(SIZE(INDEX,1).NE.SIZE(INPUT,1)) THEN
        WRITE(*,*) 'SORTR_IO: INDEX DIM DOES NOT MATCH INPUT DIM'
        RETURN
      ENDIF
      
      N = SIZE(INPUT)
      IF(.NOT.ALLOCATED(OUTPUT)) ALLOCATE(OUTPUT(N))
      OUTPUT = INPUT
      CALL DLASRT(DIR,OUTPUT,N,INDEX)

      END SUBROUTINE SORTR_IO
!=======================================================================
SUBROUTINE SORTR_IND(DIR,INPUT,INDEX,OUTPUT)
!=======================================================================
      IMPLICIT NONE
      INTEGER:: DIR
      REAL(P8),DIMENSION(:):: INPUT
      INTEGER,DIMENSION(:):: INDEX
      REAL(P8),DIMENSION(:),ALLOCATABLE,OPTIONAL:: OUTPUT
      REAL(P8),DIMENSION(:),ALLOCATABLE:: INDEX_R
      
      ALLOCATE(INDEX_R(SIZE(INDEX)))

      IF(PRESENT(OUTPUT)) THEN
        CALL SORTR_IO(DIR,INPUT,INDEX_R,OUTPUT)
      ELSE
        CALL SORTR_I(DIR,INPUT,INDEX_R)
      ENDIF

      INDEX = INT(INDEX_R)
      DEALLOCATE(INDEX_R)

      END SUBROUTINE SORTR_IND
!=======================================================================
SUBROUTINE SORTC(DIR,INPUT)
!=======================================================================
      IMPLICIT NONE
      INTEGER:: DIR,N
      COMPLEX(P8),DIMENSION(:):: INPUT
      REAL(P8),DIMENSION(:),ALLOCATABLE:: INDEX,INPUT_I

      N = SIZE(INPUT)
      ALLOCATE(INDEX(N))
      ALLOCATE(INPUT_I(N))
      INPUT_I = AIMAG(INPUT)
      CALL DLASRT(DIR,INPUT_I,N,INDEX) ! SORT BY THE IMAG PART
      INPUT = INPUT(INDEX)

      DEALLOCATE(INDEX)
      DEALLOCATE(INPUT_I)

      END SUBROUTINE SORTC
!=======================================================================
SUBROUTINE SORTC_I(DIR,INPUT,INDEX)
!=======================================================================
      IMPLICIT NONE
      INTEGER:: DIR,N
      COMPLEX(P8),DIMENSION(:):: INPUT
      REAL(P8),DIMENSION(:):: INDEX
      REAL(P8),DIMENSION(:),ALLOCATABLE:: OUTPUT_I

      IF(SIZE(INDEX,1).NE.SIZE(INPUT,1)) THEN
        WRITE(*,*) 'SORTC_I: INDEX DIM DOES NOT MATCH INPUT DIM'
        RETURN
      ENDIF

      N = SIZE(INPUT)
      ALLOCATE(OUTPUT_I(N))
      OUTPUT_I = AIMAG(INPUT)
      CALL DLASRT(DIR, OUTPUT_I,N,INDEX) ! SORT BY THE IMAG PART
      DEALLOCATE(OUTPUT_I)

      END SUBROUTINE SORTC_I
!=======================================================================
SUBROUTINE SORTC_IO(DIR,INPUT,INDEX,OUTPUT)
!=======================================================================
      IMPLICIT NONE
      INTEGER:: DIR,N
      COMPLEX(P8),DIMENSION(:):: INPUT
      COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: OUTPUT
      REAL(P8),DIMENSION(:):: INDEX
      REAL(P8),DIMENSION(:),ALLOCATABLE:: OUTPUT_I

      IF(SIZE(INDEX,1).NE.SIZE(INPUT,1)) THEN
        WRITE(*,*) 'SORTC_I: INDEX DIM DOES NOT MATCH INPUT DIM'
        RETURN
      ENDIF

      N = SIZE(INPUT)
      ALLOCATE(OUTPUT_I(N))
      OUTPUT_I = AIMAG(INPUT)
      CALL DLASRT(DIR, OUTPUT_I,N,INDEX) ! SORT BY THE IMAG PART
      DEALLOCATE(OUTPUT_I)

      IF(.NOT.ALLOCATED(OUTPUT)) ALLOCATE(OUTPUT(N))
      OUTPUT = INPUT(INDEX)    

      END SUBROUTINE SORTC_IO
!=======================================================================
SUBROUTINE SORTC_IND(DIR,INPUT,INDEX,OUTPUT)
!=======================================================================
      IMPLICIT NONE
      INTEGER:: DIR
      COMPLEX(P8),DIMENSION(:):: INPUT
      COMPLEX(P8),DIMENSION(:),ALLOCATABLE,OPTIONAL:: OUTPUT
      INTEGER,DIMENSION(:):: INDEX
      REAL(P8),DIMENSION(:),ALLOCATABLE:: INDEX_R

      ALLOCATE(INDEX_R(SIZE(INDEX)))

      IF(PRESENT(OUTPUT)) THEN
        CALL SORTC_IO(DIR,INPUT,INDEX_R,OUTPUT)
      ELSE
        CALL SORTC_I(DIR,INPUT,INDEX_R)
      ENDIF

      INDEX = INT(INDEX_R)
      DEALLOCATE(INDEX_R)

      END SUBROUTINE SORTC_IND
!=======================================================================
SUBROUTINE DLASRT(DIR,INPUT,N,INDEX)
!=======================================================================
! [USAGE]:
! SORT THE INPUT DATA; RETURN THE SORTED INPUT DATA AND THE INDEX
! [PARAMETERS]:
!    DIR: 0, Sort D in decreasing order
!    DIR: 1, Sort D in increasing order
! [UPDATES]
! -- LAPACK computational routine --
! -- LAPACK is a software package provided by Univ. of Tennessee,    --
! -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! RE-CODED BY JINGE WANG @ 4/30/2021
!=======================================================================
      IMPLICIT NONE
      INTEGER:: N, DIR
      REAL(P8),INTENT(INOUT)::  INPUT(N)
      REAL(P8),INTENT(OUT),OPTIONAL:: INDEX(N)

      INTEGER:: SELECT,INFO
      parameter( SELECT = 20 )
      INTEGER:: ENDD, I, J, START, STKPNT
      INTEGER:: STACK( 2, 32 )
      REAL(P8),DIMENSION(2):: D1, D2, D3, DMNMX, TMP
      REAL(P8):: D(N,2)

      ! Test the input parameters.
      info = 0
      IF( ( DIR.NE.0 ) .AND. (DIR.NE.1) ) THEN
          info = -1
          WRITE(*,*) 'DLASRT: DIR MUST BE EITHER O OR 1'
          RETURN
      END IF

      ! Combine DATA,INDEX into N
      INDEX = (/ (I, I = 1, N) /)
      D(:,1) = INPUT(:)
      D(:,2) = INDEX(:)

      stkpnt = 1
      stack( 1, 1 ) = 1
      stack( 2, 1 ) = N
      10 CONTINUE
      start = stack( 1, stkpnt )
      endd = stack( 2, stkpnt )
      stkpnt = stkpnt - 1

      ! For array size smaller than SELECT, do Insertion sort
      IF( endd-start.LE.SELECT .AND. endd-start.GT.0 ) THEN
      ! Do Insertion sort on D( START:ENDD )
          IF( dir.EQ.0 ) THEN
      ! Sort into decreasing order
            DO 30 i = start + 1, endd
                DO 20 j = i, start + 1, -1
                  IF( d( j,1 ).GT.d( j-1,1 ) ) THEN
                      dmnmx = d( j,: )
                      d( j,: ) = d( j-1,: )
                      d( j-1,: ) = dmnmx
                  ELSE
                      GO TO 30
                  END IF
      20          CONTINUE
      30       CONTINUE

          ELSE
      ! Sort into increasing order
            DO 50 i = start + 1, endd
                DO 40 j = i, start + 1, -1
                  IF( d( j,1 ).LT.d( j-1,1 ) ) THEN
                      dmnmx = d( j,: )
                      d( j,: ) = d( j-1,: )
                      d( j-1,: ) = dmnmx
                  ELSE
                      GO TO 50
                  END IF
      40          CONTINUE
      50       CONTINUE
          END IF

      ! For array size bigger than SELECT, do QUISKSORT
      ELSE IF( endd-start.GT.SELECT ) THEN
      ! Partition D( START:ENDD ) and stack parts, largest one first
      ! Choose partition entry as median of 3

          d1 = d( start,: )
          d2 = d( endd,: )
          i = ( start+endd ) / 2
          d3 = d( i,: )
          IF( d1(1).LT.d2(1) ) THEN
            IF( d3(1).LT.d1(1) ) THEN
                dmnmx = d1
            ELSE IF( d3(1).LT.d2(1) ) THEN
                dmnmx = d3
            ELSE
                dmnmx = d2
            END IF
          ELSE
            IF( d3(1).LT.d2(1) ) THEN
                dmnmx = d2
            ELSE IF( d3(1).LT.d1(1) ) THEN
                dmnmx = d3
            ELSE
                dmnmx = d1
            END IF
          END IF

          IF( dir.EQ.0 ) THEN

      ! Sort into decreasing order

            i = start - 1
            j = endd + 1
      60       CONTINUE
      70       CONTINUE
            j = j - 1
            IF( d( j,1 ).LT.dmnmx(1) ) GO TO 70
      80       CONTINUE
            i = i + 1
            IF( d( i,1 ).GT.dmnmx(1) ) GO TO 80
            IF( i.LT.j ) THEN
                tmp = d( i,: )
                d( i,: ) = d( j,: )
                d( j,: ) = tmp
                GO TO 60
            END IF
            IF( j-start.GT.endd-j-1 ) THEN
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = start
                stack( 2, stkpnt ) = j
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = j + 1
                stack( 2, stkpnt ) = endd
            ELSE
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = j + 1
                stack( 2, stkpnt ) = endd
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = start
                stack( 2, stkpnt ) = j
            END IF
          ELSE

      ! Sort into increasing order

            i = start - 1
            j = endd + 1
      90       CONTINUE
      100       CONTINUE
            j = j - 1
            IF( d( j,1 ).GT.dmnmx(1) ) GO TO 100
      110       CONTINUE
            i = i + 1
            IF( d( i,1 ).LT.dmnmx(1) ) GO TO 110
            IF( i.LT.j ) THEN
                tmp = d( i,: )
                d( i,: ) = d( j,: )
                d( j,: ) = tmp
                GO TO 90
            END IF
            IF( j-start.GT.endd-j-1 ) THEN
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = start
                stack( 2, stkpnt ) = j
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = j + 1
                stack( 2, stkpnt ) = endd
            ELSE
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = j + 1
                stack( 2, stkpnt ) = endd
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = start
                stack( 2, stkpnt ) = j
            END IF
          END IF
      END IF
      IF( stkpnt.GT.0 ) GO TO 10
      
      INPUT = D(:,1)
      INDEX = D(:,2)

      RETURN
      END SUBROUTINE DLASRT
!=======================================================================   
!============================ FUNCTIONS ================================
!=======================================================================
      FUNCTION ATOF(COMM)
!=======================================================================
! [USAGE]:
! CONVERT ASCII TO FLOATING NUMBER
! [INPUTS]:
! COMM >> ASCII STRING
! [OUTPUTS]:
! ATOF >> CONVERTED FLOATING NUMBER FROM COMM
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN):: COMM
      REAL(P8) :: ATOF

      CHARACTER(LEN=72) :: STR
      REAL(P8) :: VALUE

      IF(LEN(COMM)>72) THEN
        WRITE(*,*) 'ATOF: CHAR TOO LONG'
        STOP
      ENDIF

      STR=COMM(1:LEN(COMM))//REPEAT(" ",72)
      READ(UNIT=STR,FMT='(F72.0)',ERR=911) VALUE
      ATOF = VALUE

      RETURN

  911 CONTINUE
      WRITE(0,*) 'ATOF: ERROR IN CONVERSION'
      WRITE (*,*) 'COMM=',COMM
      ATOF = 0.D0

      RETURN
      END FUNCTION ATOF
!=======================================================================
      FUNCTION ITOA3(I)
!=======================================================================
! [USAGE]:
! CONVERT INTEGER TO CHARACTER(LEN=3)
! [INPUTS]:
! I >> INTEGER (0 - 999)
! [OUTPUTS]:
! ITOA3 >> CONVERTED CHARACTER FROM I
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      INTEGER:: I
      CHARACTER(LEN=3):: ITOA3

      WRITE(ITOA3,60) I

 60   FORMAT(I3)

      IF(ITOA3(1:1).EQ.' ') ITOA3(1:1)='0'
      IF(ITOA3(2:2).EQ.' ') ITOA3(2:2)='0'
      IF(ITOA3(3:3).EQ.' ') ITOA3(3:3)='0'

      IF(ITOA3(2:2).EQ.'-') THEN
        ITOA3(2:2)='0'
        ITOA3(1:1)='-'
      ENDIF

      RETURN
      END FUNCTION ITOA3
!=======================================================================
      FUNCTION ITOA4(I)
!=======================================================================
! [USAGE]:
! CONVERT INTEGER TO CHARACTER(LEN=4)
! [INPUTS]:
! I >> INTEGER (0 - 999)
! [OUTPUTS]:
! ITOA4 >> CONVERTED CHARACTER FROM I
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      INTEGER,INTENT(IN):: I
      CHARACTER(LEN=4):: ITOA4

      WRITE(ITOA4,60) I

  60   FORMAT(I4)

      IF(ITOA4(1:1).EQ.' ') ITOA4(1:1)='0'
      IF(ITOA4(2:2).EQ.' ') ITOA4(2:2)='0'
      IF(ITOA4(3:3).EQ.' ') ITOA4(3:3)='0'
      IF(ITOA4(4:4).EQ.' ') ITOA4(4:4)='0'

      RETURN
      END FUNCTION ITOA4
!=======================================================================
      FUNCTION RTOA6(I)
!=======================================================================
! [USAGE]:
! CONVERT INTEGER TO CHARACTER(LEN=4)
! [INPUTS]:
! I >> INTEGER (0 - 999)
! [OUTPUTS]:
! ITOA4 >> CONVERTED CHARACTER FROM I
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!=======================================================================
      REAL:: I
      CHARACTER(LEN=6):: RTOA6

      WRITE(RTOA6,'(F06.2)') I

      RETURN
      END FUNCTION RTOA6
!=======================================================================
SUBROUTINE LINSPACE(from,to,array)
!=======================================================================
! [USAGE]:
! GENERATES EVENLY SPACED NUMBERS FROM 'from' TO 'to'
! [INPUTS]:
! from(input) >> starting value
! to(input) >> end value
! array(output) >> dimension of arrays gives number of values       
! [OUTPUTS]:
! array >> evenly spaced numbers
! [UPDATES]:
! CODED BY JIGNE WANG @ MAR 10 2021         
!=======================================================================
      IMPLICIT NONE

      REAL(P8):: from, to
      REAL(P8),DIMENSION(:) :: array
      REAL(P8) :: RANGE
      INTEGER :: I,N

      N = SIZE(array)
      RANGE = to - from

      IF (N == 0) RETURN
      IF (N == 1) THEN
        array(1) = from
        RETURN
      END IF

      DO I = 1,N
        array(I) = from + RANGE*(I-1)/(N-1)
      END DO

    END SUBROUTINE LINSPACE
!=======================================================================
!=======================================================================
    FUNCTION ZEROS1(NI,NJ)
!=======================================================================
! [USAGE]:                                   
! CREATE AN NI X NJ ZERO MATRIX
! [INPUTS]:                                  
! NI >> # OF THE ROWS OF THE PRODUCED ZERO MATRIX
! NJ >> # OF THE COLUMNS OF THE PRODUCED ZERO MATRIX
! [OUTPUTS]:
! ZEROS1 >> NI X NJ ZERO MATRIX
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      INTEGER(P4), INTENT(IN)    :: NI,NJ
      REAL(P8), DIMENSION(NI,NJ) :: ZEROS1

      ZEROS1 = 0.D0

      RETURN
      END FUNCTION ZEROS1
!=======================================================================
      FUNCTION EYE1(NI) ! FAMILY OF EYE
!=======================================================================
! [USAGE]:                                   
! CREATE AN IDENTITY MATRIX OF ORDER NI
! [INPUTS]:                                  
! NI >> THE ORDER OF THE IDENTITY MATRIX TO BE CREATED
! [OUTPUTS]:
! EYE1 >> NI X NI IDENTITY MATRIX
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      INTEGER(P4),INTENT(IN) :: NI
      ! REAL(P8),DIMENSION(:,:),ALLOCATABLE:: EYE1
      REAL(P8),DIMENSION(NI,NI):: EYE1
      INTEGER(P4) :: I

      ! ALLOCATE(EYE1(NI,NI))
      EYE1=0.D0

      DO I=1,NI
        EYE1(I,I)=1.D0
      ENDDO

      RETURN
      END FUNCTION EYE1
!=======================================================================
      FUNCTION DIAG1R(A)
!=======================================================================
! [USAGE]:                                   
! EXTRACT THE DIAGONAL ELEMENTS OF AN INPUT MATRIX A
! [INPUTS]:                                  
! A >> M X N RECTANGULAR REAL MATRIX
! [OUTPUTS]:
! DIAG1R >> DIAGONAL ELEMENTS OF A. DIM=MIN(#COL OF A,#ROW OF A)
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),DIMENSION(:,:) :: A
      REAL(P8),DIMENSION(MIN(SIZE(A,1),SIZE(A,2))):: DIAG1R
      INTEGER(P4) :: I

      DO I=1,MIN(SIZE(A,1),SIZE(A,2))
        DIAG1R(I)=A(I,I)
      ENDDO

      RETURN
      END FUNCTION DIAG1R
!=======================================================================
      FUNCTION DIAG1C(A)
!=======================================================================
! [USAGE]:                                   
! EXTRACT THE DIAGONAL ELEMENTS OF AN INPUT MATRIX A
! [INPUTS]:                                  
! A >> M X N RECTANGULAR COMPLEX MATRIX
! [OUTPUTS]:
! DIAG1C >> DIAGONAL ELEMENTS OF A. DIM=MIN(#COL OF A,#ROW OF A)
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8),DIMENSION(:,:) :: A
      COMPLEX(P8),DIMENSION(MIN(SIZE(A,1),SIZE(A,2))):: DIAG1C
      INTEGER(P4) :: I

      DO I=1,MIN(SIZE(A,1),SIZE(A,2))
        DIAG1C(I)=A(I,I)
      ENDDO

      RETURN
      END FUNCTION DIAG1C
!=======================================================================
      FUNCTION DIAG2R(A)
!=======================================================================
! [USAGE]:                                   
! USING THE ELEMENTS OF A, CREATE THE DIAGONAL MATRIX WHERE (D_II) = A_I
! [INPUTS]:                                  
! A >> N X 1 REAL VECTOR
! [OUTPUTS]:
! DIAG2R >> DIAGONAL MATRIX WHERE ITS I-TH COMPONENT IS EQUAL TO A_I
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      REAL(P8),DIMENSION(:) :: A
      REAL(P8),DIMENSION(SIZE(A),SIZE(A)):: DIAG2R
      INTEGER(P4) :: I

      DIAG2R = 0.D0

      DO I = 1,SIZE(A)
        DIAG2R(I,I) = A(I)
      ENDDO

      RETURN
      END FUNCTION DIAG2R
!=======================================================================
      FUNCTION DIAG2C(A)
!=======================================================================
! [USAGE]:                                   
! USING THE ELEMENTS OF A, CREATE THE DIAGONAL MATRIX WHERE (D_II) = A_I
! [INPUTS]:                                  
! A >> N X 1 COMPLEX VECTOR
! [OUTPUTS]:
! DIAG2C >> DIAGONAL MATRIX WHERE ITS I-TH COMPONENT IS EQUAL TO A_I
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      COMPLEX(P8),DIMENSION(:) :: A
      COMPLEX(P8),DIMENSION(SIZE(A),SIZE(A)):: DIAG2C
      INTEGER(P4) :: I

      DIAG2C = 0.D0

      DO I = 1,SIZE(A)
        DIAG2C(I,I) = A(I)
      ENDDO

      RETURN
      END FUNCTION DIAG2C
! ======================================================================

END MODULE MOD_MISC