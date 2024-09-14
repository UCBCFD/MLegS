MODULE MOD_EIG ! LEVEL 1 MODULE
! FULL MATRIX OPERATORS AND EIGVAL RELATED OPERATIONS
      USE OMP_LIB
      USE MPI
      USE MOD_MISC, ONLY : P4,P8,PI,IU, EYE                              ! LEVEL 0
      IMPLICIT NONE
      PRIVATE
!=======================================================================
!======================== PUBLIC DECLARATION ===========================
!=======================================================================
      ! OPERATORS
      PUBLIC :: OPERATOR(.MUL.)                                          ! MATRIX/VECTOR MULTIPLICATION OPERATOR
      PUBLIC :: OPERATOR(.LR.)                                           ! LEFT-RIGHT MATRIX PASTE OPERATOR
      PUBLIC :: OPERATOR(.UD.)                                           ! UP-DOWN MATRIX PASTE OPERATOR
      ! INVERSE MATRIX CALCULATION
      PUBLIC :: INV
      ! MATRIX WHOSE COLUMNS ARE THE EIGVECS
      PUBLIC :: GENEIG
      ! GENERALIZED EIGENVALUE PROBLEM
      PUBLIC :: EIGVEC
      ! FLIP MATRIX 
      PUBLIC :: FLIPUD                                                    ! UPSIDE DOWN
      PUBLIC :: FLIPLR                                                    ! LEFT TO RIGHT
      ! LU DECOMPOSITION
      PUBLIC :: LU
      ! SOLVE AX = B VIA LU DECOMPOSITION OF A
      PUBLIC :: SOLVE
!=======================================================================
!============================ INTERFACES ===============================
!=======================================================================
      INTERFACE OPERATOR(.MUL.)
        MODULE PROCEDURE MULRR
        MODULE PROCEDURE MULCC
        MODULE PROCEDURE MULCR
        MODULE PROCEDURE MULRC
        MODULE PROCEDURE MULRR1
        MODULE PROCEDURE MULRC1
        MODULE PROCEDURE MULCR1
        MODULE PROCEDURE MULCC1
        MODULE PROCEDURE MULC1R
      END INTERFACE

      INTERFACE OPERATOR(.LR.)
        MODULE PROCEDURE PASTE_LR_RR
        MODULE PROCEDURE PASTE_LR_CC
      END INTERFACE

      INTERFACE OPERATOR(.UD.)
        MODULE PROCEDURE PASTE_UD_RR
        MODULE PROCEDURE PASTE_UD_CC
      END INTERFACE

      INTERFACE INV
        MODULE PROCEDURE INVDR
        MODULE PROCEDURE INVDC
      END INTERFACE

      INTERFACE GENEIG
        MODULE PROCEDURE GENEIG1
      END INTERFACE

      INTERFACE EIGVEC
        MODULE PROCEDURE EIGVEC1
        MODULE PROCEDURE EIGVEC2
      END INTERFACE

      INTERFACE FLIPUD
        MODULE PROCEDURE FLIPUDR
        MODULE PROCEDURE FLIPUDC
      END INTERFACE

      INTERFACE FLIPLR
        MODULE PROCEDURE FLIPLRR
        MODULE PROCEDURE FLIPLRC
      END INTERFACE

      INTERFACE LU
        MODULE PROCEDURE LU_REAL0
        MODULE PROCEDURE LU_COMPLEX0
      END INTERFACE

      INTERFACE SOLVE
        MODULE PROCEDURE SOLVE_REAL1
        MODULE PROCEDURE SOLVE_REAL2
        MODULE PROCEDURE SOLVE_REAL1C
        MODULE PROCEDURE SOLVE_REAL2C
        MODULE PROCEDURE SOLVE_COMPLEX1
        MODULE PROCEDURE SOLVE_COMPLEX2
      END INTERFACE
CONTAINS
!=======================================================================
!============================ SUBROUTINES ==============================
!=======================================================================
      SUBROUTINE LU_REAL0(A,IPIV)
!=======================================================================
! [USAGE]: 
! PERFORM THE LU FACTORIZATION WITH RESPECT TO AN INPUT M X N REAL
! MATRIX A. ON EXIT, A IS REPLACED WITH L+U-I
! [PARAMETERS]:
! A >> ON ENTRY, THE M X N MATRIX TO BE FACTORED
!      ON EXIT, THE FACTORS L AND U VIA FACTORIZATION (L+U-I)
! [DEPENDENCY]:
! 1. DGETRF(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 9 2020
! IMPLEMENTED DGETRF FOR LU FACTORIZATION @ NOV 9 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL DGETRF
      REAL(P8),DIMENSION(:,:),INTENT(INOUT):: A
      INTEGER(P4),DIMENSION(:),INTENT(OUT)::IPIV

      INTEGER(P4)::M,N,INFO

      M = SIZE(A,1)
      N = SIZE(A,2)

      IF(SIZE(IPIV) .NE. MIN(M,N)) THEN
        WRITE(*,*) 'LU_REAL0: INCONSISTENCY IN IPIV DIMENSION.'
        STOP
      ENDIF

      CALL DGETRF(M,N,A,M,IPIV,INFO)

      IF(INFO.NE.0)THEN
        WRITE(*,*) 'LU_REAL0: LU FACTORIZATION FAILURE.'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE LU_REAL0
!=======================================================================
      SUBROUTINE LU_COMPLEX0(A,IPIV)
!=======================================================================
! [USAGE]: 
! PERFORM THE LU FACTORIZATION WITH RESPECT TO AN INPUT M X N REAL
! MATRIX A. ON EXIT, A IS REPLACED WITH L+U-I
! [PARAMETERS]:
! A >> ON ENTRY, THE M X N MATRIX TO BE FACTORED
!      ON EXIT, THE FACTORS L AND U VIA FACTORIZATION (L+U-I)
! [DEPENDENCY]:
! 1. ZGETRF(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 9 2020
! IMPLEMENTED ZGETRF FOR LU FACTORIZATION @ NOV 9 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL ZGETRF
      COMPLEX(P8),DIMENSION(:,:),INTENT(INOUT):: A
      INTEGER(P4),DIMENSION(:),INTENT(OUT)::IPIV

      INTEGER(P4)::M,N,INFO

      M = SIZE(A,1)
      N = SIZE(A,2)

      IF(SIZE(IPIV) .NE. MIN(M,N)) THEN
        WRITE(*,*) 'LU_REAL0: INCONSISTENCY IN IPIV DIMENSION.'
        STOP
      ENDIF

      CALL ZGETRF(M,N,A,M,IPIV,INFO)

      IF(INFO.NE.0)THEN
        WRITE(*,*) 'LU_COMPLEX0: LU FACTORIZATION FAILURE.'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE LU_COMPLEX0
!=======================================================================
      SUBROUTINE SOLVE_REAL1(A,B)
!=======================================================================
! [USAGE]: 
! SOLVE THE SYSTEMS OF LINEAR EQUATIONS AX = B  
! MATRIX A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [PARAMETERS]:
! A >> N X N REAL LHS MATRIX AFTER LU FACTORZATION
!      IF NOT LU-FACTORIZED, THEN A WILL BE FACTORIZED ON EXIT
! B >> ON ENTRY, REAL N X 1 RHS VECTOR OF THE SYSTEM
!      ON EXIT, THE SOLUTION OF THE SYSTEM
! [DEPENDENCY]:
! 1. LU_REAL0(A) @ MOD_EIG
! 2. DGETRS(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 9 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL DGETRS
      REAL(P8),DIMENSION(:,:),INTENT(INOUT):: A
      REAL(P8),DIMENSION(:),INTENT(INOUT):: B

      INTEGER(P4),DIMENSION(SIZE(A,1)):: IPIV
      INTEGER(P4):: NI,NJ,INFO

      IF(SIZE(A,1).NE.SIZE(A,2)) THEN
        WRITE (*,*) 'SOLVE_REAL1: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF(SIZE(A,2).NE.SIZE(B)) THEN
        WRITE (*,*) 'SOLVE_REAL1: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI = SIZE(A,1)
      NJ = 1

      CALL LU_REAL0(A,IPIV)

      CALL DGETRS('N',NI,NJ,A,NI,IPIV,B,NI,INFO)

      IF(INFO.NE.0)THEN
        WRITE(*,*) 'SOLVE_REAL1: LINEAR SYSTEM SOLVING FAILURE.'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE SOLVE_REAL1
!=======================================================================
      SUBROUTINE SOLVE_REAL2(A,B)
!=======================================================================
! [USAGE]: 
! SOLVE THE SYSTEMS OF LINEAR EQUATIONS AX = B  
! MATRIX A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [PARAMETERS]:
! A >> N X N REAL LHS MATRIX AFTER LU FACTORZATION
!      IF NOT LU-FACTORIZED, THEN A WILL BE FACTORIZED ON EXIT
! B >> ON ENTRY, REAL N X K RHS VECTOR OF THE SYSTEM
!      ON EXIT, THE SOLUTION OF THE SYSTEM
! [DEPENDENCY]:
! 1. LU_REAL0(A) @ MOD_EIG
! 2. DGETRS(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 9 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL DGETRS
      REAL(P8),DIMENSION(:,:),INTENT(INOUT):: A
      REAL(P8),DIMENSION(:,:),INTENT(INOUT):: B

      INTEGER(P4),DIMENSION(SIZE(A,1)):: IPIV
      INTEGER(P4):: NI,NJ,INFO

      IF(SIZE(A,1).NE.SIZE(A,2)) THEN
        WRITE (*,*) 'SOLVE_REAL2: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF(SIZE(A,2).NE.SIZE(B,1)) THEN
        WRITE (*,*) 'SOLVE_REAL2: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI = SIZE(A,1)
      NJ = SIZE(B,2)

      CALL LU_REAL0(A,IPIV)

      CALL DGETRS('N',NI,NJ,A,NI,IPIV,B,NI,INFO)

      IF(INFO.NE.0)THEN
        WRITE(*,*) 'SOLVE_REAL2: LINEAR SYSTEM SOLVING FAILURE.'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE SOLVE_REAL2
!=======================================================================
      SUBROUTINE SOLVE_REAL1C(A,B)
!=======================================================================
! [USAGE]: 
! SOLVE THE SYSTEMS OF LINEAR EQUATIONS AX = B  
! MATRIX A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [PARAMETERS]:
! A >> N X N REAL LHS MATRIX AFTER LU FACTORZATION
!      IF NOT LU-FACTORIZED, THEN A WILL BE FACTORIZED ON EXIT
! B >> ON ENTRY, COMPLEX N X 1 RHS VECTOR OF THE SYSTEM
!      ON EXIT, THE SOLUTION OF THE SYSTEM
! [DEPENDENCY]:
! 1. LU_REAL0(A) @ MOD_EIG
! 2. DGETRS(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 9 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL DGETRS
      REAL(P8),DIMENSION(:,:),INTENT(INOUT):: A
      COMPLEX(P8),DIMENSION(:),INTENT(INOUT):: B

      REAL(P8),DIMENSION(SIZE(B)):: WR,WI

      INTEGER(P4),DIMENSION(SIZE(A,1)):: IPIV
      INTEGER(P4):: NI,NJ,INFO

      IF(SIZE(A,1).NE.SIZE(A,2)) THEN
        WRITE (*,*) 'SOLVE_REAL1C: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF(SIZE(A,2).NE.SIZE(B)) THEN
        WRITE (*,*) 'SOLVE_REAL1C: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI = SIZE(A,1)
      NJ = 1

      WR = REAL(B)
      WI = AIMAG(B)

      CALL LU_REAL0(A,IPIV)

      CALL DGETRS('N',NI,NJ,A,NI,IPIV,WR,NI,INFO)

      IF(INFO.NE.0)THEN
        WRITE(*,*) 'SOLVE_REAL1C: LINEAR SYSTEM SOLVING FAILURE.'
        STOP
      ENDIF

      CALL DGETRS('N',NI,NJ,A,NI,IPIV,WI,NI,INFO)

      IF(INFO.NE.0)THEN
        WRITE(*,*) 'SOLVE_REAL1C: LINEAR SYSTEM SOLVING FAILURE.'
        STOP
      ENDIF

      B = CMPLX(WR,WI,P8)

      RETURN
      END SUBROUTINE SOLVE_REAL1C
!=======================================================================
      SUBROUTINE SOLVE_REAL2C(A,B)
!=======================================================================
! [USAGE]: 
! SOLVE THE SYSTEMS OF LINEAR EQUATIONS AX = B  
! MATRIX A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [PARAMETERS]:
! A >> N X N REAL LHS MATRIX AFTER LU FACTORZATION
!      IF NOT LU-FACTORIZED, THEN A WILL BE FACTORIZED ON EXIT
! B >> ON ENTRY, COMPLEX N X K RHS VECTOR OF THE SYSTEM
!      ON EXIT, THE SOLUTION OF THE SYSTEM
! [DEPENDENCY]:
! 1. LU_REAL0(A) @ MOD_EIG
! 2. DGETRS(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 9 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL DGETRS
      REAL(P8),DIMENSION(:,:),INTENT(INOUT):: A
      COMPLEX(P8),DIMENSION(:,:),INTENT(INOUT):: B

      REAL(P8),DIMENSION(SIZE(B,1),SIZE(B,2)):: WR,WI

      INTEGER(P4),DIMENSION(SIZE(A,1)):: IPIV
      INTEGER(P4):: NI,NJ,INFO

      IF(SIZE(A,1).NE.SIZE(A,2)) THEN
        WRITE (*,*) 'SOLVE_REAL2C: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF(SIZE(A,2).NE.SIZE(B,1)) THEN
        WRITE (*,*) 'SOLVE_REAL2C: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI = SIZE(A,1)
      NJ = SIZE(B,2)

      WR = REAL(B)
      WI = AIMAG(B)

      CALL LU_REAL0(A,IPIV)

      CALL DGETRS('N',NI,NJ,A,NI,IPIV,WR,NI,INFO)

      IF(INFO.NE.0)THEN
        WRITE(*,*) 'SOLVE_REAL2C: LINEAR SYSTEM SOLVING FAILURE.'
        STOP
      ENDIF

      CALL DGETRS('N',NI,NJ,A,NI,IPIV,WI,NI,INFO)

      IF(INFO.NE.0)THEN
        WRITE(*,*) 'SOLVE_REAL2C: LINEAR SYSTEM SOLVING FAILURE.'
        STOP
      ENDIF

      B = CMPLX(WR,WI,P8)

      RETURN
      END SUBROUTINE SOLVE_REAL2C
!=======================================================================
      SUBROUTINE SOLVE_COMPLEX1(A,B)
!=======================================================================
! [USAGE]: 
! SOLVE THE SYSTEMS OF LINEAR EQUATIONS AX = B  
! MATRIX A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [PARAMETERS]:
! A >> N X N COMPLEX LHS MATRIX AFTER LU FACTORZATION
!      IF NOT LU-FACTORIZED, THEN A WILL BE FACTORIZED ON EXIT
! B >> ON ENTRY, COMPLEX N X 1 RHS VECTOR OF THE SYSTEM
!      ON EXIT, THE SOLUTION OF THE SYSTEM
! [DEPENDENCY]:
! 1. LU_COMPLEX0(A) @ MOD_EIG
! 2. ZGETRS(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 9 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL ZGETRS
      COMPLEX(P8),DIMENSION(:,:),INTENT(INOUT):: A
      COMPLEX(P8),DIMENSION(:),INTENT(INOUT):: B

      INTEGER(P4),DIMENSION(SIZE(A,1)):: IPIV
      INTEGER(P4):: NI,NJ,INFO

      IF(SIZE(A,1).NE.SIZE(A,2)) THEN
        WRITE (*,*) 'SOLVE_COMPLEX1: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF(SIZE(A,2).NE.SIZE(B)) THEN
        WRITE (*,*) 'SOLVE_COMPLEX1: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI = SIZE(A,1)
      NJ = 1

      CALL LU_COMPLEX0(A,IPIV)

      CALL ZGETRS('N',NI,NJ,A,NI,IPIV,B,NI,INFO)

      IF(INFO.NE.0)THEN
        WRITE(*,*) 'SOLVE_COMPLEX1: LINEAR SYSTEM SOLVING FAILURE.'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE SOLVE_COMPLEX1
!=======================================================================
      SUBROUTINE SOLVE_COMPLEX2(A,B)
!=======================================================================
! [USAGE]: 
! SOLVE THE SYSTEMS OF LINEAR EQUATIONS AX = B  
! MATRIX A NEEDS TO BE LU-FACTORIZED IN ADVANCE AS L+U-I
! [PARAMETERS]:
! A >> N X N COMPLEX LHS MATRIX AFTER LU FACTORZATION
!      IF NOT LU-FACTORIZED, THEN A WILL BE FACTORIZED ON EXIT
! B >> ON ENTRY, COMPLEX N X K RHS VECTOR OF THE SYSTEM
!      ON EXIT, THE SOLUTION OF THE SYSTEM
! [DEPENDENCY]:
! 1. LU_COMPLEX0(A) @ MOD_EIG
! 2. ZGETRS(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 9 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL ZGETRS
      COMPLEX(P8),DIMENSION(:,:),INTENT(INOUT):: A
      COMPLEX(P8),DIMENSION(:,:),INTENT(INOUT):: B

      INTEGER(P4),DIMENSION(SIZE(A,1)):: IPIV
      INTEGER(P4):: NI,NJ,INFO

      IF(SIZE(A,1).NE.SIZE(A,2)) THEN
        WRITE (*,*) 'SOLVE_COMPLEX2: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF(SIZE(A,2).NE.SIZE(B,1)) THEN
        WRITE (*,*) 'SOLVE_COMPLEX2: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI = SIZE(A,1)
      NJ = SIZE(B,2)

      CALL LU_COMPLEX0(A,IPIV)

      CALL ZGETRS('N',NI,NJ,A,NI,IPIV,B,NI,INFO)

      IF(INFO.NE.0)THEN
        WRITE(*,*) 'SOLVE_COMPLEX2: LINEAR SYSTEM SOLVING FAILURE.'
        STOP
      ENDIF

      RETURN
      END SUBROUTINE SOLVE_COMPLEX2
!=======================================================================
!============================ FUNCTIONS ================================
!=======================================================================
      FUNCTION MULRR(A,B) ! FAMILY OF OPERATOR(.MUL.)
!=======================================================================
! [USAGE]:
! REAL-REAL MATRIX MULTIPLICATION
! [INPUTS]:
! A >> M X P REAL MATRIX
! B >> P X N REAL MATRIX
! [OUTPUTS]:
! MULRR >> M X N REAL MATRIX WHICH IS EQUAL TO AB
! [DEPENDENCIES]:
! 1. DGEMM(~) FUNCTION @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
! CHG SGEMM (SINGLE PRECISION) TO DGEMM (DOUBLE PRECISION) @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL DGEMM ! FORCED TO USE DGEMM IN LAPACK
      REAL(P8), DIMENSION(:,:), INTENT(IN)    :: A,B
      REAL(P8), DIMENSION(SIZE(A,1),SIZE(B,2)):: MULRR
      
      INTEGER(P4)         :: I,J,K
      REAL(P8), PARAMETER :: ZERO = 0.D0
      REAL(P8), PARAMETER :: ONE  = 1.D0

      IF (SIZE(A,2).NE.SIZE(B,1)) THEN
        WRITE (*,*) 'MUL: INVALID SHAPE.'
        WRITE (*,*) 'ARGUMENT1 = ',SIZE(A,1),'X',SIZE(A,2)
        WRITE (*,*) 'ARGUMENT2 = ',SIZE(B,1),'X',SIZE(B,2)
        STOP
      ENDIF

      I=SIZE(A,1)
      J=SIZE(B,2)
      K=SIZE(A,2)
      !WRITE(6,*) I,J,K !For debug
      CALL DGEMM('N','N',I,J,K,ONE,A,I,B,K,ZERO,MULRR,I)

      RETURN
      END FUNCTION MULRR
!=======================================================================
      FUNCTION MULCC(A,B) ! FAMILY OF OPERATOR(.MUL.)
!=======================================================================
! [USAGE]:
! COMPLEX-COMPLEX MATRIX MULTIPLICATION
! [INPUTS]:
! A >> M X P COMPLEX MATRIX
! B >> P X N COMPLEX MATRIX
! [OUTPUTS]:
! MULRR >> M X N COMPLEX MATRIX WHICH IS EQUAL TO AB
! [DEPENDENCIES]:
! 1. ZGEMM(~) FUNCTION @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
! CHG CGEMM (SINGLE PRECISION) TO ZGEMM (DOUBLE PRECISION) @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL ZGEMM ! FORCED TO USE ZGEMM IN LAPACK
      COMPLEX(P8), DIMENSION(:,:), INTENT(IN)    :: A,B
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(B,2)):: MULCC
      
      INTEGER(P4)            :: I,J,K
      COMPLEX(P8), PARAMETER :: ZERO = 0.D0
      COMPLEX(P8), PARAMETER :: ONE  = 1.D0

      IF (SIZE(A,2).NE.SIZE(B,1)) THEN
        WRITE (*,*) 'MUL: INVALID SHAPE.'
        WRITE (*,*) 'ARGUMENT1 = ',SIZE(A,1),'X',SIZE(A,2)
        WRITE (*,*) 'ARGUMENT2 = ',SIZE(B,1),'X',SIZE(B,2)
        STOP
      ENDIF

      I=SIZE(A,1)
      J=SIZE(B,2)
      K=SIZE(A,2)
      CALL ZGEMM('N','N',I,J,K,ONE,A,I,B,K,ZERO,MULCC,I)
      
      RETURN
      END FUNCTION MULCC
!=======================================================================
      FUNCTION MULCR(A,B) ! FAMILY OF OPERATOR(.MUL.)
!=======================================================================
! [USAGE]:
! COMPLEX-REAL MATRIX MULTIPLICATION
! [INPUTS]:
! A >> M X P COMPLEX MATRIX
! B >> P X N REAL MATRIX
! [OUTPUTS]:
! MULCR >> M X N COMPLEX MATRIX WHICH IS EQUAL TO AB
! [DEPENDENCIES]:
! 1. MULRR(A,B) FUNCTION @ MOD_EIG
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8), DIMENSION(:,:), INTENT(IN) :: A
      REAL(P8),    DIMENSION(:,:), INTENT(IN) :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(B,2)) :: MULCR

      REAL(P8), DIMENSION(SIZE(A,1),SIZE(B,2)) :: W1,W2

      IF (SIZE(A,2).NE.SIZE(B,1)) THEN
        WRITE (*,*) 'MUL: INVALID SHAPE.'
        WRITE (*,*) 'ARGUMENT1 = ',SIZE(A,1),'X',SIZE(A,2)
        WRITE (*,*) 'ARGUMENT2 = ',SIZE(B,1),'X',SIZE(B,2)
        STOP
      ENDIF

      W1 = MULRR(REAL(A,P8),B)
      W2 = MULRR(AIMAG(A),B)

      MULCR = CMPLX(W1,W2,P8)
      
      RETURN
      END FUNCTION MULCR
!=======================================================================
      FUNCTION MULRC(A,B) ! FAMILY OF OPERATOR(.MUL.)
!=======================================================================
! [USAGE]:
! REAL-COMPLEX MATRIX MULTIPLICATION
! [INPUTS]:
! A >> M X P REAL MATRIX
! B >> P X N COMPLEX MATRIX
! [OUTPUTS]:
! MULRC >> M X N COMPLEX MATRIX WHICH IS EQUAL TO AB
! [DEPENDENCIES]:
! 1. MULRR(A,B) FUNCTION @ MOD_EIG
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),    DIMENSION(:,:),INTENT(IN) :: A
      COMPLEX(P8), DIMENSION(:,:),INTENT(IN) :: B
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(B,2)) :: MULRC
      
      REAL(P8),DIMENSION(SIZE(A,1),SIZE(B,2)) :: W1,W2
      
      IF (SIZE(A,2).NE.SIZE(B,1)) THEN
        WRITE (*,*) 'MUL: INVALID SHAPE.'
        WRITE (*,*) 'ARGUMENT1 = ',SIZE(A,1),'X',SIZE(A,2)
        WRITE (*,*) 'ARGUMENT2 = ',SIZE(B,1),'X',SIZE(B,2)
        STOP
      ENDIF
      
      W1 = MULRR(A,REAL(B,P8))
      W2 = MULRR(A,AIMAG(B))
      
      MULRC = CMPLX(W1,W2,P8)
      
      RETURN
      END FUNCTION MULRC
!=======================================================================
      FUNCTION MULRR1(A,B) ! FAMILY OF OPERATOR(.MUL.)
!=======================================================================
! [USAGE]:
! REAL MATRIX - REAL VECTOR MULTIPLICATION
! [INPUTS]:
! A >> M X N REAL MATRIX
! B >> N X 1 REAL VECTOR
! [OUTPUTS]:
! MULRR1 >> M X 1 REAL VECTOR WHICH IS EQUAL TO AB
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),DIMENSION(:,:),INTENT(IN):: A
      REAL(P8),DIMENSION(:),INTENT(IN):: B
      REAL(P8),DIMENSION(SIZE(A,1)):: MULRR1
      
      INTEGER :: NV
      
      NV = SIZE(B)
      
      IF(NV.NE.SIZE(A,2)) THEN
        WRITE (*,*) 'MULRR1: SIZE MISMATCH.'
        WRITE (*,*) '[',SIZE(A,1),',',SIZE(A,2),']*[',NV,']'
        STOP
      ENDIF
      
      MULRR1 = MATMUL(A(:,1:NV),B)
      
      RETURN
      END FUNCTION MULRR1
!=======================================================================
      FUNCTION MULRC1(A,B) ! FAMILY OF OPERATOR(.MUL.)
!=======================================================================
! [USAGE]:
! REAL MATRIX - COMPLEX VECTOR MULTIPLICATION
! [INPUTS]:
! A >> M X N REAL MATRIX
! B >> N X 1 COMPLEX VECTOR
! [OUTPUTS]:
! MULRC1 >> M X 1 COMPLEX VECTOR WHICH IS EQUAL TO AB
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),DIMENSION(:,:),INTENT(IN):: A
      COMPLEX(P8),DIMENSION(:),INTENT(IN):: B
      COMPLEX(P8),DIMENSION(SIZE(A,1)):: MULRC1

      INTEGER :: NV

      NV = SIZE(B)

      IF(NV.NE.SIZE(A,2)) THEN
        WRITE (*,*) 'MULRC1: SIZE MISMATCH.'
        WRITE (*,*) '[',SIZE(A,1),',',SIZE(A,2),']*[',NV,']'
        STOP
      ENDIF

      MULRC1 = MATMUL(A(:,1:NV),B)
      
      RETURN
      END FUNCTION MULRC1
!=======================================================================
      FUNCTION MULCR1(A,B) ! FAMILY OF OPERATOR(.MUL.)
!=======================================================================
! [USAGE]:
! COMPLEX MATRIX - REAL VECTOR MULTIPLICATION
! [INPUTS]:
! A >> M X N COMPLEX MATRIX
! B >> N X 1 REAL VECTOR
! [OUTPUTS]:
! MULCR1 >> M X 1 COMPLEX VECTOR WHICH IS EQUAL TO AB
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8),DIMENSION(:,:),INTENT(IN):: A
      REAL(P8),DIMENSION(:),INTENT(IN):: B
      COMPLEX(P8),DIMENSION(SIZE(A,1)):: MULCR1
      
      INTEGER :: NV
      
      NV = SIZE(B)
      
      IF(NV.NE.SIZE(A,2)) THEN
        WRITE (*,*) 'MULCR1: SIZE MISMATCH.'
        WRITE (*,*) '[',SIZE(A,1),',',SIZE(A,2),']*[',NV,']'
        STOP
      ENDIF
      
      MULCR1 = MATMUL(A(:,1:NV),B)
      
      RETURN
      END FUNCTION MULCR1
!=======================================================================
      FUNCTION MULCC1(A,B) ! FAMILY OF OPERATOR(.MUL.)
!=======================================================================
! [USAGE]:
! COMPLEX MATRIX - COMPLEX VECTOR MULTIPLICATION
! [INPUTS]:
! A >> M X N COMPLEX MATRIX
! B >> N X 1 COMPLEX VECTOR
! [OUTPUTS]:
! MULCC1 >> M X 1 COMPLEX VECTOR WHICH IS EQUAL TO AB
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8),DIMENSION(:,:),INTENT(IN):: A
      COMPLEX(P8),DIMENSION(:),INTENT(IN):: B
      COMPLEX(P8),DIMENSION(SIZE(A,1)):: MULCC1

      INTEGER :: NV

      NV = SIZE(B)

      IF(NV.GT.SIZE(A,2)) THEN
        WRITE (*,*) 'MULCC1: SIZE MISMATCH.'
        WRITE (*,*) '[',SIZE(A,1),',',SIZE(A,2),']*[',NV,']'
        STOP
      ENDIF

      MULCC1 = MATMUL(A(:,1:NV),B)

      RETURN
      END FUNCTION MULCC1
!=======================================================================
      FUNCTION MULC1R(A,B) ! FAMILY OF OPERATOR(.MUL.)
!=======================================================================
! [USAGE]:
! COMPLEX VECTOR - REAL MATRIX MULTIPLICATION
! [INPUTS]:
! A >> 1 X M COMPLEX VECTOR
! B >> M X N REAL MATRIX
! [OUTPUTS]:
! MULC1R >> 1 X N COMPLEX VECTOR WHICH IS EQUAL TO AB
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8),DIMENSION(:),INTENT(IN):: A
      REAL(P8),DIMENSION(:,:),INTENT(IN):: B
      COMPLEX(P8),DIMENSION(SIZE(A,1)):: MULC1R

      IF(SIZE(A).NE.SIZE(B,1)) THEN
        WRITE (*,*) 'MULC1R: SIZE MISMATCH.'
        WRITE (*,*) '[',SIZE(A),']*[',SIZE(B,1),',',SIZE(B,2),']'
        STOP
      ENDIF

      MULC1R = (/ MATMUL(A,B) /)

      RETURN
      END FUNCTION MULC1R
!=======================================================================
      FUNCTION PASTE_LR_RR(A,B) ! FAMILY OF OPERATOR(.LR.)
!=======================================================================
! [USAGE]:
! CONCATENATE TWO REAL MATRICES A, B LIKE [A | B]
! [INPUTS]:
! A >> M X N REAL MATRIX
! B >> M X P REAL MATRIX
! [OUTPUTS]:
! PASTE_LR_RR >> M X (N+P) REAL MATRIX
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),DIMENSION(:,:),INTENT(IN):: A,B
      REAL(P8),DIMENSION(SIZE(A,1),SIZE(A,2)+SIZE(B,2)):: PASTE_LR_RR

      IF(SIZE(A,1).NE.SIZE(B,1)) THEN
        WRITE (*,*) 'PASTE_LR_RR: SIZE MISMATCH.'
        WRITE (*,*) 'SIZEOF(A)=[',SIZE(A,1),',',SIZE(A,2),']'
        WRITE (*,*) 'SIZEOF(B)=[',SIZE(B,1),',',SIZE(B,2),']'
      ENDIF

      PASTE_LR_RR(:,:SIZE(A,2)) = A
      PASTE_LR_RR(:,SIZE(A,2)+1:) = B

      RETURN
      END FUNCTION PASTE_LR_RR
!=======================================================================
      FUNCTION PASTE_LR_CC(A,B) ! FAMILY OF OPERATOR(.LR.)
!=======================================================================
! [USAGE]:
! CONCATENATE TWO COMPLEX MATRICES A, B LIKE [A | B]
! [INPUTS]:
! A >> M X N COMPLEX MATRIX
! B >> M X P COMPLEX MATRIX
! [OUTPUTS]:
! PASTE_LR_CC >> M X (N+P) COMPLEX MATRIX
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8),DIMENSION(:,:),INTENT(IN):: A,B
      COMPLEX(P8),DIMENSION(SIZE(A,1),SIZE(A,2)+SIZE(B,2)):: PASTE_LR_CC

      IF(SIZE(A,1).NE.SIZE(B,1)) THEN
        WRITE (*,*) 'PASTE_LR_CC: SIZE MISMATCH.'
        WRITE (*,*) 'SIZEOF(A)=[',SIZE(A,1),',',SIZE(A,2),']'
        WRITE (*,*) 'SIZEOF(B)=[',SIZE(B,1),',',SIZE(B,2),']'
      ENDIF

      PASTE_LR_CC(:,:SIZE(A,2)) = A
      PASTE_LR_CC(:,SIZE(A,2)+1:) = B

      RETURN
      END FUNCTION PASTE_LR_CC
!=======================================================================
      FUNCTION PASTE_UD_RR(A,B) ! FAMILY OF OPERATOR(.UD.)
!=======================================================================
! [USAGE]:                                / A \
! CONCATENATE TWO REAL MATRICES A, B LIKE |---|
! [INPUTS]:                               \ B /
! A >> M X N REAL MATRIX
! B >> P X N REAL MATRIX
! [OUTPUTS]:
! PASTE_UD_RR >> (M+P) X N REAL MATRIX
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      REAL(P8),DIMENSION(:,:),INTENT(IN):: A,B
      REAL(P8),DIMENSION(SIZE(A,1)+SIZE(B,1),SIZE(A,2)):: PASTE_UD_RR

      IF(SIZE(A,2).NE.SIZE(B,2)) THEN
        WRITE (*,*) 'PASTE_UD_RR: SIZE MISMATCH.'
        WRITE (*,*) 'SIZEOF(A)=[',SIZE(A,1),',',SIZE(A,2),']'
        WRITE (*,*) 'SIZEOF(B)=[',SIZE(B,1),',',SIZE(B,2),']'
      ENDIF

      PASTE_UD_RR(:SIZE(A,1),:) = A
      PASTE_UD_RR(SIZE(A,1)+1:,:) = B

      RETURN
      END FUNCTION PASTE_UD_RR
!=======================================================================
      FUNCTION PASTE_UD_CC(A,B) ! FAMILY OF OPERATOR(.UD.)
!=======================================================================
! [USAGE]:                                   / A \
! CONCATENATE TWO COMPLEX MATRICES A, B LIKE |---|
! [INPUTS]:                                  \ B /
! A >> M X N COMPLEX MATRIX
! B >> P X N COMPLEX MATRIX
! [OUTPUTS]:
! PASTE_UD_CC >> (M+P) X N COMPLEX MATRIX
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      COMPLEX(P8),DIMENSION(:,:),INTENT(IN):: A,B
      COMPLEX(P8),DIMENSION(SIZE(A,1)+SIZE(B,1),SIZE(A,2)):: PASTE_UD_CC

      IF(SIZE(A,2).NE.SIZE(B,2)) THEN
        WRITE (*,*) 'PASTE_UD_CC: SIZE MISMATCH.'
        WRITE (*,*) 'SIZEOF(A)=[',SIZE(A,1),',',SIZE(A,2),']'
        WRITE (*,*) 'SIZEOF(B)=[',SIZE(B,1),',',SIZE(B,2),']'
      ENDIF

      PASTE_UD_CC(:SIZE(A,1),:) = A
      PASTE_UD_CC(SIZE(A,1)+1:,:) = B

      RETURN
      END FUNCTION PASTE_UD_CC
!=======================================================================
      FUNCTION INVDR(A) ! FAMILY OF INV
!=======================================================================
! [USAGE]:                                   
! FIND THE INVERSE OF A REAL SQUARE MATRIX A
! [INPUTS]:                                  
! A >> M X N REAL MATRIX (M SHOULD BE EQUAL TO N)
! [OUTPUTS]:
! INVDR >> INVERSE OF A (N X N)
! [DEPENDENCIES]:
! 1. DGETRF(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! 2. DGETRI(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
! IMPLEMENT LAPACK'S DGETRI FUNCTION @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL DGETRF
      EXTERNAL DGETRI
      REAL(P8),DIMENSION(:,:):: A
      REAL(P8),DIMENSION(SIZE(A,1),SIZE(A,2)):: INVDR

      REAL(P8),DIMENSION(SIZE(A,1)):: WK
      INTEGER(P4),DIMENSION(SIZE(A,1)):: IP
      INTEGER(P4):: NI,NJ, INFO

      NI=SIZE(A,1)
      NJ=SIZE(A,2)

      IF(NI.NE.NJ) THEN
        WRITE (*,*) 'INVDR: MATRIX NOT SQUARE.'
        STOP
      ENDIF

      INVDR = A

      CALL DGETRF(NI,NI,INVDR,NI,IP,INFO)

      IF (INFO .NE. 0) THEN
        WRITE (*,*) 'INVDR: SINGULAR MATRIX. CANNOT BE INVERTED.'
        STOP
      ENDIF

      CALL DGETRI(NI,INVDR,NI,IP,WK,NI,INFO)

      IF (INFO .NE. 0) THEN
        WRITE (*,*) 'INVDR: MATRIX INVERSION FAILURE.'
        STOP
      ENDIF

      RETURN
      END FUNCTION INVDR
!=======================================================================
      FUNCTION INVDC(A) ! FAMILY OF INV
!=======================================================================
! [USAGE]:                                   
! FIND THE INVERSE OF A COMPLEX SQUARE MATRIX A
! [INPUTS]:                                  
! A >> M X N COMPLEX MATRIX (M SHOULD BE EQUAL TO N)
! [OUTPUTS]:
! INVDC >> INVERSE OF A (N X N)
! [DEPENDENCIES]:
! 1. ZGETRF(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! 2. ZGETRI(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
! IMPLEMENT LAPACK'S ZGETRI FUNCTION @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL ZGETRF
      EXTERNAL ZGETRI
      COMPLEX(P8),DIMENSION(:,:):: A
      COMPLEX(P8),DIMENSION(SIZE(A,1),SIZE(A,2)):: INVDC

      COMPLEX(P8),DIMENSION(SIZE(A,1)):: WK
      INTEGER(P4),DIMENSION(SIZE(A,1)):: IP
      INTEGER(P4):: NI,NJ,INFO

      NI=SIZE(A,1)
      NJ=SIZE(A,2)

      IF(NI.NE.NJ) THEN
        WRITE (*,*) 'INVDC: MATRIX NOT SQUARE.'
        STOP
      ENDIF

      INVDC = A

      CALL ZGETRF(NI,NI,INVDC,NI,IP,INFO)

      IF (INFO .NE. 0) THEN
        WRITE (*,*) 'INVDC: SINGULAR MATRIX. CANNOT BE INVERTED.'
        STOP
      ENDIF

      CALL ZGETRI(NI,INVDC,NI,IP,WK,NI,INFO)

      IF (INFO .NE. 0) THEN
        WRITE (*,*) 'INVDC: MATRIX INVERSION FAILURE.'
        STOP
      ENDIF

      RETURN
      END FUNCTION INVDC
!=======================================================================
      FUNCTION GENEIG1(A,B) ! FAMILY OF GENEIG
!=======================================================================
! [USAGE]:                                   
! SOLVE THE GENRALIZED EIGENVALUE PROBLEM AV=BVD FOR GIVEN A AND B
! AND YIELDS EIGENVALUES (DIAG(D)).
! FOR STANDARD EIGENVALUE PROBLEMS, SET B AS AN IDENTITY MATRIX.
! [INPUTS]:                                  
! A >> N X N COMPLEX MATRIX
! B >> N X N COMPLEX MATRIX (OPTIONAL. IF NONE, THEN SOLVE AV=VD)
! [OUTPUTS]:
! GENEIG1 >> EIGENVALUES OF THE GENERALIZED EIGENVALUE PROBLEM AV = BVD
! [DEPENDENCIES]:
! 1. ZGGEV(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! 2. EYE1(A)  @ MOD_EIG
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
! IMPLEMENT LAPACK'S ZGGEV FUNCTION @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      EXTERNAL ZGGEV
      COMPLEX(P8), DIMENSION(:,:):: A
      COMPLEX(P8), DIMENSION(:,:),OPTIONAL:: B
      COMPLEX(P8), DIMENSION(SIZE(A,1)):: GENEIG1

      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(A,2)):: W1,W2

      COMPLEX(P8), DIMENSION(SIZE(A,1)):: AEIG,BEIG
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(A,2)):: VL,VR
      COMPLEX(P8), DIMENSION(2*SIZE(A,1)):: WK
      REAL(P8)   , DIMENSION(8*SIZE(A,1)):: RWK

      INTEGER(P4) :: I,NI,LW,INFO

      IF(SIZE(A,1).NE.SIZE(A,2)) THEN
        WRITE (*,*) 'GENEIG1: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF(PRESENT(B)) THEN
        IF(SIZE(A).NE.SIZE(B)) THEN
          WRITE (*,*) 'GENEIG1: MATRIX SIZE INCONSISTENT'
          RETURN
        ENDIF
      ENDIF

      NI=SIZE(A,1)

      IF(PRESENT(B)) THEN
        W1=A
        W2=B
      ELSE
        W1=A
        W2=EYE(NI) !EYE1(NI)
      ENDIF

      LW = SIZE(WK)
      CALL ZGGEV('N','N',NI,W1,NI,W2,NI,AEIG,BEIG, &
                            VL,NI,VR,NI,WK,LW,RWK,INFO)

      IF (INFO .NE. 0) THEN
        WRITE (*,*) 'GENEIG1: EIGENVALUE GENERATION FAILURE.'
        STOP
      ENDIF

      DO I = 1, NI
        IF(ABS(BEIG(I)).EQ.0.D0) THEN
          GENEIG1(I) = HUGE(0.D0)                                        ! ASSIGN INFINITY TO AVOID THE DIVISION-BY-ZERO ERROR
        ELSE
          GENEIG1(I) = AEIG(I)/BEIG(I)
        ENDIF
      ENDDO

      RETURN
      END FUNCTION GENEIG1
!=======================================================================
      FUNCTION EIGVEC1(LR,A,EIV) ! FAMILY OF EIGVEC
!=======================================================================
! [USAGE]:                                   
! CALCULATE THE EIGENVECTORS FOR A SINGLE EIGENVALUE INPUT EIV
! THE EIGENVALUE PROBLEM IS SOLVED (AV=VD)
! [INPUTS]:                
! LR >> 1. 'L' : GET LEFT EIGENVECTORS / 2. 'R' : GET RIGHT EIGENVECTORS                  
! A >> N X N COMPLEX MATRIX
! EIV >> A COMPLEX EIGVAL (OPTIONAL. IF GIVEN, CALCULATE EIV'S EIGVEC)
! [OUTPUTS]:
! EIGVEC1 >> (N X 1) EIGENVECTOR(S) CORRESPONDING TO EIV
!             FOR NON-REPEATED EIGENVALUES, 1ST COLUMN ONLY MEANINGFUL
!             FOR K-REPEATED EIGENVALUES, 1~KTH COLUMNS ARE MEANINGFUL
!             IF NOT EIGENVALUE, SIMPLY OUTPUT A N X N ZERO MATRIX
! K >> DIMENSION OF THE EIGENSPACE WITH EIV
! [DEPENDENCIES]:
! 1. ZGGEV(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! 2. EYE(A)  @ MOD_EIG
! [UPDATES]:
! CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN)                 :: LR
      COMPLEX(P8), DIMENSION(:,:), INTENT(IN) :: A
      COMPLEX(P8), INTENT(IN), OPTIONAL       :: EIV
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(A,2)) :: EIGVEC1
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(A,2)) :: W1,W2

      COMPLEX(P8), DIMENSION(SIZE(A,1)):: AEIG,BEIG
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(A,2)):: VV
      COMPLEX(P8), DIMENSION(2*SIZE(A,1)):: WK
      REAL(P8)   , DIMENSION(8*SIZE(A,1)):: RWK

      INTEGER(P4) :: I,NI,LW,INFO,COUNTER
      COMPLEX(P8), DIMENSION(SIZE(A,1)):: EIGVALS

      IF(SIZE(A,1).NE.SIZE(A,2)) THEN
        WRITE (*,*) 'EIGVEC1: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      NI=SIZE(A,1)

      W1=A
      W2=EYE(NI) !EYE1(NI)

      LW = SIZE(WK)
      EIGVEC1 = (0.D0, 0.D0)

      IF(LR.EQ.'L') THEN
        CALL ZGGEV('V','N',NI,W1,NI,W2,NI,AEIG,BEIG, &
                              EIGVEC1,NI,VV,NI,WK,LW,RWK,INFO)
      ELSEIF(LR.EQ.'R') THEN
        CALL ZGGEV('N','V',NI,W1,NI,W2,NI,AEIG,BEIG, &
                              VV,NI,EIGVEC1,NI,WK,LW,RWK,INFO)
      ELSE
        WRITE (*,*) 'EIGVEC1: INCORRECT LEFT/RIGHT EIGENVECTOR SPECIFIER.'
        STOP
      ENDIF

      IF (INFO .NE. 0) THEN
        WRITE (*,*) 'EIGVEC1: EIGENVECTOR GENERATION FAILURE.'
        STOP
      ENDIF

      IF(PRESENT(EIV)) THEN
        COUNTER = 0
        DO I = 1, NI
          IF(ABS(BEIG(I)).EQ.0.D0) THEN
            EIGVALS(I) = HUGE(0.D0)                                        ! ASSIGN INFINITY TO AVOID THE DIVISION-BY-ZERO ERROR
          ELSE
            EIGVALS(I) = AEIG(I)/BEIG(I)
          ENDIF
          
          IF(EIV .EQ. EIGVALS(I)) THEN
            COUNTER = COUNTER + 1
            EIGVEC1(:,COUNTER) = EIGVEC1(:,I)
          ENDIF
        ENDDO

        EIGVEC1(:,COUNTER+1:) = (0.D0, 0.D0)
      ENDIF

      RETURN
      END FUNCTION EIGVEC1
!=======================================================================
      FUNCTION EIGVEC2(LR,A,B,EIV) ! FAMILY OF EIGVEC
!=======================================================================
! [USAGE]:                                   
! CALCULATE THE EIGENVECTORS FOR A SINGLE EIGENVALUE INPUT EIV
! THE GENERALIZED EIGENVALUE PROBLEM IS SOLVED (AV=BVD)
! [INPUTS]:                
! LR >> 1. 'L' : GET LEFT EIGENVECTORS / 2. 'R' : GET RIGHT EIGENVECTORS                  
! A >> N X N COMPLEX MATRIX
! B >> N X N COMPLEX MATRIX
! EIV >> A COMPLEX EIGVAL (OPTIONAL. IF GIVEN, CALCULATE EIV'S EIGVEC)
! [OUTPUTS]:
! EIGVEC2 >> (N X 1) EIGENVECTOR(S) CORRESPONDING TO EIV
!             FOR NON-REPEATED EIGENVALUES, 1ST COLUMN ONLY MEANINGFUL
!             FOR K-REPEATED EIGENVALUES, 1~KTH COLUMNS ARE MEANINGFUL
!             IF NOT EIGENVALUE, SIMPLY OUTPUT A N X N ZERO MATRIX
! K >> DIMENSION OF THE EIGENSPACE WITH EIV
! [DEPENDENCIES]:
! 1. ZGGEV(~) @ LAPACK (OR MKL IF USING INTEL COMPILER)
! [UPDATES]:
! CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      IMPLICIT NONE
      CHARACTER*1, INTENT(IN)                 :: LR
      COMPLEX(P8), DIMENSION(:,:), INTENT(IN) :: A, B
      COMPLEX(P8), INTENT(IN), OPTIONAL       :: EIV
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(A,2)) :: EIGVEC2
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(A,2)) :: W1,W2

      COMPLEX(P8), DIMENSION(SIZE(A,1)):: AEIG,BEIG
      COMPLEX(P8), DIMENSION(SIZE(A,1),SIZE(A,2)):: VV
      COMPLEX(P8), DIMENSION(2*SIZE(A,1)):: WK
      REAL(P8)   , DIMENSION(8*SIZE(A,1)):: RWK

      INTEGER(P4) :: I,NI,LW,INFO,COUNTER
      COMPLEX(P8), DIMENSION(SIZE(A,1)):: EIGVALS

      IF(SIZE(A,1).NE.SIZE(A,2)) THEN
        WRITE (*,*) 'EIGVEC2: MATRIX NOT SQUARE'
        WRITE (*,*) SIZE(A,1),'X',SIZE(A,2)
        STOP
      ENDIF

      IF(SIZE(A).NE.SIZE(B)) THEN
        WRITE (*,*) 'EIGVEC2: MATRIX SIZE INCONSISTENT'
        STOP
      ENDIF

      NI=SIZE(A,1)

      W1=A
      W2=B

      LW = SIZE(WK)
      EIGVEC2 = (0.D0, 0.D0)

      IF(LR.EQ.'L') THEN
        CALL ZGGEV('V','N',NI,W1,NI,W2,NI,AEIG,BEIG, &
                              EIGVEC2,NI,VV,NI,WK,LW,RWK,INFO)
      ELSEIF(LR.EQ.'R') THEN
        CALL ZGGEV('N','V',NI,W1,NI,W2,NI,AEIG,BEIG, &
                              VV,NI,EIGVEC2,NI,WK,LW,RWK,INFO)
      ELSE
        WRITE (*,*) 'EIGVEC2: INCORRECT LEFT/RIGHT EIGENVECTOR SPECIFIER.'
        STOP
      ENDIF

      IF (INFO .NE. 0) THEN
        WRITE (*,*) 'EIGVEC2: EIGENVECTOR GENERATION FAILURE.'
        STOP
      ENDIF

      IF(PRESENT(EIV)) THEN
        COUNTER = 0
        DO I = 1, NI
          IF(ABS(BEIG(I)).EQ.0.D0) THEN
            EIGVALS(I) = HUGE(0.D0)                                        ! ASSIGN INFINITY TO AVOID THE DIVISION-BY-ZERO ERROR
          ELSE
            EIGVALS(I) = AEIG(I)/BEIG(I)
          ENDIF
          
          IF(EIV .EQ. EIGVALS(I)) THEN
            COUNTER = COUNTER + 1
            EIGVEC2(:,COUNTER) = EIGVEC2(:,I)
          ENDIF
        ENDDO

        EIGVEC2(:,COUNTER+1:) = (0.D0, 0.D0)
      ENDIF

      RETURN
      END FUNCTION EIGVEC2
!=======================================================================
      FUNCTION FLIPUDR(A)
!=======================================================================
! [USAGE]:                                   
! FLIP A REAL MATRIX UPSIDE DOWN
! [INPUTS]:                                  
! A >> M X N REAL MATRIX
! [OUTPUTS]:
! FLIPUDR >> VERTICALLY FLIPPED MATRIX A
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      REAL(P8),DIMENSION(:,:),INTENT(IN) :: A
      REAL(P8),DIMENSION(SIZE(A,1),SIZE(A,2)):: FLIPUDR
      INTEGER(P4) :: I,N

      N = SIZE(A,1)

      DO I=1,N
        FLIPUDR(I,:) = A(N-I+1,:)
      ENDDO

      RETURN
      END FUNCTION FLIPUDR
!=======================================================================
      FUNCTION FLIPUDC(A)
!=======================================================================
! [USAGE]:                                   
! FLIP A COMPLEX MATRIX UPSIDE DOWN
! [INPUTS]:                                  
! A >> M X N COMPLEX MATRIX
! [OUTPUTS]:
! FLIPUDR >> VERTICALLY FLIPPED MATRIX A
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      COMPLEX(P8),DIMENSION(:,:),INTENT(IN) :: A
      COMPLEX(P8),DIMENSION(SIZE(A,1),SIZE(A,2)):: FLIPUDC
      INTEGER(P4) :: I,N

      N = SIZE(A,1)

      DO I=1,N
        FLIPUDC(I,:) = A(N-I+1,:)
      ENDDO

      RETURN
      END FUNCTION FLIPUDC
!=======================================================================
      FUNCTION FLIPLRR(A)
!=======================================================================
! [USAGE]:                                   
! FLIP A REAL MATRIX FROM LEFT TO RIGHT
! [INPUTS]:                                  
! A >> M X N REAL MATRIX
! [OUTPUTS]:
! FLIPUDR >> HORIZONTALLY FLIPPED MATRIX A
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      REAL(P8),DIMENSION(:,:),INTENT(IN) :: A
      REAL(P8),DIMENSION(SIZE(A,1),SIZE(A,2)):: FLIPLRR
      INTEGER(P4) :: J,N

      N = SIZE(A,2)
      DO J=1,N
        FLIPLRR(:,J) = A(:,N-J+1)
      ENDDO

      RETURN
      END FUNCTION FLIPLRR
!=======================================================================
      FUNCTION FLIPLRC(A)
!=======================================================================
! [USAGE]:                                   
! FLIP A COMPLEX MATRIX FROM LEFT TO RIGHT
! [INPUTS]:                                  
! A >> M X N COMPLEX MATRIX
! [OUTPUTS]:
! FLIPUDR >> HORIZONTALLY FLIPPED MATRIX A
! [UPDATES]:
! ORIGINALLY WRITTEN BY TATSUHITO MATSUSHIMA 
! RE-CODED BY SANGJOON LEE @ NOV 8 2020
!=======================================================================
      COMPLEX(P8),DIMENSION(:,:),INTENT(IN) :: A
      COMPLEX(P8),DIMENSION(SIZE(A,1),SIZE(A,2)):: FLIPLRC
      INTEGER(P4) :: J,N

      N = SIZE(A,2)
      DO J=1,N
        FLIPLRC(:,J) = A(:,N-J+1)
      ENDDO

      RETURN
      END FUNCTION FLIPLRC
!=======================================================================
END MODULE MOD_EIG
!=======================================================================
!                      ORIGINAL COMMENTS (BY TATSU)
!=======================================================================
!     USE EIG  MAKES THE INTERFACE EXPLICIT AND ALLOWS
!     THE USE OF THE GENERIC NAMES
!=======================================================================
! OPERATORS:
!   .MUL.   MATRIX/VECTOR MULTIPLICATION
!   .LR.    LEFT-RIGHT MATRIX PASTE
!   .UD.    UP-DOWN MATRIX PASTE
!
! FUNCTIONS:
!   GENEIG(A)       RETURNS AN ONE-DIMENSIONAL COMPLEX ARRAY 
!                   CONTAINING THE EIGENVALUES OF MATRIX A(:,:).
!                   MATRIX A MUST BE SQUARE AND COMPLEX.
!   GENEIG(A,B) 
!                   SOLVES GENERALIZED EIGENVALUE PROBLEM A.X=(LAMB)B.X
!
!   EIGVEC(A,EIG)   RETURNS A MATRIX WHOSE COLUMNS ARE THE 
!                   EIGENVECTORS CORRESPONDING TO THE EIGENVALUES EIG.
!                   A(:,:) AND EIG(:) ARE COMPLEX ARRAYS. NUMBER OF 
!                   COLUMNS IN THE RESULT IS THE SAME AS THE NUMBER OF 
!                   ELEMENTS IN EIG. EIGENVECTORS ARE NORMALIZED SUCH
!                   THAT ITS NORM IS ONE SUM(CONJG(EIVEC)*EIVEC)=1.
!
!   EIGVEC(A,B,EIG) IS FOR GENERALIZED EIGENVALUE PROBLEM.
!
!   EIGVEC(A,EIG,ADJ), EIGVEC(A,B,EIG,ADJ) IN ADDITION RETURNS THE 
!             ADJOINT EIGENVECTORS IN COLUMNS OF COMPLEX(P8) ADJ(:,:).
!             ADJOINT IS NORMALIZED SUCH THAT SUM(ADJ*EIVEC) = 1
!             (NOT SUM(CONJG(ADJ)*EIVEC), NOR SUM(CONJG(ADJ)*A*ADJ).
!
!   INV(A)  RETURNS THE INVERSE MATRIX OF A. 
!           MATRIX A MUST BE SQUARE AND CAN BE EITHER 
!           COMPLEX OR REAL.  
!
!   ZEROS(NI,NJ)  RETURNS A REAL NI*NJ MATRIX WHOSE ELEMENTS ARE ZERO.
!                 NI,NJ ARE INTEGERS.
!
!   EYE(NI) RETURNS A REAL IDENTITY MATRIX OF ORDER NI. NI IS INTEGER.
!
!   DIAG(X) IF X IS ONE DIMENSIONAL, DIAG(X) IS A TWO DIMENSIONAL
!           ARRAY WHOSE DIAGONALS ARE THE ELEMENTS OF X. IF X IS TWO
!           DIMENSIONAL, DIAG(X) IS ONE DIMENSIONAL WHOSE ELEMENTS ARE
!           THE DIAGONALS OF X.
!
!   FLIPLR(X) FLIPS MATRIX X IN THE LEFT/RIGHT DIRECTION.
!
!   FLIPUD(X) FLIPS MATRIX X IN THE UP/DOWN DIRECTION.
!
! SUBROUTINES
!   LU(A,W) DECOMPOSES REAL(P8)/COMPLEX(P8) A(:,:) INTO LU. 
!           INTEGER W(N,2) IS A WORK AREA WITH N=SIZE(A,1)
!
!   SOLVE(A,W,B) SOLVES A.X=B WHERE A(:,:), W(:,:) IS THE RESULT OF LU 
!                ABOVE, REAL(P8)/COMPLEX(P8) B(:,:). ON OUTPUT, 
!                B IS REPLACED BY X.
!
!   SPY(A)  SHOWS THE STRUCTURE OF THE MATRIX A(:,:) BASED ON THE VALUES
!           OF ABS(A). MATRIX A CAN BE COMPLEX OR REAL. 
!
!   MCAT(A) PRINTS THE VALUES OF MATRIX A(:,:). IF A IS ONE DIMENSIONAL, 
!           VALUES ARE PRINTED OUT IN A GREATER(~16 DIGITS) ACCURACY.
!           MATRIX A CAN BE COMPLEX OR REAL.
!     
!   MSAVE(A,'AA') SAVES MATRIX A TO DISK UNDER THE NAME 'AA'.
!             MATRIX A CAN BE COMPLEX OR REAL AND CAN BE ONE/TWO 
!             DIMENSIONAL. IF A IS COMPLEX VALUED, THE NUMBER OF ROWS
!             ARE DOUBLED. THE FORMAT IS      
!                REAL(A(1,1)) IMAG(A(1,1)) REAL(A(1,2)) ....
!                REAL(A(2,1)) IMAG(A(2,1)) REAL(A(2,2)) ....
!                                   .
!                                   .
!
!             THE FILE CAN BE READ BY MLOAD BELOW. 
!             ALSO, MATLAB CAN READ THE FILE AS FOLLOWS. 
!
!             >> LOAD AA
!
!             IF THE MATRIX IS COMPLEX, A M-FILE
!
!             FUNCTION [X]=MLOADC(A)
!             %  MLOADC   CONVERT ASCII FILE TO COMPLEX MATRIX
!             %         EXAMPLE:   
!             %         >> LOAD ABC
!             %         >> ABC=MLOADC(ABC);
!             %
!             %         SEE ALSO MSAVEC, FORTRAN ROUTINES MSAVE, MLOAD
!             %
!             [NI,NJ]=SIZE(A);
!             X=A(:,1:2:NJ-1)+1I*A(:,2:2:NJ);
!          
!             CAN RECOVER THE MATRIX.
!
!   MLOAD('AA',A) LOADS THE MATRIX A(:,:). MATRIX A CAN BE COMPLEX OR 
!             REAL. THE SIZE OF A MUST AGREE WITH THE SIZE OF THE FILE
!             TO BE READ. MLOAD CAN LOAD A MATRIX SAVED BY MATLAB BY
!
!             >> SAVE AA AA -ASCII -DOUBLE
!
!             TO SAVE A COMPLEX MATRIX, FOLLOWING PREPROCESSING BY A 
!             M-FILE IS USEFUL.
!
!             FUNCTION [X]=MSAVEC(A)
!             %  MSAVEC   CONVERT A COMPLEX MATRIX TO A REAL MATRIX 
!             %         FOR ASCII SAVE.
!             %         EXAMPLE:   
!             %         >> A=MSAVEC(ABC);
!             %         >> SAVE A A -ASCII -DOUBLE
!             %
!             %         SEE ALSO MLOADC, FORTRAN ROUTINES MSAVE, MLOAD
!             %
!             [NI,NJ]=SIZE(A);
!             X=ZEROS(NI,NJ*2);
!             FOR J=1:NJ
!               X(:,J*2-1)=REAL(A(:,J));
!               X(:,J*2  )=IMAG(A(:,J));
!             END
!
!   MLOAD('AA',A,NI,NJ) LOADS A REAL OR COMPLEX MATRIX OF
!             UNKNOWN SIZE. OUTPUTS NI,NJ ARE THE SIZE OF THE LOADED 
!             MATRIX. ARRAY A MUST BE LARGE ENOUGH TO STORE THE NI*NJ 
!             MATRIX.
!=======================================================================