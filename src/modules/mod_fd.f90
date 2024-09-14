MODULE MOD_FD ! LEVEL 1 MODULE
! UTILITY FUNCTIONS USED BY MARCH
    USE omp_lib
    USE MPI
    USE MOD_MISC, ONLY : P4,P8,PI,IU                                   ! LEVEL 0
    IMPLICIT NONE
    PRIVATE
!=======================================================================
!============================ PARAMETERS ===============================
!=======================================================================
! CURRENTLY NONE
!=======================================================================
!======================== PUBLIC DECLARATION ===========================
!=======================================================================
    ! ! SMOOTHING BY CONVOLUTION
    ! PUBLIC:: SMOOTH
    ! NUMERICAL INTEGRATION
    PUBLIC:: CUMTRAPZ
    ! ! UNUSED FUNCTIONS
    ! PUBLIC:: DC,DDC,OPER
CONTAINS
!=======================================================================
!============================= FUNCTIONS ===============================
!=======================================================================
    FUNCTION SMOOTH(A)
!=======================================================================
! [USAGE]: 
! SMOOTHING A SEQUENCE A BY CONVOLUTION
! [INPUTS]:
! A >> A SEQUENCE CONTAINING NOISES
! [OUTPUTS]:
! SMOOTH >> SMOOTHENED SEQUENCE USING CONVOLUTION
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 20 2020
!=======================================================================
    IMPLICIT NONE
    COMPLEX(P8),DIMENSION(:):: A
    COMPLEX(P8),DIMENSION(SIZE(A)):: SMOOTH

    REAL(P8),DIMENSION(3):: FAE=(/ 0.2D0, 0.3D0, 0.5D0 /)
    REAL(P8):: F
    INTEGER:: NI,I,NB,NE

    NI=SIZE(A)

    IF(NI.LT.3) THEN
      WRITE(*,*) 'SMOOTH: SIZE TOO SMALL'
      STOP
    ENDIF

    SMOOTH=0
    SMOOTH(1)=A(1)
    SMOOTH(NI-2:NI)= SMOOTH(NI-2:NI) +A(NI)*FAE

    DO I=2,NI-1
      F = (1+((NI-I)/(NI-1.0D0)))*0.5D0
      SMOOTH(I-1)= SMOOTH(I-1)+A(I)*(1-F)/2
      SMOOTH(I  )= SMOOTH(I  )+A(I)*F
      SMOOTH(I+1)= SMOOTH(I+1)+A(I)*(1-F)/2
    ENDDO

    RETURN
    END FUNCTION SMOOTH
!=======================================================================
    FUNCTION DDC(R)
!=======================================================================
! [USAGE]: 
! SMOOTHING A SEQUENCE A BY CONVOLUTION
! [INPUTS]:
! A >> A SEQUENCE CONTAINING NOISES
! [OUTPUTS]:
! SMOOTH >> SMOOTHENED SEQUENCE USING CONVOLUTION
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 20 2020
!=======================================================================
    REAL(P8),DIMENSION(:):: R
    REAL(P8),DIMENSION(SIZE(R),3):: DDC
    REAL(P8):: X0,X1,X2
    INTEGER:: NI,I
    NI = SIZE(R)
    DDC(1,:)=0
    DDC(NI,:)=0
    DO I=2,NI-1
    X0=R(I-1)
    X1=R(I  )
    X2=R(I+1)
    DDC(I,1) = 2/(X0-X1)/(X0-X2)
    DDC(I,2) = 2/(X1-X0)/(X1-X2)
    DDC(I,3) = 2/(X2-X0)/(X2-X1)
    ENDDO

    RETURN
    END FUNCTION DDC
!=======================================================================
    FUNCTION DC(R)
!=======================================================================
    REAL(P8),DIMENSION(:):: R
    REAL(P8),DIMENSION(SIZE(R),3):: DC

    REAL(P8):: X0,X1,X2
    INTEGER:: NI,I

    NI = SIZE(R)
    DC(1,:)=0
    DC(NI,:)=0

    DO I=2,NI-1
      X0=R(I-1)
      X1=R(I  )
      X2=R(I+1)
      DC(I,1) = (X1-X2)/(X0-X1)/(X0-X2)
      DC(I,2) = (2*X1-X0-X2)/(X1-X0)/(X1-X2)
      DC(I,3) = (X1-X0)/(X2-X0)/(X2-X1)
    ENDDO

    RETURN
    END FUNCTION DC
!=======================================================================
    FUNCTION OPER(OP,F)
!=======================================================================
    REAL(P8),DIMENSION(:):: F
    REAL(P8),DIMENSION(:,:):: OP
    REAL(P8),DIMENSION(SIZE(F)):: OPER

    INTEGER:: NI,I

    NI = SIZE(F)

    OPER(1)=0
    OPER(NI)=0

    OPER(2:NI-1) = OP(2:NI-1,1)*F(1:NI-2)  &
    + OP(2:NI-1,2)*F(2:NI-1)  &
    + OP(2:NI-1,3)*F(3:NI  )

    RETURN
    END FUNCTION OPER
!=======================================================================
    PURE FUNCTION CUMTRAPZ(X,Y)
!=======================================================================
! CUMULATIVE TRAPEZOIDAL RULE FOR NUMERICAL INTEGRATION - SERIAL
! -> x * cumsum((y(1:end-1,:) + y(2:end,:)),1)]/2
! ======================================================================
    IMPLICIT NONE
    REAL(P8),DIMENSION(:),INTENT(IN):: X, Y
    REAL(P8),DIMENSION(SIZE(Y)):: CUMTRAPZ, DX

    INTEGER:: NI

    NI = SIZE(X)

    DX = X(2:NI)-X(1:NI-1)

    CUMTRAPZ = 0.D0
    CUMTRAPZ(2:NI) = CUMSUM(DX * (Y(1:NI-1)+Y(2:NI)) / 2)

    RETURN
    END FUNCTION CUMTRAPZ
!=======================================================================
    PURE FUNCTION CUMSUM(X)
!=======================================================================
! CUMULATIVE SUM - SERIAL
! ======================================================================
    IMPLICIT NONE
    REAL(P8),DIMENSION(:),INTENT(IN):: X
    REAL(P8),DIMENSION(SIZE(X)):: CUMSUM

    INTEGER:: NI,I

    NI = SIZE(X)

    CUMSUM(1) = X(1)
    DO I = 2,NI
      CUMSUM(I) = CUMSUM(I-1) + X(I)
    ENDDO

    RETURN
    END FUNCTION CUMSUM
!=======================================================================


END MODULE MOD_FD