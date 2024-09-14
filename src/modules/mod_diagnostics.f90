MODULE MOD_DIAGNOSTICS ! LEVEL 5 MODULE
   USE OMP_LIB
   USE MPI
   USE MOD_MISC, ONLY : P4,P8,PI,IU,SPY,MSAVE,MLOAD,ITOA3             ! LEVEL 0
!XUSE USE MOD_FD                                                         ! LEVEL 1
   USE MOD_EIG                                                        ! LEVEL 1
   USE MOD_LIN_LEGENDRE                                               ! LEVEL 1
!XUSE USE MOD_BANDMAT                                                    ! LEVEL 1
   USE MOD_SCALAR3                                                    ! LEVEL 2
   USE MOD_FFT                                                        ! LEVEL 2.5
   USE MOD_LEGOPS                                                     ! LEVEL 3
!XUSE USE MOD_LAYOUT                                                     ! LEVEL 3
   USE MOD_MARCH                                                      ! LEVEL 4
   IMPLICIT NONE
   PRIVATE
!=======================================================================
!============================ PARAMETERS ===============================
!=======================================================================
   ! 1> DIMENSIONS
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:):: MONITOR_MK          ! EIGENVALUE MONITOR FOR DIFFERENT AZIMUTHAL AND AXIAL WAVENUMBERS

   ! 2> HYPERV
   REAL(P8),DIMENSION(:),ALLOCATABLE :: HYPER_R, HYPER_T, HYPER_X     ! USED IN HYPSET AND HYPERV3

!=======================================================================
!======================== PUBLIC DECLARATION ===========================
!=======================================================================
   ! MONITOR EIGENVALUES
   PUBLIC:: MONITOR_EIG
   ! ENERGY GROWTH RATE
   PUBLIC:: EGROWTH
   PUBLIC:: EGROWTH_MODIFIED
   ! ENERGY SPECTRUM
   PUBLIC:: ENERGY_SPEC, ENERGY_SPEC_MODIFIED
   PUBLIC:: PRINT_ENERGY_SPECTRUM
   ! WRITE ENERGY AND MODE DATA TO FILE
   PUBLIC:: WRITE_ENERGY_DATA_TO_FILE
   PUBLIC:: WRITE_MODE_DATA_TO_FILE
   ! CALCULATE INTEGRATION OF PRODUCT*(1-X)^2
   PUBLIC:: PRODUCT_MK
   PUBLIC:: PRODUCT_MK_MODIFIED
   ! SAVE COLLOCATION INFO
   PUBLIC:: COLLOC_INFO
   ! CALCULATE MAXIMUM AZIMUTHAL VELOCITY
   PUBLIC:: MAXVELP
   ! NORMALIZATIONS:
   PUBLIC:: RENORMALIZE, REWIPE
   PUBLIC:: RICH_RE, HYPERV3
CONTAINS
!=======================================================================
!============================ SUBROUTINES ==============================
!=======================================================================
SUBROUTINE MONITOR_EIG(PSI0,CHI0,PSI1,CHI1,TIME0,TIME1)
! ======================================================================
!> COMPUTATION OF EIGENVALUE AND MEASURE OF SHAPE PRESERVATION.
!> J. BARRANCO, FEBRURARY 2000 
! ======================================================================
    IMPLICIT NONE
    TYPE(SCALAR):: PSI0,CHI0,PSI1,CHI1
    REAL(P8):: TIME0,TIME1
    INTEGER:: STATUS

    TYPE(SCALAR):: P0,P1,C0,C1
    INTEGER:: NUM, I, NN, MM, KK, J
    REAL(P8):: DT, L2NP, LINFP, L2NC, LINFC, FAC
    REAL(P8), ALLOCATABLE, DIMENSION(:):: W, STANDARDDEVIATION
    COMPLEX(P8), ALLOCATABLE, DIMENSION(:):: EIGVALS, EIGVALS1,&
         &F1, TEMP, TEMP1
    COMPLEX(P8), PARAMETER:: II=(0._P8,1._P8) 
    COMPLEX(P8), ALLOCATABLE, DIMENSION(:,:):: TEMP_TBL
    COMPLEX(P8), DIMENSION(NDIMR,NTCHOPDIM,NXCHOPDIM):: P0_GLB, P1_GLB, C0_GLB, C1_GLB

    IF (.NOT.ALLOCATED(MONITOR_MK)) RETURN

    NUM = SIZE(MONITOR_MK,1)
    ALLOCATE(EIGVALS(NUM), EIGVALS1(NUM),F1(NR),TEMP(128),TEMP1(128),&
         & STANDARDDEVIATION(NUM), W(NR))
    CALL ALLOCATE(P0)
    CALL ALLOCATE(P1)
    CALL ALLOCATE(C0)
    CALL ALLOCATE(C1)
    P0 = PSI0
    P1 = PSI1
    C0 = CHI0
    C1 = CHI1
    P0%LN=0._P8
    P1%LN=0._P8
    C0%LN=0._P8
    C1%LN=0._P8
    CALL RTRAN(P0,1) ! PFF SPACE
    CALL RTRAN(P1,1)
    CALL RTRAN(C0,1)
    CALL RTRAN(C1,1)
    
    DT = TIME1-TIME0
    
    EIGVALS = 0.D0
    EIGVALS1 = 0.D0
    STANDARDDEVIATION = 0.D0
    ALLOCATE(TEMP_TBL(8,NUM))
    TEMP_TBL = 0.D0
    DO I=1,NUM
       IF (MONITOR_MK(I,1)<0) MONITOR_MK(I,:) = -MONITOR_MK(I,:)
       MM = MONITOR_MK(I,1)+1 ! INDEX
       KK = MONITOR_MK(I,2)+1 ! INDEX
       !IF (KK<1) KK=KK+NXCHOPDIM
       IF (KK<1) KK=KK+2*NXCHOP-1
       IF (MM==1) THEN
          NN=2 ! AVOID M=0,K=0,N=0
       ELSE
          NN=1
       END IF
       
       IF ((PSI1%INTH.LT.MM).AND.(PSI1%INTH+SIZE(PSI1%E,2).GE.MM)) THEN
         IF ((PSI1%INX.LT.KK).AND.(PSI1%INX+SIZE(PSI1%E,3).GE.KK)) THEN
            ! Original:FFF
            !  TEMP1(1:8) = II*LOG(PSI1%E(NN:NN+7,MM,KK)/PSI0%E(NN:NN+7,MM,KK))/DT ! used EXP(-II*Sigma*t)
            !  TEMP(1:8) = II*LOG(CHI1%E(NN:NN+7,MM,KK)/CHI0%E(NN:NN+7,MM,KK))/DT ! used EXP(-II*Sigma*t)
            TEMP1(1:8) = LOG(PSI1%E(NN:NN+7,MM-PSI1%INTH,KK-PSI1%INX) & 
                            /PSI0%E(NN:NN+7,MM-PSI1%INTH,KK-PSI1%INX))/DT
            TEMP(1:8) = LOG(CHI1%E(NN:NN+7,MM-PSI1%INTH,KK-PSI1%INX) &
                            /CHI0%E(NN:NN+7,MM-PSI1%INTH,KK-PSI1%INX))/DT
            
            IF (SUM(CHI0%E(NN:NN+7,MM-PSI1%INTH,KK-PSI1%INX)).EQ.0) THEN
               WRITE(*,*) 'ALL-ZERO CHI: MM = ',MM,', KK = ',KK
               TEMP(1:8) = TEMP1(1:8)
            ENDIF
            ! WRITE(*,*) MM,KK
            ! CALL MCAT(PSI1%E(NN:NN+7,MM-PSI1%INTH,KK-PSI1%INX))
            ! WRITE(*,*) ''
            ! Edited:  PFF
            !TEMP1(1:8) = II*LOG(P1%E(NN:NN+7,MM,KK)/P0%E(NN:NN+7,MM,KK))/DT
            !TEMP(1:8) = II*LOG(C1%E(NN:NN+7,MM,KK)/C0%E(NN:NN+7,MM,KK))/DT
            ! Original:FFF ALL
            !EMP1 = II*LOG(PSI1%E(NN:,MM,KK)/PSI0%E(NN:,MM,KK))/DT
            !TEMP = II*LOG(CHI1%E(NN:,MM,KK)/CHI0%E(NN:,MM,KK))/DT  
            ! Edited:  PFF ALL
            !TEMP1 = II*LOG(P1%E(NN:,MM,KK)/P0%E(NN:,MM,KK))/DT
            !TEMP = II*LOG(C1%E(NN:,MM,KK)/C0%E(NN:,MM,KK))/DT        

            ! DEBUG:
            !WRITE(*,*) PSI1%E(NN:NN+7,MM,KK),';',PSI0%E(NN:NN+7,MM,KK)

            ! ! Note: below is trying to print eigenvalues of the first 8 collocation points,
            ! !       but both PSI and CHI are in FFF space right now.
            ! PRINT *,'MODE:',MM-1,KK-1
            ! PRINT *,'EIGENVALUES FOR FIRST 8 RADIAL COLLOCATION POINTS'
            ! PRINT *, TEMP(1:8)
            ! !PRINT *,'EIGENVALUES FOR ALL RADIAL COLLOCATION POINTS'
            ! !PRINT *, TEMP(:)

            TEMP_TBL(:,I) = TEMP(1:8)

            EIGVALS(I) = SUM(TEMP(1:8))/8
            EIGVALS1(I) = SUM(TEMP1(1:8))/8

            !CALL MSAVE(TEMP, 'TEMP')
            !STOP

            ! !NOW THAT EIGENVALUES ARE CALCULATED USING
            ! ! PSI AND CHI, SEE IF THEY ARE DIFFERENT.
            ! IF (ABS(EIGVALS(I)-EIGVALS1(I)).GE.1.0E-10) THEN
            !    PRINT *,'DIFFERENCE IN EIGENVALUES BASED ON CHI AND PSI!! &
            !          & FOR MODE:',MM-1,KK-1, &
            !          & EIGVALS(I)-EIGVALS1(I)
            !    PRINT*,'------------------------------------------';
            ! END IF
            
            !CALCULATE STD. DEV. BETWEEN RADIAL MODES. 
            !THE EIGENFUNCTION SHOULD GROW AT THE SAME RATE 
            !AT EACH RADIAL COLLOCATION POINT 
            STANDARDDEVIATION(I)=0.0_P8
            DO J=1,8
               STANDARDDEVIATION(I) = STANDARDDEVIATION(I) +&
                     & (TEMP(J)-EIGVALS(I))*(CONJG(TEMP(J)-EIGVALS(I)))
            END DO
            STANDARDDEVIATION(I)= (STANDARDDEVIATION(I)/8.0_P8)**0.5_P8
         ENDIF
      ENDIF
    END DO
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, EIGVALS, NUM, MPI_DOUBLE_COMPLEX, MPI_SUM, &
                     MPI_COMM_IVP, IERR)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, EIGVALS1, NUM, MPI_DOUBLE_COMPLEX, MPI_SUM, &
                     MPI_COMM_IVP, IERR)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, STANDARDDEVIATION, NUM, MPI_DOUBLE_PRECISION, MPI_SUM, &
                     MPI_COMM_IVP, IERR)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, TEMP_TBL, NUM*8, MPI_DOUBLE_COMPLEX, MPI_SUM, &
                     MPI_COMM_IVP, IERR)
    IF (MPI_RANK.EQ.0) THEN
      DO I = 1,NUM
         PRINT *,'MODE:',MONITOR_MK(I,1),MONITOR_MK(I,2)
            PRINT *,'EIGENVALUES FOR FIRST 8 RADIAL COLLOCATION POINTS'
            PRINT *, TEMP_TBL(1:8,I)
      ENDDO

      !NOW THAT EIGENVALUES ARE CALCULATED USING
      ! PSI AND CHI, SEE IF THEY ARE DIFFERENT.
      IF (ABS(EIGVALS(I)-EIGVALS1(I)).GE.1.0E-10) THEN
         PRINT *,'DIFFERENCE IN EIGENVALUES BASED ON CHI AND PSI!! &
               & FOR MODE:',MM-1,KK-1, &
               & EIGVALS(I)-EIGVALS1(I)
         PRINT*,'------------------------------------------';
      END IF

      PRINT *, '---------------------------------------------------------'
      PRINT *, 'T0,T1,DT= ',TIME0,TIME1,DT  

    ENDIF
                     
    ! PFF SPACE: P0, P1, C0, C1
    CALL MASSEMBLE(P0%E,P0_GLB,1)
    CALL MASSEMBLE(P1%E,P1_GLB,1)
    CALL MASSEMBLE(C0%E,C0_GLB,1)
    CALL MASSEMBLE(C1%E,C1_GLB,1)
    IF (MPI_RANK.EQ.0) THEN
      W = TFM%W*(ELL/(1._P8-TFM%X))**2
      W = TFM%W
      DO I=1,NUM
         MM = MONITOR_MK(I,1)+1
         KK = MONITOR_MK(I,2)+1
         IF (KK<1) KK=KK+2*NXCHOP-1 !+NX
         F1 = P0_GLB(1:NR,MM,KK)
         FAC = SQRT(DOT_PRODUCT(F1,F1*W))
         !  F1 = (P1%E(1:NR,MM,KK)*EXP(II*EIGVALS(I)*DT) &
         !       & -P0%E(1:NR,MM,KK))/FAC
         F1 = (P1_GLB(1:NR,MM,KK)*EXP(-EIGVALS(I)*DT) &
               & -P0_GLB(1:NR,MM,KK))/FAC
         L2NP = SQRT(DOT_PRODUCT(F1,F1*W))
         LINFP = MAXVAL(ABS(F1))
         F1 = C0_GLB(1:NR,MM,KK)
         FAC = SQRT(DOT_PRODUCT(F1,F1*W))
         !  F1 = (C1%E(1:NR,MM,KK)*EXP(II*EIGVALS(I)*DT) &
         !       & -C0%E(1:NR,MM,KK))/FAC
         F1 = (C1_GLB(1:NR,MM,KK)*EXP(-EIGVALS(I)*DT) &
               & -C0_GLB(1:NR,MM,KK))/FAC
         L2NC = SQRT(DOT_PRODUCT(F1,F1*W))
         LINFC = MAXVAL(ABS(F1))
         !--
         OPEN(UNIT=2000,FILE='eigenValueData.dat',STATUS='UNKNOWN',&
               & IOSTAT=STATUS, POSITION='APPEND')
   700    FORMAT(3E23.15)
         WRITE(2000, 700) TIME1, EIGVALS(I)
         CLOSE(2000)
         !--
         PRINT *, 'FOR THE MODE(M,K)   = ', MONITOR_MK(I,1), MONITOR_MK(I,2)
         PRINT *, 'AVERAGE EIGENVALUE  = ', EIGVALS(I)
         PRINT *, 'STANDARD DEVIATION  = ', STANDARDDEVIATION(I)
         PRINT *, 'PSI L_2 NORM        = ', L2NP
         PRINT *, 'PSI L_INFINITY      = ', LINFP
         PRINT *, 'CHI L_2 NORM        = ', L2NC
         PRINT *, 'CHI L_INFINITY      = ', LINFC
      END DO
      PRINT *, '---------------------------------------------------------'
      ENDIF

    DEALLOCATE(EIGVALS,F1,W,TEMP_TBL)
    CALL DEALLOCATE(P1)
    CALL DEALLOCATE(P0)
    CALL DEALLOCATE(C0)
    CALL DEALLOCATE(C1)
END SUBROUTINE MONITOR_EIG
! ======================================================================
  
SUBROUTINE EGROWTH(PSI0,CHI0,PSI,CHI,TIME0,TIME1)
! ======================================================================
   IMPLICIT NONE
   TYPE(SCALAR):: PSI, CHI, PSI0, CHI0
   REAL(P8):: TIME0,TIME1
   REAL(P8), ALLOCATABLE, DIMENSION(:,:):: ESPEC, ESPEC0
   REAL(P8):: DT
   REAL(P8), PARAMETER:: TOL=1.0E-16
   INTEGER:: MM,KK
   !> CALCULATION OF GROWTH RATES FROM ENERGY SPECTRUM.
   !> J. BARRANCO, FEBRUARY 2000
   
   DT = TIME1-TIME0
   ALLOCATE(ESPEC(NTCHOP,NXCHOP),ESPEC0(NTCHOP,NXCHOP))
   ESPEC  = ENERGY_SPEC(PSI,CHI)
   ESPEC0 = ENERGY_SPEC(PSI0,CHI0)
   
   IF (MPI_RANK.EQ.0) THEN
      PRINT *, '--------------------------------------------------'
      PRINT *, 'T0, T1, DT  =  ', TIME0, TIME1, DT
      WHERE((ESPEC0>TOL).AND.(ESPEC>TOL))
         ESPEC=LOG(ESPEC/ESPEC0)/2.0_P8/DT
      ELSEWHERE
         ESPEC=0.0_P8
      END WHERE
      PRINT *, 'GROWTH RATES G(M,K):'
      PRINT '(7E16.8)' , ((ESPEC(MM,KK),KK=1,7),MM=1,7)
      PRINT *, '--------------------------------------------------'
   ENDIF

   DEALLOCATE(ESPEC,ESPEC0)

END SUBROUTINE EGROWTH



SUBROUTINE EGROWTH_MODIFIED(PSI0,CHI0,PSI,CHI,TIME0,TIME1)
! ======================================================================
   IMPLICIT NONE
   TYPE(SCALAR):: PSI, CHI, PSI0, CHI0
   REAL(P8):: TIME0,TIME1
   REAL(P8), ALLOCATABLE, DIMENSION(:,:):: ESPEC, ESPEC0
   REAL(P8):: DT
   REAL(P8), PARAMETER:: TOL=1.0E-16
   INTEGER::MM,KK

   DT = TIME1-TIME0
   !ALLOCATE(ESPEC(NTCHOP,NX),ESPEC0(NTCHOP,NX))
   ALLOCATE(ESPEC(NTCHOP,NXCHOPDIM),ESPEC0(NTCHOP,NXCHOPDIM))
   ESPEC  = ENERGY_SPEC_MODIFIED(PSI,CHI)
   ESPEC0 = ENERGY_SPEC_MODIFIED(PSI0,CHI0)

   IF (MPI_RANK.EQ.0) THEN
      PRINT *, '--------------------------------------------------'
      PRINT *, 'T0, T1, DT  =  ', TIME0, TIME1, DT
      WHERE((ESPEC0>TOL).AND.(ESPEC>TOL))
         ESPEC=LOG(ESPEC/ESPEC0)/2.0_P8/DT
      ELSEWHERE
         ESPEC=0.0_P8
      END WHERE
      PRINT *, 'GROWTH RATES G(M,K=0,K=1,K=2...):'
      PRINT '(7E16.8)' , ((ESPEC(MM,KK),KK=1,7),MM=1,7)
      PRINT *,'---'
      PRINT *, 'GROWTH RATES G(M,K=-1,K=-2,K=-3...):'
      !PRINT '(7E16.8)' , ((ESPEC(MM,NX+1-KK),KK=1,7),MM=1,7)
      PRINT '(7E16.8)' , ((ESPEC(MM,2*NXCHOP-KK),KK=1,7),MM=1,7)
      PRINT *, '--------------------------------------------------'
   ENDIF

   DEALLOCATE(ESPEC,ESPEC0)
END SUBROUTINE EGROWTH_MODIFIED
! ======================================================================
  
SUBROUTINE PRINT_ENERGY_SPECTRUM(PSI,CHI,TAG)
! ======================================================================
    IMPLICIT NONE
    TYPE(SCALAR):: PSI, CHI
    INTEGER::TAG
    REAL(P8), ALLOCATABLE, DIMENSION(:,:):: ESPEC
    INTEGER::MM,KK,FILESTATUS
    REAL(P8)::TIME
    
    IF (TAG.EQ.0) THEN
       !PRINT SPECTRUM FOR POSITIVE M'S AND K'S ONLY
       
       ALLOCATE(ESPEC(NTCHOP,NXCHOP))
       ESPEC  = ENERGY_SPEC(PSI,CHI)
       IF (MPI_RANK.EQ.0) THEN
         PRINT *, 'ENERGY SPECTRUM  E(M,K=0,K=1,K=2...):'
         WRITE(*,FMT=111) ((ESPEC(MM,KK),KK=1,7),MM=1,7)
       ENDIF
       
    ELSE IF (TAG.EQ.1) THEN
       !PRINT SPECTRUM FOR ALL K'S BUT POSITIVE M'S ONLY

       ALLOCATE(ESPEC(NTCHOP,NXCHOPDIM)) !NX))
       ESPEC  = ENERGY_SPEC_MODIFIED(PSI,CHI)
       IF (MPI_RANK.EQ.0) THEN
         PRINT *, 'ENERGY SPECTRUM  E(M,K=0,K=1,K=2...):'
         WRITE(*,FMT=111) ((ESPEC(MM,KK),KK=1,7),MM=1,7)
         PRINT *,'---'
         PRINT *, 'ENERGY SPECTRUM  E(M,K=-1,K=-2,K=-3...):'
         WRITE(*,FMT=111) ((ESPEC(MM,NXCHOPDIM+1-KK),KK=1,7),MM=1,7)
       ENDIF
    END IF

    IF (MPI_RANK.EQ.0) THEN
      OPEN(UNIT=666,FILE='especData_M.dat',STATUS='UNKNOWN',&
      & IOSTAT=FILESTATUS, POSITION='APPEND')
      TIME = TIM%T
      !  WRITE(666, 787) TIME , (ESPEC(MM,1),MM=1,NRCHOP)        ! Only record m = 0 mode
      WRITE(666, 787) TIME , (SUM(ESPEC(MM,:)),MM=1,NTCHOP)          
      CLOSE(666)     

      OPEN(UNIT=667,FILE='especData_K.dat',STATUS='UNKNOWN',&
      & IOSTAT=FILESTATUS, POSITION='APPEND')
      TIME = TIM%T
      !  WRITE(667, 788) TIME , (ESPEC(1,KK),KK=1,NXCHOP)        ! Only record k = 0 mode
      WRITE(667, 788) TIME , (SUM(ESPEC(:,KK)),KK=1,NXCHOPDIM)              
      CLOSE(667)  
    ENDIF

    DEALLOCATE(ESPEC)
111 FORMAT(7E14.6)
787 FORMAT(100E24.16)
788 FORMAT(100E24.16)

END SUBROUTINE PRINT_ENERGY_SPECTRUM
! ======================================================================

SUBROUTINE WRITE_ENERGY_DATA_TO_FILE(PSI,CHI)
! ======================================================================
   IMPLICIT NONE
   TYPE(SCALAR):: PSI, CHI
   REAL(P8), ALLOCATABLE, DIMENSION(:,:):: ESPEC
   REAL(P8):: TIME
   INTEGER::I, MM, KK, FILESTATUS,A,B
   INTEGER, ALLOCATABLE, DIMENSION(:,:):: TEMP 

   A = SIZE(MONITOR_MK,1)
   B = SIZE(MONITOR_MK,2)

   ALLOCATE( TEMP(A,B) )
   TEMP = MONITOR_MK

   ALLOCATE(ESPEC(NTCHOP,NXCHOPDIM)) !NX))
   !WRITE(*,*) 'oh yeah'
   ESPEC  = ENERGY_SPEC_MODIFIED(PSI,CHI)
   
   IF (MPI_RANK.EQ.0) THEN
      TIME = TIM%T
      OPEN(UNIT=555,FILE='energyData.dat',STATUS='UNKNOWN',&
            & IOSTAT=FILESTATUS, POSITION='APPEND')
      
      DO I=1,SIZE(TEMP,1)
         IF (TEMP(I,1)<0) TEMP(I,:) = -TEMP(I,:)
         TEMP(I,1) = TEMP(I,1)+1
         TEMP(I,2) = TEMP(I,2)+1
         IF (TEMP(I,2)<1) TEMP(I,2)=NXCHOPDIM+TEMP(I,2) !TEMP(I,2)=NX+TEMP(I,2)
      END DO

      WRITE(555, 786) TIME , &
            & ESPEC( TEMP(1,1),TEMP(1,2) ), &
            & ESPEC( TEMP(2,1),TEMP(2,2) ), &
            & ESPEC( TEMP(3,1),TEMP(3,2) ), &
            & ESPEC( TEMP(4,1),TEMP(4,2) )         
      CLOSE(555)    

786 FORMAT(5E23.15)
   ENDIF

   DEALLOCATE(ESPEC)
   DEALLOCATE(TEMP)

END SUBROUTINE WRITE_ENERGY_DATA_TO_FILE
! ======================================================================

SUBROUTINE WRITE_MODE_DATA_TO_FILE(PSI,CHI,SPACE_FLAG)
! ======================================================================
! [USAGE]:
! WRITE THE MONITORED MODE'S DATA TO FILE. 
! OUTPUT SPACE IS GIVEN BY SPACE_FLAG
! INPUT MUST BE EITHER IN FFF_SPACE OR PFF_SPACE
! ======================================================================
   IMPLICIT NONE
   TYPE(SCALAR):: PSI, CHI
   TYPE(SCALAR):: PSI_COPY, CHI_COPY
   REAL(P8):: TIME
   INTEGER::I, MM, KK, FILESTATUS, A, B, SPACE_FLAG, R_SIZE, JJ
   INTEGER, ALLOCATABLE, DIMENSION(:,:):: TEMP
   COMPLEX(P8), DIMENSION(:,:,:),ALLOCATABLE:: PSI_GLB, CHI_GLB

   IF ((PSI%SPACE.NE.FFF_SPACE).AND.(PSI%SPACE.NE.PFF_SPACE)) STOP
   IF ((SPACE_FLAG.NE.FFF_SPACE).AND.(SPACE_FLAG.NE.PFF_SPACE)) STOP

   IF (PSI%SPACE .EQ. SPACE_FLAG) THEN
      ALLOCATE(PSI_GLB(SIZE(PSI%E,1),NTCHOPDIM,NXCHOPDIM))
      ALLOCATE(CHI_GLB(SIZE(PSI%E,1),NTCHOPDIM,NXCHOPDIM))
      CALL MASSEMBLE(PSI%E,PSI_GLB,1)
      CALL MASSEMBLE(CHI%E,CHI_GLB,1)
   ELSE
      IF (PSI%SPACE .EQ. FFF_SPACE) THEN
         CALL ALLOCATE(PSI_COPY, FFF_SPACE)
         CALL ALLOCATE(CHI_COPY, FFF_SPACE)
         PSI_COPY = PSI
         CHI_COPY = CHI
         CALL RTRAN(PSI_COPY,1)
         CALL RTRAN(CHI_COPY,1)
      ELSE
         CALL ALLOCATE(PSI_COPY, PFF_SPACE)
         CALL ALLOCATE(CHI_COPY, PFF_SPACE)
         PSI_COPY = PSI
         CHI_COPY = CHI
         CALL RTRAN(PSI_COPY,-1)
         CALL RTRAN(CHI_COPY,-1)
      ENDIF
      ALLOCATE(PSI_GLB(SIZE(PSI_COPY%E,1),NTCHOPDIM,NXCHOPDIM))
      ALLOCATE(CHI_GLB(SIZE(CHI_COPY%E,1),NTCHOPDIM,NXCHOPDIM))
      CALL MASSEMBLE(PSI_COPY%E,PSI_GLB,1)
      CALL MASSEMBLE(CHI_COPY%E,CHI_GLB,1)
      CALL DEALLOCATE(PSI_COPY)
      CALL DEALLOCATE(CHI_COPY)
   ENDIF 

   IF (MPI_RANK.EQ.0) THEN
      A = SIZE(MONITOR_MK,1)
      B = SIZE(MONITOR_MK,2)

      R_SIZE = SIZE(PSI_GLB,1)
      ! R_SIZE = R_SIZE/2 ! ONLY SAVE THE FIRST HALF

      ALLOCATE( TEMP(A,B) )
      TEMP = MONITOR_MK

      TIME = TIM%T

      DO I=1,SIZE(TEMP,1)
         IF (TEMP(I,1)<0) TEMP(I,:) = -TEMP(I,:)
         TEMP(I,1) = TEMP(I,1)+1
         TEMP(I,2) = TEMP(I,2)+1
         IF (TEMP(I,2)<1) TEMP(I,2)=2*NXCHOP+TEMP(I,2)-1
      END DO

      DO I=1,A
         OPEN(UNIT=777,FILE=TRIM(ADJUSTL(FILES%SAVEDIR))//'modeData_M_'//ITOA3(TEMP(I,1))&
               & //'_K_'//ITOA3(TEMP(I,2))//'.dat', STATUS='UNKNOWN',&
               & IOSTAT=FILESTATUS, POSITION='APPEND')
         
            !TRIM(ADJUSTL(FILES%SAVEDIR))//...
         WRITE(777, 987) TIME , (PSI_GLB(JJ,TEMP(I,1),TEMP(I,2)),JJ=1,R_SIZE), (CHI_GLB(JJ,TEMP(I,1),TEMP(I,2)),JJ=1,R_SIZE)
         ! WRITE(777, 987) TIME , (CHI_GLB(JJ,TEMP(I,1),TEMP(I,2)),JJ=1,R_SIZE)

         ! WRITE(777, 986) TIME , &
         !       & PSI_GLB(I, TEMP(1,1),TEMP(1,2) ), &
         !       & CHI_GLB(I, TEMP(1,1),TEMP(1,2) ), &
         !       & PSI_GLB(I, TEMP(2,1),TEMP(2,2) ), &
         !       & CHI_GLB(I, TEMP(2,1),TEMP(2,2) ), &
         !       & PSI_GLB(I, TEMP(3,1),TEMP(3,2) ), &
         !       & CHI_GLB(I, TEMP(3,1),TEMP(3,2) ), &
         !       & PSI_GLB(I, TEMP(4,1),TEMP(4,2) ), &  
         !       & CHI_GLB(I, TEMP(4,1),TEMP(4,2) )   
               
         CLOSE(777)       
      END DO

      ! 986 FORMAT(1E10.3,16E23.15)
      987 FORMAT(1E10.3,*(', ',S,E14.6E3,SP,E14.6E3,'i'))

      DEALLOCATE(TEMP)
   ENDIF
   DEALLOCATE(PSI_GLB,CHI_GLB)

END SUBROUTINE WRITE_MODE_DATA_TO_FILE
! ======================================================================

SUBROUTINE COLLOC_INFO()
!=======================================================================
   IMPLICIT NONE

      ! SAVE THE RADIAL COLLOCATION POINTS (R)
      CALL MSAVE(TFM%R, TRIM(ADJUSTL(FILES%SAVEDIR))//&
      'r_colloc_pts.dat')

      ! SAVE THE X COLLOCATION POINTS (X)
      CALL MSAVE(TFM%X, TRIM(ADJUSTL(FILES%SAVEDIR))//&
         'x_colloc_pts.dat')

      ! SAVE THE GAUSS_LEGENDRE WEIGHTS (W)
      CALL MSAVE(TFM%W, TRIM(ADJUSTL(FILES%SAVEDIR))//&
         'gau_leg_weights.dat')

      ! SAVE THE AZIMUTHAL COLLOCATION POINTS (THETA)
      CALL MSAVE(TFM%TH, TRIM(ADJUSTL(FILES%SAVEDIR))//&
         't_colloc_pts.dat')

      ! SAVE THE AXIAL COLLOCATION POINTS (Z)
      CALL MSAVE(TFM%Z, TRIM(ADJUSTL(FILES%SAVEDIR))//&
         'z_colloc_pts.dat')

   RETURN
END SUBROUTINE COLLOC_INFO

!=======================================================================
!============================ FUNCTIONS ================================
!=======================================================================
FUNCTION PRODUCT_MK(A,B)
! ======================================================================
!  J. BARRANCO, 01/22/1999
!  FOR EACH (M,K), COMPUTE INTEGRAL OF
!  (1-X^2)*A*B RDR D(PHI) DZ OVER THE ENTIRE DOMAIN.
!  CALL IN R-PHYSICAL, PHI/Z-FOURIER SPACE.
! ======================================================================
   IMPLICIT NONE
   TYPE(SCALAR):: A,B
   REAL(P8):: WK1(NR)
   REAL(P8):: PRODUCT_MK(NTCHOP,NXCHOP)
   INTEGER:: MM,KK,KC,COUNT
   
   IF(A%SPACE.NE.PFF_SPACE .OR. B%SPACE.NE.PFF_SPACE) THEN
   IF (MPI_RANK.EQ.0) THEN
      PRINT *,'ERROR: PRODUCT_MK() -- NOT IN PFF_SPACE.'
      PRINT *,'A%SPACE, B%SPACE = ',A%SPACE,B%SPACE
   ENDIF
   STOP
   ENDIF

   IF ((A%INTH.EQ.0).AND.(A%INX.EQ.0)) THEN
      IF(ABS(A%LN)>1.0E-8 .OR. ABS(B%LN)>1.0E-8) THEN
         PRINT *,'WARNING: PRODUCT_MK() -- LOGTERM NOT ZERO.'
         PRINT *,'A%LN, B%LN = ',A%LN,B%LN
      ENDIF
   ENDIF

   PRODUCT_MK = 0.D0
   DO MM=1,SIZE(A%E,2) !NTCHOP
      DO KK=1,SIZE(A%E,3) !NXCHOP
         KC = KK+A%INX ! GLB_COORDINATE
         IF (KC==1) THEN
            WK1 = REAL(A%E(:NR,MM,KK)*CONJG(B%E(:NR,MM,KK)) &
                    & +B%E(:NR,MM,KK)*CONJG(A%E(:NR,MM,KK)))
         ELSEIF (KC.LE.NXCHOP) THEN
            ! WK1 = REAL(A%E(:NR,MM,KK)*CONJG(B%E(:NR,MM,KK)) &
            !         & +B%E(:NR,MM,KC)*CONJG(A%E(:NR,MM,KC)))
            WK1 = REAL(A%E(:NR,MM,KK)*CONJG(B%E(:NR,MM,KK)))
         ELSE
            KC = NXCHOPDIM+2-KC
            WK1 = REAL(B%E(:NR,MM,KK)*CONJG(A%E(:NR,MM,KK)))
         ENDIF
         PRODUCT_MK(MM+A%INTH,KC) = PRODUCT_MK(MM+A%INTH,KC) &
                     +4.0_P8*PI*ZLEN0*ELL2*DOT_PRODUCT(WK1,TFM%W)
      END DO
   END DO
   COUNT = NTCHOP*NXCHOP
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, PRODUCT_MK, COUNT, &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_IVP, IERR)
   PRODUCT_MK(:,1) = 0.5_P8*PRODUCT_MK(:,1)
   PRODUCT_MK(1,:) = 0.5_P8*PRODUCT_MK(1,:)

   RETURN
END FUNCTION PRODUCT_MK
! ======================================================================
  
FUNCTION PRODUCT_MK_MODIFIED(A,B)
! ======================================================================
IMPLICIT NONE
TYPE(SCALAR):: A,B
REAL(P8):: WK1(NR)
REAL(P8):: PRODUCT_MK_MODIFIED(NTCHOP,NXCHOPDIM) !NX)
INTEGER:: MM,KK,KC,COUNT

   IF(A%SPACE.NE.PFF_SPACE .OR. B%SPACE.NE.PFF_SPACE) THEN
      IF (MPI_RANK.EQ.0) THEN
         PRINT *,'ERROR: PRODUCT_MK() -- NOT IN PFF_SPACE.'
         PRINT *,'A%SPACE, B%SPACE = ',A%SPACE,B%SPACE
      ENDIF
      STOP
   ENDIF

   IF ((A%INTH.EQ.0).AND.(A%INX.EQ.0)) THEN
      IF(ABS(A%LN)>1.0E-8 .OR. ABS(B%LN)>1.0E-8) THEN
         PRINT *,'WARNING: PRODUCT_MK() -- LOGTERM NOT ZERO.'
         PRINT *,'A%LN, B%LN = ',A%LN,B%LN
      ENDIF
   ENDIF

   PRODUCT_MK_MODIFIED = 0.D0
   DO MM=1,SIZE(A%E,2) !NTCHOP
      DO KK=1,SIZE(A%E,3) !NXCHOP
         !KC = NX-KK+2
         WK1 = REAL( A%E(:NR,MM,KK)*CONJG(B%E(:NR,MM,KK)) )

         ! {MM,KK}&{-MM,-KK} combined
         ! note: We save non-negative M's, but all K's.
         !       So, for {0,KK}, we should not add {0,-KK} during the calculation,
         !       as {0,-KK} is actually saved. Hence, all MM = 0 modes' product is
         !       off by a factor of 2.
         PRODUCT_MK_MODIFIED(MM+A%INTH,KK+A%INX) = &
                  4.0_P8*PI*ZLEN0*ELL2*DOT_PRODUCT(WK1,TFM%W)

         ! IF (KK.NE.1) THEN
         !    WK1 = REAL( A%E(:NR,MM,KC)*CONJG(B%E(:NR,MM,KC)) ) 
         !    PRODUCT_MK_MODIFIED(MM,KC) = 4.0_P8*PI*ZLEN0*ELL2* &
         !       & DOT_PRODUCT(WK1,TFM%W)
         ! END IF
      END DO
   END DO

   COUNT = NTCHOP*NXCHOPDIM
   CALL MPI_ALLREDUCE(MPI_IN_PLACE, PRODUCT_MK_MODIFIED, COUNT, &
                      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_IVP, IERR)
   PRODUCT_MK_MODIFIED(1,:) = 0.5_P8*PRODUCT_MK_MODIFIED(1,:)

   RETURN
END FUNCTION PRODUCT_MK_MODIFIED
! ======================================================================
  
FUNCTION ENERGY_SPEC(PSI,CHI)
!  ------------------------------------------------------
!  J. BARRANCO, 01/22/1999
!  COMPUTES ENERGY IN EACH MODE (M,K).
!  CALL IN R-FUNCTION, PHI/Z-FOURIER SPACE.
!  ------------------------------------------------------    
   IMPLICIT NONE 
   TYPE(SCALAR):: PSI,CHI,W1,W2
   REAL(P8):: ENERGY_SPEC(NTCHOP,NXCHOP)

   CALL ALLOCATE(W1)
   CALL ALLOCATE(W2)
   CALL CHOPSET(3)

   ! ENERGY IN TOROIDAL COMPONENT
   W1 = PSI
   IF ((PSI%INX.EQ.0).AND.(PSI%INTH.EQ.0)) W1%LN = 2.0_P8*PSI%LN
   CALL DELSQH(W1,W2)
   W1%LN = 0.0_P8
   W2%LN = 0.0_P8
   CALL RTRAN(W1,1)
   CALL RTRAN(W2,1)
   ENERGY_SPEC = -PRODUCT_MK(W1,W2)

   ! ENERGY IN POLOIDAL COMPONENT
   CALL DEL2(CHI,W1)
   CALL DELSQH(CHI,W2)
   W1%LN = 0.0_P8
   W2%LN = 0.0_P8
   CALL RTRAN(W1,1)
   CALL RTRAN(W2,1)
   ENERGY_SPEC = ENERGY_SPEC + PRODUCT_MK(W1,W2)
   CALL CHOPSET(-3)
   CALL DEALLOCATE(W1)
   CALL DEALLOCATE(W2)
   
END FUNCTION ENERGY_SPEC
! ======================================================================

FUNCTION ENERGY_SPEC_MODIFIED(PSI,CHI)
! ======================================================================
    IMPLICIT NONE
    TYPE(SCALAR):: PSI,CHI,W1,W2
    REAL(P8):: ENERGY_SPEC_MODIFIED(NTCHOP,NXCHOPDIM) !NX)
   !  TYPE(SCALAR):: RUR,RUP,UZ,RUR2,RUP2,UZ2,UP,RUPXX
    
    CALL ALLOCATE(W1)
    CALL ALLOCATE(W2)
    CALL CHOPSET(3)
    
    ! ENERGY IN TOROIDAL COMPONENT
    W1 = PSI
    IF ((PSI%INX.EQ.0).AND.(PSI%INTH.EQ.0)) W1%LN = 2.0_P8*PSI%LN
    CALL DELSQH(W1,W2)
    W1%LN = 0.0_P8
    W2%LN = 0.0_P8
    CALL RTRAN(W1,1)
    CALL RTRAN(W2,1)
    ENERGY_SPEC_MODIFIED = -PRODUCT_MK_MODIFIED(W1,W2)
    
    ! ENERGY IN POLOIDAL COMPONENT
    CALL DEL2(CHI,W1)
    CALL DELSQH(CHI,W2)
    W1%LN = 0.0_P8
    W2%LN = 0.0_P8
    CALL RTRAN(W1,1)
    CALL RTRAN(W2,1)
    ENERGY_SPEC_MODIFIED = ENERGY_SPEC_MODIFIED + PRODUCT_MK_MODIFIED(W1,W2)
    CALL CHOPSET(-3)
    CALL DEALLOCATE(W1)
    CALL DEALLOCATE(W2)
    
   ! CALL ALLOCATE(UP)
   ! CALL ALLOCATE(RUR)
   ! CALL ALLOCATE(RUP)
   ! CALL ALLOCATE(UZ)
   ! CALL PC2VEL(PSI,CHI,RUR,RUP,UZ)
   ! CALL MULXMDIVXP(RUP,UP,.TRUE.)
   ! CALL RTRAN(UP,1)
   ! CALL RTRAN(RUR,1)
   ! CALL RTRAN(RUP,1)
   ! CALL RTRAN(UZ,1)

   ! CALL ALLOCATE(RUPXX)
   ! CALL ALLOCATE(RUR2)
   ! CALL ALLOCATE(RUP2)
   ! CALL ALLOCATE(UZ2)
   ! CALL PC2VEL(PSI,CHI,RUR2,RUP2,UZ2)
   ! CALL DIVXM(RUP2,RUPXX)
   ! CALL DIVXM(RUPXX,RUP2)

   ! CALL RTRAN(RUR2,1)
   ! CALL RTRAN(RUP2,1)
   ! CALL RTRAN(UZ2,1)

   ! ! WRITE(*,*) RUP%E(:NR,1,1)/TFM%R
   ! ! WRITE(*,*) SUM(PRODUCT_MK(RUR,RUR))
   ! WRITE(*,*) SUM(PRODUCT_MK(UP,RUP2))
   ! ! WRITE(*,*) SUM(PRODUCT_MK(UZ,UZ))
   ! CALL DEALLOCATE(UP)
   ! CALL DEALLOCATE(RUR)
   ! CALL DEALLOCATE(RUP)
   ! CALL DEALLOCATE(UZ)
   ! CALL DEALLOCATE(RUR2)
   ! CALL DEALLOCATE(RUP2)
   ! CALL DEALLOCATE(UZ2)
   ! CALL DEALLOCATE(RUPXX)

END FUNCTION ENERGY_SPEC_MODIFIED
! ======================================================================

SUBROUTINE RENORMALIZE(PSI0,CHI0,PSI1,CHI1,PSIi,CHIi,veli)
! ======================================================================
! [USAGE]: 
! RENORMALIZE THE CURRENT TIME STEP:
!   1. {0,0} mode back to original (PSIi,CHIi)
!   2. {m_bar,k_bar} mode back to original
!   3. {m_1,k_1} mode refactored to have a max|UP| of 1 (no phase change)
!   4. {m_2,k_2} mode refactored according to the same factor for {m_1,k_1} 
! [PARAMETERS]:
! PSI0 >> CURRENT  TOROIDAL TERM 
! CHI0 >> CURRENT  POLOIDAL TERM 
! PSI1 >> TEMPORAR TOROIDAL TERM PRE-ALLOCATED
! CHI1 >> TEMPORAR POLOIDAL TERM PRE-ALLOCATED
! PSIi >> ORIGINAL TOROIDAL TERM (FFF) when monitor starts
! CHIi >> ORIGINAL POLOIDAL TERM (FFF)
! veli >> INITIAL MAXIMUM AZIMUTHAL VELOCITY
! [DEPENDENCIES]:
! 1. MAXVELP(~) @ MOD_MARCH
! [NOTE]:
! MONITOR_MK MUST BE CONSISTENT WITH THE PERTURBATION ADDED.
! MONITOR_MK{1,:} -> M_BAR, K_BAR
! MONITOR_MK{2,:} -> M_1, K_1
! MONITOR_MK{3,:} -> M_2, K_2
! [UPDATES]:
! CODED BY JINGE WANG @ MAR 11 2022
! ======================================================================
   IMPLICIT NONE
   TYPE(SCALAR):: PSI0,CHI0,PSI1,CHI1,PSIi,CHIi
   INTEGER:: MON_SIZE, MON_II, MM, KK
   REAL(P8):: vel0,veli

   ! WIPE PSI1&CHI1
   PSI1%ln = 0.D0 
   CHI1%ln = 0.D0
   PSI1%e = 0.D0
   CHI1%e = 0.D0

   ! OBTAIN MAXIMUM VELOCITY FOR V_{m1,k1}
   vel0 = MAXVELP(PSI0,CHI0)
   
   ! RENORMALIZE
   MON_SIZE = SIZE(MONITOR_MK,1)
   DO MON_II = 1,MIN(MON_SIZE,3)
      ! NOTE: THE LINE BELOW IS KINDA INCORRECT. FOR {-M,-K}, {M,K} SHOULD BE MONITORED. HOWEVER THE CODE CHANGES 
      !       IT TO {M,-K}. ALSO, THE LINE CHANGES MONITOR_MK, WHICH IS A GLOBAL VARIABLE -> THIS IS BAD.
      !       CURRENTLY, ALL 3 ADDED MODES ARE POSITIVE, SO THIS IS FINE FOR NOW.
      IF (MONITOR_MK(MON_II,1)<0) MONITOR_MK(MON_II,:) = -MONITOR_MK(MON_II,:)
      MM = MONITOR_MK(MON_II,1)+1 ! INDEX
      KK = MONITOR_MK(MON_II,2)+1 ! INDEX
      IF (KK<1) KK=KK+2*NXCHOP-1

      IF ((PSI1%INTH.LT.MM).AND.(PSI1%INTH+SIZE(PSI1%E,2).GE.MM)) THEN
            IF ((PSI1%INX.LT.KK).AND.(PSI1%INX+SIZE(PSI1%E,3).GE.KK)) THEN
               IF (MON_II.GT.1) THEN 
               ! RENORMALIZE THE TESTED EIGENVECTOR: alpha*V_{m1,k1}+beta*V_{m2,k2} SUCH THAT V_{m1,k1} HAS THE 
               ! CORRECT MAXIMUM VELOCITY.
                  PSI1%E(:,MM-PSI1%INTH,KK-PSI1%INX) = PSI0%E(:,MM-PSI0%INTH,KK-PSI0%INX)*veli/vel0
                  CHI1%E(:,MM-CHI1%INTH,KK-CHI1%INX) = CHI0%E(:,MM-CHI0%INTH,KK-CHI0%INX)*veli/vel0
               ELSE ! REFIX THE BASE FLOW: V_{m_bar,k_bar}
                  PSI1%E(:,MM-PSI1%INTH,KK-PSI1%INX) = PSIi%E(:,MM-PSI0%INTH,KK-PSI0%INX)
                  CHI1%E(:,MM-CHI1%INTH,KK-CHI1%INX) = CHIi%E(:,MM-CHI0%INTH,KK-CHI0%INX)
               ENDIF
            ENDIF
      ENDIF
   ENDDO
   
   IF ((PSI1%INX.EQ.0).AND.(PSI1%INTH.EQ.0)) THEN ! REFIX THE BASE FLOW: V_{0,0}
      PSI1%E(:,1,1) = PSIi%e(:,1,1) 
      CHI1%E(:,1,1) = CHIi%e(:,1,1)
   ENDIF
   
   PSI0%E = PSI1%E
   CHI0%E = CHI1%E
   PSI0%LN = PSIi%LN
   CHI0%LN = CHIi%LN

   CALL CHOPDO(PSI0)
   CALL CHOPDO(CHI0)

   ! vel0 = MAXVELP(PSI0,CHI0)
   ! WRITE(*,*) MPI_RANK,",",vel0,'normalized'

   IF (MPI_RANK.EQ.0) WRITE(*,*) 'RENORMALIZED'
   call WRITE_ENERGY_DATA_TO_FILE(PSI0,CHI0)

END SUBROUTINE RENORMALIZE
! ======================================================================

SUBROUTINE REWIPE(PSI0,CHI0,PSI1,CHI1,PSIi,CHIi)
! ======================================================================
! [USAGE]: 
! RENORMALIZE THE CURRENT TIME STEP:
!   1. {0,0} mode back to original (PSIi,CHIi)
!   2. {m_bar,k_bar} mode back to original
!   3. {m_1,k_1} & {m_2,k_2} stay
!   4. all other modes are wiped
! [PARAMETERS]:
! PSI0 >> CURRENT  TOROIDAL TERM 
! CHI0 >> CURRENT  POLOIDAL TERM 
! PSI1 >> TEMPORAR TOROIDAL TERM PRE-ALLOCATED
! CHI1 >> TEMPORAR POLOIDAL TERM PRE-ALLOCATED
! PSIi >> ORIGINAL TOROIDAL TERM (FFF) when monitor starts
! CHIi >> ORIGINAL POLOIDAL TERM (FFF)
! [NOTE]:
! MONITOR_MK MUST BE CONSISTENT WITH THE PERTURBATION ADDED.
! MONITOR_MK{1,:} -> M_BAR, K_BAR
! MONITOR_MK{2,:} -> M_1, K_1
! MONITOR_MK{3,:} -> M_2, K_2
! [UPDATES]:
! CODED BY JINGE WANG @ MAR 11 2022
! ======================================================================
   IMPLICIT NONE
   TYPE(SCALAR):: PSI0,CHI0,PSI1,CHI1,PSIi,CHIi
   INTEGER:: MON_SIZE, MON_II, MM, KK

   ! WIPE PSI1&CHI1
   PSI1%ln = 0.D0 
   CHI1%ln = 0.D0
   PSI1%e = 0.D0
   CHI1%e = 0.D0
   
   ! RETURN & WIPE
   MON_SIZE = SIZE(MONITOR_MK,1)
   DO MON_II = 1,MIN(MON_SIZE,3)
      ! NOTE: THE LINE BELOW IS KINDA INCORRECT. FOR {-M,-K}, {M,K} SHOULD BE MONITORED. HOWEVER THE CODE CHANGES 
      !       IT TO {M,-K}. ALSO, THE LINE CHANGES MONITOR_MK, WHICH IS A GLOBAL VARIABLE -> THIS IS BAD.
      !       CURRENTLY, ALL 3 ADDED MODES ARE POSITIVE, SO THIS IS FINE FOR NOW.
      IF (MONITOR_MK(MON_II,1)<0) MONITOR_MK(MON_II,:) = -MONITOR_MK(MON_II,:)
      MM = MONITOR_MK(MON_II,1)+1 ! INDEX
      KK = MONITOR_MK(MON_II,2)+1 ! INDEX
      IF (KK<1) KK=KK+2*NXCHOP-1

      IF ((PSI1%INTH.LT.MM).AND.(PSI1%INTH+SIZE(PSI1%E,2).GE.MM)) THEN
            IF ((PSI1%INX.LT.KK).AND.(PSI1%INX+SIZE(PSI1%E,3).GE.KK)) THEN
               IF (MON_II.GT.1) THEN 
               ! KEEP THE MODES
                  PSI1%E(:,MM-PSI1%INTH,KK-PSI1%INX) = PSI0%E(:,MM-PSI0%INTH,KK-PSI0%INX)
                  CHI1%E(:,MM-CHI1%INTH,KK-CHI1%INX) = CHI0%E(:,MM-CHI0%INTH,KK-CHI0%INX)
               ELSE ! REFIX THE BASE FLOW: V_{m_bar,k_bar}
                  PSI1%E(:,MM-PSI1%INTH,KK-PSI1%INX) = PSIi%E(:,MM-PSI0%INTH,KK-PSI0%INX)
                  CHI1%E(:,MM-CHI1%INTH,KK-CHI1%INX) = CHIi%E(:,MM-CHI0%INTH,KK-CHI0%INX)
               ENDIF
            ENDIF
      ENDIF
   ENDDO
   
   IF ((PSI1%INX.EQ.0).AND.(PSI1%INTH.EQ.0)) THEN ! REFIX THE BASE FLOW: V_{0,0}
      PSI1%E(:,1,1) = PSIi%e(:,1,1) 
      CHI1%E(:,1,1) = CHIi%e(:,1,1)
   ENDIF
   
   PSI0%E = PSI1%E
   CHI0%E = CHI1%E
   PSI0%LN = PSIi%LN
   CHI0%LN = CHIi%LN

   CALL CHOPDO(PSI0)
   CALL CHOPDO(CHI0)

   ! vel0 = MAXVELP(PSI0,CHI0)
   ! WRITE(*,*) MPI_RANK,",",vel0,'normalized'

   IF (MPI_RANK.EQ.0) WRITE(*,*) 'RENORMALIZED'
   call WRITE_ENERGY_DATA_TO_FILE(PSI0,CHI0)

END SUBROUTINE REWIPE
! ======================================================================

SUBROUTINE RICH_RE(PSI,CHI,PSIN,CHIN,PSI1,CHI1,PSIi,CHIi,veli)
!=======================================================================
! [USAGE]: 
! UPDATE POLOIDAL-TOROIDAL TERMS OF THE VELOCITY FIELD BY 1 TIME STEP
! USING RICHARDSON EXTRAPOLATION WITH RENORMALIZATION
! [PARAMETERS]:
! PSI >> TOROIDAL TERM IN A SCALAR-TYPE VARIABLE
! CHI >> POLOIDAL TERM IN A SCALAR-TYPE VARIABLE
! PSIN >> NONLINEAR COMPONENT OF THE TOROIDAL TERM
! CHIN >> NONLINEAR COMPONENT OF THE POLOIDAL TERM
! [DEPENDENCIES]:
! 1. (DE)ALLOCATE(~) @ MOD_SCALAR3
! 2. NONLIN(~) @ MOD_LEGOPS
! 3. VISC1(~) @ MOD_MARCH
! 4. HYPERV(~) @ MOD_MARCH
! [UPDATES]:
! CODED BY JINGE WANG @ MAR 11 2022
!=======================================================================
   IMPLICIT NONE
   TYPE(SCALAR):: PSI,CHI,PSIN,CHIN
   TYPE(SCALAR):: PSI1,CHI1,PSIi,CHIi

   TYPE(SCALAR):: PSI2,CHI2
   TYPE(SCALAR):: PN,CN
   REAL(P8):: DT, HDT, veli

   DT = TIM%DT
   HDT = DT*0.5D0
   TIM%T = TIM%T + DT
   TIM%N = TIM%N + 1
   ADV%X = ADV%X + ADV%UX*DT
   ADV%Y = ADV%Y + ADV%UY*DT

   CALL ALLOCATE(PSI2)
   CALL ALLOCATE(CHI2)

   CALL ALLOCATE( PN )
   CALL ALLOCATE( CN )

   ! FIRST HALF-STEP
   CALL NONLIN(PSI,CHI,PSIN,CHIN)
   PSI2%E =PSI%E + (HDT*PSIN%E)
   CHI2%E =CHI%E + (HDT*CHIN%E)
   CALL REWIPE(PSI2,CHI2,PSI1,CHI1,PSIi,CHIi)
   CALL VISC1(PSI2,CHI2,HDT)
   CALL HYPERV(PSI2,CHI2,HDT)

   ! SECOND HALF-STEP
   CALL NONLIN(PSI2,CHI2,PN,CN)
   PSI2%E =PSI2%E + (HDT*PN%E)
   CHI2%E =CHI2%E + (HDT*CN%E)
   CALL REWIPE(PSI2,CHI2,PSI1,CHI1,PSIi,CHIi)
   CALL VISC1(PSI2,CHI2,HDT)
   CALL HYPERV(PSI2,CHI2,HDT)

   ! FULL STEP
   PSI%E =PSI%E +(DT)*PSIN%E
   CHI%E =CHI%E +(DT)*CHIN%E

   CALL VISC1(PSI,CHI,DT)
   CALL HYPERV(PSI,CHI,DT)

   PSI%E =(2.0D0*PSI2%E)-PSI%E
   CHI%E =(2.0D0*CHI2%E)-CHI%E
   CALL RENORMALIZE(PSI,CHI,PSI1,CHI1,PSIi,CHIi,veli)

   CALL DEALLOCATE( PSI2 )
   CALL DEALLOCATE( CHI2 )
   CALL DEALLOCATE( PN )
   CALL DEALLOCATE( CN )

   RETURN
END SUBROUTINE RICH_RE
!=======================================================================

FUNCTION MAXVELP(PSIi,CHIi)
! ======================================================================
! [USAGE]: 
! OBTAIN THE MAXIMUM -AZIMUTHAL- VELOCITY OF PSI&CHI
! [INPUTS]:
! PSIi >> 
! CHIi >>
! [OUTPUTS]:
! MAXVEL >> MAXIMUM AZIMUTHAL VELOCITY
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 20 2020
! ======================================================================
   IMPLICIT NONE
   TYPE(SCALAR):: PSIi,CHIi
   TYPE(SCALAR):: RURi,RUPi,UZi
   INTEGER:: MM,KK
   REAL(P8):: MAXVELP

   CALL ALLOCATE(RURi)
   CALL ALLOCATE(RUPi)
   CALL ALLOCATE( UZi)

   CALL CHOPSET(3)
   CALL PC2VEL(PSIi,CHIi,RURi,RUPi,UZi)
   CALL CHOPSET(-3)
   CALL RTRAN(RUPi,1)

   IF (MONITOR_MK(2,1)<0) MONITOR_MK(2,:) = -MONITOR_MK(2,:)
   MM = MONITOR_MK(2,1)+1 ! INDEX
   KK = MONITOR_MK(2,2)+1 ! INDEX
   IF (KK<1) KK=KK+2*NXCHOP-1

   IF ((PSIi%INTH.LT.MM).AND.(PSIi%INTH+SIZE(PSIi%E,2).GE.MM)) THEN
      IF ((PSIi%INX.LT.KK).AND.(PSIi%INX+SIZE(PSIi%E,3).GE.KK)) THEN
            MAXVELP = MAXVAL(ABS(RUPi%E(1:NR,MM-PSIi%INTH,KK-PSIi%INX)/TFM%R))
      ELSE
            MAXVELP = 0
      ENDIF
   ELSE
      MAXVELP = 0
   ENDIF
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,MAXVELP,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_IVP,IERR)
   CALL MPI_BARRIER(MPI_COMM_IVP,IERR)

   CALL DEALLOCATE(RURi)
   CALL DEALLOCATE(RUPi)
   CALL DEALLOCATE( UZi)
END FUNCTION MAXVELP
! ======================================================================

SUBROUTINE HYPSET()
! ======================================================================
! CALCULATE THE HYPERVISCOSITY
! EXP[ - NUP * (FILTER * MODE # / FACTOR)^P * DT]
! ======================================================================
  IMPLICIT NONE

  REAL(P8) :: NU_R, NU_T, NU_X
  INTEGER :: NT_UP, NX_UP

  REAL(P8) :: FACTOR_NU
  INTEGER :: II, JJ
  INTEGER :: FACTOR_UP
  
  REAL(P8) :: FACTOR_R, FACTOR_T, FACTOR_X
  REAL(P8) :: INDEX_R(NRCHOP), INDEX_T(NTCHOPDIM), INDEX_X(NXCHOPDIM)
  REAL(P8) :: FILTER_R(NRCHOP), FILTER_T(NTCHOPDIM), FILTER_X(NXCHOPDIM)
!  REAL(P8) :: HYPER_R(NRCHOP), HYPER_T(NTCHOPDIM), HYPER_X(NXCHOPDIM)

  FACTOR_NU = 2.D0

  FACTOR_R = 1.D0
  FACTOR_T = 0.8
  FACTOR_X = 1.D0

  ! IF PERTURBATION/DISTURBANCE SIZE IS 0.1, ITS FACTOR_UP-TH HARMONICS
  ! WILL GET BELOW MACHINE ROUND-OFF ERROR
  FACTOR_UP = 16 ! DP
  NT_UP = MIN(FACTOR_UP*MAX(MAXVAL(MONITOR_MK(1:3,1)),1),NTCHOP-1)
  NX_UP = MIN(FACTOR_UP*MAX(MAXVAL(ABS(MONITOR_MK(1:3,2))),1),NXCHOP-1)

  ! CREATE TANH FILTERS
  ! R: CENTER AT 1/2*NRCHOP WIDTH NRCHOP/4
  INDEX_R = (/(II, II=1,NRCHOP)/)
  FILTER_R = 0.5*(TANH((INDEX_R/REAL(NRCHOP)-1.D0/2.D0)*4.D0)+1.D0)
  ! T: CENTER AT 2/3*NT_UP WIDTH NT_UP/4
  INDEX_T = (/(II, II=0,NTCHOPDIM-1)/)
  FILTER_T = 0.5*(TANH((INDEX_T/REAL(NT_UP)-1.D0/2.D0)*5.D0)+1.D0)
  ! X: CENTER AT 2/3*NX_UP WIDTH NX_UP/4
  DO II=1,NXCHOPDIM
    INDEX_X(II) = II-1
    IF(II.GT.NXCHOP) THEN
        INDEX_X(II)=ABS(-(NXCHOPDIM-II+1))
    ENDIF
  ENDDO
  FILTER_X = 0.5*(TANH((INDEX_X/REAL(NX_UP)-1.D0/2.D0)*8.D0)+1.D0)

  ! FORCE THE LAST MODE TO DECAY AT EXP(-1*DT) WHILE ALSO USING TANH TO
  ! PREVENT CHANGES TO THE INTERESTED MODES
  NU_R = FACTOR_NU*(FACTOR_R*FILTER_R(NRCHOP))**(-VISC%P) 
  NU_T = FACTOR_NU*(FACTOR_T*FILTER_T(NTCHOPDIM))**(-VISC%P+2)
  NU_X = FACTOR_NU*(FACTOR_X*FILTER_X(NXCHOP))**(-VISC%P+2)

  ! CREATE HYPERV VECTORS
  ALLOCATE(HYPER_R(NRCHOP), HYPER_T(NTCHOPDIM), HYPER_X(NXCHOPDIM))
  HYPER_R = EXP(-(NU_R*TIM%DT)*(FILTER_R*INDEX_R/NRCHOP)**VISC%P)
  HYPER_T = EXP(-(NU_T*TIM%DT)*(FILTER_T*INDEX_T/NTCHOP)**(VISC%P-2))
  HYPER_X = EXP(-(NU_X*TIM%DT)*(FILTER_X*INDEX_X/NXCHOP)**(VISC%P-2))

   ! SAVE HYPER VECTORS
   IF (MPI_RANK.EQ.0) THEN
      open(UNIT=888,FILE='hyperv.dat',STATUS='UNKNOWN',ACTION='WRITE')
      ! write(888,518) NU_R, FILTER_R
      ! write(888,518) NU_T, FILTER_T
      ! write(888,518) NU_X, FILTER_X
      write(888,518) NU_R, HYPER_R
      write(888,518) NU_T, HYPER_T
      write(888,518) NU_X, HYPER_X
      close(888)
   518 FORMAT((E23.16),*(' ,',E23.16))
      ! CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
   ENDIF

  RETURN
END SUBROUTINE
! ======================================================================

SUBROUTINE HYPERV3(PSI,CHI)
! ======================================================================
! APPLY HYPERVISCOSITY USING EXPONENTIAL SCALING
! ONLY ACTIVE WHEN VISC%SW = 3
! ======================================================================
  IMPLICIT NONE
  TYPE(SCALAR):: PSI,CHI
  INTEGER:: MM, KK, NN

   IF (VISC%SW.EQ.3) THEN
    
      IF (.NOT.(ALLOCATED(HYPER_R))) CALL HYPSET()

      DO KK = 1,SIZE(PSI%E,3); DO MM = 1,SIZE(PSI%E,2)

         NN = NRCHOPS(MM + PSI%INTH)
         PSI%E(:NN,MM,KK) = PSI%E(:NN,MM,KK)*HYPER_R(:NN)*HYPER_T(MM+PSI%INTH)*HYPER_X(KK+PSI%INX)
         CHI%E(:NN,MM,KK) = CHI%E(:NN,MM,KK)*HYPER_R(:NN)*HYPER_T(MM+PSI%INTH)*HYPER_X(KK+PSI%INX)

      ENDDO; ENDDO
      
   ENDIF

  RETURN
END SUBROUTINE
! ======================================================================


END MODULE MOD_DIAGNOSTICS