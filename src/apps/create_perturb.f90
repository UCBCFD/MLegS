!=======================================================================
!     
!     WRITTEN BY JINGE WANG
!     WRITTEN BY SANGJOON (JOON) LEE
!     DEPT. OF MECHANICAL ENGINEERING
!     UNIV. OF CALIFORNIA AT BERKELEY
!     EMAIL: JINGE@BERKELEY.EDU
!     EMAIL: SANGJOONLEE@BERKELEY.EDU
!     
!     UC BERKELEY CFD LAB
!     HTTPS://CFD.ME.BERKELEY.EDU/
!
!     NONCOMMERCIAL USE WITH COPYRIGHTED (C) MARK
!     UNDER DEVELOPMENT FOR RESEARCH PURPOSE 
!
!=======================================================================
PROGRAM create_perturb
!=======================================================================
! [USAGE]: 
! CREATE new_perturb.input FOR PROGRAM ADDPERTURB
! [UPDATES]:
! LAST UPDATE ON FEB 3, 2021
!=======================================================================
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
      USE MOD_DIAGNOSTICS
      USE MOD_INIT

      IMPLICIT NONE
      TYPE(SCALAR):: PSI0,CHI0,RUR0,RUP0,UZ0,ROR0,ROP0,OZ0
      TYPE(SCALAR):: PSIU,CHIU,RURU,RUPU,UZU,RORU,ROPU,OZU
      TYPE(SCALAR):: PSI1,CHI1,PSI2,CHI2,PSI3,CHI3
      TYPE(SCALAR):: W1,W2
      INTEGER     :: MIND,AKIND,MREAD,IS,NDIM,EIG_INDEX
      REAL(P8)    :: AKREAD,rp1,ip1,EIG_ENE
      COMPLEX(P8),DIMENSION(:,:),ALLOCATABLE:: H,EIG_VEC
      COMPLEX(P8),DIMENSION(:),ALLOCATABLE:: EIG_VAL,EIG_VEC_DESIRED
      COMPLEX(P8):: EIG_VAL_DESIRED

      INTEGER     :: I,J,K

      WRITE(*,*) 'PROGRAM STARTED'

      call execute_command_line('echo "%read.input" | ./bin/init_exec')

      CALL PRINT_REAL_TIME()   ! @ MOD_MISC

      ! TYPE '%read.input' TO READ THE INPUT VALUES FROM read.input
      CALL READIN(5)           ! @ MOD_INIT
      ! call execute_command_line('echo "%read.input"')

      CALL ALLOCATE(PSI0)
      CALL ALLOCATE(CHI0)
      CALL ALLOCATE(RUR0)
      CALL ALLOCATE(RUP0)
      CALL ALLOCATE( UZ0)
      CALL ALLOCATE(ROR0)
      CALL ALLOCATE(ROP0)
      CALL ALLOCATE( OZ0)

      ! LOAD INITIAL PSI (PSI0) AND CHI (CHI0) DATA
      ! PERFORM 'INIT' AS A PREREQUISITE TO MAKE PSI0 AND CHI0
      CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI0,PSI0)
      CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHI0,CHI0)
      CALL CHOPDO(PSI0)
      CALL CHOPDO(CHI0)

      ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VORTICITY VARIALBES
      ! RETURNS (R*OMEGA0_R,R*OMEGA0_THETA,OMEGA0_Z) IN FFF SPACE
      CALL PC2VOR(PSI0,CHI0,ROR0,ROP0,OZ0)

      ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VORTICITY VARIALBES
      ! RETURNS (R*U0_R,R*U0_THETA,U0_Z) IN FFF SPACE
      CALL PC2VEL(PSI0,CHI0,RUR0,RUP0,UZ0)

      WRITE(*,*) ''
      WRITE(*,102) MINC
      WRITE(*,*) 'M_IND = INPUT - 1'
      WRITE(*,*) 'CHOOSE AN AZIMUTHAL WAVENUMBER (M) INDEX; M_IND ='
      WRITE(*,101) NTCHOPDIM
      READ(5,'(I72)') MIND
      WRITE(*,*) ''
      WRITE(*,103) ZLEN
      WRITE(*,104) NXCHOP
      WRITE(*,105) NXCHOPDIM, NXCHOP
      WRITE(*,*) 'CHOOSE AN AXIAL WAVENUMBER (AK) INDEX; AK_IND ='
      WRITE(*,101) NXCHOPDIM
      READ(5,'(I72)') AKIND
  101  FORMAT(' *** POSSIBLE INDICES: 1 ~ ', I3)
  102  FORMAT(' M  = ', I3, ' *(M_IND - 1)')
  103  FORMAT(' AK = 2*PI*AK_IND/', F0.3)
  104  FORMAT(' AK_IND = INPUT - 1          (INPUT <=', I3, ')')
  105  FORMAT('          -(', I3,' - INPUT + 1) (INPUT >', I3, ')')

      MREAD = M(MIND)
      AKREAD = AK(MIND,AKIND)

      WRITE(*,*) 'AZIMUTHAL AND AXIAL WAVENUMBERS:'
      WRITE(*,106) MREAD, AKREAD
  106  FORMAT('M = ',I8,' & K = ',F8.3)

      ALLOCATE(H(2*NRCHOPS(MIND),2*NRCHOPS(MIND)))

      CALL ALLOCATE(PSIU)
      CALL ALLOCATE(CHIU)
      CALL ALLOCATE(RURU)
      CALL ALLOCATE(RUPU)
      CALL ALLOCATE( UZU)
      CALL ALLOCATE(RORU)
      CALL ALLOCATE(ROPU)
      CALL ALLOCATE( OZU)

      CALL ALLOCATE(PSI1)
      CALL ALLOCATE(CHI1)

      CALL ALLOCATE(PSI2)
      CALL ALLOCATE(CHI2)

      CALL ALLOCATE(PSI3)
      CALL ALLOCATE(CHI3)

      DO I = 1,NRCHOPS(MIND)
        PSIU%LN = 0.D0
        CHIU%LN = 0.D0
        PSIU%E = 0.D0
        CHIU%E = 0.D0
        PSIU%E(I,MIND,AKIND) = 1.D0

        CALL CHOPSET(3)

        IF (VISC%SW.EQ.1) THEN
          IF(VISC%NU.EQ.0.0) THEN
            WRITE(*,*) 'VISC: VISC%NU CANNOT BE ZERO'
            STOP
          ENDIF
          CALL DEL2(PSIU,PSI3)
          CALL DEL2(CHIU,CHI3)
          PSI3%E = PSI3%E * VISC%NU
          CHI3%E = CHI3%E * VISC%NU
        ELSE
          PSI3%E = 0.D0
          CHI3%E = 0.D0          
          PSI3%LN = 0.D0
          CHI3%LN = 0.D0
        ENDIF            

        CALL PC2VOR(PSIU,CHIU,RORU,ROPU,OZU)
        CALL PC2VEL(PSIU,CHIU,RURU,RUPU,UZU)

        CALL VPROD(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        CALL PROJECT(RURU,RUPU,UZU,PSI1,CHI1)
        
        CALL VPROD(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        CALL PROJECT(RORU,ROPU,OZU,PSI2,CHI2)

        CALL CHOPSET(-3)

        PSI1%E = -PSI1%E + PSI2%E + PSI3%E
        CHI1%E = -CHI1%E + CHI2%E + CHI3%E

        H(:NRCHOPS(MIND),  I) = PSI1%E(:NRCHOPS(MIND),MIND,AKIND)
        H(NRCHOPS(MIND)+1:,I) = CHI1%E(:NRCHOPS(MIND),MIND,AKIND)
        
        PSIU%LN = 0.D0
        CHIU%LN = 0.D0
        PSIU%E = 0.D0
        CHIU%E = 0.D0
        CHIU%E(I,MIND,AKIND) = 1.D0

        CALL CHOPSET(3)
        IF (VISC%SW.EQ.1) THEN
          CALL DEL2(PSIU,PSI3)
          CALL DEL2(CHIU,CHI3)
          PSI3%E = PSI3%E * VISC%NU
          CHI3%E = CHI3%E * VISC%NU
        ELSE
          PSI3%E = 0.D0
          CHI3%E = 0.D0          
          PSI3%LN = 0.D0
          CHI3%LN = 0.D0
        ENDIF
        CALL CHOPSET(-3)  !IDEL2 must have the original chopset.       

        CALL IDEL2(CHIU,CHIU)

        CALL CHOPSET(3)

        CALL PC2VOR(PSIU,CHIU,RORU,ROPU,OZU)
        CALL PC2VEL(PSIU,CHIU,RURU,RUPU,UZU)

        CALL VPROD(ROR0,ROP0,OZ0,RURU,RUPU,UZU)
        CALL PROJECT(RURU,RUPU,UZU,PSI1,CHI1)
        
        CALL VPROD(RUR0,RUP0,UZ0,RORU,ROPU,OZU)
        CALL PROJECT(RORU,ROPU,OZU,PSI2,CHI2)

        CALL CHOPSET(-3)

    PSI1%E = -PSI1%E + PSI2%E + PSI3%E
    CHI1%E = -CHI1%E + CHI2%E + CHI3%E

        H(:NRCHOPS(MIND),  I+NRCHOPS(MIND)) &
                            = PSI1%E(:NRCHOPS(MIND),MIND,AKIND)
        H(NRCHOPS(MIND)+1:,I+NRCHOPS(MIND)) &
                            = CHI1%E(:NRCHOPS(MIND),MIND,AKIND)

!            WRITE(*,107) I, NRCHOPS(MIND)
      ENDDO
!     107  FORMAT('PROGRESS: ',I3.3,'/',I3.3)

      
      !! CREATE PERTURBATION: new_perturb.input
      EIG_VAL = GENEIG(H) ! store eigenvalues
      EIG_VEC = EIGVEC('R', H) ! store eigenvectors [Nx1]
      CALL MCAT(EIG_VAL) ! print all eigenvalues
      WRITE(6,*) ''
      WRITE(6,*) '---CREATE PERTURBATION---'
      WRITE(6,*) 'CHOOSE EIGEN_VAL # ='
      WRITE(6,101) 2*NRCHOPS(MIND)
      READ(5,'(I72)') EIG_INDEX ! I now becomes the index of the desired eigenvalue
      EIG_VAL_DESIRED = EIG_VAL(EIG_INDEX)

      EIG_VEC_DESIRED = EIG_VEC(:,EIG_INDEX)
      PSIU%LN = 0.D0
      CHIU%LN = 0.D0
      PSIU%E = 0.D0
      CHIU%E = 0.D0
      PSIU%E(:NRCHOPS(MIND),MIND,AKIND) = EIG_VEC_DESIRED(:NRCHOPS(MIND))
      CHIU%E(:NRCHOPS(MIND),MIND,AKIND) = EIG_VEC_DESIRED(NRCHOPS(MIND)+1:)
      CALL IDEL2(CHIU,CHIU) ! Computed Eigenvec were in terms of [Psi,DEL2(Chi)]

      ! Energy Calculation with Psi, Chi
      CALL ALLOCATE(W1)
      CALL ALLOCATE(W2)
      ! Psi part:
      W1=PSIU
      CALL DELSQH(W1,W2)
      CALL RTRAN(W2,1)
      W1=PSIU
      CALL RTRAN(W1,1)
      EIG_ENE = -SUM(PRODUCT_MK_MODIFIED(W2,W1))
      ! Chi part:
      CALL DEL2(CHIU,W1)
      CALL RTRAN(W1,1)
      CALL DELSQH(CHIU,W2)
      CALL RTRAN(W2,1)
      EIG_ENE = EIG_ENE + SUM(PRODUCT_MK_MODIFIED(W1,W2))
      ! Deallocate:
      CALL DEALLOCATE(W1)
      CALL DEALLOCATE(W2)
      PSIU%E = PSIU%E / (EIG_ENE**0.5)
      CHIU%E = CHIU%E / (EIG_ENE**0.5)
      CALL PRINT_ENERGY_SPECTRUM(PSIU,CHIU,0)

      EIG_VEC_DESIRED(:NRCHOPS(MIND))   = PSIU%E(:NRCHOPS(MIND),MIND,AKIND)
      EIG_VEC_DESIRED(NRCHOPS(MIND)+1:) = CHIU%E(:NRCHOPS(MIND),MIND,AKIND)


      open(10,FILE='new_perturb.input',STATUS='unknown',ACTION='WRITE',IOSTAT=IS)
      if (IS.ne.0) then
        print *, 'ERROR: createperturb -- Could not creat new file new_perturb.dat'
        GOTO 1
      end if

      ! Number of perturbation (currently: 1 at a time)
      I = 1
      WRITE(10,'(I2)') I

      ! Scale
      WRITE(6,*) ''
      WRITE(6,*) 'SELECT SCALE(EPSILON)'
      WRITE(6,*) 'REAL(SCALE) ='
      READ(5,*) rp1
      WRITE(6,*) 'IMAG(SCALE) ='
      READ(5,*) ip1
      WRITE(10,'(F9.7,1X,F9.7)') rp1, ip1
                  
      ! Component
      WRITE(10,*) 'right eigenfunction'
      WRITE(10,108) 'ndim','m','k'
      NDIM = size(EIG_VEC_DESIRED,1)/2
      WRITE(10,109) NDIM, MREAD, AKREAD
        ! Below are not actually used:
        WRITE(10,108) 'q','h','b'
        WRITE(10,110) QPAIR%Q(1),QPAIR%H(1),QPAIR%B(1) ! MIGHT NEED TO CHANGE!!!!!! 
        WRITE(10,108) 'nu_pow','nu','lmap'
        WRITE(10,110) VISC%NUP,VISC%NU,ELL
        WRITE(10,'(9X,A9,9X,A9)') 're(sigma)','im(sigma)'
        WRITE(10,'(9X,F9.7,9X,F9.7)') REAL(EIG_VAL_DESIRED),AIMAG(EIG_VAL_DESIRED)
        WRITE(10,*) 'radial spectral coefficients'
        WRITE(10,111) 're(psi)','im(psi)','re(chi)','im(chi)'
      DO I = 1,NDIM
        WRITE(10,112) REAL(EIG_VEC_DESIRED(I)),AIMAG(EIG_VEC_DESIRED(I)),REAL(EIG_VEC_DESIRED(NDIM+I)),AIMAG(EIG_VEC_DESIRED(NDIM+I))
      ENDDO
      close(10)

  108  FORMAT(3(9X,A9))
  109  FORMAT(2(9X,I9),9X,F9.7)
  110  FORMAT(3(9X,F9.5))
  111  FORMAT(4(9X,A9))
  112  FORMAT(4(G20.12,:,''))

      ! Save to PSI0 AND CHI0
      IF (MREAD .EQ. 0) THEN 
        WRITE(*,*) "WARNING: MUST USE ADDPERTURB"
        GOTO 1 
      ENDIF
      PSI0%E(:NRCHOPS(MIND),MIND,AKIND) = PSI0%E(:NRCHOPS(MIND),MIND,AKIND) + PSIU%E(:NRCHOPS(MIND),MIND,AKIND)*cmplx(rp1,ip1,p8)
      CHI0%E(:NRCHOPS(MIND),MIND,AKIND) = CHI0%E(:NRCHOPS(MIND),MIND,AKIND) + CHIU%E(:NRCHOPS(MIND),MIND,AKIND)*cmplx(rp1,ip1,p8)
      where(ABS(PSI0%E)<1.0e-24) PSI0%E=0.0_p8
      where(ABS(CHI0%E)<1.0e-24) CHI0%E=0.0_p8

      call msave(PSI0, TRIM(ADJUSTL(FILES%SAVEDIR))//files%psi0)
      call msave(CHI0, TRIM(ADJUSTL(FILES%SAVEDIR))//files%chi0)


1        DEALLOCATE( H )
      DEALLOCATE(EIG_VAL)
      DEALLOCATE(EIG_VEC)
      DEALLOCATE(EIG_VEC_DESIRED)

      CALL DEALLOCATE(PSI1)
      CALL DEALLOCATE(CHI1)
      CALL DEALLOCATE(PSI2)
      CALL DEALLOCATE(CHI2)
      CALL DEALLOCATE(PSI3)
      CALL DEALLOCATE(CHI3)

      CALL DEALLOCATE(PSIU)
      CALL DEALLOCATE(CHIU)
      CALL DEALLOCATE(RURU)
      CALL DEALLOCATE(RUPU)
      CALL DEALLOCATE( UZU)
      CALL DEALLOCATE(RORU)
      CALL DEALLOCATE(ROPU)
      CALL DEALLOCATE( OZU)

      CALL DEALLOCATE(PSI0)
      CALL DEALLOCATE(CHI0)
      CALL DEALLOCATE(RUR0)
      CALL DEALLOCATE(RUP0)
      CALL DEALLOCATE( UZ0)
      CALL DEALLOCATE(ROR0)
      CALL DEALLOCATE(ROP0)
      CALL DEALLOCATE( OZ0)

      WRITE(*,*) ''
      WRITE(*,*) 'PROGRAM FINISHED'
      CALL PRINT_REAL_TIME()  ! @ MOD_MISC

      ! Run shell script directly
      ! call execute_command_line('./vort.sh')

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================
! N/A
!=======================================================================
END PROGRAM create_perturb
!=======================================================================