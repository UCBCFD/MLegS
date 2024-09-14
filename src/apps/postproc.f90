!=======================================================================
!
!     WRITTEN BY SANGJOON (JOON) LEE
!     DEPT. OF MECHANICAL ENGINEERING
!     UNIV. OF CALIFORNIA AT BERKELEY
!     EMAIL: SANGJOONLEE@BERKELEY.EDU
!     
!     UC BERKELEY CFD LAB
!     HTTPS://CFD.ME.BERKELEY.EDU/
!
!     NONCOMMERCIAL USE WITH COPYRIGHTED (C) MARK
!     UNDER DEVELOPMENT FOR RESEARCH PURPOSE 
!
!=======================================================================
PROGRAM POSTPROC
!=======================================================================
! [USAGE]: 
! READ TOROIDAL (PSII) AND POLODIAL (CHII) TERMS OF THE VELOCITY FIELD
! AND CREATE THE 3-D OR SECTIONAL VELOCITY OR VORTICITY FIELD
! OUTPUT DATA MAY BE COMPATIBLE WITH MATLAB
! [NOTE]:
! POSTPROCESS%JOB = 0: SAVE OMEGA_Z IN THE R-THETA PLANE
! POSTPROCESS%JOB = 1: SAVE R*U_THETA IN THE R-THETA PLANE
! POSTPROCESS%SLICEINT = 999: 3-DIMENSIONAL FIELD
! [UPDATES]:
! BASED ON TATSU'S ORIGINAL CODE POSTPROCESS
! LAST UPDATE ON NOV 23, 2020
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
      TYPE(SCALAR):: PSI,CHI,RUR,RUP,UZ,ROR,ROP,OZ 
      CHARACTER(LEN=72) :: PSIFILE 
      CHARACTER(LEN=72) :: CHIFILE 
      CHARACTER(LEN=72) :: FILENAME, WRITENAME, OMEGAZFILE, UTHETAFILE
      CHARACTER(LEN=8)  :: TEXT
      CHARACTER(LEN=3)  :: NUM
      INTEGER :: DATATYPE, RANGE, START, FINISH, STATUS, I

      WRITE(*,*) 'PROGRAM STARTED'
      CALL PRINT_REAL_TIME()   ! @ MOD_MISC

      ! TYPE '%read.input' TO READ THE INPUT VALUES FROM read.input
      CALL READIN(5)           ! @ MOD_INIT

      CALL ALLOCATE(PSI)
      CALL ALLOCATE(CHI)
      CALL ALLOCATE(RUR)
      CALL ALLOCATE(RUP)
      CALL ALLOCATE(UZ)
      CALL ALLOCATE(ROR)
      CALL ALLOCATE(ROP)
      CALL ALLOCATE(OZ)

      ! SAVE OMEGA_Z IN THE R-THETA PLANE
      IF(POSTPROCESS%JOB.EQ.0) THEN
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSII,PSI)
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHII,CHI)
          
          ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VORTICITY VARIALBES
          ! RETURNS (R*OMEGA_R,R*OMEGA_THETA,OMEGA_Z) IN FFF SPACE
          CALL PC2VOR(PSI,CHI,ROR,ROP,OZ)

          ! R*OMEGA_R (FFF) --> R*OMEGA_R (PPP)
          CALL TOFP(ROR)

          ! R*OMEGA_THETA (FFF) --> R*OMEGA_THETA (PPP)
          CALL TOFP(ROP)

          ! OMEGA_Z (FFF) --> OMEGA_Z (PPP)
          CALL TOFP(OZ)

        DO I=POSTPROCESS%START,POSTPROCESS%FINISH
          WRITE(NUM, 21) I

          IF (POSTPROCESS%SLICEINT.EQ.999) THEN
            OMEGAZFILE = 'omegaZ_3D_'//NUM//'.dat'
          ELSE
            OMEGAZFILE = 'omegaZ_RTplane_'//NUM//'.dat'  
          ENDIF

          CALL MSAVE_SLICES_IN_RTHETA_PLANE(OZ,TRIM(ADJUSTL(FILES%SAVEDIR))//OMEGAZFILE,&
                                            POSTPROCESS%SLICEINT)
        ENDDO

      ! SAVE R*(UTHETA) IN THE R-THETA PLANE
      ELSE IF(POSTPROCESS%JOB.EQ.1) THEN
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSII,PSI)
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHII,CHI)

          ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VELOCITY VARIALBES
          ! RETURNS (R*V_R,R*V_THETA,V_Z) IN FFF SPACE
          CALL PC2VEL(PSI,CHI,RUR,RUP,UZ)

          ! R*V_R (FFF) --> R*V_R (PPP)
          CALL TOFP(RUR)

          ! R*V_THETA (FFF) --> R*V_THETA (PPP)
          CALL TOFP(RUP)

          ! V_Z (FFF) --> V_Z (PPP)
          CALL TOFP(UZ)

        DO I=POSTPROCESS%START,POSTPROCESS%FINISH
          WRITE(NUM, 21) I

          UTHETAFILE = 'ru_Theta_RTplane_'//NUM//'.dat'
          CALL MSAVE_SLICES_IN_RTHETA_PLANE(RUP, TRIM(ADJUSTL(FILES%SAVEDIR))//UTHETAFILE,&
                                            POSTPROCESS%SLICEINT)
        ENDDO
      ELSE

        WRITE(*,*) 'POSTPROCESS: WRONG JOB # INPUT'
      ENDIF

      CALL DEALLOCATE(PSI)
      CALL DEALLOCATE(CHI)
      CALL DEALLOCATE(RUR)
      CALL DEALLOCATE(RUP)
      CALL DEALLOCATE(UZ)
      CALL DEALLOCATE(ROR)
      CALL DEALLOCATE(ROP)
      CALL DEALLOCATE(OZ)

 21   FORMAT (I3.3)

      WRITE(*,*) ''
      WRITE(*,*) 'PROGRAM FINISHED'
      CALL PRINT_REAL_TIME()  ! @ MOD_MISC

CONTAINS
!=======================================================================
!=================== PROGRAM-DEPENDENT SUBROUTINES =====================
!=======================================================================
      SUBROUTINE MSAVE_SLICES_IN_RTHETA_PLANE(A,FILENAME,ZPLANE)
!=======================================================================
      TYPE(SCALAR):: A
      CHARACTER(LEN=*):: FILENAME
      INTEGER:: ZPLANE

      INTEGER:: STATUS,MM,KK,I

      IF(A%SPACE.EQ.FFF_SPACE) THEN
        WRITE(*,*) 'MSAVE_SLICES_IN_RTHETA_PLANE: INPUT IN FFF SPACE'
        STOP
      ENDIF

      IF (ZPLANE.EQ.999) THEN
        CALL MSAVE(A%E(:NR,:NTH,:NX), FILENAME)
      ELSE
        CALL MSAVE(A%E(:NR,:NTH,ZPLANE), FILENAME)
      ENDIF

      RETURN
      END SUBROUTINE MSAVE_SLICES_IN_RTHETA_PLANE
!=======================================================================
      SUBROUTINE MSAVE_SLICES_IN_RZ_PLANE(A,FILENAME,THETAPLANE)
!=======================================================================
      TYPE(SCALAR):: A
      CHARACTER(LEN=*):: FILENAME
      INTEGER:: THETAPLANE

      INTEGER:: STATUS,MM,KK,I

      IF(A%SPACE.EQ.FFF_SPACE) THEN
        WRITE(*,*) 'MSAVE_SLICES_IN_RZ_PLANE: INPUT IN FFF SPACE'
        STOP
      ENDIF

      CALL MSAVE(A%E(:NR,THETAPLANE,:NX), FILENAME)

      RETURN
      END SUBROUTINE MSAVE_SLICES_IN_RZ_PLANE
!=======================================================================
END PROGRAM POSTPROC
!=======================================================================