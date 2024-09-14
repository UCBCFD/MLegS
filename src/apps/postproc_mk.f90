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
PROGRAM POSTPROC_MK
!=======================================================================
! [USAGE]: 
! READ TOROIDAL (PSII) AND POLODIAL (CHII) TERMS OF THE VELOCITY FIELD
! WIPE OUT THE BASE FLOW
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
      CHARACTER(LEN=72) :: WRITENAME, OMEGAZFILE, UTHETAFILE
      CHARACTER(LEN=8)  :: TEXT
      CHARACTER(LEN=3)  :: NUM
      INTEGER :: DATATYPE, RANGE, START, FINISH, I, MIND, KIND

      WRITE(*,*) 'PROGRAM STARTED'
      CALL PRINT_REAL_TIME()   ! @ MOD_MISC

      ! TYPE '%read.input' TO READ THE INPUT VALUES FROM read.input
      CALL READIN(5)           ! @ MOD_INIT
      
      WRITE(*,*) 'M#：'
      READ(5,'(I72)') MIND
      MIND = MIND + 1
      WRITE(*,*) 'K#：'
      READ(5,'(I72)') KIND
      KIND = KIND + 1

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
        IF (POSTPROCESS%FINISH .EQ. 1) THEN
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSII,PSI)
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHII,CHI)

          PSI%E(:,1,1) = 0.D0;
          CHI%E(:,1,1) = 0.D0;
          
          ! Edited: avoid XXDX warning
          !PSI%E(NRCHOPS(2),2,1) = 0.
          !CHI%E(NRCHOPS(2),2,1) = 0.

          ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VORTICITY VARIALBES
          ! RETURNS (R*OMEGA_R,R*OMEGA_THETA,OMEGA_Z) IN FFF SPACE
          CALL PC2VOR(PSI,CHI,ROR,ROP,OZ)

          ! R*OMEGA_R (FFF) --> R*OMEGA_R (PFF)
          CALL RTRAN(ROR,1)

          ! R*OMEGA_THETA (FFF) --> R*OMEGA_THETA (PFF)
          CALL RTRAN(ROP,1)

          ! OMEGA_Z (FFF) --> OMEGA_Z (PFF)
          CALL RTRAN(OZ,1)

        
          WRITE(NUM, 21) I

          OMEGAZFILE = 'omegaZ_MK_'//NUM//'.dat'    
          CALL MSAVE_SLICES_IN_MK(OZ, TRIM(ADJUSTL(FILES%SAVEDIR))//OMEGAZFILE,&
                                            MIND,KIND)

        ELSE
          DO I=POSTPROCESS%START,POSTPROCESS%FINISH

            CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI(I),PSI)
            CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHI(I),CHI)

            PSI%E(:,1,1) = 0.D0;
            CHI%E(:,1,1) = 0.D0;

            ! Edited: avoid XXDX warning
            ! PSI%E(NRCHOPS(2),2,1) = 0.
            ! CHI%E(NRCHOPS(2),2,1) = 0.

            ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VORTICITY VARIALBES
            ! RETURNS (R*OMEGA_R,R*OMEGA_THETA,OMEGA_Z) IN FFF SPACE
            CALL PC2VOR(PSI,CHI,ROR,ROP,OZ)

            ! R*OMEGA_R (FFF) --> R*OMEGA_R (PFF)
            CALL RTRAN(ROR,1)

            ! R*OMEGA_THETA (FFF) --> R*OMEGA_THETA (PFF)
            CALL RTRAN(ROP,1)

            ! OMEGA_Z (FFF) --> OMEGA_Z (PFF)
            CALL RTRAN(OZ,1)


            WRITE(NUM, 21) I

            OMEGAZFILE = 'omegaZ_MK_'//NUM//'.dat'  
            CALL MSAVE_SLICES_IN_MK(OZ, TRIM(ADJUSTL(FILES%SAVEDIR))//OMEGAZFILE,&
                                            MIND,KIND)
          ENDDO
        ENDIF
      ! SAVE R*(UTHETA) IN THE R-THETA PLANE
      ELSE IF(POSTPROCESS%JOB.EQ.1) THEN
        IF(POSTPROCESS%FINISH.EQ.1) THEN
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSII,PSI)
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHII,CHI)

          PSI%E(:,1,1) = 0.D0;
          CHI%E(:,1,1) = 0.D0;

          ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VELOCITY VARIALBES
          ! RETURNS (R*V_R,R*V_THETA,V_Z) IN FFF SPACE
          CALL PC2VEL(PSI,CHI,RUR,RUP,UZ)

          ! R*V_R (FFF) --> R*V_R (PFF)
          CALL RTRAN(RUR,1)

          ! R*V_THETA (FFF) --> R*V_THETA (PFF)
          CALL RTRAN(RUP,1)

          ! V_Z (FFF) --> V_Z (PFF)
          CALL RTRAN(UZ,1)

          WRITE(NUM, 21) I

          UTHETAFILE = 'ru_Theta_MK_'//NUM//'.dat'
          CALL MSAVE_SLICES_IN_MK(RUP, TRIM(ADJUSTL(FILES%SAVEDIR))//UTHETAFILE,&
                                            MIND,KIND)
        ELSE
          DO I=POSTPROCESS%START,POSTPROCESS%FINISH

            CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI(I),PSI)
            CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHI(I),CHI)

            PSI%E(:,1,1) = 0.D0;
            CHI%E(:,1,1) = 0.D0;

            ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VELOCITY VARIALBES
            ! RETURNS (R*V_R,R*V_THETA,V_Z) IN FFF SPACE
            CALL PC2VEL(PSI,CHI,RUR,RUP,UZ)

            ! R*V_R (FFF) --> R*V_R (PPP)
            CALL RTRAN(RUR,1)

            ! R*V_THETA (FFF) --> R*V_THETA (PPP)
            CALL RTRAN(RUP,1)

            ! V_Z (FFF) --> V_Z (PPP)
            CALL RTRAN(UZ,1)


            WRITE(NUM, 21) I

            UTHETAFILE = 'ru_Theta_MK_'//NUM//'.dat'
            CALL MSAVE_SLICES_IN_MK(RUP, TRIM(ADJUSTL(FILES%SAVEDIR))//UTHETAFILE,&
                                            MIND,KIND)
          ENDDO
        ENDIF

      ! SAVE R*(UR) IN THE R-THETA PLANE
      ELSE IF(POSTPROCESS%JOB.EQ.2) THEN
        IF(POSTPROCESS%FINISH.EQ.1) THEN
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSII,PSI)
          CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHII,CHI)

          PSI%E(:,1,1) = 0.D0;
          CHI%E(:,1,1) = 0.D0;

          ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VELOCITY VARIALBES
          ! RETURNS (R*V_R,R*V_THETA,V_Z) IN FFF SPACE
          CALL PC2VEL(PSI,CHI,RUR,RUP,UZ)

          ! R*V_R (FFF) --> R*V_R (PFF)
          CALL RTRAN(RUR,1)

          ! R*V_THETA (FFF) --> R*V_THETA (PFF)
          CALL RTRAN(RUP,1)

          ! V_Z (FFF) --> V_Z (PFF)
          CALL RTRAN(UZ,1)

          WRITE(NUM, 21) I

          UTHETAFILE = 'ru_r_MK_'//NUM//'.dat'
          CALL MSAVE_SLICES_IN_MK(RUR, TRIM(ADJUSTL(FILES%SAVEDIR))//UTHETAFILE,&
                                            MIND,KIND)
        ELSE
          DO I=POSTPROCESS%START,POSTPROCESS%FINISH

            CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%PSI(I),PSI)
            CALL MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%CHI(I),CHI)

            PSI%E(:,1,1) = 0.D0;
            CHI%E(:,1,1) = 0.D0;

            ! CONVERT THE POLOIDAL-TOROIDAL VARIALBES TO VELOCITY VARIALBES
            ! RETURNS (R*V_R,R*V_THETA,V_Z) IN FFF SPACE
            CALL PC2VEL(PSI,CHI,RUR,RUP,UZ)

            ! R*V_R (FFF) --> R*V_R (PPP)
            CALL RTRAN(RUR,1)

            ! R*V_THETA (FFF) --> R*V_THETA (PPP)
            CALL RTRAN(RUP,1)

            ! V_Z (FFF) --> V_Z (PPP)
            CALL RTRAN(UZ,1)


            WRITE(NUM, 21) I

            UTHETAFILE = 'ru_r_MK_'//NUM//'.dat'
            CALL MSAVE_SLICES_IN_MK(RUR, TRIM(ADJUSTL(FILES%SAVEDIR))//UTHETAFILE,&
                                            MIND,KIND)
          ENDDO
        ENDIF

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
      SUBROUTINE MSAVE_SLICES_IN_MK(A,FILENAME,M_IND,K_IND)
!=======================================================================
      TYPE(SCALAR):: A
      CHARACTER(LEN=*):: FILENAME
      INTEGER:: M_IND,K_IND

      IF(A%SPACE.NE.PFF_SPACE) THEN
        WRITE(*,*) 'MSAVE_SLICES_IN_RZ_PLANE: INPUT NOT IN PFF SPACE'
        STOP
      ENDIF

      CALL MSAVE(A%E(:NR,M_IND,K_IND), FILENAME)

      RETURN
      END SUBROUTINE MSAVE_SLICES_IN_MK
!=======================================================================
END PROGRAM POSTPROC_MK
!=======================================================================