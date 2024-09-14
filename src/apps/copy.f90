!=======================================================================
!
!     WRITTEN BY JINGE WANG
!     DEPARTMENT OF MECHANICAL ENGINEERING
!     UNIVERSITY OF CALIFORNIA AT BERKELEY
!     E-MAIL: JINGE@BERKELEY.EDU
!     
!     UC BERKELEY CFD LAB
!     URL: HTTPS://CFD.ME.BERKELEY.EDU/
!
!     NONCOMMERCIAL USE WITH COPYRIGHTED (C) MARK
!     UNDER DEVELOPMENT FOR RESEARCH PURPOSE 
!     LAST UPDATE ON FEB 3, 2021
!
!=======================================================================
PROGRAM COPY
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
          TYPE(SCALAR):: psi,chi

    
          WRITE(*,*) 'PROGRAM STARTED'
          CALL PRINT_REAL_TIME()  ! @ MOD_MISC
          call readin(5)
          
          call allocate(psi)
          call allocate(chi)

          call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//files%psi0,psi)
          call mload(TRIM(ADJUSTL(FILES%SAVEDIR))//files%chi0,chi)
          CALL MSAVE(psi, TRIM(ADJUSTL(FILES%SAVEDIR))// files%psii)
          CALL MSAVE(chi, TRIM(ADJUSTL(FILES%SAVEDIR))// files%chii)

          CALL DEALLOCATE(psi)
          CALL DEALLOCATE(chi)
    
          WRITE(*,*) ''
          WRITE(*,*) 'PROGRAM FINISHED'
          CALL PRINT_REAL_TIME()  ! @ MOD_MISC
    
    END PROGRAM COPY
    !=======================================================================