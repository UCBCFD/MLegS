MODULE MOD_INIT ! LEVEL 6 MODULE
      USE OMP_LIB
      USE MPI
      USE MOD_MISC, ONLY : P4,P8,PI,IU,SPY,MSAVE,MLOAD,ATOF,ITOA3        ! LEVEL 0
!XUSE USE MOD_FD                                                         ! LEVEL 1
!XUSE USE MOD_EIG                                                        ! LEVEL 1
!XUSE USE MOD_LIN_LEGENDRE                                               ! LEVEL 1
!XUSE USE MOD_BANDMAT                                                    ! LEVEL 1
      USE MOD_SCALAR3                                                    ! LEVEL 2
      USE MOD_FFT                                                        ! LEVEL 2.5
      USE MOD_LEGOPS                                                     ! LEVEL 3
      USE MOD_LAYOUT                                                     ! LEVEL 3
      USE MOD_MARCH                                                      ! LEVEL 4
      USE MOD_DIAGNOSTICS                                                ! LEVEL 5
      IMPLICIT NONE
      PRIVATE
!=======================================================================
!============================ PARAMETERS ===============================
!=======================================================================
      INTEGER:: LNUMBER = 0                                              ! USED IN READCOM
      INTEGER:: IECHO = 0                                                ! USED IN READCOM
      INTEGER:: IEXEC = 0                                                ! USED IN READCOM
      INTEGER:: IESAV,IREAD                                              ! USED IN READCOM
!=======================================================================
!======================== PUBLIC DECLARATION ===========================
!=======================================================================
      ! READ COMMANDS
      PUBLIC :: READCOM
      ! SHIFT THE CONTENTS OF AN INPUT ARRAY BY OPEN
      PUBLIC :: SHIFTCOM
      ! READ INITIAL VALS FROM KEYBOARD OR INPUT FILE (ENTER %FILENAME)
      PUBLIC :: READIN
      ! SYNC READIN INFO FROM MASTER MODE
      PUBLIC :: READ_SYNC
!=======================================================================
!============================ INTERFACES ===============================
!=======================================================================
      INTERFACE READCOM
        MODULE PROCEDURE READCOM0
        MODULE PROCEDURE READCOM1
        MODULE PROCEDURE ECHO
      END INTERFACE

      INTERFACE READIN
        MODULE PROCEDURE READIN0
        MODULE PROCEDURE READIN_MATLAB
        MODULE PROCEDURE READIN_almost
      END INTERFACE 
CONTAINS
!=======================================================================
!============================ SUBROUTINES ==============================
!=======================================================================
      SUBROUTINE ECHO(C) ! FAMILY OF READCOM
!=======================================================================
! [USAGE]: 
! TOGGLES THE ECHO MODE (ON/OFF).
! [PARAMETERS]:
! C >> (OPTIONAL) 'ECHO' TO TURN ON THE ECHO MODE
!      IF BLANK, TOGGLES THE MODE BETWEEN ON AND OFF 
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 20 2020
!=======================================================================
      CHARACTER(LEN=*),OPTIONAL:: C

      IF(PRESENT(C)) THEN
        IF((C .EQ. 'echo') .OR. (C .EQ. 'ECHO')) THEN 
          IECHO = 1
        ELSE
          IECHO = 0
        ENDIF
      ELSE
        IECHO = 1-IECHO
      ENDIF

      RETURN
      END SUBROUTINE ECHO
!=======================================================================
      SUBROUTINE READCOM0(I) ! FAMILY OF READCOM
!=======================================================================
! [USAGE]: 
! RESET THE LINE NUMBER TO READ TO I
! [PARAMETERS]:
! I >> RESET LINE NUMBER
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 20 2020
!=======================================================================
      INTEGER :: I

      LNUMBER = I

      RETURN
      END SUBROUTINE READCOM0
!=======================================================================
      SUBROUTINE READCOM1(IR,CMD,LNO) ! FAMILY OF REACOM
!=======================================================================
! [USAGE]: 
! READ COMMANDS FROM UNIT IR
! [PARAMETERS]:
! IR >> UNIT NUMBER
! CMD >> COMMAND STRINGS
! LNO >> (OPTIONAL) ON EXIT, RETURNS THE CURRENT LINE NUMBER
! [NOTES]:
! READCOM(IR,COM) READS COMMANDS FROM UNIT IR.  IR=5 MEANS KEYBOARD.
!               IR: INTEGER, COM(:): CHARACTERS.
!               IF THE LINE READ FROM IR IS
!
!                    A 'B CD' E,F 
!
!                THEN THE RESULT IS
!
!                COM(1) = 'A'
!                COM(2) = 'B CD'
!                COM(3) = 'E'
!                COM(4) = 'F'
!                COM(5) = ' ' ...
!
!                - IF THE LINE BEGINS WITH '#', THEN THE LINE IS SKIPPED.
!                - IF THE LINE IS, FOR EXAMPLE, '%AUTO', THEN 
!                  THE CONTENT OF FILE AUTO IS READ AS IF IT IS INPUT FROM 
!                  UNIT IR. THIS FACILITY CAN'T BE NESTED.
!                - A LINE WHICH BEGINS WITH '!' IS A SHELL ESCAPE.
!
! READCOM(IR,COM,LNUMBER), WHERE INTEGER LNUMBER, RETURNS THE CURRENT LINE
!                  NUMBER.
!
! READCOM(I)     RESETS THE LINE NUMBER TO I.
!
! READCOM('ECHO')    SETS THE ECHO MODE ON.
! READCOM('NOECHO')  RESETS THE ECHO MODE ON.
! CALL READCOM       TOGGLES THE ECHO MODE.
!
! READING FORMAT  
!  (2)COMMAND [LNO] / COMMENT
!    CHARACTERS (',) ARE REPLACED WITH A SPACE BEFORE ANY ANALYZING
! [UPDATES]:
! PROC#0 WILL TAKE KEYBOARD INPUT AND BROADCAST IT TO ALL PROCS
!=======================================================================
      IMPLICIT NONE
      EXTERNAL UNIXCOM, TABTOSP
      INTEGER :: IR
      CHARACTER(LEN=*),DIMENSION(:):: CMD
      INTEGER,OPTIONAL:: LNO

      INTEGER,PARAMETER:: NLEN=200!200
      CHARACTER(LEN=NLEN):: C
      LOGICAL:: OPEN
      INTEGER:: I,III,N,NDIM

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)

      NDIM = SIZE(CMD,1)
      
      DO  I=1,NDIM
        CMD(I)=' '
      ENDDO
        
      IF(IEXEC.EQ.0) THEN
        IREAD = 18
        DO 
          INQUIRE(UNIT=IREAD,OPENED=OPEN)
          IF(.NOT.OPEN) EXIT
          IREAD=IREAD+1
        ENDDO
        OPEN(UNIT=IREAD,FILE='read.input',STATUS='OLD')
        IEXEC=1
        IESAV=IECHO
        IECHO=0
      ENDIF


 1    CONTINUE

      IF(IEXEC.NE.1) THEN
        IF (MPI_RANK.EQ.0) THEN
          WRITE(6,*) 'START READIN'
          READ(IR,10) C ! Read input
        ENDIF
        CALL MPI_BCAST(C,NLEN,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR) ! MUST BE MUTED IF USING READ_SYNC
        
      ELSE
        READ(IREAD,10,ERR=9) C
        IF(C.EQ.'%END') GOTO 9 ! Reach the end
        GOTO 19 ! Move to next line

 9    CONTINUE

        CLOSE(IREAD)
        IEXEC=0
        IECHO=IESAV
        GOTO 1
      ENDIF

      IF(C(1:1).EQ.'%') THEN ! read filename
        IREAD=18
        DO
          INQUIRE(UNIT=IREAD,OPENED=OPEN)
          IF(.NOT.OPEN) EXIT
          IREAD=IREAD+1
        ENDDO

        OPEN(UNIT=IREAD,FILE=C(2:NLEN),STATUS='OLD',ACTION='READ')
        IEXEC=1
        IESAV=IECHO
        ! IECHO=1
        IF (MPI_RANK.EQ.0) THEN
          IECHO=1
        ELSE
          IECHO=0
        ENDIF
        GOTO 1
      ENDIF

 19   CONTINUE

      LNUMBER=LNUMBER + 1

      IF(IECHO.EQ.1) THEN
        I=201

 11   CONTINUE
      
        I=I-1
        IF((C(I:I).EQ.' '.OR.C(I:I).EQ.'\0').AND.I.GT.1) GOTO 11
        WRITE(*,*) C(1:I)
      ENDIF
      
 10   FORMAT(A200)

      IF(C(1:1).EQ.'#') THEN
        GOTO 1
      ENDIF

      IF(C(1:1).EQ.'!') THEN
        CALL EXECUTE_COMMAND_LINE(C(2:72))
        GOTO 1
      ENDIF
      
 3    CONTINUE

      III=INDEX(C,',')

      IF(III.NE.0) THEN
        C=C(1:III-1)//' '//C(III+1:NLEN)
        GOTO 3
      ENDIF

 5    CONTINUE

      III=INDEX(C,CHAR(9))

      IF(III.NE.0) THEN
        C=C(1:III-1)//' '//C(III+1:NLEN)
        GOTO 5
      ENDIF

      N = 1

 2    CONTINUE

      IF(C.NE.' ') THEN
        IF(C(1:1).EQ.' '.OR.C(1:1).EQ.',') THEN
          C=C(2:NLEN)
          GOTO 2
        ENDIF
      ENDIF

      IF(C(1:1).EQ.'"') THEN
        C=C(2:NLEN)
        CMD(N)=C(1:INDEX(C,'"')-1)
        C=C(INDEX(C,'"')+1:NLEN)
      ELSEIF(C(1:1).EQ.'''') THEN
        C=C(2:NLEN)
        CMD(N)=C(1:INDEX(C,'''')-1)
        C=C(INDEX(C,'''')+1:NLEN)
      ELSEIF(C(1:1).EQ.'(') THEN
        C=C(2:NLEN)
        CMD(N)=C(1:INDEX(C,')')-1)
        C=C(INDEX(C,')')+1:NLEN)
      ELSE
        CMD(N)=C(1:INDEX(C,' '))
        C=C(INDEX(C,' ')+1:NLEN)
      ENDIF

      N=N+1

      IF(N.GT.NDIM) THEN
        WRITE(*,*) 'READCOM: TOO MANY FIELDS'
        WRITE(*,*) 'UNPROCESSED FIELD=',C
        GOTO 999
      ENDIF

      IF(C.NE.' ') GOTO 2
      
 999  CONTINUE

      IF(PRESENT(LNO)) THEN
        LNO=LNUMBER
      ENDIF

      RETURN
      END SUBROUTINE READCOM1
!=======================================================================
      SUBROUTINE SHIFTCOM(COM)
!=======================================================================
! [USAGE]:
! SHIFTS THE CONTENTS OF ARRAY BY ONE. THE COM(1) IS 
! REPLACED BY COM(2), COM(2) BY COM(3) AND SO ON. 
! THE LAST ELEMENT IS FILLED BY ' '.
! [PARAMETERS]:
! CMD >> COMMAND STRINGS
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!======================================================================= 
      IMPLICIT NONE
      CHARACTER(LEN=*),DIMENSION(:)::  COM

      INTEGER :: I,NDIM

      NDIM=SIZE(COM,1)

      DO I=1,NDIM-1
        COM(I)=COM(I+1)
      ENDDO

      COM(NDIM)=' '

      RETURN
      END SUBROUTINE SHIFTCOM
!======================================================================= 
      SUBROUTINE READIN0(IR)
!=======================================================================
! [USAGE]:
! READ COMMANDS AND ASSIGN INITIAL VALUES INTO THE VARIABLES
! [PARAMETERS]:
! IR >> IF 5, READ COMMAND FROM KEYBOARD
!       IF 5 AND ENTER '%FILENAME', START READING COMMANDS FROM FILE
! [DEPENDENCIES]:
! 1. LEGINIT(~) @ MOD_LIN_SCALAR3
! 2. READCOM(~) @ MOD_INIT
! [UPDATES]:
! RE-CODED BY SANGJOON LEE @ NOV 10 2020
!======================================================================= 
      IMPLICIT NONE
      INTEGER,INTENT(IN):: IR

      CHARACTER(LEN=72),DIMENSION(20):: COM
      INTEGER :: NRIN, NTHIN, NXIN, I, NUM
      INTEGER :: NRCHOPIN, NTCHOPIN, NXCHOPIN
      REAL(P8):: ZLENIN,ELLIN

      INTEGER :: MINCIN = 1
      INTEGER :: MKLINKIN = 0
      INTEGER :: IDEBUG = 0
      INTEGER :: CHECKER = 0

      !INTEGER :: LNUMBER = 0 !Initilization
      LNUMBER = 0

      ! mod_march:
      DIAGV%IT1=40
      DIAGV%IT2=40

      ADV%SW = 0
      ADV%INT = 10000
      ADV%UX = 0
      ADV%UY = 0
      ADV%UZ = 0

      NADD%ANG = 0
      NADD%U = 0
      NADD%UZ = 0
      NADD%STRAIN=0

      RMV%SW=0
      RMV%INT=10000

      VISC%ADJSW=0
      ! end mod_march

      CALL READCOM(LNUMBER)
      ! WRITE(6,*) 'START READIN'
      DO WHILE (CHECKER .LT. 1000)
        CHECKER = CHECKER + 1
        CALL READCOM(IR,COM,LNUMBER)
        IF(IDEBUG.EQ.1) WRITE(*,'(6A16)') 'DEBUG MODE: ',(COM(I),I=1,5)

        ! 1> SCALAR3 - DIMENSIONS
        IF(COM(1).EQ.'DIMENSION') THEN
          CALL READCOM(IR,COM)
          NRIN = ATOF(COM(1)) 
          NTHIN= ATOF(COM(2))
          NXIN = ATOF(COM(3))
        ELSEIF(COM(1).EQ.'CHOPS') THEN
          CALL READCOM(IR,COM)
          NRCHOPIN = ATOF(COM(1)) 
          NTCHOPIN = ATOF(COM(2)) 
          NXCHOPIN = ATOF(COM(3))
        ELSEIF(COM(1).EQ.'ELL') THEN
          CALL READCOM(IR,COM)
          ELLIN = ATOF(COM(1))
          ELL = ELLIN
        ELSEIF(COM(1).EQ.'ZLEN') THEN
          CALL READCOM(IR,COM)
          IF (ATOF(COM(1)).GT.1.) THEN
            ZLENIN = ATOF(COM(1))
          ELSE ! If ZLENIN <= 1, assumes it's the wavenumber
            ZLENIN = 2.*PI/ATOF(COM(1))
          ENDIF
        ELSEIF(COM(1).EQ.'MINC') THEN
          CALL READCOM(IR,COM)
          MINCIN = ATOF(COM(1))
        ELSEIF(COM(1).EQ.'MKLINK') THEN
          CALL READCOM(IR,COM)
          MKLINKIN = ATOF(COM(1))

        ! 2> LAYOUT - VORTEX PAIR(QPAIR)
        ELSEIF(COM(1).EQ.'VORTEXPAIR') THEN
          CALL READCOM(IR,COM)
          QPAIR%N = ATOF(COM(1))
          DO I=1,QPAIR%N
            CALL READCOM(IR,COM)
            QPAIR%X(I) = ATOF(COM(1))
            QPAIR%Y(I) = ATOF(COM(2))
            QPAIR%Q(I) = ATOF(COM(3))
            QPAIR%H(I) = ATOF(COM(4))
            QPAIR%B(I) = ATOF(COM(5))
          ENDDO
          QPAIR%NOISE = 0
          QPAIR%DX = 0
          QPAIR%DY = 0
          QPAIR%DP = 0
          QPAIR%K  = 1
          QPAIR%R  = 1
          QPAIR%DR = 0
          QPAIR%RP = 0
          QPAIR%RK = 0
        ELSEIF(COM(1).EQ.'VORTEXPAIR_NOISE') THEN
          CALL READCOM(IR,COM)
          QPAIR%NOISE=ATOF(COM(1))
        ELSEIF(COM(1).EQ.'VORTEXPAIR_CROW') THEN
          DO I=1,QPAIR%N
            CALL READCOM(IR,COM)
            QPAIR%DX(I)=ATOF(COM(1))
            QPAIR%DY(I)=ATOF(COM(2))
            QPAIR%DP(I)=ATOF(COM(3))
            QPAIR%K(I) =ATOF(COM(4))
          ENDDO
        ELSEIF(COM(1).EQ.'VORTEXPAIR_SECTION') THEN
          DO I=1,QPAIR%N
            CALL READCOM(IR,COM)
            QPAIR%R(I) =ATOF(COM(1))
            QPAIR%DR(I)=ATOF(COM(2))
            QPAIR%RK(I)=ATOF(COM(3))
            QPAIR%RP(I)=ATOF(COM(4))
          ENDDO
  
        ! 3> MARCH - TIME DATA
        ELSEIF(COM(1).EQ.'DT') THEN
          CALL READCOM(IR,COM)
          TIM%DT   =ATOF(COM(1))
          TIM%LIMIT=ATOF(COM(2))
          TIM%N    = 0
          TIM%T    = 0
        ELSEIF(COM(1).EQ.'TIME_INTEG') THEN
          CALL READCOM(IR,COM)
          TIM%SCHEME = ATOF(COM(1))
  
        ! 4> MARCH - VISCOSITY DATA
        ELSEIF(COM(1).EQ.'VISCOSITY') THEN
          CALL READCOM(IR,COM)
          VISC%SW  = ATOF(COM(1))
          VISC%NU  = ATOF(COM(2))
          IF (VISC%NU .GE. 1.D0) THEN
            VISC%NU = VISC%NU**-1.D0
          ENDIF
          VISC%NUP = ATOF(COM(3))
          IF (VISC%NUP .GE. 1.D0) THEN
            VISC%NUP = VISC%NUP**-1.D0
          ENDIF
          VISC%P   = ATOF(COM(4))
        ELSEIF(COM(1).EQ.'NUP_ADJUST') THEN
          CALL READCOM(IR,COM)
          VISC%ADJSW  = ATOF(COM(1))
          VISC%ADJINT = ATOF(COM(2))
  
        ! 4> MARCH - FREESTREAM DATA & LEGOPS - NONLINADD
        ELSEIF(COM(1).EQ.'FREESTREAM') THEN
          CALL READCOM(IR,COM)
          ADV%SW = ATOF(COM(1))
          ADV%UX = ATOF(COM(2))
          ADV%UY = ATOF(COM(3))
          ADV%UZ = ATOF(COM(4))
          ADV%INT = ATOF(COM(5))
          NADD%U = SQRT(ADV%UX**2+ADV%UY**2)
          NADD%UZ= ADV%UZ
          IF(NADD%U.EQ.0.0) THEN
            NADD%ANG = 0.0
          ELSEIF(ADV%UY >= 0) THEN
            NADD%ANG = ACOS(ADV%UX/NADD%U)
          ELSE 
            NADD%ANG = -1.0*ACOS(ADV%UX/NADD%U)
          ENDIF
        ELSEIF(COM(1).EQ.'STRAIN') THEN
          CALL READCOM(IR,COM)
          NADD%STRAIN = ATOF(COM(1))

        ! 5> MARCH - REMOVAL (RMV)
        ELSEIF(COM(1).EQ.'REMOVE') THEN
          CALL READCOM(IR,COM)
          RMV%SW = ATOF(COM(1))
          RMV%INT = ATOF(COM(2))
  
        !6> MARCH - FILE INOUT
        ELSEIF(COM(1).EQ.'SAVEDIR') THEN
          CALL READCOM(IR,COM)
          FILES%SAVEDIR=COM(1)
        ELSEIF(COM(1).EQ.'PSI0') THEN
          CALL READCOM(IR,COM)
          FILES%PSI0=COM(1)
        ELSEIF(COM(1).EQ.'CHI0') THEN
          CALL READCOM(IR,COM)
          FILES%CHI0=COM(1)
        ELSEIF(COM(1).EQ.'PSII') THEN
          CALL READCOM(IR,COM)
          FILES%PSII=COM(1)
        ELSEIF(COM(1).EQ.'CHII') THEN
          CALL READCOM(IR,COM)
          FILES%CHII=COM(1)
        ELSEIF(COM(1).EQ.'OUTPUTS') THEN
          CALL READCOM(IR,COM)
          FILES%NE = ATOF(COM(1))
          ALLOCATE( FILES%T(FILES%NE) )
          ALLOCATE( FILES%PSI(FILES%NE) )
          ALLOCATE( FILES%CHI(FILES%NE) )
          DO I=1,FILES%NE
            CALL READCOM(IR,COM)
            FILES%T(I)   = ATOF(COM(1))
            FILES%PSI(I) = COM(2)
            FILES%CHI(I) = COM(3)
          ENDDO
  
        !7> MARCH - DIAGNOST & VELMON
        ELSEIF(COM(1).EQ.'DIAGNOST') THEN
          CALL READCOM(IR,COM)
          DIAGV%IT1  = ATOF(COM(1))
          DIAGV%IT2  = ATOF(COM(2))
          VELMON%INT = DIAGV%IT1
          VELMON%N   = 0
  
        !8> MARCH - LINEAR
        ELSEIF(COM(1).EQ.'LINEAR') THEN
          CALL READCOM(IR,COM)
          LIN%MM=ATOF(COM(1))
          LIN%KK=ATOF(COM(2))
  
        !9> MARCH - EIGENVALUE MONITOR KIT & DIAGNOSTICS - MONITOR_MK
        ELSEIF(COM(1).EQ.'MONITOR_MK') THEN
          CALL READCOM(IR,COM)
          NUM = ATOF(COM(1))
          ALLOCATE(MONITOR_MK(NUM,2))
          DO I=1,NUM
            CALL READCOM(IR,COM)
            MONITOR_MK(I,1) = ATOF(COM(1))
            MONITOR_MK(I,2) = ATOF(COM(2))
          ENDDO
          CALL READCOM(IR,COM)
          MONITORDATA%START = ATOF(COM(1))
          MONITORDATA%INTERVAL = ATOF(COM(2))
  
        !10> MARCH - POSTPROCESS KIT
        ELSEIF(COM(1).EQ.'POSTPROCESS') THEN
          CALL READCOM(IR,COM)
          POSTPROCESS%START  = ATOF(COM(1))
          POSTPROCESS%FINISH  = ATOF(COM(2))
          POSTPROCESS%JOB = ATOF(COM(3))
          POSTPROCESS%SLICEINT = ATOF(COM(4))
  
        !11> MISCELLANOUS
        ELSEIF(COM(1).EQ.'DEBUG') THEN
          CALL READCOM(IR,COM)
          IDEBUG = ATOF(COM(1))
        ELSEIF(COM(1).EQ.' ') THEN
          CONTINUE
        ELSEIF(COM(1).EQ.'END') THEN
          EXIT
        ELSE 
          WRITE(*,*) 'READIN: UNKNOWN COMMAND IN LINE ',LNUMBER
        ENDIF
      ENDDO       

      IF (CHECKER .GE. 1000) THEN
        WRITE(*,*) 'READIN: COMMAND TOO LONG (MAX. OF 1000 LINES' //&
                   ' EXPECTED)'
        STOP
      ENDIF

      NR = NRIN
      NTH = NTHIN
      NX = NXIN
      NRCHOP = NRCHOPIN
      NTCHOP = NTCHOPIN
      NXCHOP = NXCHOPIN
      ZLEN = ZLENIN
      ZLEN0 = ZLEN
      ELL = ELLIN
      MKLINK = MKLINKIN
      MINC = MINCIN

      ! CALL LEGINIT(NRIN,NTHIN,NXIN,NRCHOPIN,NTCHOPIN,NXCHOPIN,&
      ! ZLENIN,ELLIN,MKLINKIN,MINCIN)

      !IF(LIN%KK.LT.0) LIN%KK = NXCHOPDIM+LIN%KK+1
      IF(LIN%KK.LT.0) LIN%KK = 2*NXCHOP+LIN%KK

      RETURN
      END SUBROUTINE READIN0
!======================================================================= 
SUBROUTINE READ_SYNC()
!=======================================================================
! [USAGE]:
! SYNC ALL READIN INFO THROUGH ALL PROCESSORS. MUST MUTE MPI_BCAST IN RE
! ADCOM
!======================================================================= 
      IMPLICIT NONE

      INTEGER :: NUM, I

      INTEGER :: MINCIN = 1
      INTEGER :: MKLINKIN = 0
      INTEGER :: IDEBUG = 0
      INTEGER :: CHECKER = 0

      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)
      ! IF (MPI_RANK.NE.0) THEN
      !   WRITE(*,*) 'NOT IN MASTER PROC. ABORTED'
      !   RETURN
      ! ENDIF

      CALL MPI_BCAST(NR,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NTH,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NX,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

      CALL MPI_BCAST(NRCHOP,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NTCHOP,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NXCHOP,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

      CALL MPI_BCAST(ELL,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ZLEN,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

      CALL MPI_BCAST(MINC,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(MKLINK,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

      CALL MPI_BCAST(QPAIR%N,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(QPAIR%NOISE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      DO I=1,QPAIR%N
        CALL MPI_BCAST(QPAIR%X(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(QPAIR%Y(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(QPAIR%Q(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(QPAIR%H(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(QPAIR%B(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(QPAIR%DX(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(QPAIR%DY(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(QPAIR%DP(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(QPAIR%K(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(QPAIR%R(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(QPAIR%DR(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(QPAIR%RK(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(QPAIR%RP(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      ENDDO
      
      CALL MPI_BCAST(TIM%DT,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(TIM%LIMIT,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(TIM%N,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(TIM%T,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(TIM%SCHEME,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  
      CALL MPI_BCAST(VISC%SW,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(VISC%NU,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(VISC%NUP,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(VISC%P,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(VISC%ADJSW,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(VISC%ADJINT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  
      CALL MPI_BCAST(ADV%SW,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ADV%UX,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ADV%UY,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ADV%UZ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ADV%INT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

      CALL MPI_BCAST(NADD%U,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NADD%UZ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NADD%ANG,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NADD%STRAIN,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  
      CALL MPI_BCAST(RMV%SW,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(RMV%INT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

      CALL MPI_BCAST(FILES%SAVEDIR,72,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(FILES%PSI0,72,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(FILES%CHI0,72,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(FILES%PSII,72,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(FILES%CHII,72,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(FILES%NE,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      IF (MPI_RANK.NE.0) THEN
      ALLOCATE( FILES%T(FILES%NE) )
      ALLOCATE( FILES%PSI(FILES%NE) )
      ALLOCATE( FILES%CHI(FILES%NE) )
      ENDIF
      DO I=1,FILES%NE
        CALL MPI_BCAST(FILES%T(I),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(FILES%PSI(I),72,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(FILES%CHI(I),72,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
      ENDDO
  
      CALL MPI_BCAST(DIAGV%IT1,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DIAGV%IT2,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(VELMON%INT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(VELMON%N,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  
      CALL MPI_BCAST(LIN%MM,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(LIN%KK,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  
      CALL MPI_BCAST(DIAGV%IT1,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  
      IF (MPI_RANK.EQ.0) THEN
        NUM = SIZE(MONITOR_MK,1)
      ENDIF
      CALL MPI_BCAST(NUM,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      IF (MPI_RANK.NE.0) ALLOCATE(MONITOR_MK(NUM,2))
      DO I=1,NUM
        CALL MPI_BCAST(MONITOR_MK(I,1),1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(MONITOR_MK(I,2),1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      ENDDO

      CALL MPI_BCAST(MONITORDATA%START,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(MONITORDATA%INTERVAL,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

      CALL MPI_BCAST(POSTPROCESS%START,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(POSTPROCESS%FINISH,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(POSTPROCESS%JOB,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(POSTPROCESS%SLICEINT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  
      CALL MPI_BCAST(POSTPROCESS%START,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

      RETURN
      END SUBROUTINE READ_SYNC
!=======================================================================
subroutine READIN_MATLAB(FILENAME,M_IN,AK_IN,SWITCH)
! ======================================================================
      implicit NONE
      CHARACTER(LEN=*),INTENT(IN):: FILENAME
      INTEGER,INTENT(OUT):: M_IN, SWITCH
      REAL(P8),INTENT(OUT):: AK_IN

      LOGICAL:: OPEN
      CHARACTER(LEN=72),DIMENSION(20):: COM

      IREAD = 18
      DO 
        INQUIRE(UNIT=IREAD,OPENED=OPEN)
        IF(.NOT.OPEN) EXIT
        IREAD=IREAD+1
      ENDDO
      OPEN(UNIT=IREAD,FILE=FILENAME,STATUS='OLD')
      IEXEC=1
      IESAV=IECHO
      IECHO=0

      LNUMBER = 0
      CALL READCOM(LNUMBER)

      CALL READCOM(IREAD,COM)
      M_IN = ATOF(COM(1))
      CALL READCOM(IREAD,COM)
      AK_IN = ATOF(COM(1))
      CALL READCOM(IREAD,COM)
      SWITCH = ATOF(COM(1))
      CALL READCOM(IREAD,COM)

      END subroutine READIN_MATLAB
! ======================================================================
      subroutine READIN_almost(FILENAME,M_IN,AK_IN,EIG_IND,M_IN2,AK_IN2)
! ======================================================================
      implicit NONE
      CHARACTER(LEN=*),INTENT(IN):: FILENAME
      INTEGER,INTENT(OUT):: M_IN,EIG_IND,M_IN2
      REAL(P8),INTENT(OUT):: AK_IN,AK_IN2

      LOGICAL:: OPEN
      CHARACTER(LEN=72),DIMENSION(20):: COM

      IREAD = 18
      DO 
        INQUIRE(UNIT=IREAD,OPENED=OPEN)
        IF(.NOT.OPEN) EXIT
        IREAD=IREAD+1
      ENDDO
      OPEN(UNIT=IREAD,FILE=FILENAME,STATUS='OLD')
      IEXEC=1
      IESAV=IECHO
      IECHO=0

      LNUMBER = 0
      CALL READCOM(LNUMBER)

      CALL READCOM(IREAD,COM)
      M_IN = ATOF(COM(1))
      CALL READCOM(IREAD,COM)
      AK_IN = ATOF(COM(1))
      CALL READCOM(IREAD,COM)
      EIG_IND = ATOF(COM(1))
      CALL READCOM(IREAD,COM)
      M_IN2 = ATOF(COM(1))
      CALL READCOM(IREAD,COM)
      AK_IN2 = ATOF(COM(1))
      CALL READCOM(IREAD,COM)

      END subroutine READIN_almost
! ======================================================================

END MODULE MOD_INIT
