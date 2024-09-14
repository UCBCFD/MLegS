program vort5
! -------------------------
    USE OMP_LIB
    USE MPI
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
! -------------------------
implicit none
! -------------------------
type(scalar):: psi_p,chi_p,psi_pn,chi_pn,del2chi_pn
type(scalar):: psi_m,chi_m,psi_mn,chi_mn,del2chi_mn
type(scalar):: psi_t,chi_t,del2chi_t
integer:: ii,MM,KK,MBAR,KBAR,MEND1,KEND1,MEND2,KEND2
logical:: MM_true,KK_true,MEND1_true,KEND1_true,MEND2_true,KEND2_true,MEND2_neg
REAL(P8):: time_start,time_end
complex(P8),dimension(:,:),ALLOCATABLE:: PERTURB_M1, PERTURB_M2
real(p8),dimension(1,1):: status

character(LEN=200) :: FileName_p, FileName_m, FileName_C, num
! complex(p8),dimension(:,:),allocatable:: prod_result
! complex(p8),dimension(:),allocatable:: eig_1_psi,eig_2_psi
! complex(p8),dimension(:),allocatable:: eig_1_chi,eig_2_chi

CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
    WRITE(*,*) 'The threading support is less than that demanded.'
    CALL MPI_ABORT(MPI_COMM_WORLD,1,IERR)
ENDIF
CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, IERR)

CALL READCOM('NOECHO')
CALL READIN(5)
CALL LEGINIT()

status(1,1)=0
if (MPI_RANK.eq.0) then
    WRITE(*,*) 'PROGRAM STARTED'
    call msave(status, 'status.dat')
    CALL PRINT_REAL_TIME()
    WRITE(*,*) 'psi/chi(+) filename: '
    READ(*,10) FileName_p
    WRITE(*,*) 'psi/chi(-) filename: '
    READ(*,10) FileName_m
    WRITE(*,*) 'MM(index): '
    READ(*,*) MM
    WRITE(*,*) 'KK(index): '
    READ(*,*) KK
    WRITE(*,*) 'MBAR(index): '
    READ(*,*) MBAR
    WRITE(*,*) 'KBAR(index): '
    READ(*,*) KBAR
    WRITE(*,*) 'Uz: ', NADD%UZ
endif
CALL MPI_BCAST(FileName_p,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(FileName_m,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(MM,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(KK,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(MBAR,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(KBAR,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
10 FORMAT(A200)

!> initial values
call allocate(psi_p); call allocate(chi_p)
call allocate(psi_m); call allocate(chi_m)

call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_p))//"_psi.dat",psi_p)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_p))//"_chi.dat",chi_p)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_m))//"_psi.dat",psi_m)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_m))//"_chi.dat",chi_m)

call allocate(psi_t); call allocate(chi_t); call allocate(del2chi_t)
call allocate(psi_pn); call allocate(chi_pn); call allocate(del2chi_pn)
call allocate(psi_mn); call allocate(chi_mn); call allocate(del2chi_mn)
allocate(PERTURB_M1(2*NRCHOPDIM,2*NRCHOPDIM)); PERTURB_M1 = 0.D0
allocate(PERTURB_M2(2*NRCHOPDIM,2*NRCHOPDIM)); PERTURB_M2 = 0.D0

MEND1 = MM + MBAR - 1; KEND1 = KK + KBAR - 1
MEND2 = MM - MBAR; KEND2 = KK - KBAR

IF (MEND2.LT.0) THEN
    MEND2 = -MEND2
    MEND2_neg = .TRUE.
ELSE
    MEND2_neg = .FALSE.
ENDIF
MEND2 = MEND2 + 1
KEND2 = KEND2 + 1
IF (KEND2.LE.0) THEN
    KEND2 = KEND2 + 2*NXCHOP-1
ENDIF

MM = MM - psi_t%INTH; MM_true = ((MM.GT.0).AND.(MM.LE.SIZE(psi_t%E,2)))
KK = KK - psi_t%INX ; KK_true = ((KK.GT.0).AND.(KK.LE.SIZE(psi_t%E,3)))
MEND1 = MEND1 - psi_t%INTH; MEND1_true = ((MEND1.GT.0).AND.(MEND1.LE.SIZE(psi_t%E,2)))
KEND1 = KEND1 - psi_t%INX ; KEND1_true = ((KEND1.GT.0).AND.(KEND1.LE.SIZE(psi_t%E,3)))
MEND2 = MEND2 - psi_t%INTH; MEND2_true = ((MEND2.GT.0).AND.(MEND2.LE.SIZE(psi_t%E,2)))
KEND2 = KEND2 - psi_t%INX ; KEND2_true = ((KEND2.GT.0).AND.(KEND2.LE.SIZE(psi_t%E,3)))

IF ((MM_true).AND.(KK_true)) THEN
    WRITE(*,*) '{MM,KK}',M(MM+psi_t%INTH),AK(MM+psi_t%INTH,KK+psi_t%INX)
ENDIF
IF ((MEND1_true).AND.(KEND1_true)) THEN
    WRITE(*,*) '{MBAR,KBAR}+{MM,KK}',M(MEND1+psi_t%INTH),AK(MEND1+psi_t%INTH,KEND1+psi_t%INX)
ENDIF
IF ((MEND2_true).AND.(KEND2_true)) THEN
    WRITE(*,*) '-{MBAR,KBAR}+{MM,KK}',M(MEND2+psi_t%INTH),AK(MEND2+psi_t%INTH,KEND2+psi_t%INX)
ENDIF

write(*,*) '{m1,k1}',MEND1,KEND1,MPI_RANK
write(*,*) '{m2,k2}',MEND2,KEND2,MPI_RANK
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
! CALL MPI_ABORT(MPI_COMM_IVP,1,IERR)

DO II = 1,NRCHOPDIM

    ! PSI
    ! ======================================================================
    psi_t%e = 0.d0
    chi_t%e = 0.d0
    del2chi_t%e = 0.d0

    IF ((MM_true).AND.(KK_true)) THEN
        psi_t%E(II,MM,KK) = 1.D0
        psi_p%E(:,MM,KK) = psi_t%E(:,MM,KK)
        psi_m%E(:,MM,KK) = psi_t%E(:,MM,KK)
        chi_p%E(:,MM,KK) = chi_t%E(:,MM,KK)
        chi_m%E(:,MM,KK) = chi_t%E(:,MM,KK)
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    CALL NONLIN(psi_p,chi_p,psi_pn,chi_pn)
    CALL NONLIN(psi_m,chi_m,psi_mn,chi_mn)
    CALL DEL2(chi_pn,del2chi_pn)
    CALL DEL2(chi_mn,del2chi_mn)

    IF ((MEND1_true).AND.(KEND1_true)) THEN
        PERTURB_M1(1:NRCHOPDIM,II) = psi_pn%E(:,MEND1,KEND1)-psi_mn%E(:,MEND1,KEND1)
        PERTURB_M1(NRCHOPDIM+1:,II) = del2chi_pn%E(:,MEND1,KEND1)-del2chi_mn%E(:,MEND1,KEND1)
    ENDIF
    IF ((MEND2_true).AND.(KEND2_true)) THEN
        PERTURB_M2(1:NRCHOPDIM,II) = psi_pn%E(:,MEND2,KEND2)-psi_mn%E(:,MEND2,KEND2)
        PERTURB_M2(NRCHOPDIM+1:,II) = del2chi_pn%E(:,MEND2,KEND2)-del2chi_mn%E(:,MEND2,KEND2)
    ENDIF

    ! CHI
    ! ======================================================================
    psi_t%e = 0.d0
    chi_t%e = 0.d0
    del2chi_t%e = 0.d0

    IF ((MM_true).AND.(KK_true)) THEN
        del2chi_t%E(II,MM,KK) = 1.D0
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    CALL IDEL2(del2chi_t,chi_t)
    IF ((MM_true).AND.(KK_true)) THEN
        psi_p%E(:,MM,KK) = psi_t%E(:,MM,KK)
        psi_m%E(:,MM,KK) = psi_t%E(:,MM,KK)
        chi_p%E(:,MM,KK) = chi_t%E(:,MM,KK)
        chi_m%E(:,MM,KK) = chi_t%E(:,MM,KK)
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

    CALL NONLIN(psi_p,chi_p,psi_pn,chi_pn)
    CALL NONLIN(psi_m,chi_m,psi_mn,chi_mn)
    CALL DEL2(chi_pn,del2chi_pn)
    CALL DEL2(chi_mn,del2chi_mn)

    IF ((MEND1_true).AND.(KEND1_true)) THEN
        PERTURB_M1(1:NRCHOPDIM,NRCHOPDIM+II) = psi_pn%E(:,MEND1,KEND1)-psi_mn%E(:,MEND1,KEND1)
        PERTURB_M1(NRCHOPDIM+1:,NRCHOPDIM+II) = del2chi_pn%E(:,MEND1,KEND1)-del2chi_mn%E(:,MEND1,KEND1)
    ENDIF
    IF ((MEND2_true).AND.(KEND2_true)) THEN
        PERTURB_M2(1:NRCHOPDIM,NRCHOPDIM+II) = psi_pn%E(:,MEND2,KEND2)-psi_mn%E(:,MEND2,KEND2)
        PERTURB_M2(NRCHOPDIM+1:,NRCHOPDIM+II) = del2chi_pn%E(:,MEND2,KEND2)-del2chi_mn%E(:,MEND2,KEND2)
    ENDIF

ENDDO

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

IF ((MEND1_true).AND.(KEND1_true)) THEN
    call save_perturbm(PERTURB_M1,MEND1+psi_t%INTH-1,KEND1+psi_t%INX-1)
ENDIF
IF ((MEND2_true).AND.(KEND2_true)) THEN
    if (MEND2_neg) PERTURB_M2 = CONJG(PERTURB_M2)
    call save_perturbm(PERTURB_M2,MEND2+psi_t%INTH-1,KEND2+psi_t%INX-1)
ENDIF

CALL DEALLOCATE(psi_t); CALL DEALLOCATE(chi_t); CALL DEALLOCATE(del2chi_t)
CALL DEALLOCATE(psi_p); CALL DEALLOCATE(chi_p)
CALL DEALLOCATE(psi_pn); CALL DEALLOCATE(chi_pn); CALL DEALLOCATE(del2chi_pn)
CALL DEALLOCATE(psi_m); CALL DEALLOCATE(chi_m)
CALL DEALLOCATE(psi_mn); CALL DEALLOCATE(chi_mn); CALL DEALLOCATE(del2chi_mn)
DEALLOCATE(PERTURB_M1,PERTURB_M2)

time_end = mpi_wtime()
IF (MPI_RANK.eq.0) THEN
    WRITE(*,*) 'PROGRAM COMPLETED'
    CALL PRINT_REAL_TIME()
    WRITE(*,*) 'EXECUTION TIME: ',time_end-time_start,'seconds'
ENDIF
call MPI_BARRIER(MPI_COMM_IVP,IERR)
call MPI_FINALIZE(IERR)

! ======================================================================
contains
subroutine save_perturbm(perturb_m,m_1,k_1)
complex(p8),DIMENSION(:,:):: perturb_m
integer:: m_1,k_1
integer:: nn

open(UNIT=777,FILE='IVP_PERTURBM_'//ITOA3(m_1)//'_'//ITOA3(k_1)//'.dat',&
&STATUS='UNKNOWN')

do nn = 1,2*NRCHOPDIM
    WRITE(777,320) perturb_m(nn,:)
enddo

close(777)
320 FORMAT((S,E14.6E3,SP,E14.6E3,'i'),*(','S,E14.6E3,SP,E14.6E3,'i'))

end subroutine save_perturbm


end program vort5
