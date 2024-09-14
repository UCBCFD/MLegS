program vort4
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
type(scalar):: psi,chi,del2chi
type(scalar):: ror, rop, oz
integer:: it,pp,ii,jj,kk
REAL(P8):: time_start,time_end
real(p8),dimension(1,1):: status
real(p8), allocatable, dimension(:,:):: tempMonitor_mk, EMK
real(p8):: time0, tmp, vel_o
character(LEN=200) :: FileName_R, FileName_L, FileName_C, num
complex(p8),dimension(:,:),allocatable:: prod_result
complex(p8),dimension(:),allocatable:: eig_1_psi,eig_2_psi
complex(p8),dimension(:),allocatable:: eig_1_chi,eig_2_chi


CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,MPI_THREAD_MODE,IERR)
IF (MPI_THREAD_MODE.LT.MPI_THREAD_SERIALIZED) THEN
    WRITE(*,*) 'The threading support is lesser than that demanded.'
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
    WRITE(*,*) 'psi/chi(R) filename: '
    READ(*,10) FileName_R
    WRITE(*,*) 'Uz: ', NADD%UZ
endif
CALL MPI_BCAST(FileName_R,200,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
10 FORMAT(A200)

! degen pair
call allocate(psi)
call allocate(chi)
call allocate(del2chi)

!> initial values
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"_psi.dat",psi)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//TRIM(ADJUSTL(FileName_R))//"_chi.dat",chi)

call allocate(ror)
call allocate(rop)
call allocate( oz)

call CHOPSET(3)

CALL PC2VOR(psi,chi,ror,rop,oz)
CALL RTRAN(ror,1)
CALL RTRAN(rop,1)
CALL RTRAN( oz,1)
CALL PROJECT(ror,rop,oz,psi,del2chi)
CALL IDEL2(del2chi,chi)
call CHOPSET(-3)

allocate(EMK(NTCHOP,NXCHOP))
write(*,*) MPI_RANK
EMK = ENERGY_SPEC(psi,chi)
WRITE(*,*) EMK(5,33)

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL DEALLOCATE(psi)
CALL DEALLOCATE(chi)
CALL DEALLOCATE(del2chi)

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
subroutine save_compare(psi_new,chi_new,psi_old,chi_old,m_i,k_i)
complex(p8),DIMENSION(:):: psi_new,chi_new,psi_old,chi_old
integer:: m_i,k_i
integer:: nn

open(UNIT=777,FILE='mode_compare_'//ITOA3(m_i)//'_'//ITOA3(k_i)//'.dat',&
&STATUS='UNKNOWN')

do nn = 1,size(psi_new)
    WRITE(777,320) psi_new(nn),psi_old(nn),chi_new(nn),chi_old(nn)
enddo

close(777)
320 FORMAT((S,E14.6E3,SP,E14.6E3,'i'),*(','S,E14.6E3,SP,E14.6E3,'i'))

end subroutine save_compare


end program vort4
