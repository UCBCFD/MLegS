program vort
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
type(scalar):: psi,chi,psio,chio,dpsi,dchi,psi_o,chi_o
integer:: iii,it,mm,kk,pp
REAL(P8):: time_start,time_end
real(p8),dimension(1,1):: status
real(p8), allocatable, dimension(:,:):: tempMonitor_mk
real(p8):: time0, tmp
character(len=72) :: psiFileName, chiFileName, num  

!    ! debug:
!   COMPLEX(P8),DIMENSION(:,:,:),ALLOCATABLE:: GLB_ARRAY


! CALL MPI_INIT(IERR)
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
endif

call allocate(psi)
call allocate(chi)
call allocate(psio)
call allocate(chio)
call allocate(psi_o)
call allocate(chi_o)
call allocate(dpsi)
call allocate(dchi)

!> initial values
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%psii,psi)
call MLOAD(TRIM(ADJUSTL(FILES%SAVEDIR))//FILES%chii,chi)

time_start = mpi_wtime()
call PRINT_ENERGY_SPECTRUM(psi,chi,1)

!> startup : very first time step
call diagnost(psi,chi)

!dpsi and dchi are initially empty, then they are assigned 
!the nonlinear part of the first step
call rich(psi,chi,dpsi,dchi)
call diagnost(psi,chi)

!  psi_o = psi
!  chi_o = chi
!  time0 = tim%t

iii = tim%limit/tim%dt
files%n = 1

do it=1,iii

   if(monitorData%start.eq.tim%n) then 
      ! use this to calculate eigenvalues. this is the psi/chi at t1
      psi_o = psi
      chi_o = chi
      time0 = tim%t
   end if
   
   !> integrate
   call adamsb(psi,chi,dpsi,dchi)
   call diagnost(psi,chi)
   
   !> monitor eigenmodes
   call WRITE_ENERGY_DATA_TO_FILE(psi,chi)
!     call WRITE_MODE_DATA_TO_FILE(psi,chi)
   tmp = (real(tim%n)-real(monitorData%start))/real(monitorData%interval)
   if (tmp.ge.1.0_p8 .and. abs(tmp-int(tmp)).lt.1e-08) then
      call egrowth_modified(psi_o,chi_o,psi,chi,time0,tim%t)
      call monitor_eig(psi_o,chi_o,psi,chi,time0,tim%t)
   !   call PRINT_ENERGY_SPECTRUM(psi,chi,1)
      call PRINT_ENERGY_SPECTRUM(psi,chi,0)      
   end if
   
   !> freestream adjust
   if(adv%sw.eq.1 .and. mod(it,adv%int).eq.0) then
      call freeadj
      call rich(psi,chi,dpsi,dchi)
      call diagnost(psi,chi)
   endif
   
   !> vorticity removal
   if(rmv%sw.ne.0 .and. mod(it,rmv%int).eq.0) then
      call remove(psi,chi)
      call rich(psi,chi,dpsi,dchi)
      call diagnost(psi,chi)
   endif
   
   !> hyperviscosity adjust
   if(visc%adjsw.ne.0 .and. mod(it,visc%adjint).eq.0) then
      call hypadj(psi,chi,visc%nup)
      call diagnost(psi,chi)
   endif
   
   !> output
   if(files%t(files%n).le.tim%t) then
      call msave(psi, TRIM(ADJUSTL(FILES%SAVEDIR))//files%psi(files%n))
      call msave(chi, TRIM(ADJUSTL(FILES%SAVEDIR))//files%chi(files%n))
      files%n = files%n + 1
      if(files%n > files%ne) goto 999
      
      ! 'status.dat' is currently un-used
      ! Setting status.dat's stored value to 1 stops the program after 
      ! storing the results at the previous savepoint
      if (MPI_RANK.eq.0) then 
      call mload('status.dat',status)
      if(status(1,1).eq.1.0) goto 999
      endif

   endif
enddo

999 continue
time_end = mpi_wtime()
!  ALLOCATE(GLB_ARRAY(NRCHOPDIM,NTCHOPDIM,NXCHOPDIM))
!  call MASSEMBLE(PSI%E,GLB_ARRAY,1)
if (MPI_RANK.eq.0) THEN
print *,tim%n,' steps'
WRITE(*,*) 'PROGRAM STARTED'
CALL PRINT_REAL_TIME()
WRITE(*,*) 'EXECUTION TIME: ',time_end-time_start,'seconds'
! CALL MCAT(GLB_ARRAY(:,1,1))
ENDIF
call MPI_BARRIER(MPI_COMM_IVP,IERR)
!stop
!CALL MPI_ABORT(MPI_COMM_IVP,1,IERR)
call MPI_FINALIZE(IERR)

end program vort
