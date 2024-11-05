program test_assemble
! ======================== TEST SUBROUTINE #3 ==========================
! [DESCRIPTION]
! This subroutine tests the assembly and exchange of a scalar    
! ======================================================================
USE MPI
USE mlegs_bndmat
USE mlegs_envir
USE mlegs_genmat
USE mlegs_misc
USE mlegs_scalar
! -------------------------
implicit none
! -------------------------
integer:: ii, jj, kk
real(p8), dimension(1,1):: status
real(P8):: time_start, time_end, time0
type(scalar):: A_loc, B_loc, A_loc_new
integer:: glb_sz(3), axis_comm(3)
complex(p8), dimension(:,:,:),allocatable:: A_glb

! Initialize MPI
call MPI_init(MPI_err)
if (MPI_err.ne.MPI_success) then
write(*,*) 'MPI initialization failed.'
call MPI_abort(MPI_comm_world, error_flag_comm, MPI_err)
endif

! Get the global communicator and start timing
time_start = MPI_wtime()
comm_glb = MPI_comm_world
call MPI_comm_rank(comm_glb, rank_glb, MPI_err)
if (rank_glb.eq.0) then
    write(*,*) 'program started'
    call print_real_time()
endif

! Set 2D processor grid
call set_comm_grps(comm_glb) ! default: 2D grid / pencil
! call set_comm_grps(comm_glb,(/ 4, 1 /)) ! 1D grid / slab

! Initialize scalars with different axis communications
! Create a scalar with global size 10x12x10
glb_sz = (/ 6, 8, 10 /)
! A: dim 1 is distributed by comm_grp(1), dim 2 by comm_grp(2), dim 3 is complete
call A_loc%init(glb_sz, (/ 1, 2, 0 /))
! B: dim 1 is complete, dim 2 is distributed by comm_grp(1), dim 3 by comm_grp(2)
call B_loc%init(glb_sz, (/ 0, 1, 2 /))

! Set the scalar values
! A and B both correspond to the same global array:
do ii = 1, A_loc%loc_sz(1)
    do jj = 1, size(A_loc%e, 2)
        do kk = 1, size(A_loc%e, 3)
            A_loc%e(ii, jj, kk) = (ii + A_loc%loc_st(1)) * 10**4 + (jj + A_loc%loc_st(2)) * 10**2 + (kk + A_loc%loc_st(3))
        enddo
    enddo
enddo
do ii = 1, B_loc%loc_sz(1)
    do jj = 1, B_loc%loc_sz(2)
        do kk = 1, B_loc%loc_sz(3)
            B_loc%e(ii, jj, kk) = (ii + B_loc%loc_st(1)) * 10**4 + (jj + B_loc%loc_st(2)) * 10**2 + (kk + B_loc%loc_st(3))
        enddo
    enddo
enddo
! Copy A_loc to A_loc_new
A_loc_new = A_loc

! ================ TEST 1: EXCHANGE THE SCALAR VALUES ==================
! Exchange the scalar values
! A_loc: (1,2,0) -> (1,0,2) -> (0,1,2)
! after exchange, A_loc should be the same as B_loc
! 1. printout the local size of A before exchange
write(*,*) 'rank = ', rank_glb, 'A_loc%loc_sz = ', size(A_loc%e, 1), size(A_loc%e, 2), size(A_loc%e, 3)
call MPI_barrier(comm_glb, MPI_err)
! 2. exchange the scalar values
call A_loc%exchange(3,2) ! (1,2,0) -> (1,0,2)
! 3. printout the local size of A after 1st exchange
write(*,*) 'rank = ', rank_glb, 'A_loc%loc_sz = ', size(A_loc%e, 1), size(A_loc%e, 2), size(A_loc%e, 3)
call MPI_barrier(comm_glb, MPI_err)
! 4. exchange the scalar values
call A_loc%exchange(2,1) ! (1,0,2) -> (0,1,2)
! 5. printout the local size of A after 2nd exchange
write(*,*) 'rank = ', rank_glb, 'A_loc%loc_sz = ', size(A_loc%e, 1), size(A_loc%e, 2), size(A_loc%e, 3)
call MPI_barrier(comm_glb, MPI_err)

! Assemble local scalars into global arrays and save them to files
! note: currently, assemble only supports arrays with either 1st or 2nd axis complete
! call save_glb(A_loc%assemble(1), "A_glb_012")
call save_glb(A_loc%assemble(), "A_glb_012")
! call save_glb(B_loc%assemble(1), "B_glb_012")
call save_glb(B_loc%assemble(), "B_glb_012")

! ================ TEST 2: ASSEMBLE THE SCALAR VALUES ==================
! Exchange the scalar values
! 1. printout the local size of A before exchange
write(*,*) 'rank = ', rank_glb, 'A_loc_new%loc_sz = ', size(A_loc_new%e, 1), size(A_loc_new%e, 2), size(A_loc_new%e, 3)
call MPI_barrier(comm_glb, MPI_err)
! 2. assemble the scalar values -> assemble automatically exchanges the scalar values first before assembling
! call save_glb(A_loc_new%assemble(3), "A_glb_120")
call save_glb(A_loc_new%assemble(), "A_glb_120")

! 2. exchange the scalar values
call A_loc_new%exchange(3,2) ! (1,2,0) -> (1,0,2)
! 3. printout the local size of A after exchange
write(*,*) 'rank = ', rank_glb, 'A_loc_new%loc_sz = ', size(A_loc_new%e, 1), size(A_loc_new%e, 2), size(A_loc_new%e, 3)
call MPI_barrier(comm_glb, MPI_err)
! if (rank_glb.eq.0) write(*,*) int(real(A_loc_new%e))
! 4. assemble the scalar values
! call save_glb(A_loc_new%assemble(2), "A_glb_102")
call save_glb(A_loc_new%assemble(), "A_glb_102")

! 2. exchange the scalar values
call A_loc_new%exchange(2,3) ! (1,0,2) -> (1,2,0)
call A_loc_new%exchange(3,1) ! (1,2,0) -> (0,2,1)
! 3. printout the local size of A after exchange
write(*,*) 'rank = ', rank_glb, 'A_loc_new%loc_sz = ', size(A_loc_new%e, 1), size(A_loc_new%e, 2), size(A_loc_new%e, 3)
call MPI_barrier(comm_glb, MPI_err)
! 4. assemble the scalar values
! call save_glb(A_loc_new%assemble(1), "A_glb_021")
call save_glb(A_loc_new%assemble(), "A_glb_021")

! 2. exchange the scalar values
call A_loc_new%exchange(1,2) ! (0,2,1) -> (2,0,1)
! 3. printout the local size of A after exchange
write(*,*) 'rank = ', rank_glb, 'A_loc_new%loc_sz = ', size(A_loc_new%e, 1), size(A_loc_new%e, 2), size(A_loc_new%e, 3)
call MPI_barrier(comm_glb, MPI_err)
! if (rank_glb.eq.0) write(*,*) int(real(A_loc_new%e))
! 4. assemble the scalar values
! call save_glb(A_loc_new%assemble(2), "A_glb_201")
call save_glb(A_loc_new%assemble(), "A_glb_201")

! =============== TEST 3: DISASSEMBLE THE SCALAR VALUES ================
allocate(A_glb(glb_sz(1), glb_sz(2), glb_sz(3)))

! A_glb = A_loc_new%assemble(2) ! (2,0,1)
! A_loc = A_loc_new
! call A_loc_new%dealloc()
! call A_loc_new%init(glb_sz, (/ 2, 0, 1 /))
! call A_loc_new%disassemble(2, A_glb)

call A_loc_new%exchange(2,1) ! (2,0,1) -> (0,2,1)
! A_glb = A_loc_new%assemble(1)
A_glb = A_loc_new%assemble()
A_loc = A_loc_new
! A_loc%e = A_loc%e * 2
call A_loc_new%dealloc()

! Disassemble the global array into local scalars
call A_loc_new%init(glb_sz, (/ 0, 2, 1 /))
! call A_loc_new%disassemble(A_glb, 1)
call A_loc_new%disassemble(A_glb)

do ii = 1, A_loc_new%loc_sz(1)
    do jj = 1, size(A_loc_new%e, 2)
        do kk = 1, size(A_loc_new%e, 3)
            if (A_loc_new%e(ii, jj, kk) .ne. A_loc%e(ii, jj, kk) ) then
                write(*,*) 'rank = ', rank_glb, 'A_loc_new%e(ii, jj, kk) = ', A_loc_new%e(ii, jj, kk), &
                'A_loc%e(ii, jj, kk) = ', A_loc%e(ii, jj, kk)
            endif
        enddo
    enddo
enddo

! =========================== END OF TESTS =============================
! Deallocate local scalars and global arrays
call A_loc%dealloc()
call B_loc%dealloc()
deallocate(A_glb)

! Final printout
call MPI_barrier(comm_glb, MPI_err)
time_end = MPI_wtime()
if (rank_glb.eq.0) then
    write(*,*) 'program started'
    call print_real_time()
    write(*,*) 'execution time: ', time_end - time_start, 'seconds'
endif

! Synchronize and finalize MPI
call MPI_barrier(comm_glb, MPI_err)
call MPI_finalize(MPI_err)

! ======================================================================
contains
! ======================================================================

! Subroutine to save mode data to a file
subroutine save_mode(t, psi_new, chi_new, m_i, k_i)
! ======================================================================
    complex(p8), dimension(:):: psi_new, chi_new
    real:: t
    integer:: m_i, k_i
    integer:: nn

    open(unit=777, file='mode_track_' // itoa(m_i) // '_' // itoa(k_i) // '.dat', &
            status='unknown', action='write', position='append')

    write(777, 320) t, psi_new(1:size(psi_new)), chi_new(1:size(psi_new))

    close(777)
320   format(f10.3, ',', (s, e14.6e3, sp, e14.6e3, 'i'), *(',',s, e14.6e3, sp, e14.6e3, 'i'))
end subroutine save_mode
! ======================================================================

! Subroutine to save global data to a file
subroutine save_glb(glb_data, filename)
! ======================================================================
    complex(p8), dimension(:,:,:):: glb_data
    character(*):: filename
    integer:: im, ik

    if (rank_glb.eq.0) then
        open(unit=777, file=adjustl(trim(filename)) // '.dat', &
                status='unknown', action='write')

        do im = 1, size(glb_data, 2)
            do ik = 1, size(glb_data, 3)
                write(777, 256) real(glb_data(:, im, ik))
            enddo
        enddo

        close(777)
    endif

    call MPI_barrier(comm_glb, MPI_err)
256   format((f10.0), *(',', f10.0))
end subroutine save_glb
! ======================================================================

end program test_assemble
