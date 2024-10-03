submodule (mlegs_scalar) mlegs_scalar_mpi
  implicit none

contains

    !> initialise a scalar
    module procedure scalar_initialise
        integer :: i, comm_grp_idx, nproc_i, rank_i, loc_sz_i, loc_st_i

        this%glb_sz = glb_sz

        ! multi-dimension support
        select case (size(axis_comm))
            ! 1D: r
            case (1)
                if ((glb_sz(2).ne.1).or.(glb_sz(3).ne.1)) then
                    if (rank_glb.eq.0) write(*,*) "ERROR: scalar_initialise is expecting 1D data"
                    call mpi_abort(comm_glb,error_flag_comm,mpi_ierr)
                endif
                if ((is_warning).and.(rank_glb.eq.0)) write(*,*) "WARNING: scalar_initialise detects 1D (radial) configuration"
                this%axis_comm = 0
                this%axis_comm(1) = axis_comm(1)
            ! 2D: r-theta or r-z
            case (2)
                ! find the dimension that is not present: theta or z
                i = findloc(glb_sz,1)
                if (i.eq.1) then
                    if (rank_glb.eq.0) write(*,*) "ERROR: scalar_initialise does not support 2D (theta-z) configuration"
                    call mpi_abort(comm_glb,error_flag_comm,mpi_ierr)
                endif
                ! set axis_comm
                this%axis_comm(i) = 0
                this%axis_comm(3-i) = axis_comm(2)
                this%axis_comm(1) = axis_comm(1)
            ! 3D: r-theta-z
            case (3)
                this%axis_comm = axis_comm
        end select

        !> along each dimension
        do i = 1,3
            comm_grp_idx = this%axis_comm(i)
            !> dimension that is not distributed
            if (comm_grp_idx.eq.0) then
                this%loc_sz(i) = glb_sz(i)
                this%loc_st(i) = 0
            !> dimension that gets distributed
            else
                if (comm_grp_idx.gt.count(comm_grps.ne.0)) then
                    if (rank_glb.eq.0) write(*,*) "ERROR: scalar_initialise does not have enough communicator groups available"
                    call mpi_abort(comm_glb,error_flag_comm,mpi_ierr)
                endif
                rank_i = rank_grps(comm_grp_idx)
                nproc_i = nprocs_grps(comm_grp_idx)
                call decompose(glb_sz(i),nproc_i,rank_i,loc_sz_i,loc_st_i)
                this%loc_sz(i) = loc_sz_i
                this%loc_st(i) = loc_st_i
            endif
        enddo

        !> allocate
        call this%alloc()

    end procedure

    !> allocate a scalar
    module procedure scalar_alloc
        allocate(this%e(this%loc_sz(1),this%loc_sz(2),this%loc_sz(3)))
    end procedure

    !> deallocate a scalar
    module procedure scalar_dealloc
        if (.not.associated(this%e)) then
            if ((is_warning).and.(rank_glb.eq.0)) write(*,*) "WARNING: scalar_dealloc tries to deallocate memory that has not been allocated"
            return
        endif
        deallocate(this%e)
        nullify(this%e)
    end procedure

    !> re-orient distributed data
    module procedure scalar_exchange !(this,axis_old,axis_new)
        integer :: this_comm, loc_sz_new(3), i
        integer, dimension(:), allocatable :: subarray_type_old, subarray_type_new
        complex(p8), dimension(:,:,:), target, allocatable :: e_new

        !> check whether the old dimension is distributed or not
        if (this%axis_comm(axis_old).ne.0) then
            if (rank_glb.eq.0) write(*,*) "ERROR: scalar_exchange requires the data to be non-distributed along the old dimension: ",itoa(axis_old)
            call mpi_abort(comm_glb,error_flag_comm,mpi_ierr)
        endif

        !> obtain communicator group
        this_comm = comm_grps(this%axis_comm(axis_new))

        !> create subarray datatype for before and after swap configuration
        subarray_type_old = create_new_type3D(this_comm,this%loc_sz,axis_old,scalar_element_type)
        subarray_type_new = create_new_type3D(this_comm,loc_sz_new,axis_new,scalar_element_type)

        !> create after-swap new local size and data
        loc_sz_new = this%loc_sz
        loc_sz_new(axis_new) = this%glb_sz(new)
        loc_sz_new(axis_old) = local_size(this%glb_sz(axis_old),this_comm)
        allocate(e_new(loc_sz_new(1),loc_sz_new(2),loc_sz_new(3)))

        !> perform swap
        call exchange_3dcomplex_fast(this_comm,this%e,subarray_type_old,e_new,subarray_type_new)

        !> reassign pointer, loc_sz, loc_st, and axis_comm
        call this%dealloc()
        this%e => e_new
        this%loc_sz = loc_sz_new
        this%loc_st(axis_old) = local_index(this%glb_sz(axis_old),this_comm)
        this%loc_st(axis_new) = 0
        this%axis_comm(axis_old) = this%axis_comm(axis_new)
        this%axis_comm(axis_new) = 0

        !> clean-up
        do i = 1,size(subarray_type_old)
            call mpi_type_free(subarray_type_old(i),mpi_ierr)
            call mpi_type_free(subarray_type_new(i),mpi_ierr)
        enddo
        deallocate(subarray_type_new,subarray_type_old)
    end procedure

    !> assemble distributed data into a single proc
    module procedure scalar_assemble
        integer :: n1_glb, n2_glb, n3_glb, n1_loc, n2_loc, n3_loc
        integer :: element_size, subcomm_l, subcomm_r
        integer(kind = MPI_ADDRESS_KIND) :: extend_size
        integer :: subarray_type_temp, subarray_type, displacement_loc, recvcount_loc
        integer, dimension(:), allocatable :: displacement, recvcount
        complex(p8), dimension(:,:,:), allocatable :: local_array_temp, subglobal_array, global_array_temp

        if (cp8_size.eq.0) then
            call mpi_type_size(scalar_element_type,cp8_size,mpi_ierr)
        endif
        element_size = cp8_size
        
        n1_glb = this%glb_sz(1); n2_glb = this%glb_sz(2); n3_glb = this%glb_sz(3)
        n1_loc = this%loc_sz(1); n2_loc = this%loc_sz(2); n3_loc = this%loc_sz(3)

        allocate(array_glb(n1_glb,n2_glb,n3_glb))

        if (axis.eq.1) then
            ! axis,    SUBCOMMS_L  ,   SUBCOMMS_R
            !  N1 , N2/axis_comm(2), N3/axis_comm(3)
            subcomm_l = comm_grps(this%axis_comm(2))
            subcomm_r = comm_grps(this%axis_comm(3))

            allocate(subglobal_array(n1_glb, n2_glb, n3_loc))
            allocate(local_array_temp(n1_glb, n3_loc, n2_loc))
            !      #1   #2   #3
            ! old: N1 x N2 x N3
            ! new: N1 x N3 x N2
            local_array_temp = reshape(this%e, shape(local_array_temp), order = [1,3,2])

            !> create custom MPI datatype
            call mpi_type_vector(n3_loc, n1_glb, n2_glb*n1_glb, scalar_element_type, subarray_type_temp, mpi_ierr)
            extend_size = element_size * n1_glb
            recvcount_loc = n2_loc
            displacement_loc = local_index(n2_glb, subcomm_l)

        elseif (axis.eq.2) then
            !   SUBCOMMS_L   , axis,   SUBCOMMS_R
            ! N1/axis_comm(1),  N2 , N3/axis_comm(3)
            subcomm_l = comm_grps(this%axis_comm(1))
            subcomm_r = comm_grps(this%axis_comm(3))

            allocate(subglobal_array(n2_glb, n1_glb, n3_loc))
            allocate(local_array_temp(n2_glb, n3_loc, n1_loc))
            !      #1   #2   #3
            ! old: N1 x N2 x N3
            ! new: N2 x N3 x N1
            local_array_temp = reshape(this%e, shape(local_array_temp), order = [2,3,1])

            call mpi_type_vector(n3_loc, n2_glb, n2_glb*n1_glb, scalar_element_type, subarray_type_temp, mpi_ierr)
            extend_size = element_size * n2_glb
            recvcount_loc = n1_loc
            displacement_loc = local_index(n1_glb, subcomm_l)

        else 
            if (rank_glb.eq.0) write(*,*) "ERROR: scalar_assemble does not support the current axis"
            call mpi_abort(comm_glb,error_flag_comm,mpi_ierr)
        endif

        !> create elemental subarray of local array
        call mpi_type_create_resized(subarray_type_temp, 0, extend_size, subarray_type, mpi_ierr)
        call mpi_type_commit(subarray_type, mpi_ierr)

        !> in ea sub_comm_l, gather local array into subglobal array in subcomm_l's proc#0
        allocate(recvcount(0:count_proc(subcomm_l)-1))
        allocate(displacement(0:count_proc(subcomm_l)-1))
        call mpi_allgather(recvcount_loc,1,MPI_INTEGER,recvcount,1,MPI_INTEGER,subcomm_l,mpi_ierr)
        call mpi_allgather(displacement_loc,1,MPI_INTEGER,displacement,1,MPI_INTEGER,subcomm_l,mpi_ierr)
        call mpi_gatherv(local_array_temp, n1_loc*n2_loc*n3_loc, scalar_element_type, &
            subglobal_array, recvcount, displacement, subarray_type, 0, subcomm_l, mpi_ierr)

        !> in the subcomm_r that contains subcomm_l's proc#0:
        !> gather subglobal array into global array that is stored in subcomm_r's proc#0 (=> global proc#0)
        if (local_proc(subcomm_l).eq.0) then

            deallocate(recvcount,displacement)
            allocate(recvcount(0:count_proc(subcomm_r)-1))
            allocate(displacement(0:count_proc(subcomm_r)-1))
            recvcount_loc = n1_glb*n2_glb*n3_loc
            displacement_loc = n1_glb*n2_glb*local_index(n3_glb,subcomm_r)

            call mpi_allgather(recvcount_loc,1,MPI_INTEGER,recvcount,1,MPI_INTEGER,subcomm_r,mpi_ierr)
            call mpi_allgather(displacement_loc,1,MPI_INTEGER,displacement,1,MPI_INTEGER,subcomm_r,mpi_ierr)

            if (axis.eq.1) then
                call mpi_gatherv(subglobal_array, n1_glb*n2_glb*n3_loc, scalar_element_type, &
                    array_glb, recvcount, displacement, scalar_element_type, 0, subcomm_r, mpi_ierr)
            elseif (axis.eq.2) then
                allocate(global_array_temp(n2_glb,n1_glb,n3_glb))
                call mpi_gatherv(subglobal_array, n1_glb*n2_glb*n3_loc, scalar_element_type, &
                    global_array_temp, recvcount, displacement, scalar_element_type, 0, subcomm_r, mpi_ierr)
                array_glb = reshape(global_array_temp, shape(array_glb), order = [2,1,3])
                deallocate(global_array_temp)
            endif

        endif
             
        !> clean-up
        call mpi_type_free(subarray_type, mpi_ierr)
        call mpi_type_free(subarray_type_temp, mpi_ierr)
        deallocate(local_array_temp, subglobal_array, recvcount, displacement)

    end procedure

    !> disassemble data into distributed procs
    module procedure scalar_disassemble
        integer :: n1_glb, n2_glb, n3_glb, n1_loc, n2_loc, n3_loc
        integer :: element_size, subcomm_l, subcomm_r
        integer(kind = MPI_ADDRESS_KIND) :: extend_size
        integer :: subarray_type_temp, subarray_type, displacement_loc, recvcount_loc
        integer, dimension(:), allocatable :: displacement, recvcount
        complex(p8), dimension(:,:,:), allocatable :: local_array_temp, subglobal_array, global_array_temp

        if (cp8_size.eq.0) then
            call mpi_type_size(scalar_element_type,cp8_size,mpi_ierr)
        endif
        element_size = cp8_size
        
        n1_glb = this%glb_sz(1); n2_glb = this%glb_sz(2); n3_glb = this%glb_sz(3)
        n1_loc = this%loc_sz(1); n2_loc = this%loc_sz(2); n3_loc = this%loc_sz(3)

        !> check global array dimensions
        if ((size(array_glb,1).ne.n1_glb).or.(size(array_glb,2).ne.n2_glb).or.(size(array_glb,3).ne.n3_glb)) then
            if (rank_glb.eq.0) write(*,*) "ERROR: scalar_disassemble has inconsistent dimension between array_glb and this%glb_sz"
            return
        endif

        !> define left and right subcommunicator groups
        if (axis.eq.1) then
            ! axis,    SUBCOMMS_L  ,   SUBCOMMS_R
            !  N1 , N2/axis_comm(2), N3/axis_comm(3)
            subcomm_l = comm_grps(this%axis_comm(2))
            subcomm_r = comm_grps(this%axis_comm(3))
        elseif (axis.eq.2) then
            !   SUBCOMMS_L   , axis,   SUBCOMMS_R
            ! N1/axis_comm(1),  N2 , N3/axis_comm(3)
            subcomm_l = comm_grps(this%axis_comm(1))
            subcomm_r = comm_grps(this%axis_comm(3))
        else
            if (rank_glb.eq.0) write(*,*) "ERROR: scalar_disassemble does not support the current axis"
            return
        endif

        !> unpack global array in global proc#0 to subglobal arrays
        !> subglobal arrays are slab decomposed along the 3rd dim and stored in prod#0 of each subcomm_l
        if (local_proc(subcomm_l).eq.0) then

            allocate(recvcount(0:count_proc(subcomm_r)-1))
            allocate(displacement(0:count_proc(subcomm_r)-1))

            recvcount_loc = n1_glb * n2_glb * n3_loc
            displacement_loc = n1_glb * n2_glb * local_index(n3_glb,subcomm_r)
            call mpi_allgather(recvcount_loc,1,MPI_INTEGER,recvcount,1,MPI_INTEGER,subcomm_r,mpi_ierr)
            call mpi_allgather(displacement_loc,1,MPI_INTEGER,displacement,1,MPI_INTEGER,subcomm_r,mpi_ierr)

            if (axis.eq.1) then
                allocate(subglobal_array(n1_glb,n2_glb,n3_loc))
                call mpi_scatterv(array_glb, recvcount, displacement, scalar_element_type, &
                    subglobal_array, n1_glb*n2_glb*n3_loc, scalar_element_type, 0, subcomm_r, mpi_ierr)
            elseif (axis.eq.2) then
                allocate(global_array_temp(n2_glb,n1_glb,n3_glb))
                global_array_temp = reshape(array_glb, shape(global_array_temp), order = [2,1,3])
                call mpi_scatterv(global_array_temp, recvcount, displacement, scalar_element_type, &
                    subglobal_array, n1_glb*n2_glb*n3_loc, scalar_element_type, 0, subcomm_r, mpi_ierr)
                deallocate(global_array_temp)
            endif

            deallocate(recvcount,displacement)

        endif

        if (axis.eq.1) then

            ! CUSTOM VECTOR TYPE:
            ! N2_loc x Vector {1:N1_glb,n2,1:N3_loc} per PROC
            call mpi_type_vector(n3_loc, n1_glb, n2_glb*n1_glb, scalar_element_type, subarray_type_temp, mpi_ierr)
            extend_size = element_size * n1_glb
            recvcount_loc = n2_loc
            displacement_loc = local_index(n2_glb,subcomm_l)

            ! RECEIVED DATA ORDER:
            ! {1:N1_glb,1,1:N3_loc}, ... , {1:N1_glb,N2_loc,1:N3_loc}
            allocate(local_array_temp(n1_glb,n3_loc,n2_loc))

        else

            ! CUSTOM VECTOR TYPE:
            ! N3_loc x Vector {1:N2_glb,n1,1:N3_loc} per PROC
            call mpi_type_vector(n3_loc, n2_glb, n2_glb*n1_glb, scalar_element_type, subarray_type_temp, mpi_ierr)
            extend_size = element_size * N2_glb
            recvcount_loc = n1_loc
            displacement_loc = local_index(n1_glb,subcomm_l)

            ! RECEIVED DATA ORDER:
            ! {1:N2_glb,1,1:N3_loc}, ... , {1:N2_glb,N1_loc,1:N3_loc}
            allocate(local_array_temp(n2_glb,n3_loc,n1_loc))

        ENDIF

        !> create elemental subarray of local array
        call mpi_type_create_resized(subarray_type_temp, 0, extend_size, subarray_type, mpi_ierr)
        call mpi_type_commit(subarray_type, mpi_ierr)

        !> unpack subglobal array that is slab decomposed along the 3rd dim and stored in prod#0 of each subcomm_l
        !> to local arrays (pencils) stored in individual procs in that subcomm_l
        allocate(recvcount(0:count_proc(subcomm_l)-1))
        allocate(displacement(0:count_proc(subcomm_l)-1))
        call mpi_allgather(recvcount_loc,1,MPI_INTEGER,recvcount,1,MPI_INTEGER,subcomm_l,mpi_ierr)
        CALL mpi_allgather(displacement_loc,1,MPI_INTEGER,displacement,1,MPI_INTEGER,subcomm_l,mpi_ierr)
        CALL mpi_scatterv(subglobal_array, recvcount, displacement, subarray_type, &
                        local_array_temp, n1_loc*n2_loc*n3_loc, scalar_element_type, 0, subcomm_l, mpi_ierr)

        !> copy value to this%e
        if (.not.associated(this%e)) call this%alloc()
        if (axis.eq.1) then
            !      #1   #2   #3
            ! old: N1 X N3 X N2
            ! new: N1 X N2 X N3
            this%e = reshape(local_array_temp, shape(this%e), order = [1,3,2])

        elseif (axis.EQ.2) THEN
            !      #1   #2   #3
            ! old: N2 X N3 X N1
            ! new: N1 X N2 X N3
            this%e = reshape(local_array_temp, shape(this%e), order = [3,1,2])
        endif

        call mpi_type_free(subarray_type,mpi_ierr)
        call mpi_type_free(subarray_type_temp,mpi_ierr)
        if (allocated(subglobal_array)) deallocate(subglobal_array)
        deallocate(local_array_temp,recvcount,displacement)

    end procedure

    !> create subcommunicator groups
    module procedure subcomm_cart_2d
        integer :: i
        integer, allocatable :: subcomms(:)

        call subcomm_cart(comm, 2, subcomms)

        do i = 1,2
            comm_grps(i) = subcomms(i)
            rank_grps(i) = local_proc(subcomms(i))
            nprocs_grps(i) = count_proc(subcomms(i))
        enddo

        deallocate(subcomms)
        
    end procedure

! ======================================================================================================== !
! VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV INTERNAL (PRIVATE) SUBROUTINES/FUNCTIONS VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV !
! ======================================================================================================== !
    subroutine DECOMPOSE(NSIZE,NPROCS,PROC_NUM,NSIZE_PROC,INDEX_PROC)
! ======================================================================
! [USAGE]:
! DECOMPOSE 1D DOMAIN INTO ALL PROCESSORS
! [PARAMETERS]:
! NSIZE >> TOT NUM OF ELEMENTS ALONG THAT DIMENSION 
! NPROCS >> TOT NUM OF PROCS ALONG THAT DIMENSION
! PROC_NUM >> RANK OF THAT PROCESSOR
! NSIZE_PROC >> LOC NUM OF ELEMENTS ALONG THAT DIMENSION
! INDEX_PROC >> START INDEX OF THAT PROCESSOR ALONG THAT DIMENSION
! [NOTE]:
! INDEX STARTS FROM 0
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: NSIZE, NPROCS, PROC_NUM, NSIZE_PROC, INDEX_PROC
    INTEGER:: Q,R

    Q = NSIZE/NPROCS
    R = MOD(NSIZE,NPROCS)
    IF (R.GT.PROC_NUM) THEN
        NSIZE_PROC = Q+1
        INDEX_PROC = NSIZE_PROC*PROC_NUM
    ELSE
        NSIZE_PROC = Q
        INDEX_PROC = NSIZE_PROC*PROC_NUM + R
    ENDIF

    ! NSIZE_PROC = Q + (R.GT.PROC_NUM) ! DOESN'T WORK IN FORTRAN
    ! INDEX_PROC = Q * PROC_NUM + MIN(R,PROC_NUM)

    RETURN
    end subroutine DECOMPOSE
! ======================================================================
    subroutine SUBARRAY(ELEMENT_DATA_TYPE,NDIM,DATASIZE_PROC,DIM,NPROCS,NEW_DATA_TYPE)
! ======================================================================
! [USAGE]: CREATE SUBARRAY DATATYPES FOR EACH PROC. ORIGINAL ARRAY IS DE
! COMPOSED ALONG DIM INTO SUBARRAYS.
! [PARAMETERS]: 
! OLD_DATA_TYPE >> ORIGINAL MPI ELEMENTARY DATATYPE OF THE ARRAY
! NDIM >> NUMBER OF DIMENSIONS OF THE LOCAL ARRAY
! DATASIZE_PROC >> ORIGINAL DATASIZE OF THE LOCAL ARRAY
! DIM >> DIMENSION ALONG WHICH THE LOCAL ARRAY IS FURTHER DECOMPOSED
! NPROCS >> NUMBER OF PROCS ALONG THAT DIM
! NEW_DATA_TYPE >> MPI DERIVED DATATYPE FOR ALL PROCS (DIM = NPROCS)
! [NOTE]:
! 1. PROC INDEX STARTS FROM 0
! 2. DIMN INDEX STARTS FROM 1
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: ELEMENT_DATA_TYPE, NDIM, DIM, NPROCS
    INTEGER,DIMENSION(NDIM):: DATASIZE_PROC, SUBDSIZE_PROC, SUBSTARTS_PROC
    INTEGER,DIMENSION(0:NPROCS-1):: NEW_DATA_TYPE

    INTEGER:: I 

    DO I = 1,NDIM
        SUBDSIZE_PROC(I) = DATASIZE_PROC(I)
        SUBSTARTS_PROC(I) = 0
    ENDDO

    ! ALONG THAT DIM
    DO I = 0,NPROCS-1
        ! DECOMPOSE A DOMAIN OF DATASIZE_PROC(DIM) INTO NPROCS
        ! -> DECOMPOSED DOMAIN SIZE ALONG THAT DIM: NSIZE_PROC(DIM)
        ! -> OTHER DIMS KEEP THE ORIGINAL DOMAIN SIZE
        ! -> DECOMPOSED DOMAIN START INDEX ALONG THAT DIM: SUBSTARTS_PROC(DIM)
        ! -> OTHER DIMS KEEP THE ORIGINAL START INDEX
        CALL DECOMPOSE(DATASIZE_PROC(DIM),NPROCS,I,SUBDSIZE_PROC(DIM),SUBSTARTS_PROC(DIM))
        ! CREATE A NDIM-DIMENSION SUBARRAY DATATYPE FOR EACH PROC
        CALL MPI_TYPE_CREATE_SUBARRAY(NDIM,DATASIZE_PROC,SUBDSIZE_PROC,SUBSTARTS_PROC, & 
            MPI_ORDER_FORTRAN,ELEMENT_DATA_TYPE,NEW_DATA_TYPE(I),mpi_ierr)
        CALL MPI_TYPE_COMMIT(NEW_DATA_TYPE(I),mpi_ierr)
    ENDDO

    RETURN
    END subroutine SUBARRAY
! ======================================================================
    function create_new_type3D(COMM,DATA_SIZE_PROC_OLD,DIM_OLD,DATA_TYPE)
! ======================================================================
! [USAGE]:
! CREATE SUBARRAY DATATYPE FOR 3D ARRAY OF DATA_SIZE_PROC_OLD STORED IN
! COMMUNICATOR GROUP: COMM
! [NOTE]:
! COMPANION FUNC FOR EXCHANGE_3DCOMPLEX_FAST
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: NDIM = 3
    INTEGER:: COMM, DIM_OLD, DATA_TYPE
    INTEGER,DIMENSION(3),INTENT(IN):: DATA_SIZE_PROC_OLD
    INTEGER,DIMENSION(:),ALLOCATABLE:: create_new_type3D

    INTEGER:: I , NPROC_COMM

    ! DETERMINE THE NUMBER OF PROCS IN PROCESSOR GROUP: COMM
    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,mpi_ierr)
    ALLOCATE(create_new_type3D(0:NPROC_COMM-1))

    ! CREATE SUBARRAY DATATYPES
    ! OLD LOCAL ARRAY: ARRAY_OLD OF SIZE_OLD IS DECOMPOSED INTO SUBARRAY ALONG DIM_OLD
    CALL SUBARRAY(DATA_TYPE,NDIM,DATA_SIZE_PROC_OLD,DIM_OLD,NPROC_COMM,create_new_type3D)

    end function create_new_type3D
! ======================================================================
    subroutine EXCHANGE_3DCOMPLEX_FAST(COMM,ARRAY_PROC_OLD,DATA_TYPE_OLD &
                                           ,ARRAY_PROC_NEW,DATA_TYPE_NEW)
! ======================================================================
! [USAGE]:
! FOR PROCESSOR GROUP 'COMM', SWAP axis_old AND axis_new.
! [PARAMETERS]:
! COMM >> PROCESSOR GROUP
! DATA_SIZE_PROC >> SIZE OF THE LOCAL ARRAY
! ARRAY_PROC >> LOCAL ARRAY TO BE SWAPPED
! DIM >> DIMENSIONS TO BE SWAPPED
! [NOTE]
! 1. USE THIS ONLY IF THE DATATYPE ONLY NEEDS TO BE USED ONCE
! 2. PROC INDEX STARTS FROM 0
! 3. DIMN INDEX STARTS FROM 1
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: NDIM
    INTEGER:: COMM
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(IN):: ARRAY_PROC_OLD
    COMPLEX(P8),DIMENSION(:,:,:),INTENT(INOUT):: ARRAY_PROC_NEW
    INTEGER,DIMENSION(:),INTENT(IN):: DATA_TYPE_OLD, DATA_TYPE_NEW

    INTEGER:: I , NPROC_COMM
    INTEGER,DIMENSION(:),ALLOCATABLE:: counts, displs !, DATA_TYPE_OLD_COPY

    ! CHECK
    NDIM = 3
    IF (RANK(ARRAY_PROC_OLD).NE.NDIM) THEN
    ! ARRAY_PROC_OLD must have NDIM ranks
    ! i.e NDIM = 3, ARRAY_PROC_OLD(:,:,:)
        WRITE(*,*) 'ERROR: EXCHANGE_3DCOMPLEX_FAST - GLOBAL PROC#'//ITOA(RANK_GLB)//' - INCOMPATIBLE INPUT ARRAY DIMS'
    ENDIF

    ! SWAP SUBARRAYS
    ! EACH SUBARRAY IS COUNTED AS ONE UNIT
    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,mpi_ierr)
    ALLOCATE(counts(0:NPROC_COMM-1)); counts = 1
    ALLOCATE(displs(0:NPROC_COMM-1)); displs = 0
    ! DO I = 0,NPROC_COMM-1
    ! counts(I) = 1; ! SWAP ONE SUBARRAY A TIME
    ! displs(I) = 0; ! DIRECTLY START FROM THE FIRST MEMORY LOC OF THE ARRAY
    ! enddo

    CALL MPI_ALLTOALLW(ARRAY_PROC_OLD, counts, displs, DATA_TYPE_OLD, &
    ARRAY_PROC_NEW, counts, displs, DATA_TYPE_NEW, COMM, mpi_ierr)

    DEALLOCATE(counts, displs)

    end subroutine EXCHANGE_3DCOMPLEX_FAST
! ======================================================================
!                           UTILITY FUNCTIONS                           
! ======================================================================
    function local_size(NSIZE,COMM)
! ======================================================================
! [USAGE]:
! CALCULATE THE LOCAL SIZE
! [NOTE]:
! 1. ASSUMES THE PROCESSOR GROUP COMM IS ALIGNED 1D.
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM, NSIZE
    INTEGER:: local_size

    INTEGER:: NPROC_COMM, RANK_COMM, INDEX_COMM 

    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,mpi_ierr)
    CALL MPI_COMM_RANK(COMM,RANK_COMM,mpi_ierr)
    CALL DECOMPOSE(NSIZE,NPROC_COMM,RANK_COMM,local_size,INDEX_COMM)

    end function local_size
! ======================================================================
function local_index(NSIZE,COMM)
! ======================================================================
! [USAGE]:
! CALCULATE THE LOCAL INDEX
! [NOTE]:
! 1. ASSUMES THE PROCESSOR GROUP COMM IS ALIGNED 1D.
! 2. INDEX COUNTS FROM 0
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM, NSIZE
    INTEGER:: local_index

    INTEGER:: NPROC_COMM, RANK_COMM, SIZE_COMM 

    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,mpi_ierr)
    CALL MPI_COMM_RANK(COMM,RANK_COMM,mpi_ierr)
    CALL DECOMPOSE(NSIZE,NPROC_COMM,RANK_COMM,SIZE_COMM,local_index)

    end function local_index
! ======================================================================
function local_proc(COMM)
! ======================================================================
! [USAGE]:
! CALCULATE THE LOCAL MPI RANK IN THAT COMM GROUP
! [NOTE]:
! 1. ASSUMES THE PROCESSOR GROUP COMM IS ALIGNED 1D.
! 2. RANK COUNTS FROM 0
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM
    INTEGER:: local_proc

    CALL MPI_COMM_RANK(COMM, local_proc, mpi_ierr)

    end function local_proc
! ======================================================================
    function count_proc(COMM)
! ======================================================================
! [USAGE]:
! COUNT THE NUMBER OF PROCS IN THAT PROC_GROUP
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM, count_proc 

    CALL MPI_COMM_SIZE(COMM, count_proc, mpi_ierr)

    end function count_proc
! ======================================================================
    SUBROUTINE SUBCOMM_CART(COMM,NDIM,SUBCOMMS)
! ======================================================================
! [USAGE]:
! CREATE COMMUNICATORS (SUBCOMMS(I)) FOR EACH DIMENSION(I) THAT HAS CART
! ESIAN TOPOLOGY. 
! [EXAMPLE]: 
! FOR N1xN2xN3 = NPROC_COMM (NDIM = 3) PROCESSORS,
! SUBCOMMS(1) >> N2XN3 GROUPS EACH HAS N1 PROCS
! SUBCOMMS(2) >> N1XN3 GROUPS EACH HAS N2 PROCS
! [NOTES]:
! EA PROC CALCULATES ITS CORRESPONDING SUBCOMMS. PROCS ON THE SAME LINE
! IN I-TH DIM WILL SHARE THE SAME SUBCOMMS(I)
! WRITTEN BY JINGE WANG @ SEP 29 2021
! ======================================================================
    INTEGER:: COMM
    INTEGER:: NDIM
    INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT)::SUBCOMMS

    INTEGER:: COMM_CART, NPROC_COMM, I 
    INTEGER,DIMENSION(1:NDIM):: DIMS
    LOGICAL,DIMENSION(1:NDIM):: PERIODS, REMDIMS

    DIMS = 0
    PERIODS = .FALSE.
    REMDIMS = .FALSE.
    ALLOCATE(SUBCOMMS(NDIM))

    CALL MPI_COMM_SIZE(COMM,NPROC_COMM,mpi_ierr)

    ! CREATES A DIVISION OF PROCESSORS IN A CARTESIAN GRID
    ! DIMS: NUMBER OF PROCESSORS IN EACH DIM
    CALL MPI_DIMS_CREATE(NPROC_COMM,NDIM,DIMS,mpi_ierr)
    ! write(*,*) MPI_RANK,'-',dims(1),'x',dims(2)

    ! MAKES PROCESSOR GROUP THAT ATTACHES THE CARTESIAN TOPOLOGY
    ! COMM: ORIGINAL PROCESSOR GROUP
    ! NDIM: NUMBER OF DIMS IN CARTESIAN GRID
    ! DIMS: NUMBER OF PROCESSORS IN EACH DIM
    ! PERI: WHETHER EA DIM IS PERIODIC
    ! REOR: REORDER OR NOT
    ! COMM_CART: NEW PROCESSOR GROUP
    CALL MPI_CART_CREATE(COMM,NDIM,DIMS,PERIODS,.TRUE.,COMM_CART,mpi_ierr)

    ! CREATES SUB-PROCESSOR GROUP
    ! REMDIMS: TELLS WHICH DIM(S) IS KEPT IN THE SUBGRID
    ! E.X: FOR 3D, REMDIMS = .FALSE.,.TRUE.,.FALSE.
    !      THEN CREATES DIM1xDIM3 SUBGROUPS, EA WITH DIM2 PROCESSORS
    DO I = 1,NDIM
        REMDIMS(I) = .TRUE.
        CALL MPI_CART_SUB(COMM_CART,REMDIMS,SUBCOMMS(I),mpi_ierr)
        REMDIMS(I) = .FALSE.
    ENDDO

    ! FOR 2D GRID ONLY:
    ! MAKE SURE SUBCOMMS(1) IS ALWAYS BIGGER THAN SUBCOMMS(2)
    IF (NDIM.EQ.2) THEN
        IF (count_proc(SUBCOMMS(1)).LT.count_proc(SUBCOMMS(2))) THEN
            mpi_ierr = SUBCOMMS(1)
            SUBCOMMS(1) = SUBCOMMS(2)
            SUBCOMMS(2) = mpi_ierr
        ENDIF
    ENDIF

    comm_glb = COMM_CART
    CALL MPI_COMM_RANK(comm_glb,rank_glb,mpi_ierr)

    END SUBROUTINE
! ======================================================================

end submodule