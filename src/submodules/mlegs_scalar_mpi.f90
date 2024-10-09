submodule (mlegs_scalar) mlegs_scalar_mpi
  implicit none

contains
    !> initialise a scalar
    module procedure scalar_initialise
        integer(i4) :: i, comm_grp_idx, nproc_i, rank_i, loc_sz_i, loc_st_i

        this%glb_sz = glb_sz

        ! multi-dimension support
        select case (size(axis_comm))
            ! 1D: r
            case (1)
                if ((glb_sz(2).ne.1).or.(glb_sz(3).ne.1)) then
                    if (rank_glb.eq.0) write(*,*) &
                    "ERROR: scalar_initialise is expecting 1D data"
                    call mpi_abort(comm_glb,error_flag_comm,mpi_ierr)
                endif
                if ((is_warning).and.(rank_glb.eq.0)) write(*,*) &
                "WARNING: scalar_initialise detects 1D (radial) configuration"
                this%axis_comm = 0
                this%axis_comm(1) = axis_comm(1)
            ! 2D: r-theta or r-z
            case (2)
                ! find the dimension that is not present: theta or z
                i = findloc(glb_sz,1,dim=1)
                if (i.eq.1) then
                    if (rank_glb.eq.0) write(*,*) &
                    "ERROR: scalar_initialise does not support 2D (theta-z) configuration"
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
                    if (rank_glb.eq.0) write(*,*) &
                    "ERROR: scalar_initialise does not have enough communicator groups available"
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
            if ((is_warning).and.(rank_glb.eq.0)) write(*,*) &
            "WARNING: scalar_dealloc tries to deallocate memory that has not been allocated"
            return
        endif
        deallocate(this%e)
        nullify(this%e)
    end procedure

    !> re-orient distributed data
    module procedure scalar_exchange !(this,axis_old,axis_new)
        integer(i4) :: this_comm, loc_sz_new(3), i
        integer, dimension(:), allocatable :: subarray_type_old, subarray_type_new
        complex(p8), dimension(:,:,:), target, allocatable :: e_new

        !> check whether the old dimension is distributed or not
        if (this%axis_comm(axis_old).ne.0) then
            if (rank_glb.eq.0) write(*,*) &
            "ERROR: scalar_exchange requires the data to be non-distributed along the old dimension: ",itoa(axis_old)
            call mpi_abort(comm_glb,error_flag_comm,mpi_ierr)
        endif

        !> obtain communicator group
        this_comm = comm_grps(this%axis_comm(axis_new))

        !> create subarray datatype for before and after swap configuration
        subarray_type_old = create_new_type3D(this_comm,this%loc_sz,axis_old,scalar_element_type)
        subarray_type_new = create_new_type3D(this_comm,loc_sz_new,axis_new,scalar_element_type)

        !> create after-swap new local size and data
        loc_sz_new = this%loc_sz
        loc_sz_new(axis_new) = this%glb_sz(axis_new)
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
        integer(i4) :: n1_glb, n2_glb, n3_glb, n1_loc, n2_loc, n3_loc
        integer(i4) :: element_size, subcomm_l, subcomm_r
        integer(kind = MPI_ADDRESS_KIND) :: extend_size
        integer(i4) :: subarray_type_temp, subarray_type, displacement_loc, recvcount_loc
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
        call mpi_type_create_resized(subarray_type_temp, 0_mpi_address_kind, extend_size, subarray_type, mpi_ierr)
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
        integer(i4) :: n1_glb, n2_glb, n3_glb, n1_loc, n2_loc, n3_loc
        integer(i4) :: element_size, subcomm_l, subcomm_r
        integer(kind = mpi_address_kind) :: extend_size
        integer(i4) :: subarray_type_temp, subarray_type, displacement_loc, recvcount_loc
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
        call mpi_type_create_resized(subarray_type_temp, 0_mpi_address_kind, extend_size, subarray_type, mpi_ierr)
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
        integer(i4) :: i
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
! ======================================================================
    subroutine decompose(nsize,nprocs,proc_num,nsize_proc,index_proc)
! ======================================================================
! [usage]:
! decompose 1d domain into all processors
! [parameters]:
! nsize >> tot num of elements along that dimension 
! nprocs >> tot num of procs along that dimension
! proc_num >> rank of that processor
! nsize_proc >> loc num of elements along that dimension
! index_proc >> start index of that processor along that dimension
! [note]:
! index starts from 0
! written by jinge wang @ sep 29 2021
! ======================================================================
    implicit none
    integer:: nsize, nprocs, proc_num, nsize_proc, index_proc
    integer:: q,r

    q = nsize/nprocs
    r = mod(nsize,nprocs)
    if (r.gt.proc_num) then
        nsize_proc = q+1
        index_proc = nsize_proc*proc_num
    else
        nsize_proc = q
        index_proc = nsize_proc*proc_num + r
    endif

    ! nsize_proc = q + (r.gt.proc_num) ! doesn't work in fortran
    ! index_proc = q * proc_num + min(r,proc_num)

    return
    end subroutine decompose
! ======================================================================
    subroutine subarray(element_data_type,ndim,datasize_proc,dim,nprocs,new_data_type)
! ======================================================================
! [usage]: create subarray datatypes for each proc. original array is de
! composed along dim into subarrays.
! [parameters]: 
! old_data_type >> original mpi elementary datatype of the array
! ndim >> number of dimensions of the local array
! datasize_proc >> original datasize of the local array
! dim >> dimension along which the local array is further decomposed
! nprocs >> number of procs along that dim
! new_data_type >> mpi derived datatype for all procs (dim = nprocs)
! [note]:
! 1. proc index starts from 0
! 2. dimn index starts from 1
! written by jinge wang @ sep 29 2021
! ======================================================================
    implicit none
    integer:: element_data_type, ndim, dim, nprocs
    integer,dimension(ndim):: datasize_proc, subdsize_proc, substarts_proc
    integer,dimension(0:nprocs-1):: new_data_type

    integer:: i 

    do i = 1,ndim
        subdsize_proc(i) = datasize_proc(i)
        substarts_proc(i) = 0
    enddo

    ! along that dim
    do i = 0,nprocs-1
        ! decompose a domain of datasize_proc(dim) into nprocs
        ! -> decomposed domain size along that dim: nsize_proc(dim)
        ! -> other dims keep the original domain size
        ! -> decomposed domain start index along that dim: substarts_proc(dim)
        ! -> other dims keep the original start index
        call decompose(datasize_proc(dim),nprocs,i,subdsize_proc(dim),substarts_proc(dim))
        ! create a ndim-dimension subarray datatype for each proc
        call mpi_type_create_subarray(ndim,datasize_proc,subdsize_proc,substarts_proc, & 
            mpi_order_fortran,element_data_type,new_data_type(i),mpi_ierr)
        call mpi_type_commit(new_data_type(i),mpi_ierr)
    enddo

    return
    end subroutine subarray
! ======================================================================
    function create_new_type3d(comm,data_size_proc_old,dim_old,data_type)
! ======================================================================
! [usage]:
! create subarray datatype for 3d array of data_size_proc_old stored in
! communicator group: comm
! [note]:
! companion func for exchange_3dcomplex_fast
! written by jinge wang @ sep 29 2021
! ======================================================================
    implicit none
    integer:: ndim = 3
    integer:: comm, dim_old, data_type
    integer,dimension(3),intent(in):: data_size_proc_old
    integer,dimension(:),allocatable:: create_new_type3d

    integer:: i , nproc_comm

    ! determine the number of procs in processor group: comm
    call mpi_comm_size(comm,nproc_comm,mpi_ierr)
    allocate(create_new_type3d(0:nproc_comm-1))

    ! create subarray datatypes
    ! old local array: array_old of size_old is decomposed into subarray along dim_old
    call subarray(data_type,ndim,data_size_proc_old,dim_old,nproc_comm,create_new_type3d)

    end function create_new_type3d
! ======================================================================
    subroutine exchange_3dcomplex_fast(comm,array_proc_old,data_type_old &
                                           ,array_proc_new,data_type_new)
! ======================================================================
! [usage]:
! for processor group 'comm', swap axis_old and axis_new.
! [parameters]:
! comm >> processor group
! data_size_proc >> size of the local array
! array_proc >> local array to be swapped
! dim >> dimensions to be swapped
! [note]
! 1. use this only if the datatype only needs to be used once
! 2. proc index starts from 0
! 3. dimn index starts from 1
! written by jinge wang @ sep 29 2021
! ======================================================================
    implicit none
    integer :: ndim
    integer :: comm
    integer(i4) :: array_rank
    complex(p8),dimension(:,:,:),intent(in):: array_proc_old
    complex(p8),dimension(:,:,:),intent(inout):: array_proc_new
    integer(i4) ,dimension(:),intent(in):: data_type_old, data_type_new

    integer:: i , nproc_comm
    integer,dimension(:),allocatable:: counts, displs !, data_type_old_copy

    ! swap subarrays
    ! each subarray is counted as one unit
    call mpi_comm_size(comm,nproc_comm,mpi_ierr)
    allocate(counts(0:nproc_comm-1)); counts = 1
    allocate(displs(0:nproc_comm-1)); displs = 0
    ! do i = 0,nproc_comm-1
    ! counts(i) = 1; ! swap one subarray a time
    ! displs(i) = 0; ! directly start from the first memory loc of the array
    ! enddo

    call mpi_alltoallw(array_proc_old, counts, displs, data_type_old, &
    array_proc_new, counts, displs, data_type_new, comm, mpi_ierr)

    deallocate(counts, displs)

    end subroutine exchange_3dcomplex_fast
! ======================================================================
!                           utility functions                           
! ======================================================================
    function local_size(nsize,comm)
! ======================================================================
! [usage]:
! calculate the local size
! [note]:
! 1. assumes the processor group comm is aligned 1d.
! written by jinge wang @ sep 29 2021
! ======================================================================
    integer:: comm, nsize
    integer:: local_size

    integer:: nproc_comm, rank_comm, index_comm 

    call mpi_comm_size(comm,nproc_comm,mpi_ierr)
    call mpi_comm_rank(comm,rank_comm,mpi_ierr)
    call decompose(nsize,nproc_comm,rank_comm,local_size,index_comm)

    end function local_size
! ======================================================================
function local_index(nsize,comm)
! ======================================================================
! [usage]:
! calculate the local index
! [note]:
! 1. assumes the processor group comm is aligned 1d.
! 2. index counts from 0
! written by jinge wang @ sep 29 2021
! ======================================================================
    integer:: comm, nsize
    integer:: local_index

    integer:: nproc_comm, rank_comm, size_comm 

    call mpi_comm_size(comm,nproc_comm,mpi_ierr)
    call mpi_comm_rank(comm,rank_comm,mpi_ierr)
    call decompose(nsize,nproc_comm,rank_comm,size_comm,local_index)

    end function local_index
! ======================================================================
function local_proc(comm)
! ======================================================================
! [usage]:
! calculate the local mpi rank in that comm group
! [note]:
! 1. assumes the processor group comm is aligned 1d.
! 2. rank counts from 0
! written by jinge wang @ sep 29 2021
! ======================================================================
    integer:: comm
    integer:: local_proc

    call mpi_comm_rank(comm, local_proc, mpi_ierr)

    end function local_proc
! ======================================================================
    function count_proc(comm)
! ======================================================================
! [usage]:
! count the number of procs in that proc_group
! written by jinge wang @ sep 29 2021
! ======================================================================
    integer:: comm, count_proc 

    call mpi_comm_size(comm, count_proc, mpi_ierr)

    end function count_proc
! ======================================================================
    subroutine subcomm_cart(comm,ndim,subcomms)
! ======================================================================
! [usage]:
! create communicators (subcomms(i)) for each dimension(i) that has cart
! esian topology. 
! [example]: 
! for n1xn2xn3 = nproc_comm (ndim = 3) processors,
! subcomms(1) >> n2xn3 groups each has n1 procs
! subcomms(2) >> n1xn3 groups each has n2 procs
! [notes]:
! ea proc calculates its corresponding subcomms. procs on the same line
! in i-th dim will share the same subcomms(i)
! written by jinge wang @ sep 29 2021
! ======================================================================
    integer:: comm
    integer:: ndim
    integer,dimension(:),allocatable,intent(out)::subcomms

    integer:: comm_cart, nproc_comm, i 
    integer,dimension(1:ndim):: dims
    logical,dimension(1:ndim):: periods, remdims

    dims = 0
    periods = .false.
    remdims = .false.
    allocate(subcomms(ndim))

    call mpi_comm_size(comm,nproc_comm,mpi_ierr)

    ! creates a division of processors in a cartesian grid
    ! dims: number of processors in each dim
    call mpi_dims_create(nproc_comm,ndim,dims,mpi_ierr)
    ! write(*,*) mpi_rank,'-',dims(1),'x',dims(2)

    ! makes processor group that attaches the cartesian topology
    ! comm: original processor group
    ! ndim: number of dims in cartesian grid
    ! dims: number of processors in each dim
    ! peri: whether ea dim is periodic
    ! reor: reorder or not
    ! comm_cart: new processor group
    call mpi_cart_create(comm,ndim,dims,periods,.true.,comm_cart,mpi_ierr)

    ! creates sub-processor group
    ! remdims: tells which dim(s) is kept in the subgrid
    ! e.x: for 3d, remdims = .false.,.true.,.false.
    !      then creates dim1xdim3 subgroups, ea with dim2 processors
    do i = 1,ndim
        remdims(i) = .true.
        call mpi_cart_sub(comm_cart,remdims,subcomms(i),mpi_ierr)
        remdims(i) = .false.
    enddo

    ! for 2d grid only:
    ! make sure subcomms(1) is always bigger than subcomms(2)
    if (ndim.eq.2) then
        if (count_proc(subcomms(1)).lt.count_proc(subcomms(2))) then
            mpi_ierr = subcomms(1)
            subcomms(1) = subcomms(2)
            subcomms(2) = mpi_ierr
        endif
    endif

    comm_glb = comm_cart
    call mpi_comm_rank(comm_glb,rank_glb,mpi_ierr)

    end subroutine
! ======================================================================

end submodule