submodule (mlegs_scalar) mlegs_scalar_dist
  implicit none

contains
  !> re-orient distributed data
  module procedure scalar_exchange !(this,axis_old,axis_new)
    integer(i4) :: this_comm, loc_sz_new(3), i
    integer, dimension(:), allocatable :: subarray_type_old, subarray_type_new
    complex(p8), dimension(:,:,:), pointer :: e_new

    real(p8) :: ln_
    integer(i4) :: nrchop_offset_, npchop_offset_, nzchop_offset_
    character(len=3) :: space_

    ln_ = this%ln
    nrchop_offset_ = this%nrchop_offset
    npchop_offset_ = this%npchop_offset
    nzchop_offset_ = this%nzchop_offset
    space_ = this%space

    !> check whether the old dimension is distributed or not
    if (this%axis_comm(axis_old).ne.0) then
      if (rank_glb.eq.0) write(*,*) &
      "ERROR: scalar_exchange requires the data to be non-distributed along the old dimension: ",ntoa(axis_old)
      call MPI_abort(comm_glb,error_flag_comm,MPI_err)
    endif

    !> obtain communicator group
    this_comm = comm_grps(this%axis_comm(axis_new))

    !> create after-swap new local size and data
    loc_sz_new = this%loc_sz
    loc_sz_new(axis_new) = this%glb_sz(axis_new)
    loc_sz_new(axis_old) = local_size(this%glb_sz(axis_old),this_comm)
    allocate(e_new(loc_sz_new(1),loc_sz_new(2),loc_sz_new(3)))

    !> create subarray datatype for before and after swap configuration
    subarray_type_old = create_new_type3D(this_comm,this%loc_sz,axis_old,scalar_element_type)
    subarray_type_new = create_new_type3D(this_comm,loc_sz_new,axis_new,scalar_element_type)

    !> perform swap
    call exchange_3dcomplex_fast(this_comm,this%e,subarray_type_old,e_new,subarray_type_new)

    !> reassign pointer, loc_sz, loc_st, and axis_comm
    call this%dealloc()
    this%e => e_new
    nullify(e_new)
    this%loc_sz = loc_sz_new
    this%loc_st(axis_old) = local_index(this%glb_sz(axis_old),this_comm)
    this%loc_st(axis_new) = 0
    this%axis_comm(axis_old) = this%axis_comm(axis_new)
    this%axis_comm(axis_new) = 0

    !> clean-up
    do i = 1,size(subarray_type_old)
      call MPI_type_free(subarray_type_old(i),MPI_err)
      call MPI_type_free(subarray_type_new(i),MPI_err)
    enddo
    deallocate(subarray_type_new,subarray_type_old)

    this%ln = ln_
    this%nrchop_offset = nrchop_offset_
    this%npchop_offset = npchop_offset_
    this%nzchop_offset = nzchop_offset_
    this%space = space_

  end procedure

  !> assemble distributed data into a single proc
  module procedure scalar_assemble
    integer(i4) :: axis, n1_glb, n2_glb, n3_glb, n1_loc, n2_loc, n3_loc
    integer(i4) :: element_size, subcomm_l, subcomm_r
    integer(kind = MPI_address_kind) :: extend_size
    integer(i4) :: subarray_type_temp, subarray_type, displacement_loc, recvcount_loc
    integer, dimension(:), allocatable :: displacement, recvcount
    complex(p8), dimension(:,:,:), allocatable :: local_array_temp, subglobal_array, global_array_temp
    type(scalar) :: this_copy

    if (cp8_size.eq.0) then
      call MPI_type_size(scalar_element_type,cp8_size,MPI_err)
    endif
    element_size = cp8_size
    
    n1_glb = this%glb_sz(1); n2_glb = this%glb_sz(2); n3_glb = this%glb_sz(3)
    n1_loc = this%loc_sz(1); n2_loc = this%loc_sz(2); n3_loc = this%loc_sz(3)

    allocate(array_glb(n1_glb,n2_glb,n3_glb))

    if (present(axis_input)) then
      axis = axis_input
    else
      axis = findloc(this%axis_comm,0,dim=1)
    end if            

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
      ! order: #1 of old (N1) becomes #1 of new, #2 of old (N2) becomes #3 of new, 
      !        #3 of old (N3) becomes #2 of newa
      local_array_temp = reshape(this%e, shape(local_array_temp), order = [1,3,2])

      !> create custom MPI datatype
      call MPI_type_vector(n3_loc, n1_glb, n2_glb*n1_glb, scalar_element_type, subarray_type_temp, MPI_err)
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
      ! order: #1 of old (N1) becomes #3 of new, #2 of old (N2) becomes #1 of new, 
      !        #3 of old (N3) becomes #2 of new
      local_array_temp = reshape(this%e, shape(local_array_temp), order = [3,1,2])

      call MPI_type_vector(n3_loc, n2_glb, n2_glb*n1_glb, scalar_element_type, subarray_type_temp, MPI_err)
      extend_size = element_size * n2_glb
      recvcount_loc = n1_loc
      displacement_loc = local_index(n1_glb, subcomm_l)

    elseif (axis.eq.3) then ! third axis is complete - use exchange then assemble

      if ((is_warning).and.(rank_glb.eq.0)) write(*,*) &
      "WARNING: scalar_assemble does not have native support for the 3rd axis, but it will proceed by exchanging axis first."
      ! prevent modifying the original data
      this_copy = this 
      ! exchange axis 3 to axis 1
      call this_copy%exchange(3,1) 
      ! assemble with axis 1 now complete
      array_glb(:,:,:) = this_copy%assemble(1) ! note: must have (:,:,:) after array_glb
      call this_copy%dealloc()
      return

    else
      if (rank_glb.eq.0) write(*,*) "ERROR: scalar_assemble does not support the current axis"
      call MPI_abort(comm_glb,error_flag_comm,MPI_err)
    endif

    !> create elemental subarray of local array
    call MPI_type_create_resized(subarray_type_temp, 0_MPI_address_kind, extend_size, subarray_type, MPI_err)
    call MPI_type_commit(subarray_type, MPI_err)

    !> in ea sub_comm_l, gather local array into subglobal array in subcomm_l's proc#0
    allocate(recvcount(0:count_proc(subcomm_l)-1))
    allocate(displacement(0:count_proc(subcomm_l)-1))
    call MPI_allgather(recvcount_loc,1,MPI_integer,recvcount,1,MPI_integer,subcomm_l,MPI_err)
    call MPI_allgather(displacement_loc,1,MPI_integer,displacement,1,MPI_integer,subcomm_l,MPI_err)
    call MPI_gatherv(local_array_temp, n1_loc*n2_loc*n3_loc, scalar_element_type, &
      subglobal_array, recvcount, displacement, subarray_type, 0, subcomm_l, MPI_err)

    !> in the subcomm_r that contains subcomm_l's proc#0:
    !> gather subglobal array into global array that is stored in subcomm_r's proc#0 (=> global proc#0)
    if (local_proc(subcomm_l).eq.0) then

      deallocate(recvcount,displacement)
      allocate(recvcount(0:count_proc(subcomm_r)-1))
      allocate(displacement(0:count_proc(subcomm_r)-1))
      recvcount_loc = n1_glb*n2_glb*n3_loc
      displacement_loc = n1_glb*n2_glb*local_index(n3_glb,subcomm_r)

      call MPI_allgather(recvcount_loc,1,MPI_integer,recvcount,1,MPI_integer,subcomm_r,MPI_err)
      call MPI_allgather(displacement_loc,1,MPI_integer,displacement,1,MPI_integer,subcomm_r,MPI_err)

      if (axis.eq.1) then
        call MPI_gatherv(subglobal_array, n1_glb*n2_glb*n3_loc, scalar_element_type, &
          array_glb, recvcount, displacement, scalar_element_type, 0, subcomm_r, MPI_err)
      elseif (axis.eq.2) then
        allocate(global_array_temp(n2_glb,n1_glb,n3_glb))
        call MPI_gatherv(subglobal_array, n1_glb*n2_glb*n3_loc, scalar_element_type, &
          global_array_temp, recvcount, displacement, scalar_element_type, 0, subcomm_r, MPI_err)
        !      #1   #2   #3
        ! old: N2 x N1 x N3
        ! new: N1 x N2 x N3 
        ! order: #1 of old (N2) becomes #2 of new, #2 of old (N1) becomes #1 of new, 
        !        #3 of old (N3) becomes #3 of new
        array_glb = reshape(global_array_temp, shape(array_glb), order = [2,1,3])
        deallocate(global_array_temp)
      endif

    endif
       
    !> clean-up
    call MPI_type_free(subarray_type, MPI_err)
    call MPI_type_free(subarray_type_temp, MPI_err)
    deallocate(local_array_temp, subglobal_array, recvcount, displacement)

  end procedure

  !> disassemble data into distributed procs
  module procedure scalar_disassemble
    integer(i4) :: axis, n1_glb, n2_glb, n3_glb, n1_loc, n2_loc, n3_loc, axis_comm_copy(3)
    integer(i4) :: element_size, subcomm_l, subcomm_r
    integer(kind = MPI_address_kind) :: extend_size
    integer(i4) :: subarray_type_temp, subarray_type, displacement_loc, recvcount_loc
    integer, dimension(:), allocatable :: displacement, recvcount
    complex(p8), dimension(:,:,:), allocatable :: local_array_temp, subglobal_array, global_array_temp

    if (cp8_size.eq.0) then
      call MPI_type_size(scalar_element_type,cp8_size,MPI_err)
    endif
    element_size = cp8_size
    
    n1_glb = this%glb_sz(1); n2_glb = this%glb_sz(2); n3_glb = this%glb_sz(3)
    n1_loc = this%loc_sz(1); n2_loc = this%loc_sz(2); n3_loc = this%loc_sz(3)

    !> check global array dimensions
    if (rank_glb.eq.0) then
      if ((size(array_glb,1).ne.n1_glb).or.(size(array_glb,2).ne.n2_glb).or.(size(array_glb,3).ne.n3_glb)) then
        write(*,*) "ERROR: scalar_disassemble has inconsistent dimension between array_glb and this%glb_sz"
        call MPI_abort(comm_glb,error_flag_comm,MPI_err)
      endif
    endif

    !> define left and right subcommunicator groups
    if (present(axis_input)) then
      axis = axis_input
    else
      axis = findloc(this%axis_comm,0,dim=1)
    end if
    if (axis.eq.1) then
      ! axis,    SUBCOMMS_L  ,   SUBCOMMS_R
      !  N1 , N2/axis_comm(2), N3/axis_comm(3)
      subcomm_l = comm_grps(this%axis_comm(2))
      subcomm_r = comm_grps(this%axis_comm(3))
      allocate(subglobal_array(n1_glb,n2_glb,n3_loc))
    elseif (axis.eq.2) then
      !   SUBCOMMS_L   , axis,   SUBCOMMS_R
      ! N1/axis_comm(1),  N2 , N3/axis_comm(3)
      subcomm_l = comm_grps(this%axis_comm(1))
      subcomm_r = comm_grps(this%axis_comm(3))
      allocate(subglobal_array(n2_glb, n1_glb, n3_loc))
    else
      ! if (rank_glb.eq.0) write(*,*) "ERROR: scalar_disassemble does not support the current axis"
      ! call MPI_abort(comm_glb,error_flag_comm,MPI_err)

      if ((is_warning).and.(rank_glb.eq.0)) write(*,*) &
      "WARNING: scalar_disassemble does not have native support for the 3rd axis, &", &
      "but it will proceed by exchanging axis first."

      ! Copy the current axis communicator configuration
      axis_comm_copy = this%axis_comm
      ! Deallocate the current scalar data
      call this%dealloc()
      ! Reinitialize the scalar with the new swapped axis communicator configuration
      call this%init(this%glb_sz,(/ axis_comm_copy(1), axis_comm_copy(3), axis_comm_copy(2) /))
      ! Disassemble the scalar data with the second axis complete
      call this%disassemble(array_glb,2)
      ! Swap the axis back to the original configuration
      call this%exchange(2, 3)

    endif

    !> unpack global array in global proc#0 to subglobal arrays
    !> subglobal arrays are slab decomposed along the 3rd dim and stored in prod#0 of each subcomm_l
    if (local_proc(subcomm_l).eq.0) then

      allocate(recvcount(0:count_proc(subcomm_r)-1))
      allocate(displacement(0:count_proc(subcomm_r)-1))

      recvcount_loc = n1_glb * n2_glb * n3_loc
      displacement_loc = n1_glb * n2_glb * local_index(n3_glb,subcomm_r)
      call MPI_allgather(recvcount_loc,1,MPI_integer,recvcount,1,MPI_integer,subcomm_r,MPI_err)
      call MPI_allgather(displacement_loc,1,MPI_integer,displacement,1,MPI_integer,subcomm_r,MPI_err)

      if (axis.eq.1) then
        call MPI_scatterv(array_glb, recvcount, displacement, scalar_element_type, &
          subglobal_array, n1_glb*n2_glb*n3_loc, scalar_element_type, 0, subcomm_r, MPI_err)
      elseif (axis.eq.2) then
        allocate(global_array_temp(n2_glb,n1_glb,n3_glb))
        !      #1   #2   #3
        ! old: N1 x N2 x N3
        ! new: N2 x N1 x N3 
        ! order: #1 of old (N1) becomes #2 of new, #2 of old (N2) becomes #1 of new, 
        !        #3 of old (N3) becomes #3 of new
        global_array_temp = reshape(array_glb, shape(global_array_temp), order = [2,1,3])
        call MPI_scatterv(global_array_temp, recvcount, displacement, scalar_element_type, &
          subglobal_array, n1_glb*n2_glb*n3_loc, scalar_element_type, 0, subcomm_r, MPI_err)
        deallocate(global_array_temp)
      endif

      deallocate(recvcount,displacement)

    endif

    if (axis.eq.1) then

      ! CUSTOM VECTOR TYPE:
      ! N2_loc x Vector {1:N1_glb,n2,1:N3_loc} per PROC
      call MPI_type_vector(n3_loc, n1_glb, n2_glb*n1_glb, scalar_element_type, subarray_type_temp, MPI_err)
      extend_size = element_size * n1_glb
      recvcount_loc = n2_loc
      displacement_loc = local_index(n2_glb,subcomm_l)

      ! RECEIVED DATA ORDER:
      ! {1:N1_glb,1,1:N3_loc}, ... , {1:N1_glb,N2_loc,1:N3_loc}
      allocate(local_array_temp(n1_glb,n3_loc,n2_loc))

    else

      ! CUSTOM VECTOR TYPE:
      ! N3_loc x Vector {1:N2_glb,n1,1:N3_loc} per PROC
      call MPI_type_vector(n3_loc, n2_glb, n2_glb*n1_glb, scalar_element_type, subarray_type_temp, MPI_err)
      extend_size = element_size * N2_glb
      recvcount_loc = n1_loc
      displacement_loc = local_index(n1_glb,subcomm_l)

      ! RECEIVED DATA ORDER:
      ! {1:N2_glb,1,1:N3_loc}, ... , {1:N2_glb,N1_loc,1:N3_loc}
      allocate(local_array_temp(n2_glb,n3_loc,n1_loc))

    ENDIF

    !> create elemental subarray of local array
    call MPI_type_create_resized(subarray_type_temp, 0_MPI_address_kind, extend_size, subarray_type, MPI_err)
    call MPI_type_commit(subarray_type, MPI_err)

    !> unpack subglobal array that is slab decomposed along the 3rd dim and stored in prod#0 of each subcomm_l
    !> to local arrays (pencils) stored in individual procs in that subcomm_l
    allocate(recvcount(0:count_proc(subcomm_l)-1))
    allocate(displacement(0:count_proc(subcomm_l)-1))
    call MPI_allgather(recvcount_loc,1,MPI_integer,recvcount,1,MPI_integer,subcomm_l,MPI_err)
    CALL MPI_allgather(displacement_loc,1,MPI_integer,displacement,1,MPI_integer,subcomm_l,MPI_err)
    CALL MPI_scatterv(subglobal_array, recvcount, displacement, subarray_type, &
            local_array_temp, n1_loc*n2_loc*n3_loc, scalar_element_type, 0, subcomm_l, MPI_err)

    !> copy value to this%e
    if (.not.associated(this%e)) then
      allocate(this%e(this%loc_sz(1),this%loc_sz(2),this%loc_sz(3)))
    endif
    if (axis.eq.1) then
      !      #1   #2   #3
      ! old: N1 X N3 X N2
      ! new: N1 X N2 X N3 
      ! order: #1 of old (N1) becomes #1 of new, #2 of old (N3) becomes #3 of new, 
      !        #3 of old (N2) becomes #2 of new
      this%e = reshape(local_array_temp, shape(this%e), order = [1,3,2])

    elseif (axis.EQ.2) THEN
      !      #1   #2   #3
      ! old: N2 X N3 X N1
      ! new: N1 X N2 X N3 
      ! order: #1 of old (N2) becomes #2 of new, #2 of old (N3) becomes #3 of new, 
      !        #3 of old (N1) becomes #1 of new
      this%e = reshape(local_array_temp, shape(this%e), order = [2,3,1])
    endif

    call MPI_type_free(subarray_type,MPI_err)
    call MPI_type_free(subarray_type_temp,MPI_err)
    if (allocated(subglobal_array)) deallocate(subglobal_array)
    deallocate(local_array_temp,recvcount,displacement)

  end procedure

  !> create subcommunicator groups
  module procedure subcomm_cart_2d
    integer(i4) :: i
    integer, allocatable :: subcomms(:)

    if (present(dims)) then
      call subcomm_cart(comm, 2, subcomms, dims)
    else
      call subcomm_cart(comm, 2, subcomms)
    endif

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

  subroutine subarray(element_data_type,ndim,datasize_proc,dim,nprocs,new_data_type)
  !> [usage]: create subarray datatypes for each proc. original array is de
  !> composed along dim into subarrays.
  !> [parameters]: 
  !> old_data_type >> original MPI elementary datatype of the array
  !> ndim >> number of dimensions of the local array
  !> datasize_proc >> original datasize of the local array
  !> dim >> dimension along which the local array is further decomposed
  !> nprocs >> number of procs along that dim
  !> new_data_type >> MPI derived datatype for all procs (dim = nprocs)
  !> [note]:
  !> 1. proc index starts from 0
  !> 2. dimn index starts from 1
  !> written by Jinge Wang @ sep 29 2021
    implicit none
    integer(i4) :: element_data_type, ndim, dim, nprocs
    integer(i4), dimension(ndim) :: datasize_proc, subdsize_proc, substarts_proc
    integer(i4), dimension(0:nprocs-1):: new_data_type

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
      call MPI_type_create_subarray(ndim,datasize_proc,subdsize_proc,substarts_proc, & 
        MPI_order_fortran,element_data_type,new_data_type(i),MPI_err)
      call MPI_type_commit(new_data_type(i),MPI_err)
    enddo

    return
  end subroutine

! ======================================================================================================== !

  function create_new_type3d(comm,data_size_proc_old,dim_old,data_type)
  !> [usage]:
  !> create subarray datatype for 3d array of data_size_proc_old stored in
  !> communicator group: comm
  !> [note]:
  !> companion func for exchange_3dcomplex_fast
  !> written by Jinge Wang @ sep 29 2021
    implicit none
    integer(i4) :: ndim = 3
    integer(i4) :: comm, dim_old, data_type
    integer(i4), dimension(3), intent(in):: data_size_proc_old
    integer(i4), dimension(:), allocatable:: create_new_type3d

    integer:: i , nproc_comm

    ! determine the number of procs in processor group: comm
    call MPI_comm_size(comm,nproc_comm,MPI_err)
    allocate(create_new_type3d(0:nproc_comm-1))

    ! create subarray datatypes
    ! old local array: array_old of size_old is decomposed into subarray along dim_old
    call subarray(data_type,ndim,data_size_proc_old,dim_old,nproc_comm,create_new_type3d)

    return
  end function

! ======================================================================================================== !

  subroutine exchange_3dcomplex_fast(comm,array_proc_old,data_type_old,array_proc_new,data_type_new)
  !> [usage]:
  !> for processor group 'comm', swap axis_old and axis_new.
  !> [parameters]:
  !> comm >> processor group
  !> data_size_proc >> size of the local array
  !> array_proc >> local array to be swapped
  !> dim >> dimensions to be swapped
  !> [note]
  !> 1. use this only if the datatype only needs to be used once
  !> 2. proc index starts from 0
  !> 3. dimn index starts from 1
  !> written by Jinge Wang @ sep 29 2021
    implicit none
    integer(i4) :: ndim
    integer(i4) :: comm
    integer(i4) :: array_rank
    complex(p8), dimension(:,:,:), intent(in):: array_proc_old
    complex(p8), dimension(:,:,:), intent(inout):: array_proc_new
    integer(i4), dimension(:), intent(in):: data_type_old, data_type_new

    integer(i4) :: i , nproc_comm
    integer(i4),dimension(:), allocatable:: counts, displs !, data_type_old_copy

    ! swap subarrays
    ! each subarray is counted as one unit
    call MPI_comm_size(comm,nproc_comm,MPI_err)
    allocate(counts(0:nproc_comm-1)); counts = 1
    allocate(displs(0:nproc_comm-1)); displs = 0

    call MPI_alltoallw(array_proc_old, counts, displs, data_type_old, &
    array_proc_new, counts, displs, data_type_new, comm, MPI_err)

    deallocate(counts, displs)

    return
  end subroutine

! ======================================================================================================== !

  subroutine subcomm_cart(comm,ndim,subcomms,dims_in)
!> [usage]:
!> create communicators (subcomms(i)) for each dimension(i) that has cart
!> esian topology. 
!> [example]: 
!> for n1xn2xn3 = nproc_comm (ndim = 3) processors,
!> subcomms(1) >> n2xn3 groups each has n1 procs
!> subcomms(2) >> n1xn3 groups each has n2 procs
!> [notes]:
!> ea proc calculates its corresponding subcomms. procs on the same line
!> in i-th dim will share the same subcomms(i)
!> written by Jinge Wang @ sep 29 2021
    integer(i4) :: comm
    integer(i4) :: ndim
    integer(i4), dimension(:), allocatable, intent(out) ::subcomms
    integer(i4), dimension(1:ndim), optional :: dims_in

    integer(i4) :: comm_cart, nproc_comm, i 
    integer(i4), dimension(1:ndim) :: dims
    logical, dimension(1:ndim):: periods, remdims

    dims = 0
    periods = .false.
    remdims = .false.
    allocate(subcomms(ndim))

    call MPI_comm_size(comm,nproc_comm,MPI_err)

    ! creates a division of processors in a cartesian grid
    ! dims: number of processors in each dim
    if (present(dims_in)) then
      dims = dims_in
    else
      call MPI_dims_create(nproc_comm,ndim,dims,MPI_err)
    endif
    if ((is_warning).and.(rank_glb.eq.0)) write(*,*) 'comm dims = ',dims
    ! write(*,*) MPI_rank,'-',dims(1),'x',dims(2)

    ! makes processor group that attaches the cartesian topology
    ! comm: original processor group
    ! ndim: number of dims in cartesian grid
    ! dims: number of processors in each dim
    ! peri: whether ea dim is periodic
    ! reor: reorder or not
    ! comm_cart: new processor group
    call MPI_cart_create(comm,ndim,dims,periods,.true.,comm_cart,MPI_err)

    ! creates sub-processor group
    ! remdims: tells which dim(s) is kept in the subgrid
    ! e.x: for 3d, remdims = .false.,.true.,.false.
    !      then creates dim1xdim3 subgroups, ea with dim2 processors
    do i = 1,ndim
      remdims(i) = .true.
      call MPI_cart_sub(comm_cart,remdims,subcomms(i),MPI_err)
      remdims(i) = .false.
    enddo

    ! for 2d grid only:
    ! make sure subcomms(1) is always bigger than subcomms(2)
    if (ndim.eq.2) then
      if (count_proc(subcomms(1)).lt.count_proc(subcomms(2))) then
        MPI_err = subcomms(1)
        subcomms(1) = subcomms(2)
        subcomms(2) = MPI_err
      endif
    endif

    comm_glb = comm_cart
    call MPI_comm_rank(comm_glb,rank_glb,MPI_err)

  end subroutine

end submodule