submodule (mlegs_envir) mlegs_envir_mpi
  implicit none

contains

  module procedure decompose
  !> [usage]:
  !> decompose 1d domain into all processors
  !> [parameters]:
  !> nsize >> tot num of elements along that dimension 
  !> nprocs >> tot num of procs along that dimension
  !> proc_num >> rank of that processor
  !> nsize_proc >> loc num of elements along that dimension
  !> index_proc >> start index of that processor along that dimension
  !> [note]:
  !> index starts from 0
  !> written by Jinge Wang @ sep 29 2021
    integer :: q,r

    q = nsize/nprocs
    r = mod(nsize,nprocs)
    if (r.gt.proc_num) then
      nsize_proc = q+1
      index_proc = nsize_proc*proc_num
    else
      nsize_proc = q
      index_proc = nsize_proc*proc_num + r
    endif

    return
  end procedure

  module procedure local_size
  !> [usage]:
  !> calculate the local size
  !> [note]:
  !> 1. assumes the processor group comm is aligned 1d.
  !> written by Jinge Wang @ sep 29 2021
    integer:: nproc_comm, rank_comm, index_comm 

    call MPI_comm_size(comm,nproc_comm,MPI_err)
    call MPI_comm_rank(comm,rank_comm,MPI_err)
    call decompose(nsize,nproc_comm,rank_comm,l_size,index_comm)

    return
  end procedure

  module procedure local_index
  !> [usage]:
  !> calculate the local index
  !> [note]:
  !> 1. assumes the processor group comm is aligned 1d.
  !> 2. index counts from 0
  !> written by Jinge Wang @ sep 29 2021
    integer:: nproc_comm, rank_comm, size_comm 

    call MPI_comm_size(comm,nproc_comm,MPI_err)
    call MPI_comm_rank(comm,rank_comm,MPI_err)
    call decompose(nsize,nproc_comm,rank_comm,size_comm,l_index)

    return
  end procedure

  module procedure local_proc
  !> [usage]:
  !> calculate the local MPI rank in that comm group
  !> [note]:
  !> 1. assumes the processor group comm is aligned 1d.
  !> 2. rank counts from 0
  !> written by Jinge Wang @ sep 29 2021
    call MPI_comm_rank(comm, l_proc, MPI_err)

    return
  end procedure

  module procedure count_proc
  !> [usage]:
  !> count the number of procs in that proc_group
  !> written by Jinge Wang @ sep 29 2021
    call MPI_comm_size(comm, nprocs, MPI_err)

    return
  end procedure

end submodule