submodule (mlegs_mpi) mlegs_mpi_utils
  implicit none

contains

  module procedure m_decompose
    if (rank .gt. nprocs-1) then
      write(*,*) 'm_decompose: MPI rank cannot exceed the total number of MPI processors'
      stop
    endif

    if (size .ge. nprocs) then
      chunk = size / nprocs
      start_idx = rank * chunk + 1
      end_idx = (rank + 1) * chunk
      if (rank .eq. nprocs - 1) end_idx = size
    else
      if (rank .lt. size) then
        chunk = 1
        start_idx = rank + 1
        end_idx = start_idx
      else
        chunk = 0
        start_idx = 0
        end_idx = -1
      endif
    endif
  end procedure

end submodule