module mlegs_mpi
  use mlegs_base
  implicit none
  private

  !> MPI -- global variables
  integer(i4), public :: m_err, m_rank, m_nprocs

  !> print real-time clock on display
  interface m_decompose
    module subroutine m_decompose(size, nprocs, rank, chunk, start_idx, end_idx)
      implicit none
      integer(i4), intent(in) :: size, nprocs, rank
      integer(i4), intent(inout) :: chunk, start_idx, end_idx
    end subroutine
  end interface
  public :: m_decompose

end module