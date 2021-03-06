!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief New fortran functions for compiler versions that do not support them
!--------------------------------------------------------------------------------------------------
module future
  use prec

  implicit none
  public

contains

#if defined(__GFORTRAN__) &&  __GNUC__<9 || defined(__INTEL_COMPILER) && INTEL_COMPILER<1800
!--------------------------------------------------------------------------------------------------
!> @brief substitute for the findloc intrinsic (only for integer, dimension(:) at the moment)
!--------------------------------------------------------------------------------------------------
function findloc(a,v)

  integer, intent(in), dimension(:) :: a
  integer, intent(in) :: v
  integer :: i,j
  integer, allocatable, dimension(:) ::  findloc

  allocate(findloc(count(a==v)))
  j = 1
  do i = 1, size(a)
    if (a(i)==v) then
      findloc(j) = i
      j = j + 1
    endif
  enddo
end function findloc
#endif

end module future
