!------------------------------------------------------------------------------
SUBROUTINE MiniPoisson( dt,TransientSimulation )
!------------------------------------------------------------------------------
USE Types
#ifdef _OPENMP
  USE OMP_LIB
#endif

  implicit none
!------------------------------------------------------------------------------
!  TYPE(Solver_t) :: Solver
!  TYPE(Model_t) :: Model
!  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12) 
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

  !TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, active
  INTEGER :: iter, maxiter, nColours, col, totelem, nthr, state, MaxNumNodes
  LOGICAL :: Found, VecAsm, InitHandles
  integer, allocatable :: n_active_in_col(:)
  LOGICAL :: initial_device


  !$omp target map(from: initial_device)
  initial_device = omp_is_initial_device()
  !$omp end target

  print *, 'initial_device:', initial_device



stop

!------------------------------------------------------------------------------
END SUBROUTINE MiniPoisson
!------------------------------------------------------------------------------
