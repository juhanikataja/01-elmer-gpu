module mylinearforms
  USE Types, only: dp, VECTOR_BLOCK_LENGTH, VECTOR_SMALL_THRESH, Element_t, Nodes_t

  !INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12) 
  !INTEGER, PARAMETER :: VECTOR_BLOCK_LENGTH = 128                                                                                                                                                                                                                                                                                                      
  !INTEGER, PARAMETER :: VECTOR_SMALL_THRESH = 9                                                                                                                                                                                                                                                                                                        

  PUBLIC :: LinearForms_GradUdotGradU
  

   !$omp declare target to(LinearForms_GradUdotGradU)
contains

   FUNCTION ElementMetricGPU( Elm, Nodes, nc, ndof, DetJ, nbmax, dLBasisdx, LtoGMap) RESULT(AllSuccess)
!------------------------------------------------------------------------------
    !$OMP DECLARE TARGET
     TYPE(Element_t)  :: Elm                                 !< Element structure
     TYPE(Nodes_t)    :: Nodes                               !< element nodal coordinates
     INTEGER, INTENT(IN) :: nc                               !< Number of points to map
     INTEGER :: ndof                                         !< Number of active nodes in element
     REAL(KIND=dp) :: DetJ(VECTOR_BLOCK_LENGTH)              !< SQRT of determinant of element coordinate metric at each point
     INTEGER, INTENT(IN) :: nbmax                            !< Maximum total number of basis functions in local basis
     REAL(KIND=dp) :: dLBasisdx(VECTOR_BLOCK_LENGTH,nbmax,3) !< Derivatives of element basis function with 
                                                             !<  respect to local coordinates at each point
     REAL(KIND=dp) :: LtoGMap(VECTOR_BLOCK_LENGTH,3,3)       !< Mapping between local and global coordinates
     LOGICAL :: AllSuccess                  !< Returns .FALSE. if some point in element is degenerate
!------------------------------------------------------------------------------
!       Local variables
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: dx(VECTOR_BLOCK_LENGTH,3,3)
     REAL(KIND=dp) :: Metric(VECTOR_BLOCK_LENGTH,6), &
             G(VECTOR_BLOCK_LENGTH,6)       ! Symmetric Metric(nc,3,3) and G(nc,3,3)

     REAL(KIND=dp) :: s
     INTEGER :: cdim,dim,i,j,k,l,n,ip, jj, kk
     INTEGER :: ldbasis, ldxyz, utind
     !------------------------------------------------------------------------------
     AllSuccess = .TRUE.

     ! Coordinates (single array)
     n = MIN( SIZE(Nodes % x, 1), ndof )

     ! Dimensions (coordinate system and element)
     cdim = CoordinateSystemDimension()
     dim  = elm % TYPE % DIMENSION

     ! Leading dimensions for local basis and coordinate arrays
     ldbasis = SIZE(dLBasisdx, 1)
     ldxyz = SIZE(Nodes % xyz, 1)

     ! For linear, extruded and otherwise regular elements mapping has to be computed
     ! only once, the problem is to identify these cases...
     !------------------------------------------------------------------------------
     !       Partial derivatives of global coordinates with respect to local coordinates
     !------------------------------------------------------------------------------
     ! Avoid DGEMM calls for nc small
     !IF (nc < VECTOR_SMALL_THRESH) THEN
       DO l=1,dim
         DO j=1,3
           dx(1:nc,j,l)=REAL(0,dp)
           DO k=1,n
!DIR$ UNROLL
             DO i=1,nc
               dx(i,j,l)=dx(i,j,l)+dLBasisdx(i,k,l)*Nodes % xyz(k,j)
             END DO
           END DO
         END DO
       END DO
       
     !ELSE ! Here would be some call to gpu aware DGEMM
       !DO i=1,dim
         !CALL DGEMM('N','N',nc, 3, n, &
                 !REAL(1,dp), dLbasisdx(1,1,i), ldbasis, &
                 !Nodes % xyz, ldxyz, REAL(0, dp), dx(1,1,i), VECTOR_BLOCK_LENGTH)
       !END DO
     !END IF

     !------------------------------------------------------------------------------
     !       Compute the covariant metric tensor of the element coordinate system (symmetric)
     !------------------------------------------------------------------------------
     ! Linearized upper triangular indices for accesses to G
     ! | (1,1) (1,2) (1,3) | = | 1 2 4 |
     ! |       (2,2) (2,3) |   |   3 5 |
     ! |             (3,3) |   |     6 |
     ! G is symmetric, compute only the upper triangular part of G=dx^Tdx
     DO j=1,dim
       DO i=1,j
         utind = GetSymmetricIndex(i,j)
         SELECT CASE (cdim)
         CASE(1)
           !_ELMER_OMP_SIMD
           DO l=1,nc
             G(l,utind)=dx(l,1,i)*dx(l,1,j)
           END DO
         CASE(2)
           !_ELMER_OMP_SIMD
           DO l=1,nc
             G(l,utind)=dx(l,1,i)*dx(l,1,j)+dx(l,2,i)*dx(l,2,j)
           END DO
         CASE(3)
           !_ELMER_OMP_SIMD
           DO l=1,nc
             G(l,utind)=dx(l,1,i)*dx(l,1,j)+dx(l,2,i)*dx(l,2,j)+dx(l,3,i)*dx(l,3,j)
           END DO
         END SELECT
       END DO
     END DO

     !------------------------------------------------------------------------------
     !       Convert the metric to contravariant base, and compute the SQRT(DetG)
     !------------------------------------------------------------------------------
     SELECT CASE( dim )
       !------------------------------------------------------------------------------
       !       Line elements
       !------------------------------------------------------------------------------
     CASE (1)
       ! Determinants
       ! DetJ(1:nc)  = G(1:nc,1,1)
       DetJ(1:nc)  = G(1:nc,1)

       DO i=1,nc
         IF (DetJ(i) <= TINY(REAL(1,dp))) THEN
           AllSuccess = .FALSE.
           EXIT
         END IF
       END DO

       IF (AllSuccess) THEN
         !_ELMER_OMP_SIMD
         DO i=1,nc
           ! Metric(i,1,1) = REAL(1,dp)/DetJ(i)
           Metric(i,1) = REAL(1,dp)/DetJ(i)
         END DO
         !_ELMER_OMP_SIMD
         DO i=1,nc
           DetJ(i) = SQRT( DetJ(i))
         END DO
       END IF


       !------------------------------------------------------------------------------
       !       Surface elements
       !------------------------------------------------------------------------------
     CASE (2)
       ! Determinants
       !_ELMER_OMP_SIMD
       DO i=1,nc
         ! DetJ(i) = ( G(i,1,1)*G(i,2,2) - G(i,1,2)*G(i,2,1) )
         ! G is symmetric
         DetJ(i) = G(i,1)*G(i,3)-G(i,2)*G(i,2)
       END DO

       DO i=1,nc
         IF (DetJ(i) <= TINY(REAL(1,dp))) THEN
           AllSuccess = .FALSE.
           EXIT
         END IF
       END DO

       IF (AllSuccess) THEN
         ! Since G=G^T, it holds G^{-1}=(G^T)^{-1}
         !_ELMER_OMP_SIMD
         DO i=1,nc
           s = REAL(1,dp)/DetJ(i)
           ! G is symmetric
           ! All in one go, with redundancies eliminated
           Metric(i,1) =  s*G(i,3)
           Metric(i,2) = -s*G(i,2)
           Metric(i,3) =  s*G(i,1)
         END DO
         !_ELMER_OMP_SIMD
         DO i=1,nc
           DetJ(i) = SQRT(DetJ(i))
         END DO

       END IF
       !------------------------------------------------------------------------------
       !       Volume elements
       !------------------------------------------------------------------------------
     CASE (3)
       ! Determinants
       !_ELMER_OMP_SIMD
       DO i=1,nc
         ! DetJ(i) = G(i,1,1) * ( G(i,2,2)*G(i,3,3) - G(i,2,3)*G(i,3,2) ) + &
         !           G(i,1,2) * ( G(i,2,3)*G(i,3,1) - G(i,2,1)*G(i,3,3) ) + &
         !           G(i,1,3) * ( G(i,2,1)*G(i,3,2) - G(i,2,2)*G(i,3,1) )
         ! G is symmetric
         DetJ(i) = G(i,1)*(G(i,3)*G(i,6)-G(i,5)*G(i,5)) + &
                 G(i,2)*(G(i,5)*G(i,4)-G(i,2)*G(i,6)) + &
                 G(i,4)*(G(i,2)*G(i,5)-G(i,3)*G(i,4))
       END DO

       DO i=1,nc
         IF (DetJ(i) <= TINY(REAL(1,dp))) THEN
           AllSuccess = .FALSE.
           EXIT
         END IF
       END DO

       IF (AllSuccess) THEN
         ! Since G=G^T, it holds G^{-1}=(G^T)^{-1}
         !_ELMER_OMP_SIMD
         DO i=1,nc
           s = REAL(1,dp) / DetJ(i)
           ! Metric(i,1,1) =  s * (G(i,2,2)*G(i,3,3) - G(i,3,2)*G(i,2,3))
           ! Metric(i,2,1) = -s * (G(i,2,1)*G(i,3,3) - G(i,3,1)*G(i,2,3))
           ! Metric(i,3,1) =  s * (G(i,2,1)*G(i,3,2) - G(i,3,1)*G(i,2,2))
           ! G is symmetric

           ! All in one go, with redundancies eliminated
           Metric(i,1)= s*(G(i,3)*G(i,6)-G(i,5)*G(i,5))
           Metric(i,2)=-s*(G(i,2)*G(i,6)-G(i,4)*G(i,5))
           Metric(i,3)= s*(G(i,1)*G(i,6)-G(i,4)*G(i,4))
           Metric(i,4)= s*(G(i,2)*G(i,5)-G(i,3)*G(i,4))
           Metric(i,5)=-s*(G(i,1)*G(i,5)-G(i,2)*G(i,4))
           Metric(i,6)= s*(G(i,1)*G(i,3)-G(i,2)*G(i,2))
         END DO

         !_ELMER_OMP_SIMD
         DO i=1,nc
           DetJ(i) = SQRT(DetJ(i))
         END DO

       END IF
     END SELECT

     IF (AllSuccess) THEN
       SELECT CASE(dim)
       CASE(1)
!DIR$ LOOP COUNT MAX=3
         DO i=1,cdim
           !_ELMER_OMP_SIMD
           DO l=1,nc
             LtoGMap(l,i,1) = dx(l,i,1)*Metric(l,1)
           END DO
         END DO
       CASE(2)
!DIR$ LOOP COUNT MAX=3
         DO i=1,cdim
           !_ELMER_OMP_SIMD
           DO l=1,nc
             LtoGMap(l,i,1) = dx(l,i,1)*Metric(l,1) + dx(l,i,2)*Metric(l,2)
             LtoGMap(l,i,2) = dx(l,i,1)*Metric(l,2) + dx(l,i,2)*Metric(l,3)
           END DO
         END DO
       CASE(3)
!DIR$ LOOP COUNT MAX=3
         DO i=1,cdim
           !_ELMER_OMP_SIMD
           DO l=1,nc
             LtoGMap(l,i,1) = dx(l,i,1)*Metric(l,1) + dx(l,i,2)*Metric(l,2) + dx(l,i,3)*Metric(l,4)
             LtoGMap(l,i,2) = dx(l,i,1)*Metric(l,2) + dx(l,i,2)*Metric(l,3) + dx(l,i,3)*Metric(l,5)
             LtoGMap(l,i,3) = dx(l,i,1)*Metric(l,4) + dx(l,i,2)*Metric(l,5) + dx(l,i,3)*Metric(l,6)
           END DO
         END DO
       END SELECT
     ELSE

       ! Degenerate element!
       ! return some error
       ! TODO: assuming non-degenerate elements! 
     END IF

   CONTAINS

     PURE FUNCTION GetSymmetricIndex(i,j) RESULT(utind)
       IMPLICIT NONE
       !$OMP DECLARE TARGET
       INTEGER, INTENT(IN) :: i, j
       INTEGER :: utind

       IF (i>j) THEN
         utind = i*(i-1)/2+j
       ELSE
         utind = j*(j-1)/2+i
       END IF
     END FUNCTION GetSymmetricIndex
!------------------------------------------------------------------------------
   END FUNCTION ElementMetricGPU


  SUBROUTINE LinearForms_GradUdotGradU(m, n, dim, GradU, weight, G, alpha)
    implicit none
    INTEGER, INTENT(IN) :: m, n, dim
    REAL(KIND=dp) CONTIG, INTENT(IN) :: GradU(:,:,:), weight(:)
    REAL(KIND=dp) CONTIG, INTENT(INOUT) :: G(:,:)
    REAL(KIND=dp) CONTIG, INTENT(IN), OPTIONAL :: alpha(:)

    REAL(KIND=dp) :: wrk(VECTOR_BLOCK_LENGTH,n)
    INTEGER :: i, ii, iin, j, l, k, kk, ldbasis, ldwrk, ldk, blklen
    LOGICAL :: noAlphaWeight

    ldbasis = SIZE(GradU,1)
    ldwrk = SIZE(wrk,1)
    ldk = SIZE(G,1)

    noAlphaWeight = .TRUE.
    IF (PRESENT(alpha)) noAlphaWeight = .FALSE.

    DO ii=1,m,VECTOR_BLOCK_LENGTH
      iin=MIN(ii+VECTOR_BLOCK_LENGTH-1,m)
      blklen=iin-ii+1
      
      IF (blklen < VECTOR_SMALL_THRESH) THEN
        ! Do not attempt to call BLAS for small cases to avoid preprocessing overhead
        IF (noAlphaWeight) THEN
          DO j=1,n
            DO i=1,n
              DO k=1,dim
                DO l=ii,iin
                  G(i,j) = G(i,j) + GradU(l,i,k)*GradU(l,j,k)*weight(l)
                END DO
              END DO
            END DO
          END DO
        ELSE
          DO j=1,n
            DO i=1,n
              DO k=1,dim
                DO l=ii,iin
                  G(i,j) = G(i,j) + GradU(l,i,k)*GradU(l,j,k)*weight(l)*alpha(l)
                END DO
              END DO
            END DO
          END DO
        END IF
      ELSE
        DO k=1, dim
          IF (noAlphaWeight) THEN
            DO j=1,n
              DO i=ii,iin
                wrk(i-ii+1,j)=weight(i)*GradU(i,j,k)
              END DO
            END DO
          ELSE
            DO j=1,n
              DO i=ii,iin
                wrk(i-ii+1,j)=weight(i)*alpha(i)*GradU(i,j,k)
              END DO
            END DO
          END IF
        END DO
      END IF
    END DO ! Vector blocks
  END SUBROUTINE LinearForms_GradUdotGradU

end module mylinearforms

!-----------------------------------------------------------------------------
!> A prototype solver for advection-diffusion-reaction equation.
!> This equation is generic and intended for education purposes
!> but may also serve as a starting point for more complex solvers.
!> Version supporting multithreading and SIMD friendly ElmerSolver
!> kernels. 
!------------------------------------------------------------------------------
SUBROUTINE AdvDiffSolver_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'AdvDiffSolver_init'
  
  IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
    CALL Fatal(Caller,'Implemented only in cartesian coordinates')
  END IF

END SUBROUTINE AdvDiffSolver_Init


!------------------------------------------------------------------------------
SUBROUTINE AdvDiffSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
#ifdef _OPENMP
  USE OMP_LIB
#endif

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, active
  INTEGER :: iter, maxiter, nColours, col, totelem, nthr, state, MaxNumNodes
  LOGICAL :: Found, VecAsm, InitHandles
  integer, allocatable :: n_active_in_col(:)
  logical :: initial_device


  type :: elem_ptr
    type(Element_t), pointer :: p
    type(Nodes_t), pointer :: nodes
    integer :: n, nd, nb                        ! nof nodes, nof dofs, nofbdofs
  end type

  type :: elem_list_t
    type(elem_ptr), allocatable :: elements(:)
  end type

  TYPE(elem_list_t), allocatable :: elem_lists(:)

  CHARACTER(*), PARAMETER :: Caller = 'AdvDiffSolver'
!------------------------------------------------------------------------------

  CALL Info(Caller,'------------------------------------------------')
  CALL Info(Caller,'Solving generic advection-diffusion-reaction PDE')

  CALL DefaultStart()

  
  
  maxiter = ListGetInteger( GetSolverParams(),&
      'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  nthr = 1
  !$ nthr = omp_get_max_threads()

  ! Nonlinear iteration loop:
  !--------------------------
  ! DO iter=1,maxiter

  ! System assembly:
  !----------------
  CALL DefaultInitialize()

  totelem = 0

  VecAsm = (nColours > 1) .OR. (nthr == 1)

  CALL ResetTimer( Caller//'BulkAssembly' )

  !$omp target
  initial_device = omp_is_initial_device()
  !$omp end target

  print *, 'initial_device:', initial_device

  nColours = GetNOFColours(Solver)

  allocate(elem_lists(nColours))
  allocate(n_active_in_col(nColours))

  do col = 1, nColours
    active = GetNOFActive(Solver)
    allocate(elem_lists(col) % elements(active))
  end do

  !! Tabulate elements and their ndofs/nnodes/nb 
  nColours = GetNOFColours(Solver)
  !$OMP PARALLEL &
  !$OMP SHARED(Active, Solver, nColours, VecAsm, elem_lists) &
  !$OMP PRIVATE(t, Element, n, nd, nb, col, InitHandles) & 
  !$OMP REDUCTION(max:MaxNumNodes) DEFAULT(none)
  do col=1,ncolours
    !$OMP SINGLE
    active = GetNOFActive(Solver)
    !$OMP END SINGLE

    !$OMP DO
    do t=1,active
      Element => GetActiveElement(t)
      elem_lists(col) % elements(t) % p => Element
      elem_lists(col) % elements(t) % n = GetElementNOFNodes(Element)
      elem_lists(col) % elements(t) % nd = GetElementNOFDOFs(Element)
      elem_lists(col) % elements(t) % nb = GetElementNOFBDOFs(Element)
      CALL GetElementNodesVec( elem_lists(col) % elements(t) % nodes,  elem_lists(col) % elements(t) % p)
      MaxNumNodes = max(MaxNumNodes,elem_lists(col) % elements(t) % n)
    end do
    !$OMP END DO
  end do
  !$OMP END PARALLEL

  CALL CheckTimer(Caller//'BulkAssembly', Delete=.TRUE.)

  nColours = GetNOFColours(Solver)

  DO col=1,nColours

    CALL Info( Caller,'Assembly of colour: '//I2S(col),Level=1)
    Active = GetNOFActive(Solver)

    DO t=1,Active
      totelem = totelem + 1
      Element => elem_lists(col) % elements(t) % p
      n = elem_lists(col) % elements(t) % n
      nd = elem_lists(col) % elements(t) % nd
      nb = elem_lists(col) % elements(t) % nb
      !!$omp target
      CALL LocalMatrixVec(  Element, n, nd+nb, nb, VecAsm, elem_lists(col)% elements(t) % nodes)
      !!$omp end target
    END DO
  END DO

  !return
  totelem = 0

  CALL DefaultFinishBulkAssembly()

  nColours = GetNOFBoundaryColours(Solver)
  VecAsm = (nColours > 1) .OR. (nthr == 1)

  CALL ResetTimer(Caller//'BCAssembly')

  !! don't touch boundary stuff yet
  !!$OMP PARALLEL &
  !!$OMP SHARED(Active, Solver, nColours, VecAsm) &
  !!$OMP PRIVATE(t, Element, n, nd, nb, col, InitHandles) & 
  !!$OMP REDUCTION(+:totelem) DEFAULT(NONE)
  DO col=1,nColours
    !!$OMP SINGLE
    CALL Info('ModelPDEthreaded','Assembly of boundary colour: '//I2S(col),Level=10)
    Active = GetNOFBoundaryActive(Solver)
    !!$OMP END SINGLE

       InitHandles = .TRUE. 
       !!$OMP DO
       DO t=1,Active
          Element => GetBoundaryElement(t)
          ! WRITE (*,*) Element % ElementIndex
          totelem = totelem + 1
          IF(ActiveBoundaryElement(Element)) THEN
             n  = GetElementNOFNodes(Element)
             nd = GetElementNOFDOFs(Element)
             nb = GetElementNOFBDOFs(Element)
             CALL LocalMatrixBC(  Element, n, nd+nb, nb, VecAsm, InitHandles )
          END IF
       END DO
       !!$OMP END DO
    END DO
    !!$OMP END PARALLEL

    CALL CheckTimer(Caller//'BCAssembly',Delete=.TRUE.)
        
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    ! IF( Solver % Variable % NonlinConverged == 1 ) EXIT

  ! END DO

  CALL DefaultFinish()

  
CONTAINS

! Assembly of the matrix entries arising from the bulk elements. SIMD version.
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixVec( Element, n, nd, nb, VecAsm, Nodes)
!------------------------------------------------------------------------------
    !USE LinearForms
    use MyLinearforms
    USE Integration
    use iso_c_binding
    use Types
    IMPLICIT NONE
  !$omp declare target to(LocalMatrixVec)


    INTEGER, INTENT(IN) :: n, nd, nb
    TYPE(Element_t), POINTER:: Element
    LOGICAL, INTENT(IN) :: VecAsm
    TYPE(Nodes_t), POINTER :: Nodes
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:,:),dBasisdx(:,:,:), DetJ(:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: MASS(:,:), STIFF(:,:), FORCE(:)
    REAL(KIND=dp), SAVE, ALLOCATABLE  :: DiffCoeff(:), ConvCoeff(:), ReactCoeff(:), &
         TimeCoeff(:), SourceCoeff(:), Velo1Coeff(:), Velo2Coeff(:), Velo3Coeff(:)
    REAL(KIND=dp), SAVE, ALLOCATABLE  :: VeloCoeff(:,:)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim,ngp,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    LOGICAL, SAVE :: FirstTime=.TRUE.
    TYPE(element_t) :: concrete_element

    !!$OMP THREADPRIVATE(Basis, dBasisdx, DetJ, &
    !!$OMP               MASS, STIFF, FORCE, Nodes, &
    !!$OMP               SourceCoeff, DiffCoeff, ReactCoeff, TimeCoeff, &
    !!$OMP               ConvCoeff, Velo1Coeff, Velo2Coeff, Velo3Coeff, VeloCoeff )
    !DIR$ ATTRIBUTES ALIGN:64 :: Basis, dBasisdx, DetJ
    !DIR$ ATTRIBUTES ALIGN:64 :: MASS, STIFF, FORCE
!------------------------------------------------------------------------------




    dim = CoordinateSystemDimension()
    IP = GaussPoints( Element )
    ngp = IP % n

    ! THIS IS UGLY AND DIRTY - assuming all elements are same!
    IF (FirstTime) THEN
      ALLOCATE(DiffCoeff(ngp), ConvCoeff(ngp), ReactCoeff(ngp), &
           TimeCoeff(ngp), SourceCoeff(ngp), Velo1Coeff(ngp), Velo2Coeff(ngp),&
           Velo3Coeff(ngp), VeloCoeff(ngp,3), MASS(nd,nd), STIFF(nd,nd), FORCE(nd),&
           Basis(ngp,nd), dBasisdx(ngp,nd,3), DetJ(ngp), &
           STAT=allocstat)
      DiffCoeff = 1.0_dp
      ConvCoeff=0.0_dp
      ReactCoeff=0.0_dp
      TimeCoeff=0.0_dp
      SourceCoeff=1.0_dp
      Velo1Coeff=0.0_dp
      Velo2Coeff=0.0_dp
      Velo3Coeff=0.0_dp
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
      FirstTime=.FALSE.
    END IF


    ! Do this on cpu
    !CALL GetElementNodesVec( Nodes, UElement=Element )

    ! Initialize
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp

    ! Numerical integration:
    ! Compute basis function values and derivatives at integration points
    !--------------------------------------------------------------

    print *, element%bodyid
    !!$omp target map(concrete_element, concrete_element % type, concrete_element % type % dimension)
    !$omp target map(to: element, element%bodyid, element%type, element%type%dimension)
    !print *, element
    print *, element%type%dimension
    !$omp end target


    !!$omp target data map(to:element, element%type)
    !!$omp target
    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, detJ, &
         SIZE(Basis,2), Basis, dBasisdx )
    !!$omp end target
    !!$omp end target data

    ! Compute actual integration weights (recycle the memory space of DetJ)
    DO t=1,ngp
      DetJ(t) = IP % s(t) * Detj(t)
    END DO


    !CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % Dimension , dBasisdx, DetJ, STIFF, DiffCoeff )

    !$omp target data map(to: t, ngp, nd, dBasisdx, detj, DiffCoeff, element, element% type) map(tofrom:stiff)
    !$omp target
    !print *, omp_is_initial_device()
    !CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % Dimension , dBasisdx, DetJ, STIFF, DiffCoeff )
    CALL LinearForms_GradUdotGradU(ngp, nd, Element % type % dimension, dBasisdx, DetJ, STIFF, DiffCoeff )
    !CALL LinearForms_UdotF(ngp, nd, Basis, DetJ, SourceCoeff, FORCE)
    !$omp end target
    !$omp end target data

    stop

    ! DEBUG
    !IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
    CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element, VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixVec
!------------------------------------------------------------------------------


! Assembly of the matrix entries arising from the Neumann and Robin conditions
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixBC( Element, n, nd, nb, VecAsm, InitHandles )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: n, nd, nb
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: VecAsm
    LOGICAL, INTENT(INOUT) :: InitHandles
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: F,C,Ext, Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd), LOAD(n)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q,dim
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(ValueList_t), POINTER :: BC
       
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueHandle_t), SAVE :: Flux_h, Robin_h, Ext_h

    SAVE Nodes
    !!$OMP THREADPRIVATE(Nodes,Flux_h,Robin_h,Ext_h)
!------------------------------------------------------------------------------
    BC => GetBC(Element)
    IF (.NOT.ASSOCIATED(BC) ) RETURN

    IF( InitHandles ) THEN
      CALL ListInitElementKeyword( Flux_h,'Boundary Condition','Field Flux')
      CALL ListInitElementKeyword( Robin_h,'Boundary Condition','Robin Coefficient')
      CALL ListInitElementKeyword( Ext_h,'Boundary Condition','External Field')
      InitHandles = .FALSE.
    END IF
    
    dim = CoordinateSystemDimension()

    CALL GetElementNodes( Nodes, UElement=Element )
    STIFF = 0._dp
    FORCE = 0._dp
    LOAD = 0._dp
           

    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      Weight = IP % s(t) * DetJ

      ! Evaluate terms at the integration point:
      !------------------------------------------

      ! Given flux:
      ! -----------
      F = ListGetElementReal( Flux_h, Basis, Element, Found )
      IF( Found ) THEN
        FORCE(1:nd) = FORCE(1:nd) + Weight * F * Basis(1:nd)
      END IF

      ! Robin condition (C*(u-u_0)):
      ! ---------------------------
      C = ListGetElementReal( Robin_h, Basis, Element, Found )

      IF( Found ) THEN
        Ext = ListGetElementReal( Ext_h, Basis, Element, Found )
        DO p=1,nd
          DO q=1,nd
            STIFF(p,q) = STIFF(p,q) + Weight * C * Basis(q) * Basis(p)
          END DO
        END DO
        FORCE(1:nd) = FORCE(1:nd) + Weight * C * Ext * Basis(1:nd)
      END IF
    END DO
    
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element,VecAssembly=VecAsm)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrixBC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE AdvDiffSolver
!------------------------------------------------------------------------------
