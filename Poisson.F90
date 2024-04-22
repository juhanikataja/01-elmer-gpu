
function MakeReferenceElement() result(element)
use Types
implicit none

type(element_t) :: element
type(nodes_t) :: nodes

end function

!call myElementMetric(nd,Nodes,dim,Metric,DetG,dBasisdx(i,:,:),LtoGMap)
subroutine myElementMetric(nDOFs,Nodes,dim,Metric,DetG,dLBasisdx,LtoGMap)
!$omp declare target
use Types
use DefUtils
implicit none
!------------------------------------------------------------------------------
INTEGER :: nDOFs, dim           !< Number of active nodes in element, dimension of space
TYPE(Nodes_t)    :: Nodes       !< Element nodal coordinates
REAL(KIND=dp), intent(out) :: Metric(3,3)    !< Contravariant metric tensor
REAL(KIND=dp), intent(in) :: dLBasisdx(nDOFs,3) !< Derivatives of element basis function with respect to local coordinates
REAL(KIND=dp), intent(out) :: DetG           !< SQRT of determinant of metric tensor
REAL(KIND=dp), intent(out) :: LtoGMap(3,3)   !< Transformation to obtain the referential description of the spatial gradient
!LOGICAL :: Success              !< Returns .FALSE. if element is degenerate
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
REAL(KIND=dp) :: dx(3,3),G(3,3),GI(3,3),s
REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z
INTEGER :: GeomId     
INTEGER :: cdim,i,j,k,n,imin,jmin
!------------------------------------------------------------------------------
x => Nodes % x
y => Nodes % y
z => Nodes % z

cdim = CoordinateSystemDimension()
n = MIN( SIZE(x), nDOFs )
!dim  = elm % TYPE % DIMENSION


!------------------------------------------------------------------------------
!    Partial derivatives of global coordinates with respect to local coordinates
!------------------------------------------------------------------------------
DO i=1,dim
  dx(1,i) = SUM( x(1:n) * dLBasisdx(1:n,i) )
  dx(2,i) = SUM( y(1:n) * dLBasisdx(1:n,i) )
  dx(3,i) = SUM( z(1:n) * dLBasisdx(1:n,i) )
END DO
!------------------------------------------------------------------------------
!    Compute the covariant metric tensor of the element coordinate system
!------------------------------------------------------------------------------
DO i=1,dim
  DO j=1,dim
    s = 0.0_dp
    DO k=1,cdim
      s = s + dx(k,i)*dx(k,j)
    END DO
    G(i,j) = s
  END DO
END DO
!------------------------------------------------------------------------------
!    Convert the metric to contravariant base, and compute the SQRT(DetG)
!------------------------------------------------------------------------------
SELECT CASE( dim )
!------------------------------------------------------------------------------
!      Line elements
!------------------------------------------------------------------------------
CASE (1)
  DetG  = G(1,1)


  Metric(1,1) = 1.0d0 / DetG
  DetG  = SQRT( DetG )

  !------------------------------------------------------------------------------
  !      Surface elements
  !------------------------------------------------------------------------------
CASE (2)
  DetG = ( G(1,1)*G(2,2) - G(1,2)*G(2,1) )


  Metric(1,1) =  G(2,2) / DetG
  Metric(1,2) = -G(1,2) / DetG
  Metric(2,1) = -G(2,1) / DetG
  Metric(2,2) =  G(1,1) / DetG
  DetG = SQRT(DetG)

  !------------------------------------------------------------------------------
  !      Volume elements
  !------------------------------------------------------------------------------
CASE (3)
  DetG = G(1,1) * ( G(2,2)*G(3,3) - G(2,3)*G(3,2) ) + &
    G(1,2) * ( G(2,3)*G(3,1) - G(2,1)*G(3,3) ) + &
    G(1,3) * ( G(2,1)*G(3,2) - G(2,2)*G(3,1) )


  CALL InvertMatrix3x3( G,GI,detG )
  Metric = GI
  DetG = SQRT(DetG)
END SELECT

!--------------------------------------------------------------------------------------
!    Construct a transformation X = LtoGMap such that (grad B)(f(p)) = X(p) Grad b(p),
!    with Grad the gradient with respect to the reference element coordinates p and 
!    the referential description of the spatial field B(x) satisfying B(f(p)) = b(p).
!    If cdim > dim (e.g. a surface embedded in the 3-dimensional space), X is
!    the transpose of the pseudo-inverse of Grad f.
!-------------------------------------------------------------------------------
DO i=1,cdim
  DO j=1,dim
    s = 0.0d0
    DO k=1,dim
      s = s + dx(i,k) * Metric(k,j)
    END DO
    LtoGMap(i,j) = s
  END DO
END DO

end subroutine  myElementMetric


!SUBROUTINE AdvDiffSolver_init( Model,Solver,dt,TransientSimulation )
!!------------------------------------------------------------------------------
  !USE DefUtils
  !IMPLICIT NONE
!!------------------------------------------------------------------------------
  !TYPE(Solver_t) :: Solver
  !TYPE(Model_t) :: Model
  !REAL(KIND=dp) :: dt
  !LOGICAL :: TransientSimulation
!!------------------------------------------------------------------------------
  !CHARACTER(*), PARAMETER :: Caller = 'AdvDiffSolver_init'
  !
  !IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
    !CALL Fatal(Caller,'Implemented only in cartesian coordinates')
  !END IF
!
!END SUBROUTINE AdvDiffSolver_Init

!------------------------------------------------------------------------------
SUBROUTINE AdvDiffSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Integration

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
  INTEGER :: ngp, i, dim
  LOGICAL :: Found, VecAsm, InitHandles
  integer, allocatable :: n_active_in_col(:)

  REAL(KIND=dp), ALLOCATABLE :: refBasis(:,:), refdBasisdx(:,:,:)
  TYPE(GaussIntegrationPoints_t) :: refIP

  type :: elem_ptr
    type(Element_t), pointer :: p
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


  nColours = GetNOFColours(Solver)

  allocate(elem_lists(nColours)) ! TODO: deallocate too
  allocate(n_active_in_col(nColours)) ! TODO: deallocate too

  do col = 1, nColours
    active = GetNOFActive(Solver)
    allocate(elem_lists(col) % elements(active))
  end do

  !!$omp target enter data map(to:solver%mesh%elements)
  !!$omp target enter data map(to:elem_lists)

  !! Tabulate elements and their ndofs/nnodes/nb 
  nColours = GetNOFColours(Solver)
  !$OMP PARALLEL &
  !$OMP SHARED(Active, Solver, nColours, VecAsm, elem_lists) &
  !$OMP PRIVATE(t, Element, n, nd, nb, col, InitHandles) & 
  !$OMP REDUCTION(max:MaxNumNodes) DEFAULT(none)
  do col=1,ncolours
    !$OMP SINGLE
    active = GetNOFActive(Solver)
    !!$omp target enter data map(to:elem_lists(col) % elements)
    !$OMP END SINGLE

    !$OMP DO
    do t=1,active
      Element => GetActiveElement(t)
      elem_lists(col) % elements(t) % p => Element
      elem_lists(col) % elements(t) % n = GetElementNOFNodes(Element)
      elem_lists(col) % elements(t) % nd = GetElementNOFDOFs(Element)
      elem_lists(col) % elements(t) % nb = GetElementNOFBDOFs(Element)
      MaxNumNodes = max(MaxNumNodes,elem_lists(col) % elements(t) % n)
    end do
    !$OMP END DO
  end do
  !$OMP END PARALLEL

  CALL CheckTimer(Caller//'BulkAssembly', Delete=.TRUE.)

  nColours = GetNOFColours(Solver)

  print *, '==BASIS FUNCTION VALUES========================='

  Element => elem_lists(1) % elements(1) % p
  dim = CoordinateSystemDimension()
  refIP = GaussPoints( Element )
  ngp = refIP % n
  nd = GetElementNOFDOFs(Element)
  print *, 'NGP:', ngp
  allocate(refbasis(ngp, nd), refdbasisdx(ngp, nd, 3))
  do i=1,ngp
    call NodalBasisFunctions(nd, refBasis(i,:), element, refIP%u(i), refIP%v(i), refIP%w(i))
    call NodalFirstDerivatives(nd, refdBasisdx(i,:,:), element, refIP%u(i), refIP%v(i), refIP%w(i))
    write (*,'(12F7.3)') refbasis(i,:)
  end do

  print *, '==dbasis function values========================'

  DO col=1,nColours

    CALL Info( Caller,'Assembly of colour: '//I2S(col),Level=1) ! TODO: this goes away
    Active = GetNOFActive(Solver) ! TODO: this goes away

    DO t=1,Active
      totelem = totelem + 1
      Element => elem_lists(col) % elements(t) % p
      n = elem_lists(col) % elements(t) % n
      nd = elem_lists(col) % elements(t) % nd
      nb = elem_lists(col) % elements(t) % nb
      CALL LocalMatrixVec(  Element, n, nd+nb, nb, VecAsm )
    END DO
  END DO

  !return
  totelem = 0

  CALL DefaultFinishBulkAssembly()

  nColours = GetNOFBoundaryColours(Solver)
  VecAsm = (nColours > 1) .OR. (nthr == 1)

  CALL ResetTimer(Caller//'BCAssembly')

  !! don't touch boundary stuff yet
  !$OMP PARALLEL &
  !$OMP SHARED(Active, Solver, nColours, VecAsm) &
  !$OMP PRIVATE(t, Element, n, nd, nb, col, InitHandles) & 
  !$OMP REDUCTION(+:totelem) DEFAULT(NONE)
  DO col=1,nColours
    !$OMP SINGLE
    CALL Info('ModelPDEthreaded','Assembly of boundary colour: '//I2S(col),Level=10)
    Active = GetNOFBoundaryActive(Solver)
    !$OMP END SINGLE

       InitHandles = .TRUE. 
       !$OMP DO
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
       !$OMP END DO
    END DO
    !$OMP END PARALLEL

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

! Assembly of the matrix entries arising from the bulk elements. Offload compatible version
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrixVec( Element, n, nd, nb, VecAsm )
!------------------------------------------------------------------------------
    USE LinearForms
    USE Integration
    use iso_c_binding
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, nd, nb
    TYPE(Element_t), POINTER:: Element
    LOGICAL, INTENT(IN) :: VecAsm
    TYPE(element_t) :: refElement
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:,:),dBasisdx(:,:,:), DetJ(:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: MASS(:,:), STIFF(:,:), FORCE(:)
    REAL(KIND=dp), SAVE, ALLOCATABLE  :: DiffCoeff(:), ConvCoeff(:), ReactCoeff(:), &
         TimeCoeff(:), SourceCoeff(:), Velo1Coeff(:), Velo2Coeff(:), Velo3Coeff(:)
    REAL(KIND=dp), SAVE, ALLOCATABLE  :: VeloCoeff(:,:)
    REAL(KIND=dp), target :: LtoGMap(3,3), metric(3,3), detg
    LOGICAL :: Stat,Found
    INTEGER :: j,k,m,i,t,p,q,dim,ngp,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    type(Nodes_t), SAVE :: refNodes
    type(ElementType_t) :: refElementType
    LOGICAL, SAVE :: FirstTime=.TRUE.
    real(KIND=dp) :: dLBasisdx(nd, 3)
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
      allocate(refnodes %x(maxnumnodes), &
	      refnodes %y(maxnumnodes), &
	      refnodes %z(maxnumnodes), stat=allocstat)
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
    END IF

      
    CALL GetElementNodesVec( Nodes, UElement=Element )
    print *, '==Element nodes================================='
    write (*,'(12F7.3)') nodes % x(1:nd)
    write (*,'(12F7.3)') nodes % y(1:nd)
    write (*,'(12F7.3)') nodes % z(1:nd)
    print *, '================================================'

    refElement = element

    !call GetRefPElementNodes(Element%type, refnodes % x, refnodes % y, refnodes % z)
    !INTEGER :: nDOFs                !< Number of active nodes in element
    !TYPE(Nodes_t)    :: Nodes       !< Element nodal coordinates
    !REAL(KIND=dp) :: Metric(:,:)    !< Contravariant metric tensor
    !REAL(KIND=dp) :: dLBasisdx(:,:) !< Derivatives of element basis function with respect to local coordinates
    !REAL(KIND=dp) :: DetG           !< SQRT of determinant of metric tensor
    !REAL(KIND=dp) :: LtoGMap(3,3)   !< Transformation to obtain the referential description of the spatial gradient

    ! TODO: wip loop over integration points and calculate myElementMetric
    
    ! Move this outside target region
    print *, '==basis function values========================='
    do i=1,ngp
      call NodalBasisFunctions(nd, Basis(i,:), element, IP%u(i), IP%v(i), IP%w(i))
      call NodalFirstDerivatives(nd, dBasisdx(i,:,:), element, IP%u(i), IP%v(i), IP%w(i))
      write (*,'(12F7.3)') basis(i,:)
    end do
    print *, '==dbasis function values========================'

    do i = 1,ngp
      write (*, '(I3)') i
      write (*,'(A)', advance='no') 'dx1 '
      write (*,'(12F8.3)') dBasisdx(i,:,1)
      write (*,'(A)', advance='no') 'dx2 '
      write (*,'(12F8.3)') dBasisdx(i,:,2)
      write (*,'(A)', advance='no') 'dx3 '
      write (*,'(12F8.3)') dBasisdx(i,:,3)
    end do
    print *, '================================================'

    do i=1,ngp
      do j = 1,dim
        do m = 1,nd
          dLBasisdx(m,j) = dBasisdx(i,m,j)
        end do
      end do
      call myElementMetric(nd,Nodes,dim,Metric,DetG,dbasisdx(i,:,:),LtoGMap)
    write (*,'(e10.3)') sqrt(detG)

      dbasisdx(i,:,:) = 0_dp
      do m = 1,nd
        do j=1,dim
          do k=1,dim
            dbasisdx(i,m,j) = dbasisdx(i,m,j) + dLbasisdx(m,k)*LtoGMap(j,k)
          end do 
        end do
      end do
      write (*, '(I3)') i
      write (*,'(A)', advance='no') 'dx1 '
      write (*,'(12F8.3)') dBasisdx(i,:,1)
      write (*,'(A)', advance='no') 'dx2 '
      write (*,'(12F8.3)') dBasisdx(i,:,2)
      write (*,'(A)', advance='no') 'dx3 '
      write (*,'(12F8.3)') dBasisdx(i,:,3)
    end do


    ! stop

    !SUBROUTINE GetRefPElementNodes(Element, U, V, W)
    !TYPE(ElementType_t) :: Element
    !REAL(KIND=dp) :: U(:), V(:), W(:)

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
    !    print *, element
    !print *, element%type%dimension
    !$omp end target

    !$omp target data map(to:element, element%type)
    !$omp target
    dbasisdx = 0._dp
    stat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, IP % W, detJ, SIZE(Basis,2), Basis, dBasisdx )
    print *, '==dbasis function values the old way============'
    do i = 1,ngp
      write (*, '(I3)') i
      write (*,'(A)', advance='no') 'dx1 '
      write (*,'(12F8.3)') dBasisdx(i,:,1)
      write (*,'(A)', advance='no') 'dx2 '
      write (*,'(12F8.3)') dBasisdx(i,:,2)
      write (*,'(A)', advance='no') 'dx3 '
      write (*,'(12F8.3)') dBasisdx(i,:,3)
    end do
    print *, '================================================'
    !$omp end target
    !$omp end target data
    stop

    ! Compute actual integration weights (recycle the memory space of DetJ)
    DO t=1,ngp
      DetJ(t) = IP % s(t) * Detj(t)
    END DO

    !!$omp target data map(to: ngp, nd, dBasisdx, detj, DiffCoeff) map(tofrom:stiff)
    !$omp target
    ! CALL LinearForms_GradUdotGradU(ngp, nd, Element % TYPE % Dimension , dBasisdx, DetJ, STIFF, DiffCoeff )
    CALL LinearForms_GradUdotGradU(ngp, nd, 3, dBasisdx, DetJ, STIFF, DiffCoeff )
    CALL LinearForms_UdotF(ngp, nd, Basis, DetJ, SourceCoeff, FORCE)
    !$omp end target
    !!$omp end target data
 
    
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
    !$OMP THREADPRIVATE(Nodes,Flux_h,Robin_h,Ext_h)
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
