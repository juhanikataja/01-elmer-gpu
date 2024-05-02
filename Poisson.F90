module LocalMV

CONTAINS
!#ifdef BUILDLOCALMV

SUBROUTINE ModuleLocalMatrixVecSO( n, nd, nb, VecAsm, x, y, z, dim, refbasis, refdBasisdx, ip, ngp)
!------------------------------------------------------------------------------
    !USE LinearForms
    USE Integration
    !use iso_c_binding
#ifdef ASSEMBLE
    use DefUtils ! TODO: defaultupdateequations is here but defutils may not be used due to threadprivate module variables if
                 ! declare target is set
#endif
    IMPLICIT NONE
!$omp declare target


    INTEGER, INTENT(IN) :: n, nd, nb
    !TYPE(Element_t), POINTER:: Element
    LOGICAL, INTENT(IN) :: VecAsm
    !TYPE(Nodes_t), intent(in) :: Nodes
    real(kind=dp), intent(in) :: x(:), y(:), z(:)
    INTEGER :: dim
    real(kind=dp), intent(in) :: refbasis(:,:), refdbasisdx(:,:,:)
    TYPE(GaussIntegrationPoints_t), intent(in) :: IP
!------------------------------------------------------------------------------
#define nd_ 4
#define ngp_ 4
    REAL(KIND=dp) :: Basis(ngp,nd)
    real(kind=dp) :: dBasisdx(ngp,nd,3), DetJ(ngp)
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd)
    REAL(KIND=dp) :: DiffCoeff(ngp), SourceCoeff(ngp)
    REAL(KIND=dp) :: LtoGMap(3,3), metric(3,3), detg
    LOGICAL :: Stat,Found
    INTEGER :: j,k,m,i,t,p,q,ngp,allocstat

    
#ifdef DEBUGPRINT
    INTEGER :: round = 1 ! TODO
#endif
    !type(ElementType_t) :: refElementType
    !LOGICAL, SAVE :: FirstTime=.TRUE.
    real(KIND=dp) :: dLBasisdx(nd, 3)
!------------------------------------------------------------------------------

    
    !dim = CoordinateSystemDimension()
    !IP = GaussPoints( Element )
    !ngp = IP % n

    DiffCoeff = 1._dp ! TODO: Material parameters must be communicated somehow to the solver

!#ifdef DOIT
    do i=1,ngp
      do j = 1,dim
        do m = 1,nd
          dLBasisdx(m,j) = refdBasisdx(i,m,j)
        end do
      end do

      call myElementMetric(nd,n,x,y,z,dim,Metric,DetG,dlbasisdx(:,:),LtoGMap)

      detj(i) = detg

      dbasisdx(i,:,:) = 0_dp
      do m = 1,nd
        do j=1,dim
          do k=1,dim
            dbasisdx(i,m,j) = dbasisdx(i,m,j) + dLbasisdx(m,k)*LtoGMap(j,k)
          end do 
        end do
      end do
    end do


    !basis(:,:) = refbasis(:,:)


    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp

    ! PART ONE: Collect local stiffness to STIFF using experimental method

    DO t=1,ngp
      DetJ(t) = IP % s(t) * Detj(t)
    END DO

    CALL LinearForms_GradUdotGradU(ngp, nd, dim, dBasisdx, DetJ, STIFF, DiffCoeff )
    CALL LinearForms_UdotF(ngp, nd, refBasis, DetJ, SourceCoeff, FORCE)
!#endif

#ifdef DEBUGPRINT

  if (round < 3) then
    print *, '' 
    DO t=1,nd 
      write (*, '(12F9.4)') stiff(t,:)
    end do
  end if
#endif

#ifdef DEBUGPRINT
  if (round < 3) then
    round = round + 1
  end if
#endif
 
    
#ifdef ASSEMBLE
    ! DEBUG
    !IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element)
    !CALL CondensateP( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element, VecAssembly=VecAsm)
#endif

!------------------------------------------------------------------------------
END SUBROUTINE ModuleLocalMatrixVecSO

!call myElementMetric(nd,Nodes,dim,Metric,DetG,dBasisdx(i,:,:),LtoGMap)

!#endif  ! BUILDLOCALMV

subroutine myElementMetric(nDOFs,nnodes,x,y,z,dim,Metric,DetG,dLBasisdx,LtoGMap)

use Types, only: dp, nodes_t
!use DefUtils
implicit none
!$omp declare target
!------------------------------------------------------------------------------
INTEGER :: nDOFs, dim, nnodes   !< Number of active nodes in element, dimension of space, 
! TYPE(Nodes_t), intent(in)  :: Nodes       !< Element nodal coordinates
real(kind=dp), intent(in)  :: x(nnodes),y(nnodes),z(nnodes)       !< Element nodal coordinates
REAL(KIND=dp), intent(out) :: Metric(3,3)    !< Contravariant metric tensor
REAL(KIND=dp), intent(in)  :: dLBasisdx(nDOFs,3) !< Derivatives of element basis function with respect to local coordinates
REAL(KIND=dp), intent(out) :: DetG           !< SQRT of determinant of metric tensor
REAL(KIND=dp), intent(out) :: LtoGMap(3,3)   !< Transformation to obtain the referential description of the spatial gradient
!LOGICAL :: Success              !< Returns .FALSE. if element is degenerate
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
REAL(KIND=dp) :: dx(3,3),G(3,3),GI(3,3),s
!REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z
INTEGER :: GeomId     
INTEGER :: cdim,i,j,k,n,imin,jmin
!------------------------------------------------------------------------------
!associate(x=> nodes %x, y => nodes % y, z => nodes % z)
! x => Nodes % x
! y => Nodes % y
! z => Nodes % z

cdim = dim !CoordinateSystemDimension()
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
    !end associate
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
#ifdef CHECKDIMMETRIC
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
#endif

  DetG = G(1,1) * ( G(2,2)*G(3,3) - G(2,3)*G(3,2) ) + &
    G(1,2) * ( G(2,3)*G(3,1) - G(2,1)*G(3,3) ) + &
    G(1,3) * ( G(2,1)*G(3,2) - G(2,2)*G(3,1) )


  CALL InvertMatrix3x3( G,GI,detG )
  Metric = GI
  DetG = SQRT(DetG)

#ifdef CHECKDIMMETRIC
END SELECT
#endif

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

  PURE SUBROUTINE InvertMatrix3x3( G,GI,detG )
    use Types, only: dp
    !$omp declare target
!------------------------------------------------------------------------------
    REAL(KIND=dp), INTENT(IN) :: G(3,3)
    REAL(KIND=dp), INTENT(OUT) :: GI(3,3)
    REAL(KIND=dp), intent(in) :: detG
    REAL(KIND=dp) :: s
!------------------------------------------------------------------------------
    s = 1.0 / DetG
    
    GI(1,1) =  s * (G(2,2)*G(3,3) - G(3,2)*G(2,3));
    GI(2,1) = -s * (G(2,1)*G(3,3) - G(3,1)*G(2,3));
    GI(3,1) =  s * (G(2,1)*G(3,2) - G(3,1)*G(2,2));
    
    GI(1,2) = -s * (G(1,2)*G(3,3) - G(3,2)*G(1,3));
    GI(2,2) =  s * (G(1,1)*G(3,3) - G(3,1)*G(1,3));
    GI(3,2) = -s * (G(1,1)*G(3,2) - G(3,1)*G(1,2));

    GI(1,3) =  s * (G(1,2)*G(2,3) - G(2,2)*G(1,3));
    GI(2,3) = -s * (G(1,1)*G(2,3) - G(2,1)*G(1,3));
    GI(3,3) =  s * (G(1,1)*G(2,2) - G(2,1)*G(1,2));
!------------------------------------------------------------------------------
  END SUBROUTINE InvertMatrix3x3

  SUBROUTINE LinearForms_GradUdotGradU(m, n, dim, GradU, weight, G, alpha)
    USE Types, ONLY: dp, VECTOR_BLOCK_LENGTH, VECTOR_SMALL_THRESH
    IMPLICIT NONE
    !$omp declare target 
    INTEGER, INTENT(IN) :: m, n, dim
    REAL(KIND=dp) CONTIG, INTENT(IN) :: GradU(:,:,:), weight(:)
    REAL(KIND=dp) CONTIG, INTENT(INOUT) :: G(:,:)
    REAL(KIND=dp) CONTIG, INTENT(IN)  :: alpha(:)

    INTEGER :: i, ii, iin, j, l, k, kk, ldbasis, ldk, blklen

    ldbasis = SIZE(GradU,1)
    ldk = SIZE(G,1)

    DO ii=1,m,VECTOR_BLOCK_LENGTH
      iin=MIN(ii+VECTOR_BLOCK_LENGTH-1,m)
      blklen=iin-ii+1

      DO j=1,n
        DO i=1,n
          DO k=1,dim
            DO l=ii,iin
              G(i,j) = G(i,j) + GradU(l,i,k)*GradU(l,j,k)*weight(l)*alpha(l)
            END DO
          END DO
        END DO
      END DO
        
    END DO ! Vector blocks
  END SUBROUTINE LinearForms_GradUdotGradU

  SUBROUTINE LinearForms_UdotF(m, n, U, weight, F, UdotF)
  USE Types, ONLY: dp, VECTOR_BLOCK_LENGTH, VECTOR_SMALL_THRESH
    IMPLICIT NONE
    !$omp declare target 

    INTEGER, INTENT(IN) :: m, n
    REAL(KIND=dp) CONTIG, INTENT(IN) :: U(:,:), F(:), weight(:)
    REAL(KIND=dp) CONTIG, INTENT(INOUT) :: UdotF(:)

    INTEGER :: i, ii, iin, j, blklen, l

    DO ii=1,m,VECTOR_BLOCK_LENGTH
      iin = MIN(ii+VECTOR_BLOCK_LENGTH-1,m)
      blklen= iin-ii+1
      ! Project local F to global basis
      DO i=1,n
        DO l=ii,iin
          UdotF(i) = UdotF(i) + U(l,i)*F(l)*weight(l)
        END DO
      END DO
    END DO
  END SUBROUTINE LinearForms_UdotF

end module

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
  USE Integration
  USE DefUtils
  USE Types
  use LocalMV
  use omp_lib

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
  INTEGER :: iter, maxiter, nColours, col, totelem, nthr, state, MaxNumNodes, MinNumNodes
  INTEGER :: ngp, i, dim, coorddim
  LOGICAL :: Found, VecAsm, InitHandles
  integer, allocatable :: n_active_in_col(:)
  real(KIND=dp), allocatable :: color_timing(:)
  type(nodes_t) :: Nodes

#ifdef _OPENMP
  LOGICAL :: initial_device
#endif
  REAL(KIND=dp), ALLOCATABLE :: refBasis(:,:), refdBasisdx(:,:,:)
  TYPE(GaussIntegrationPoints_t) :: refIP

  type :: elem_ptr
    type(Element_t), pointer :: p
    !type(Nodes_t) :: nodes
    integer :: n, nd, nb                        ! nof nodes, nof dofs, nofbdofs
  end type

  type :: elem_list_t
    type(elem_ptr), allocatable :: elements(:)
    real(kind=dp), allocatable :: x(:,:), y(:,:), z(:,:)
    !real(kind=dp), allocatable, target :: x(:,:), y(:,:), z(:,:)
  end type

  TYPE(elem_list_t), allocatable :: elem_lists(:)

  CHARACTER(*), PARAMETER :: Caller = 'AdvDiffSolver'


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

#ifdef _OPENMP
  initial_device = .true.
#endif

  coorddim = CoordinateSystemDimension()
  !$omp target
#ifdef _OPENMP
  initial_device = omp_is_initial_device()
#endif
  !$omp end target
#ifdef _OPENMP
  print *, 'initial_device:', initial_device
#endif


  nColours = GetNOFColours(Solver)

  allocate(elem_lists(nColours)) ! TODO: deallocate too
  allocate(n_active_in_col(nColours)) ! TODO: deallocate too

  allocate(color_timing(nColours)) ! TODO: deallocate too

  do col = 1, nColours
    active = GetNOFActive(Solver)
    allocate(elem_lists(col) % elements(active))
  end do

  MaxNumNodes = 0
  MinNumNodes = 10000
  !!$omp target enter data map(to:solver%mesh%elements)
  !!$omp target enter data map(to:elem_lists(1:nColours))

  !! Tabulate elements and their ndofs/nnodes/nb 
  nColours = GetNOFColours(Solver)
  !!$OMP PARALLEL &
  !!$OMP SHARED(Active, Solver, nColours, VecAsm, elem_lists) &
  !!$OMP PRIVATE(t, Element, n, nd, nb, col, InitHandles) & 
  !!$OMP REDUCTION(max:MaxNumNodes) DEFAULT(none)
  do col=1,ncolours
    !!$OMP SINGLE
    !!$omp target enter data map(to:elem_lists(col)) nowait
    active = GetNOFActive(Solver)
    !!$OMP END SINGLE

    !!$OMP DO
    do state=1,2
      do t=1,active
        Element => GetActiveElement(t)
        CALL GetElementNodesVec( nodes,  Element)
        if (state == 1) then
          elem_lists(col) % elements(t) % p => Element
          elem_lists(col) % elements(t) % n = GetElementNOFNodes(Element)
          elem_lists(col) % elements(t) % nd = GetElementNOFDOFs(Element)
          elem_lists(col) % elements(t) % nb = GetElementNOFBDOFs(Element)

          ! CALL GetElementNodesVec( elem_lists(col) % elements(t) % nodes,  elem_lists(col) % elements(t) % p)

          MaxNumNodes = max(MaxNumNodes, elem_lists(col) % elements(t) % n)
          MinNumNodes = min(MinNumNodes, elem_lists(col) % elements(t) % n)
        else
          elem_lists(col) % x(:, t) = nodes % x(:)
          elem_lists(col) % y(:, t) = nodes % y(:)
          elem_lists(col) % z(:, t) = nodes % z(:)
        end if
      end do

      if (state == 1) then
        allocate( elem_lists(col) % x(maxnumnodes, active), &
        elem_lists(col) % y(maxnumnodes, active), &
        elem_lists(col) % z(maxnumnodes, active))
      end if
    end do

  !!$OMP END DO
  end do
  !!$OMP END PARALLEL


  print *, 'Max/Min NumNodes:', MaxNumNodes, MinNumNodes

  nColours = GetNOFColours(Solver)

  print *, '==REFERENCE BASIS FUNCTION VALUES==============='

  Element => elem_lists(1) % elements(1) % p
  dim = CoordinateSystemDimension()
  refIP = GaussPoints( Element )
  ngp = refIP % n
  nd = GetElementNOFDOFs(Element)
  print *, 'NGP:', ngp
  allocate(refbasis(ngp, nd), refdbasisdx(ngp, nd, 3))
  refbasis = 0_dp
  refdBasisdx = 0_dp
  do i=1,ngp
    call NodalBasisFunctions(nd, refBasis(i,:), element, refIP%u(i), refIP%v(i), refIP%w(i))
    call NodalFirstDerivatives(nd, refdBasisdx(i,:,:), element, refIP%u(i), refIP%v(i), refIP%w(i))
    write (*,'(12F7.3)') refbasis(i,:)
  end do

  print *, '================================================'
#ifdef PROFILING
  nColours = min(20, nColours)
#endif

  CALL ResetTimer( Caller//'BulkAssembly' )
#ifndef NOGPU
  DO col=1,nColours
    active = size(elem_lists(col) % elements, 1)
    !$omp target enter data map(to: elem_lists(col) % elements(1:active))
  end do
  !$omp target enter data map(to: refbasis(1:ngp,1:nd), refdBasisdx(1:ngp,1:nd,1:3))
  !$omp target enter data map(to: color_timing(1:nColours))
#endif

  DO col=1,nColours

    active = size(elem_lists(col) % elements, 1)

    !CALL Info( Caller,'Assembly of colour: '//I2S(col),Level=1) ! TODO: this goes away

#ifdef NOGPU
    call ResetTimer( Caller//'ColorLoop')
#endif

#ifndef NOGPU
    color_timing(col) = omp_get_wtime() 
#endif

    !$omp target teams distribute parallel do
    DO t=1,Active
      associate( &
        n => elem_lists(col) % elements(t) %n, &
        nd => elem_lists(col) % elements(t) % nd, &
        nb => elem_lists(col) % elements(t) % nb)

      CALL ModuleLocalMatrixVecSO(  n, nd+nb, nb, VecAsm, &
                                    elem_lists(col) % x(1:n, t), &
                                    elem_lists(col) % y(1:n, t), &
                                    elem_lists(col) % z(1:n, t), &
                                    coorddim, &
                                    refbasis, refdBasisdx, refip, ngp)
      ! CALL ModuleLocalMatrixVecSO(  n, nd+nb, nb, VecAsm, &
      !                               elem_lists(col)% elements(t) % nodes, &
      !                               coorddim, &
      !                               refbasis, refdBasisdx, refip, ngp)
      end associate
    END DO
    !$omp end target teams distribute parallel do

#ifdef NOGPU
    write (*, '(A, I4, I9, A)', advance='no') 'Color, #elems (', col, active,') '
    CALL CheckTimer(Caller//'ColorLoop',Delete=.TRUE.)
#endif

#ifndef NOGPU
    color_timing(col) = omp_get_wtime() - color_timing(col)
#endif
  END DO

#ifndef NOGPU
    !$omp target exit data map(from: color_timing(1:nColours)
#endif

#ifndef NOGPU
do col = 1,nColours
   write (*, '(A, I4, A, F8.6, I9, E12.3)') 'Color ', col, ' time, #elems, quotient: ', &
     color_timing(col), size(elem_lists(col) % elements, 1), color_timing(col)/size(elem_lists(col) % elements, 1)
 end do
#endif



    CALL CheckTimer(Caller//'BulkAssembly',Delete=.TRUE.)
  stop
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

