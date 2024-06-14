module LocalMV

CONTAINS
!#ifdef BUILDLOCALMV

! This subroutine will accumulate "stiffs" array from integration weights in parts
SUBROUTINE AccumulateStiff( n, nd, x, y, z, dim, refbasis, refdBasisdx, & 
    ip, elem, stiffs, forces)
!------------------------------------------------------------------------------
    !USE LinearForms
    use types, only: dp
    !use iso_c_binding
#ifdef ASSEMBLE
    use DefUtils ! TODO: defaultupdateequations is here but defutils may not be used due to threadprivate module variables if
                 ! declare target is set
#endif
    IMPLICIT NONE
!$omp declare target

    INTEGER, INTENT(IN) :: n, nd
    real(kind=dp), intent(in) :: x(:,:), y(:,:), z(:,:)
    INTEGER, intent(in) :: dim, elem
    real(kind=dp), intent(in) :: refbasis(:), refdbasisdx(:,:), ip
    real(kind=dp), intent(inout) :: stiffs(:,:), forces(:,:)
!------------------------------------------------------------------------------
#define ngp_ 4
#define nd_ 4
    REAL(KIND=dp) :: Basis(nd)
    real(kind=dp) :: dBasisdx(nd,3)
    !REAL(KIND=dp) :: MASS(nd,nd), 
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd)
    REAL(KIND=dp) :: DiffCoeff, SourceCoeff
    REAL(KIND=dp) :: LtoGMap(3,3), detg
    INTEGER :: j,k,m,i,l,t,p,q,allocstat
    integer :: colind
    
#ifdef DEBUGPRINT
    INTEGER :: round = 1 ! TODO
#endif

!------------------------------------------------------------------------------

    
    !dim = CoordinateSystemDimension()


!#ifdef DOIT

    dbasisdx(:,:) = 0_dp
#if 1
    !LtoGMap(:,:) = 1_dp
    !detg = 1.0_dp
      
      ! TODO: testing effect of commenting this out
      call myElementMetric(nd,n,x,y,z, dim, DetG, refdbasisdx(:,:), LtoGMap, elem)

      ! detj(l) = detg ! When using array of integration points

      do k=1,dim
        do j=1,dim
          do m = 1,nd
            dbasisdx(m,j) = dbasisdx(m,j) + refdbasisdx(m,k)*LtoGMap(j,k)
          end do 
        end do
      end do
#endif
    ! end do



    ! PART ONE: Collect local stiffness to STIFF using experimental method

    ! DO t=1,ngp
    !   DetJ(t) = IP % s(t) * Detj(t)
    ! END DO


    !MASS  = 0._dp
    STIFF = 0._dp
    DiffCoeff = 1._dp ! TODO: Material parameters must be communicated somehow to the solver
#if 0
    do k = 1, dim
      do j= 1, nd
        do i = 1, nd
          !do l = 1,ngp
          !   stiff(i,j) = stiff(i,j) + dbasisdx(l,i,k)*dbasisdx(l,j,k)*diffcoeff(l)*detJ(l)*ip
          !end do

          stiff(i,j) = stiff(i,j) + dbasisdx(i,k)*dbasisdx(j,k)

          ! stiffs(elem,(i-1)*nd+j) = stiff(i,j)
          !stiffs(elem,(i-1)*nd+j) = stiffs(elem,(i-1)*nd+j) + dbasisdx(i,k)*dbasisdx(j,k)*diffcoeff*detg*ip
        end do
      end do
    end do
#endif

#if 0
    diffcoeff = diffcoeff*detg*ip
    do j = 1, nd
      k = (j-1)*nd
      do i = 1, nd
        stiffs(elem,k+i) = stiffs(elem, k+i) + diffcoeff*stiff(i,j)
      end do 
    end do
#endif


#if 1
FORCE = 1._dp
sourcecoeff = 1._dp
    do i = 1, nd
      !do l = 1, ngp
        ! force(i) = force(i) + refbasis(l)*sourcecoeff(l)*detJ(l)*ip
      !end do
      ! forces(elem,i) = force(i) ! TODO: add forces
      forces(elem,i) = forces(elem,i) + refbasis(i)*sourcecoeff*detg*ip + dbasisdx(i,1)
    end do
#endif

!------------------------------------------------------------------------------
END SUBROUTINE AccumulateStiff
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE AccumulateForce(forces, elem, nd, refbasis, sourcecoeff, detg, ip)
!------------------------------------------------------------------------------
  use types, only: dp
  IMPLICIT NONE 
  !$omp declare target
!------------------------------------------------------------------------------
  REAL(kind=dp), INTENT(OUT) :: forces(:,:)
  INTEGER, INTENT(IN) :: elem, nd
  REAL(kind=dp), INTENT(IN) :: refbasis(nd), sourcecoeff, detg, ip
!------------------------------------------------------------------------------
  INTEGER :: i

  do i = 1, nd
      forces(elem,i) = forces(elem,i) + refbasis(i)*sourcecoeff*detg*ip
  end do

END SUBROUTINE AccumulateForce

!------------------------------------------------------------------------------
SUBROUTINE get_crs_inds(val_inds, rows, cols, l2g, nd, elem)
    IMPLICIT NONE
!$omp declare target
!------------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: nd, rows(:), cols(:), l2g(:,:), elem
    INTEGER, INTENT(OUT) :: val_inds(:,:)
!------------------------------------------------------------------------------
    INTEGER :: colind, i, j, k

    do i = 1,nd
      do j = 1,nd
        colind = 0
        do k = rows(l2g(elem,i)), (rows(l2g(elem,i)+1)-1)
          colind = colind + merge(k, 0, cols(k) == l2g(elem,j))
        end do
        val_inds(elem,(i-1)*nd+j) = colind
      end do
    end do

END SUBROUTINE get_crs_inds

SUBROUTINE ModuleLocalMatrixVecSO( n, nd, nb, x, y, z, dim, refbasis, refdBasisdx, & 
    ip, ngp, elem, l2g, values, cols, rows, rhs, stiffs, forces, val_inds)
!------------------------------------------------------------------------------
    !USE LinearForms
    use types, only: dp
    USE Integration, only: GaussIntegrationPoints_t
    !use iso_c_binding
#ifdef ASSEMBLE
    use DefUtils ! TODO: defaultupdateequations is here but defutils may not be used due to threadprivate module variables if
                 ! declare target is set
#endif
    IMPLICIT NONE
!$omp declare target


    INTEGER, INTENT(IN) :: n, nd, nb, l2g(:,:)
    real(kind=dp), intent(in) :: x(:,:), y(:,:), z(:,:)
    INTEGER, intent(in) :: dim, elem
    real(kind=dp), intent(in) :: refbasis(:,:), refdbasisdx(:,:,:)
    TYPE(GaussIntegrationPoints_t), intent(in) :: IP
    integer, intent(in) :: cols(:), rows(:)
    integer, intent(out) :: val_inds(:,:)
    real(kind=dp), intent(inout) :: values(:), rhs(:)
    real(kind=dp), intent(inout) :: stiffs(:,:), forces(:,:)
!------------------------------------------------------------------------------
#define ngp_ 4
#define nd_ 4
    REAL(KIND=dp) :: Basis(ngp,nd)
    real(kind=dp) :: dBasisdx(ngp,nd,3), DetJ(ngp), dbasisdx_i(nd,3)
    !REAL(KIND=dp) :: MASS(nd,nd), 
    REAL(KIND=dp) :: STIFF(nd,nd), FORCE(nd)
    REAL(KIND=dp) :: DiffCoeff(ngp), SourceCoeff(ngp)
    REAL(KIND=dp) :: LtoGMap(3,3), detg
    INTEGER :: j,k,m,i,l,t,p,q,ngp,allocstat
    INTEGER :: l_to_val_ind(nd, nd)
    integer :: colind
    
#ifdef DEBUGPRINT
    INTEGER :: round = 1 ! TODO
#endif

    real(KIND=dp) :: dLBasisdx(nd, 3)
!------------------------------------------------------------------------------

    
    !dim = CoordinateSystemDimension()
    !IP = GaussPoints( Element )
    !ngp = IP % n


!#ifdef DOIT

    dbasisdx(:,:,:) = 0_dp

    !LtoGMap(:,:) = 1_dp
    !detg = 1.0_dp
    do l=1,ngp
      ! do j = 1,dim
      !   do m = 1,nd
      !     dLBasisdx(m,j) = refdBasisdx(m,j,l)
      !   end do
      ! end do
      
      call myElementMetric(nd,n,x,y,z, dim, DetG, refdbasisdx(:,:,l), LtoGMap, elem)
      ! call myElementMetric(nd,n,x,y,z, dim, Metric, DetG, dLBasisdx, LtoGMap, elem)

      detj(l) = detg

      do m = 1,nd
        do j=1,dim
          do k=1,dim
            dbasisdx(l,m,j) = dbasisdx(l,m,j) + refdbasisdx(m,k,l)*LtoGMap(j,k)
            ! dbasisdx(l,m,j) = dbasisdx(l,m,j) + dLBasisdx(m,k)*LtoGMap(j,k)
          end do 
        end do
      end do
    end do



    ! PART ONE: Collect local stiffness to STIFF using experimental method

    ! DO t=1,ngp
    !   DetJ(t) = IP % s(t) * Detj(t)
    ! END DO


    !MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    DiffCoeff = 1._dp ! TODO: Material parameters must be communicated somehow to the solver
    sourcecoeff = 1._dp
#if 1
    do j=1,nd
      do i = 1,nd
        do k = 1,dim
          do l = 1,ngp
            stiff(i,j) = stiff(i,j) + dbasisdx(l,i,k)*dbasisdx(l,j,k)*diffcoeff(l)*detJ(l)*ip%s(l)
          end do
          stiffs(elem,(i-1)*nd+j) = stiff(i,j)
        end do
      end do
    end do

    do i = 1,nd
      do l = 1, ngp
        force(i) = force(i) + refbasis(l,i)*sourcecoeff(l)*detJ(l)*ip%s(l)
      end do
      forces(elem,i) = force(i) ! TODO: add forces
    end do
#endif

#if 0
    CALL LinearForms_GradUdotGradU(ngp, nd, dim, dBasisdx, DetJ, STIFF, DiffCoeff )
    CALL LinearForms_UdotF(ngp, nd, refBasis, DetJ, SourceCoeff, FORCE)
#endif 

#ifdef DEBUGPRINT

  if (round < 3) then
    print *, 'stiff' , round
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
    CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element, VecAssembly=.true.)
#endif

#ifdef DEBUGPRINT
    if (round<3) print *, l2g(elem, :)
#endif
    do i = 1,nd
      do j = 1,nd
        colind = 0
        do k = rows(l2g(elem,i)), (rows(l2g(elem,i)+1)-1)
          colind = colind + merge(k, 0, cols(k) == l2g(elem,j))
        end do

#ifdef DEBUGPRINT
        if (round < 3) then 
          print *, colind, stiff(i,j), l2g(elem,i), l2g(elem,j)
        end if
#endif
        !values(colind) = values(colind) + stiff(i,j)
        val_inds(elem,(i-1)*nd+j) = colind
      end do
      !rhs(l2g(elem,i)) = rhs(l2g(elem,i)) + force(i)
    end do

!------------------------------------------------------------------------------
END SUBROUTINE ModuleLocalMatrixVecSO


SUBROUTINE loop_over_active_2(active, elemdofs, max_nd, x, y, z, dim, refbasis, refdBasisdx, &
    refip, ngp, l2g, values, cols, rows, rhs)

  use Types, only: dp
  USE Integration, only: GaussIntegrationPoints_t
  USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT 

  implicit none

  integer, intent(in) :: active, elemdofs(:,:), ngp, dim, &
    l2g(:,:), cols(:), rows(:), max_nd

  real(kind=dp), intent(in) :: x(:,:), y(:,:), z(:,:), &
    refBasis(:,:), refdBasisdx(:,:,:)

  real(kind=dp), intent(inout) :: values(:), rhs(:)

  type(GaussIntegrationPoints_t), intent(in) :: refip

  
  real(kind=dp) :: stiffs(active, max_nd*max_nd), forces(active,max_nd)
  integer :: val_inds(active, max_nd*max_nd)
  integer :: elem, n, nd, nb, i, j, ndsq

  
  write(ERROR_UNIT,'(A)') '=== TARGET DEBUG START ==='

  ndsq = max_nd*max_nd

  !$omp target data map(from: stiffs(1:active,1:ndsq), val_inds(1:active,1:ndsq), forces(1:active,1:max_nd))
  !$omp target 
  !$omp teams distribute parallel do
    do elem=1, active
      stiffs(elem,:) = 0_dp
      forces(elem,:) = 0_dp
    end do
  !$omp end teams distribute parallel do
  !$omp end target

  write(ERROR_UNIT,'(A)') '=== TARGET INIT 1 ==='

  !$omp target 
  !$omp teams distribute parallel do
    do elem=1, active
      n=elemdofs(elem,1)
      nd=elemdofs(elem,2)
      nb=elemdofs(elem,3)
      call get_crs_inds(val_inds, rows, cols, l2g, nd, elem)
    end do
  !$omp end teams distribute parallel do
  !$omp end target

  write(ERROR_UNIT,'(A)') '=== TARGET INIT 2 ==='

  !$omp target
  !$omp teams distribute parallel do 
  do elem=1, active
    n=elemdofs(elem,1)
    nd=elemdofs(elem,2)
    !nb=elemdofs(elem,3)
      call AccumulateStiff(n, nd, &
      x, &
      y, &
      z, &
      dim, &
      refbasis(:,1), refdBasisdx(:,:,1), refip%s(1), elem, &
      stiffs, forces) ! Here be stiffs and inds
  end do
  !$omp end teams distribute parallel do
  !$omp end target

  !$omp end target data
  write(ERROR_UNIT,'(A)') '=== TARGET DEBUG END ==='

#if 0
  !No data races due to coloring
  nd = max_nd ! TODO: masking!
    !nd=elemdofs(elem,2)
    do i=1, nd
      do j = 1, nd
!$omp parallel do
        do elem=1, active
          associate(colind => val_inds(elem, (i-1)*nd+j) )
            values(colind) = values(colind) + stiffs(elem, (i-1)*nd+j)
          end associate
        end do
!$omp end parallel do
      end do

      !$omp parallel do
      do elem=1, active
        rhs(l2g(elem, i)) = rhs(l2g(elem, i)) + forces(elem, i)
      end do
      !$omp end parallel do
    end do
#endif

END SUBROUTINE loop_over_active_2

    
SUBROUTINE loop_over_active(active, elemdofs, max_nd, x, y, z, dim, refbasis, refdBasisdx, &
    refip, ngp, l2g, values, cols, rows, rhs)

  use Types, only: dp
  USE Integration, only: GaussIntegrationPoints_t
  USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT 

  implicit none

  integer, intent(in) :: active, elemdofs(:,:), ngp, dim, &
    l2g(:,:), cols(:), rows(:), max_nd

  real(kind=dp), intent(in) :: x(:,:), y(:,:), z(:,:), &
    refBasis(:,:), refdBasisdx(:,:,:)

  real(kind=dp), intent(inout) :: values(:), rhs(:)

  type(GaussIntegrationPoints_t), intent(in) :: refip

  
  real(kind=dp) :: stiffs(active, max_nd*max_nd), forces(active,max_nd)
  integer :: val_inds(active, max_nd*max_nd)
  integer :: elem, n, nd, nb, i, j 

  
  write(ERROR_UNIT,'(A)') '=== TARGET DEBUG START ==='
  !$omp target data map(from: stiffs(:,:), val_inds(:,:))
  !$omp target 
  !$omp teams distribute parallel do 
  do elem=1, active
    n=elemdofs(elem,1)
    nd=elemdofs(elem,2)
    nb=elemdofs(elem,3)
    call ModuleLocalMatrixVecSO(n, nd+nb, nb, &
                                x, &
                                y, &
                                z, &
                                dim, &
                                refbasis, refdBasisdx, refip, ngp, elem, &
                                l2g, values, cols, rows, rhs, stiffs, forces, val_inds) ! Here be stiffs and inds
  end do
  !$omp end teams distribute parallel do
  !$omp end target
  !$omp end target data
  write(ERROR_UNIT,'(A)') '=== TARGET DEBUG END ==='

  !No data races due to coloring
  nd = max_nd ! TODO: masking!
    !nd=elemdofs(elem,2)
    do i=1, nd
      do j = 1, nd
!$omp parallel do
        do elem=1, active
          associate(colind => val_inds(elem, (i-1)*nd+j) )
            values(colind) = values(colind) + stiffs(elem, (i-1)*nd+j)
          end associate
        end do
!$omp end parallel do
      end do

      !$omp parallel do
      do elem=1, active
        rhs(l2g(elem, i)) = rhs(l2g(elem, i)) + forces(elem, i)
      end do
      !$omp end parallel do
    end do

END SUBROUTINE loop_over_active

!#endif  ! BUILDLOCALMV

subroutine myElementMetric(nDOFs,nnodes,xx,yy,zz,dim,DetG,dLBasisdx,LtoGMap, elem)

use Types, only: dp ! , nodes_t
!use DefUtils
implicit none
!$omp declare target
!------------------------------------------------------------------------------
INTEGER, intent(in) :: nDOFs, dim, nnodes, elem   !< Number of active nodes in element, dimension of space, 
! TYPE(Nodes_t), intent(in)  :: Nodes             !< Element nodal coordinates
real(kind=dp), intent(in)  :: xx(:,:),yy(:,:),zz(:,:)       !< Element nodal coordinates
REAL(KIND=dp) :: Metric(3,3)         !< Contravariant metric tensor
REAL(KIND=dp), intent(in)  :: dLBasisdx(nDOFs,3)  !< Derivatives of element basis function with respect to local coordinates
REAL(KIND=dp), intent(out) :: DetG                !< SQRT of determinant of metric tensor
REAL(KIND=dp), intent(out) :: LtoGMap(3,3)        !< Transformation to obtain the referential description of the spatial gradient
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
REAL(KIND=dp) :: dx(3,3),G(3,3),GI(3,3), s
!REAL(KIND=dp), DIMENSION(:), POINTER :: x,y,z
INTEGER :: GeomId     
INTEGER :: cdim,i,j,k,n,imin,jmin
!------------------------------------------------------------------------------
associate(x=>xx(elem, :), y=>yy(elem, :), z=>zz(elem, :))

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
end associate
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

#if 1
  DetG = G(1,1) * ( G(2,2)*G(3,3) - G(2,3)*G(3,2) ) + &
    G(1,2) * ( G(2,3)*G(3,1) - G(2,1)*G(3,3) ) + &
    G(1,3) * ( G(2,1)*G(3,2) - G(2,2)*G(3,1) )


  CALL InvertMatrix3x3( G,GI,detG )
  Metric = GI
  DetG = SQRT(DetG)
#endif
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

#if 1
DO i=1,cdim
  DO j=1,dim
    s = 0.0d0
    DO k=1,dim
      s = s + dx(i,k) * Metric(k,j)
    END DO
    LtoGMap(i,j) = s
  END DO
END DO
#endif 

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

    !DO ii=1,m,VECTOR_BLOCK_LENGTH
    !iin=MIN(ii+VECTOR_BLOCK_LENGTH-1,m)
    !blklen=iin-ii+1

      DO j=1,n
        DO i=1,n
          DO k=1,dim
            DO l=i,m
              G(i,j) = G(i,j) + GradU(l,i,k)*GradU(l,j,k)*weight(l)*alpha(l)
            END DO
          END DO
        END DO
      END DO
        
    !END DO ! Vector blocks
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
  USE Integration, only: GaussIntegrationPoints_t, GaussPoints
  USE DefUtils
  USE Types
  use LocalMV
  use omp_lib
  USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT 
  !use opena

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
  INTEGER :: iter, maxiter, nColours, col, totelem, nthr, state, MaxNumNodes, MinNumNodes, MaxNumDOFs
  INTEGER :: ngp, i, dim, coorddim, k
  LOGICAL :: Found, VecAsm, InitHandles
  integer, allocatable :: n_active_in_col(:)
  real(KIND=dp), allocatable :: color_timing(:)
  type(nodes_t) :: Nodes
  type(Matrix_t), POINTER :: global_stiff
  integer, allocatable :: indexes(:)

#ifdef _OPENMP
  LOGICAL :: initial_device
#endif
  REAL(KIND=dp), ALLOCATABLE :: refBasis(:,:), refdBasisdx(:,:,:), prefBasis(:,:), prefdBasisdx(:,:,:), trefbasis(:,:)
  TYPE(GaussIntegrationPoints_t) :: refIP

  type :: elem_ptr
    !type(Element_t), pointer :: p
    !type(Nodes_t) :: nodes
    integer :: n, nd, nb                        ! nof nodes, nof dofs, nofbdofs
  end type

  type :: elem_list_t
    integer, allocatable :: elemdofs(:,:)
    real(kind=dp), allocatable :: x(:,:), y(:,:), z(:,:)
    integer, allocatable :: l2g(:,:)
    !real(kind=dp), allocatable, target :: x(:,:), y(:,:), z(:,:)
  end type

  TYPE(elem_list_t), allocatable :: elem_lists(:)

  CHARACTER(*), PARAMETER :: Caller = 'AdvDiffSolver'

  INTEGER  :: SharedComm = -1, ierr


  CALL Info(Caller,'------------------------------------------------')
  CALL Info(Caller,'Solving generic advection-diffusion-reaction PDE')

  CALL DefaultStart()
  
  
  maxiter = ListGetInteger( GetSolverParams(),&
      'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  global_stiff => Solver % Matrix

  if ( global_stiff % format /= MATRIX_CRS) then
    call fatal('AdvDiffSolver', 'stiffness matrix is not CRS matrix')
  end if

#ifndef NOGPU
  !$omp target enter data map(to:global_stiff % values(:), global_stiff % rows(:), global_stiff % cols(:), global_stiff % rhs(:))
#endif


  !$ nthr = omp_get_max_threads()

  ! Nonlinear iteration loop:
  !--------------------------
  ! DO iter=1,maxiter

  ! System assembly:
  !----------------
  CALL DefaultInitialize()

  totelem = 0


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
    allocate(elem_lists(col) % elemdofs(active,3)) ! n, nd, nb per column
  end do

  !!$omp target enter data map(to:solver%mesh%elemdofs)
  !!$omp target enter data map(to:elem_lists(1:nColours))

!! Tabulate elemdofs and their ndofs/nnodes/nb 
  nColours = GetNOFColours(Solver)
  !!$OMP PARALLEL &
  !!$OMP SHARED(Active, Solver, nColours, VecAsm, elem_lists) &
  !!$OMP PRIVATE(t, Element, n, nd, nb, col, InitHandles) & 
  !!$OMP REDUCTION(max:MaxNumNodes) DEFAULT(none)
  do col=1,ncolours
    !!$OMP SINGLE
    !!$omp target enter data map(to:elem_lists(col)) nowait
    active = GetNOFActive(Solver)
  MaxNumNodes = 0
  MinNumNodes = 10000
  MaxNumDOFs = 0
    !!$OMP END SINGLE

    !!$OMP DO
    do state=1,2
      do t=1,active
        Element => GetActiveElement(t)
        CALL GetElementNodesVec( nodes,  Element)
        if (state == 1) then
          !elem_lists(col) % elemdofs(t) % p => Element
          ! elem_lists(col) % elemdofs(t) % n = GetElementNOFNodes(Element)
          ! elem_lists(col) % elemdofs(t) % nd = GetElementNOFDOFs(Element)
          ! elem_lists(col) % elemdofs(t) % nb = GetElementNOFBDOFs(Element)

          elem_lists(col) % elemdofs(t,1) = GetElementNOFNodes(Element)
          elem_lists(col) % elemdofs(t,2) = GetElementNOFDOFs(Element)
          elem_lists(col) % elemdofs(t,3) = GetElementNOFBDOFs(Element)

          ! CALL GetElementNodesVec( elem_lists(col) % elemdofs(t) % nodes,  elem_lists(col) % elemdofs(t) % p)

          MaxNumNodes = max(MaxNumNodes, elem_lists(col) % elemdofs(t,1))
          MinNumNodes = min(MinNumNodes, elem_lists(col) % elemdofs(t,1))
          MaxNumDOFs = max(MaxNumDOFs, elem_lists(col) % elemdofs(t, 2))
        else

          associate (n => elem_lists(col) % elemdofs(t,1), nd => elem_lists(col) % elemdofs(t,2))
            elem_lists(col) % x(t,1:n) = nodes % x(1:n)
            elem_lists(col) % y(t,1:n) = nodes % y(1:n)
            elem_lists(col) % z(t,1:n) = nodes % z(1:n)
            nd = GetElementDOFs(indexes, Element, Solver)
            !nd = GetElementDOFs(elem_lists(col) % l2g(t,:), Element, Solver)
            elem_lists(col)%l2g(t,1:nd) = solver % variable % perm(indexes)
          end associate
        end if
      end do

      if (state == 1) then
        allocate( elem_lists(col) % x(active, MaxNumNodes), &
                  elem_lists(col) % y(active, MaxNumNodes), &
                  elem_lists(col) % z(active, MaxNumNodes))
        if(allocated(indexes) .and. size(indexes,1) < MaxNumDOFs) then
          deallocate(indexes)
          allocate(indexes(MaxNumDOFs))
        end if
        
        if (.not. allocated(indexes)) allocate(indexes(MaxNumDOFs))

        allocate(elem_lists(col) % l2g(active, MaxNumDOFs))
      end if
    end do

  !!$OMP END DO
  end do
  !!$OMP END PARALLEL


  write (*,'(A, I3, I3)') 'Max/Min NumNodes:', MaxNumNodes, MinNumNodes

  nColours = GetNOFColours(Solver)

  Element => GetActiveElement(1)
  print *, '==REFERENCE BASIS FUNCTION VALUES==============='

  !Element => elem_lists(1) % elements(1) % p
  ! TODO: here we assume that all elements are the same
  dim = CoordinateSystemDimension()
  refIP = GaussPoints( Element )
  ngp = refIP % n
  nd = GetElementNOFDOFs(Element)
  print *, 'NGP:', ngp

  allocate(refbasis(ngp, nd), refdbasisdx(ngp, nd, 3))
  allocate(prefdBasisdx(nd,3,ngp), trefbasis(nd, ngp))

  refbasis = 0_dp
  trefbasis = 0_dp
  refdBasisdx = 0_dp
  do i=1,ngp
    call NodalBasisFunctions(nd, refBasis(i,:), element, refIP%u(i), refIP%v(i), refIP%w(i))
    call NodalFirstDerivatives(nd, refdBasisdx(i,:,:), element, refIP%u(i), refIP%v(i), refIP%w(i))
    prefdbasisdx(:,:,i) = refdBasisdx(i,:,:)
    write (*,'(12F7.3)') refbasis(i,1:nd)
  end do
  do i=1,ngp
    trefbasis(:,i) = refbasis(i,:)
  end do

  print *, '================================================'
#ifdef PROFILING
  nColours = min(20, nColours)
#endif

  CALL ResetTimer( Caller//'BulkAssembly' )
#ifndef NOGPU

!$omp target enter data map(to:elem_lists, elem_lists(1:ncolours), refbasis, prefdBasisdx, trefbasis(:,:))
!$omp target enter data map(to: refbasis(1:ngp,1:nd), prefdBasisdx(1:nd,1:3,1:ngp)) 
!$omp target enter data map(to: refip%u(1:ngp), refip%v(1:ngp), refip%w(1:ngp), refip, refip%u, refip%v, refip%w, refip%s) 
!$omp target enter data map(to: refip%s(1:ngp))

  DO col=1,nColours
    active = size(elem_lists(col) % elemdofs, 1)

!$omp target enter data map(to: elem_lists(col) % elemdofs) 
!$omp target enter data map(to: elem_lists(col) % elemdofs(:,:)) 
!$omp target enter data map(to: elem_lists(col) % x(:, :)) 
!$omp target enter data map(to: elem_lists(col) % y(:, :))
!$omp target enter data map(to: elem_lists(col) % z(:, :))
!$omp target enter data map(to: elem_lists(col) % l2g(:,:))
  end do

#endif
global_stiff % values(:) = 0_dp
  DO col=1,nColours

#ifdef NOGPU
    color_timing(col) = CPUTime()
#endif

#ifndef NOGPU
    color_timing(col) = omp_get_wtime() 
#endif

    active = size(elem_lists(col) % elemdofs, 1)
  call loop_over_active_2(active, elem_lists(col) % elemdofs, MaxNumDOFs, &
    elem_lists(col) % x, &
    elem_lists(col) % y, &
    elem_lists(col) % z, &
    dim, trefBasis, prefdBasisdx,refip, ngp, elem_lists(col) % l2g, &
    global_stiff % values, global_stiff % cols, global_stiff % rows, global_stiff % rhs)

#ifdef NOGPU
    color_timing(col) = CPUTime() - color_timing(col)
#endif

#ifndef NOGPU
    color_timing(col) = omp_get_wtime() - color_timing(col)
#endif

write (*, '(A, I4, A, F8.6, I9, E12.3)') 'Color ', col, ' time, #elems, quotient: ', &
  color_timing(col), size(elem_lists(col) % elemdofs, 1), color_timing(col)/size(elem_lists(col) % elemdofs, 1)

  END DO
! call cray_acc_set_debug_global_level(0)
#ifndef NOGPU
!do col = 1,nColours
   !write (*, '(A, I4, A, F8.6, I9, E12.3)') 'Color ', col, ' time, #elems, quotient: ', &
     !color_timing(col), size(elem_lists(col) % elemdofs, 2), color_timing(col)/size(elem_lists(col) % elemdofs, 2)
 !end do
#endif



    CALL CheckTimer(Caller//'BulkAssembly',Delete=.TRUE.)
  !stop
  totelem = 0

  CALL DefaultFinishBulkAssembly()
#if 0
  open(newunit=t, file="cols.csv")
  write(t, '(I5)') global_stiff % cols
  close(t)

  open(newunit=t, file="rows.csv")
  write(t, '(I5)') global_stiff % rows
  close(t)

  open(newunit=t, file="values.csv")
  write(t, '(F9.4)') global_stiff % values
  close(t)
#endif
   !STOP


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

