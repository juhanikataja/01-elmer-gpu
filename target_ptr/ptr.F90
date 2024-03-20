program ptr
implicit none

type p_t
  integer :: x
end type

type ph_t
  type(p_t), pointer :: p
end type

type(ph_t), pointer :: a(:)
type(ph_t), target :: ta 
type(p_t), target :: tp

integer :: n

!tp%x = 1
!ta%p=>tp
!a=>ta

allocate(a(10))
do n=1,10
  allocate(a(n)%p)
end do
a(1)%p%x=2
!$omp target teams
a(1)%p%x = 1
!$omp end target teams

!$omp target teams
!$omp parallel 
print *, a(1)%p%x
!$omp end parallel 
!$omp end target teams

end
