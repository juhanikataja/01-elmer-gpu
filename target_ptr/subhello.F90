module hellomodule
!use Types
#ifdef _OPENMP
  use omp_lib
#endif

contains

  subroutine hello()
  implicit none
  logical :: initial_device
  integer :: num_devices


  !$omp target map(from: initial_device)
  call subsubhello(initial_device, num_devices)
  !$omp end target

  num_devices = omp_get_num_devices()
  print *, "Number of available devices", num_devices
  if (initial_device) then
    write(*,*) "Running on host"
  else
    write(*,*) "Running on device"
  end if

  contains

    subroutine subsubhello(idev, ndev)
  use omp_lib
    !$omp declare target to(subsubhello)
    integer, intent(out) :: ndev
    logical, intent(out) :: idev

!    ndev = omp_get_num_devices()

    idev = omp_is_initial_device()

    end subroutine subsubhello

  end subroutine

end module hellomodule
