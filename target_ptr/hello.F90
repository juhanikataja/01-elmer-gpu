program mainhello
#ifdef _OPENMP
  use omp_lib
#endif
use hellomodule
  implicit none

  call hello()

end program
