* use declare variant
  (https://www.openmp.org/wp-content/uploads/openmp-examples-5.2.1.pdf page
  494) for info/warn etc errors which do nothing
* Declare variant candidates:
  - Messages::Warn (and other printers)
  - DefUtils::GetCurrentElement (threading focused branching and threadprivate variable)
  - DefUtils::GetElementNodes / DefUtils::GetElementNodesVec
* Pessimistic about GPU-ability:
  - DefUtils::GetElementNOFDOFS (mostly branching, very little arithmetic computations). If this is kept on cpu then there is no point GPUing GetCurrentElement
