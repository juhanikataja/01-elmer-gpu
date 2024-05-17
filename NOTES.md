* use declare variant
  (https://www.openmp.org/wp-content/uploads/openmp-examples-5.2.1.pdf page
  494) for info/warn etc errors which do nothing
* Declare variant candidates:
  - Messages::Warn (and other printers)
  - DefUtils::GetCurrentElement (threading focused branching and threadprivate variable)
  - DefUtils::GetElementNodes / DefUtils::GetElementNodesVec
* Pessimistic about GPU-ability:
  - DefUtils::GetElementNOFDOFS (mostly branching, very little arithmetic computations). If this is kept on cpu then there is no point GPUing GetCurrentElement

# 2024-04-26

* Adding some `omp target enter data directives will speed execution lot
  - Without any `data` mapping directives it is SLOW
  - With the `enter data` we get 4-0.5 execution time of cpu

# 2024-05-08

* TODO: 
  - [ ] Test cpu perf on 035, 035x2 and 015 meshes and compare to gpu perf. Make note here: 
    - [ ] pretty ok when looking at profiler (~30us per color for 015 mesh)
  - [ ] Glueing to global crs matrix
  - [ ] p-elements
