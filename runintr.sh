export CRAY_ACC_USE_UNIFIED_MEM=1
export HSA_XNACK=1
srun -t 00:05:00 -J cheese_test -p dev-g -A project_462000007 -N 1 -n 1 -G 1 ElmerSolver ./case.sif
