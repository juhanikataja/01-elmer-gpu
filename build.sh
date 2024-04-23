#ELMERSOLVER_LINK_LIBRARY=-lelmersolverstatic elmerf90 minipoisson.F90 -o minipoisson.so
ELMERSOLVER_LINK_LIBRARY=-lelmersolverstatic elmerf90 Poisson.F90 -o Poisson.so -DBUILDLOCALMV -DASSEMBLE
export NVCOMPILER_ACC_NOTIFY=3
