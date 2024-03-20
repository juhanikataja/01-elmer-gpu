# ELMERSOLVER_LINK_LIBRARY=-lelmersolverstatic 
elmerf90 Poisson.F90 -o Poisson.so -Minfo=mp
export NVCOMPILER_ACC_NOTIFY=3
