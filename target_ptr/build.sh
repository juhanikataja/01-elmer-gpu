rm hello
elmerf90 subhello.F90 -o libsubhello.so
nvfortran -lelmersolver -L. -lsubhello hello.F90 -mp=gpu -L$ELMER_HOME/lib/elmersolver -o hello

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ELMER_HOME/lib/elmersolver:. ./hello
