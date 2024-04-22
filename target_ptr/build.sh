rm hello
elmerf90 subhello.F90 -o libsubhello.so
#nvfortran -lelmersolver -L. -lsubhello hello.F90 -mp=gpu -L$ELMER_HOME/lib/elmersolver -o hello
ftn -lelmersolver -L. -lsubhello hello.F90 -fopenmp -L$ELMER_HOME/lib/elmersolver -o hello

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ELMER_HOME/lib/elmersolver:. srun -t 00:05:00 -J cheese_test -p dev-g -A project_462000007 -N 1 -n 1 -G 1 ./hello
