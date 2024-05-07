clean:
	rm ChEESE_test_*.out ChEESE_test_*.err

sbatch:
	sbatch run_ChEESE.sh

sbatchprof:
	sbatch run_ChEESE_profile.sh

interactive:
	srun -t 00:15:00 -J cheese_test -p dev-g -A project_462000007 -N 1 -n 1 -G 1 bash

buildlumi:
	ELMERSOLVER_LINK_LIBRARY=-lelmersolverstatic elmerf90 -O3 Poisson.F90 -o Poisson.so -DBUILDLOCALMV 

buildprof:
	ELMERSOLVER_LINK_LIBRARY=-lelmersolverstatic elmerf90 -O3 Poisson.F90 -o Poisson.so -DBUILDLOCALMV -DPROFILING

buildwasms:
	ELMERSOLVER_LINK_LIBRARY=-lelmersolverstatic elmerf90 -O3 Poisson.F90 -o Poisson.so -DBUILDLOCALMV -DPROFILING -hkeepfiles

buildcpu: 
	ELMERSOLVER_LINK_LIBRARY=-lelmersolverstatic elmerf90 -O3 Poisson.F90 -o Poisson.so -DBUILDLOCALMV -DNOGPU -DPROFILING

tailout: 
	tail -f `ls *.out|tail -n 1`

tailerr:
	tail -f `ls *.err|tail -n 1`

lessout:
	less `ls *.out|tail -n 1`

lesserr:
	less `ls *.err|tail -n 1`
