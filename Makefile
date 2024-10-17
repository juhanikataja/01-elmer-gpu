clean:
	rm ChEESE_test_*.out ChEESE_test_*.err

sbatch:
	sbatch run_ChEESE.sh

sbatchprof:
	sbatch run_ChEESE_profile.sh

sbatchomniperf:
	sbatch run_ChEESE_omniperf.sh

interactive:
	srun -t 00:15:00 -J cheese_test -p dev-g -A project_462000007 -N 1 -n 1 -G 1 bash

buildlumi:
	elmerf90 -O2 Poisson.F90 -o poisson.so -DBUILDLOCALMV -hlist=a

buildprof:
	elmerf90 -O2 Poisson.F90 -o poisson.so -DBUILDLOCALMV -DPROFILING -hkeepfiles

buildwasms:
	elmerf90 -O2 Poisson.F90 -o poisson.so -DBUILDLOCALMV -DPROFILING -hkeepfiles

buildcpuprof: 
	bash -c 'module load craype-accel-host; elmerf90 -g -O2 Poisson.F90 -o poisson_cpu.so -DPROFILING -DBUILDLOCALMV -DNOGPU; module load craype-accel-amd-gfx90a'

buildcpu: 
	bash -c 'module load craype-accel-host; elmerf90 -g -O2 Poisson.F90 -o poisson_cpu.so -DBUILDLOCALMV -DNOGPU; module load craype-accel-amd-gfx90a'

tailout: 
	tail -f  `ls *.out|tail -n 1`
	# watch -n 2 tail -n 25 `ls *.out|tail -n 1`

tailerr:
	watch -n 2 tail -n 25 `ls *.err|tail -n 1`

lessout:
	less `ls *.out|tail -n 1`

lesserr:
	less `ls *.err|tail -n 1`

