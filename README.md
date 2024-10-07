Poisson OpenMP offload code for Elmer LUMI hackathon
====================================================

## Introduction

This is a janky work in progress piece of code that has been functional on
cce/16+rocm/5.x.

## Running

- Have elmer installation set up (`PATH` variable, and `ELMER_HOME` variable possibly).
- Look at `run_ChEESE.sh` and `run_ChEESE_profile.sh`. (TODO: rename these so
  they are easier to write on keyboard).

## Compiling

There are useful some targets in `Makefile`. Most notably

``` 
buildlumi:
	elmerf90 -O3 Poisson.F90 -o Poisson.so -DBUILDLOCALMV 

buildwasms:
	elmerf90 -O3 Poisson.F90 -o Poisson.so -DBUILDLOCALMV -DPROFILING -hkeepfiles

buildcpu: 
	elmerf90 -g -O1 Poisson.F90 -o Poisson.so -DBUILDLOCALMV -DNOGPU
```

Here `NOGPU` is important for cpu builds. Defining it disables some `target
enter data` directives.



## Meshing

The mesh is generated with `gmsh`. This is available on most linux
distributions and also for windows. It is not installed in lumi.

The script `mesh` has two lines that 
1. generates the mesh with `gmsh`, and
2. converts it to elmer format with `ElmerGrid` command.
