#!/bin/bash
################################
## SLURM batchjob script for
## Elmer on LUMI
##
## copyleft 2023-06-21
##    CSC-IT Center for Sciencce
##
################################
 
#SBATCH --time=00:15:00
#SBATCH --job-name=ChEESE_test
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=dev-g

####### change to your project #######
#SBATCH --account=project_462000007

####### change to numbers of nodes and MPI tasks ###
####### NB: we provide meshes for 128,256,512 and 1024 partitions #####
#######     do the math by matching the product of next entries   #####
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --gpus-per-node=1
#SBATCH --gpus-per-task=1
#SBATCH --mem=250G
#SBATCH --exclusive

################## OpenMP Stuff ##########
## use only if you undersubscribe
## the MPI tasks
##########################################
#SBATCH --cpus-per-task=56
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "running OpenMP on $SLURM_CPUS_PER_TASK"
#export KMP_AFFINITY=compact
#export KMP_DETERMINISTIC_REDUCTION=yes

#export CRAY_ACC_NO_ASYNC=1

## These control USM behaviour
export CRAY_ACC_USE_UNIFIED_MEM=0
export HSA_XNACK=0
export CRAY_ACC_DEBUG=2
export LIBOMPTARGET_KERNEL_TRACE=2


###### enable CSC provided modules #########
#ml use /appl/local/csc/modulefiles
#module use ~/.modulefiles

#ml LUMI/23.09
#ml PrgEnv-cray
#ml craype-accel-amd-gfx90a
#ml rocm
#ml cray-libsci
#ml cray-hdf5/1.12.2.7
#ml cray-netcdf/4.9.0.7
#source ../cpe/setenv.sh
#ml myelmer/offload
ml 
export LIBOMPTARGET_KERNEL_TRACE=2
export LIBOMPTARGET_INFO=-1

# this loads the spack-PrgEnv-gnu cray-libsci (BLAS, LAPACK) version
###### best to use this for audits! #########
#ml use /appl/lumi/spack/23.03/0.19.2/share/spack/modules/linux-sles15-zen2
#ml use /appl/local/csc/soft/eng/elmer/spack/23.03/0.19.2/modules/tcl/linux-sles15-zen2/
###### this loads the container version of Elmer
# module load elmer/latest
###### this loads the PregEnv/gnu version ###
###### (using cray-libsci (BLAS, LAPACK)
#module load elmer/gcc-cray
###### make it so! ######### 
#
srun ElmerSolver case.sif
