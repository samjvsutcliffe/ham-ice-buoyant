#!/bin/bash

# Request resources:
#SBATCH --time=20:00:0  # 6 hours (hours:minutes:seconds)
#SBATCH -p shared
#SBATCH -n 4                    # number of MPI ranks
#SBATCH --cpus-per-task=16   # number of MPI ranks per CPU socket
#SBATCH --mem-per-cpu=1G
#SBATCH -N 1-4                    # number of compute nodes. 

module load gcc
module load aocl
module load mvapich2
export MV2_ENABLE_AFFINITY=0

#sbcl --dynamic-space-size 16000  --disable-debugger --load "build_step.lisp" --quit
mpirun ./mpi-worker --dynamic-space-size 16000 --disable-debugger 
#rm mpi-worker
