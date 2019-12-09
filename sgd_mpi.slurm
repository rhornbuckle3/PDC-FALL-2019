#!/bin/bash

#SBATCH --job-name=SGD_MPI
#SBATCH --output=res.txt
#SBATCH -e test.err

#SBATCH --ntasks-per-node 1
#SBATCH --time=1:00
#SBATCH --nodes 11
#SBATCH -p qCDER
#SBATCH--mem-per-cpu=100
#environment
module load Compilers/mvapich2_ACoRE

#srun hostname srun sleep 60pwd

#Execute
srun --mpi=pmi2 ./a.out