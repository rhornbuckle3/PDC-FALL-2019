#!/bin/bash

#SBATCH --job-name=addition_mpi
#SBATCH --output=res.txt
#SBATCH -e test.err

#SBATCH --ntasks-per-node 1
#SBATCH --time=10:00
#SBATCH --nodes 9
#SBATCH -p qCDER
#SBATCH--mem-per-cpu=100
#environment
module load Compilers/mvapich2_ACoRE

#srun hostname srun sleep 60pwd

#Execute
srun --mpi=pmi2 ./a.out