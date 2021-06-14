#!/bin/bash
#SBATCH --job-name=kmat_mul_mpi
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=20
#SBATCH --partition=regular1
#SBATCH --error=mat_mul.err
#SBATCH --output=mat_mul.out
#SBATCH --time=04:00:00


module load intel
module load openmpi3

cd /home/kchandra/ParallelProgramming/Matrix_Multiplication/mpi_mat_mul

for N in 64 128 256 512 1024 2048 4096 8192 16384
 do 
  for procs in 2 4 8 16 20 40 60 80
   do
    echo ${N} ${procs}
    (mpirun -np ${procs} ./mpi_mat_mul_naive.x ${N})> output_naive/${procs}_${N}_mpi_naive.txt
    (mpirun -np ${procs} ./mpi_mat_mul.x ${N})> output_mkl/${procs}_${N}_mpi_mkl.txt
   done
 done
