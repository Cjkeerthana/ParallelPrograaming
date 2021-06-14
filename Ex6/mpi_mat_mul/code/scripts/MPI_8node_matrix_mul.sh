#!/bin/bash
#SBATCH --job-name=kmpi_8node
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=20
#SBATCH --partition=regular1
#SBATCH --hint=nomultithread
#SBATCH --mem=40GB
#SBATCH --error=mpi_8node_mat_mul.err
#SBATCH --output=mpi_8node_mat_mul.out
#SBATCH --time=04:00:00


module load intel
module load openmpi3

cd /home/kchandra/ParallelProgramming/Matrix_Multiplication/mpi_mat_mul

for N in 16384
 do 
   procs=160 
   echo ${N} ${procs}
   #(mpirun -np ${procs} ./mpi_mat_mul_naive.x ${N})> output_naive/${N}/${procs}_mpi_naive.txt
   (mpirun -np ${procs} ./mpi_mat_mul.x ${N})> output_mkl/${N}/${procs}_mpi_mkl.txt
 done
