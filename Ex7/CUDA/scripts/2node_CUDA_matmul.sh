#!/bin/bash
#SBATCH --job-name=1nodecuda
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --partition=gpu2
#SBATCH --hint=nomultithread
#SBATCH --mem=60GB
#SBATCH --error=1nodeCuda_matmul.err
#SBATCH --output=mpi_1nodeCuda_matmul.out
#SBATCH --time=04:00:00


module load intel
module load openmpi3
module load cuda

cd /home/kchandra/ParallelProgramming/CUDA_Programming/simple_prg/Matrix_Multiplication/CUDA

procs=4
nodes=2
for N in 2048 4096 8192 16384 32768
 do
  echo ${N}
  (mpirun -np ${procs} ./mat_mul ${N})> output/${nodes}/${N}_cuda_matmul.txt
 done
