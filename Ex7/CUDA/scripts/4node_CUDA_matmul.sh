#!/bin/bash
#SBATCH --job-name=4nodecuda
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --partition=gpu2
#SBATCH --hint=nomultithread
#SBATCH --mem=60GB
#SBATCH --error=4nodeCuda_matmul.err
#SBATCH --output=4nodeCuda_matmul.out
#SBATCH --time=04:00:00


module load intel
module load openmpi3
module load cuda

cd /home/kchandra/ParallelProgramming/CUDA_Programming/simple_prg/Matrix_Multiplication/CUDA

procs=8
nodes=4
for N in 2048 4096 8192 16384 32768
 do
  echo ${N}
  (mpirun -np ${procs} ./mat_mul ${N})> output/2_proc_per_node/${nodes}/${N}_cuda_matmul.txt
 done
