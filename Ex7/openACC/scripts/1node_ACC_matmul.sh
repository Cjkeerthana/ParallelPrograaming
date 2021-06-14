#!/bin/bash
#SBATCH --job-name=1nodeAcc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --partition=gpu2
#SBATCH --hint=nomultithread
#SBATCH --mem=60GB
#SBATCH --error=1nodeAcc_matmul.err
#SBATCH --output=1nodeAcc_matmul.out
#SBATCH --time=04:00:00


module load intel
module load openmpi3
module load pgi
module load cuda

cd /home/kchandra/ParallelProgramming/CUDA_Programming/simple_prg/Matrix_Multiplication/openACC

procs=2
nodes=1
for N in 2048 4096 8192 16384 32768
 do
  echo ${N}
  (mpirun -np ${procs} ./mat_mul_acc.x ${N})> output/${nodes}/${N}_acc_matmul.txt
 done
