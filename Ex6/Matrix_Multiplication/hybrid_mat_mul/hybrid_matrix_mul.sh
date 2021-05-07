#!/bin/bash
#SBATCH --job-name=kmat_mul_hybrid
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=20
#SBATCH --partition=regular1
#SBATCH --hint=nomultithread
#SBATCH --error=mat_mul.err
#SBATCH --output=mat_mul.out
#SBATCH --time=04:00:00


module load intel
module load openmpi3

cd /home/kchandra/ParallelProgramming/Matrix_Multiplication/hybrid_mat_mul

for N in 64 128 256 512 1024 2048 4096 8192 16384
 do 
  for procs in 2 4 6 8 16
   do
     threads=5
     tot_procs=$(($threads*$procs))
     export OMP_NUM_THREADS=${threads}
     echo ${N} ${procs} ${threads}
     (mpirun --map-by ppr:2:socket:pe=${threads} -np ${procs} ./hybrid_mat_mul_naive.x ${N})> output_naive/${tot_procs}_${N}_hybrid_naive.txt
     (mpirun --map-by ppr:2:socket:pe=${threads} -np ${procs} ./hybrid_mat_mul.x ${N})> output_mkl/${tot_procs}_${N}_hybrid_mkl.txt
   done
 done
