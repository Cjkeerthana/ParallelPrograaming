#!/bin/bash
#SBATCH --job-name=khybrid_1node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=regular1
#SBATCH --mem=40GB
#SBATCH --hint=nomultithread
#SBATCH --error=hyb_1node_mat_mul.err
#SBATCH --output=hyb_1node_mat_mul.out
#SBATCH --time=04:00:00


module load intel
module load openmpi3

cd /home/kchandra/ParallelProgramming/Matrix_Multiplication/hybrid_mat_mul

N=32768 
for procs in 1 2
 do
  threads=10
  tot_procs=$(($threads*$procs))
  export OMP_NUM_THREADS=${threads}
  echo ${N} ${procs} ${threads}
  #(mpirun --map-by ppr:1:socket:pe=${threads} -np ${procs} ./hybrid_mat_mul_naive.x ${N})> output_naive/${N}/${tot_procs}_${N}_hybrid_naive.txt
  (mpirun --map-by ppr:1:socket:pe=${threads} -np ${procs} ./hybrid_mat_mul.x ${N})> output_mkl/${N}/${tot_procs}_${N}_hybrid_mkl.txt
 done
