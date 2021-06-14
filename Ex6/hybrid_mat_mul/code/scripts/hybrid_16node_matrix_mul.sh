#!/bin/bash
#SBATCH --job-name=khybrid_16node
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=20
#SBATCH --partition=regular1
#SBATCH --hint=nomultithread
#SBATCH --mem=40GB
#SBATCH --error=hyb_16node_mat_mul.err
#SBATCH --output=hyb_16node_mat_mul.out
#SBATCH --time=04:00:00


module load intel
module load openmpi3

cd /home/kchandra/ParallelProgramming/Matrix_Multiplication/hybrid_mat_mul

for N in 16384 32768
procs=32
threads=10
tot_procs=$(($threads*$procs))
export OMP_NUM_THREADS=${threads}
echo ${N} ${procs} ${threads}
#(mpirun --map-by ppr:1:socket:pe=${threads} -np ${procs} ./hybrid_mat_mul_naive.x ${N})> output_naive/${N}/${tot_procs}_${N}_hybrid_naive.txt
(mpirun --map-by ppr:1:socket:pe=${threads} -np ${procs} ./hybrid_mat_mul.x ${N})> output_mkl/${N}/${tot_procs}_${N}_hybrid_mkl.txt
