#!/bin/bash
#SBATCH --job-name=hjacobi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=regular1
#SBATCH --hint=nomultithread
#SBATCH --error=jacobi.err
#SBATCH --output=jacobi.out
#SBATCH --time=04:00:00


module load intel
module load openmpi3

cd /home/kchandra/ParallelProgramming/Jacobi_hybrid
 
it=10
rpeek=120
cpeek=395
write_int=10
procs=2

export OMP_NUM_THREADS=10

for N in 1200 12000
do
   echo ${N}
   (mpirun -np ${procs} ./hybrid_jacobi.x ${N} ${it} ${rpeek} ${cpeek} ${write_int})> output/${N}/${procs}_hybrid.txt
   rm -r *.bin
done
