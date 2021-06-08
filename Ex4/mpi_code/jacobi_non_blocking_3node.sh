#!/bin/bash
#SBATCH --job-name=keerthana_jacobi
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=20
#SBATCH --partition=regular1
#SBATCH --error=jacobi.err
#SBATCH --output=jacobi.out
#SBATCH --time=04:00:00


module load intel/18.0.3.222
module load openmpi/2.1.3

cd /home/kchandra/ParallelProgramming/Jacobi/mpi_code
 
N=1200
it=10
rpeek=102
cpeek=393

for procs in 2 4 8 16 20 40 60
 do
    echo ${procs}
    (mpirun -np ${procs} ./mpi_jacobi_NB.x ${N} ${it} ${rpeek} ${cpeek})> ../output_non_blocking/N1200_3node/NB_${procs}.txt
 done
