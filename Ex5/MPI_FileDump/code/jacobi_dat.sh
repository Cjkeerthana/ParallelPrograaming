#!/bin/bash
#SBATCH --job-name=keerthana_jacobi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=regular1
#SBATCH --error=jacobi.err
#SBATCH --output=jacobi.out
#SBATCH --time=04:00:00


module load intel/18.0.3.222
module load openmpi/2.1.3

cd /home/kchandra/ParallelProgramming/Jacobi_IO/code
 
N=1200
it=10000
rpeek=102
cpeek=393

for procs in 2 4 8 16 20
 do
  for write_int in 2000 4000 6000 8000 10000
  do
     echo ${procs} ${write_int}
     (mpirun -np ${procs} ./mpi_jacobi_dat.x ${N} ${it} ${rpeek} ${cpeek} ${write_int})> ../out_dat/${write_int}/${procs}_${N}_dat.txt
  done
done
