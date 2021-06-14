#!/bin/bash
#SBATCH --job-name=jacobiIO
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=20
#SBATCH --partition=regular1
#SBATCH --error=jacobi.err
#SBATCH --output=jacobi.out
#SBATCH --hint=nomultithread
#SBATCH --time=04:00:00


module load intel/19.1.3.304
module load openmpi3/3.1.4
module load hdf5/1.10.5

cd /home/kchandra/ParallelProgramming/Jacobi_IO/code
 
procs=320
it=10
rpeek=10
cpeek=39
write_int=10

for N in 1200 12000
 do
  echo ${N}
  (mpirun -np ${procs} ./mpi_jacobi_bin.x ${N} ${it} ${rpeek} ${cpeek} ${write_int})> ../output/${N}/bin/${procs}_bin.txt
  (mpirun -np ${procs} ./mpi_jacobi_dat.x ${N} ${it} ${rpeek} ${cpeek} ${write_int})> ../output/${N}/dat/${procs}_dat.txt
  (mpirun -np ${procs} ./mpi_jacobi_hdf5.x ${N} ${it} ${rpeek} ${cpeek} ${write_int})> ../output/${N}/hdf5/${procs}_hdf5.txt
  rm -r *.bin *.dat *.h5
done
