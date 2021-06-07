Exercise 1 is about initializing a square identity matrix using MPI. 
The program mpi_IdMatInit.c is the first version developed. mpi_IdMatInit_v2.c is the improved version of the program.

The compilation was done using 
  -- gcc 7.5.0
  -- Open MPI 2.1.1
  
 Compilation command: mpicc -O3 mpi_IdMatInit.c -o mpi_IdMatInit.x
 Run command: mpirun -np ${procs} ./mpi_IdMatInit.x ${size_of_matrix}
