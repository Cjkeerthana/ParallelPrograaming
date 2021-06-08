# Exercise 2 - MPI Ring Exchange

## The exercise is a study of MPI ring exchange which is the process perfomed in MPI_Allgather

A description of the program and the study is described inside the code.

The compilation was done using
  * gcc 7.5.0 
  * Open MPI 2.1.1

Compile command: `mpicc -O3 mpi_ring_exchange.c -o mpi_ring_exchange.x`
Run command: ` mpirun -np ${procs} ./mpi_ring_exchange.x`
