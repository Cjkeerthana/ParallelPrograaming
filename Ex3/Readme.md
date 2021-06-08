# Exercise 3 - MPI Boundary Exchange


The compilation was done using
  * gcc 7.5.0 
  * Open MPI 2.1.1
  

Compilation command: `mpicc -O3 mpi_boundary_exchange.c -o mpi_boundary_exchange.x`

Run command:  `mpirun -np ${procs} ./mpi_boundary_exchange.x ${size_of_matrix}`   (default size: 2 * no.of procs)

This program represents a boundary exchange problem which is predominantlyused for solving Finite Differences/Finite Element where the boundary data ineach processor has to be exchanged across all processors.

To understand the boundary problem we construct a simple version here by inserting a small block in one side of the domain and allowing it to propagate to the other side of the domain in one direction.

For this we take sqaure matrix of size n with zeros and a square block of dimension 2 is inserted in the middle of the matrix. The matrix is decomposed only in 1D across the processors. The block has to propagate from the bottom edge to the top edge of the matrix.
  
Each process has a local matrix of size (nlocal + 2 boundaries) * n  The boundaries are ghost cells which are used for the purpose of exchanging the borders of the domain. The upper boundary has a copy of the upper border and the lower boundary has a copy of the lower border
  
  eg.For a domain of dimension (4+2) * 10
  ```
  0 0 0 0 0 0 0 0 0 0 --> upper boundary
  0 0 0 0 0 0 0 0 0 0 --> upper border
  0 0 0 0 0 0 0 0 0 0 -|
  0 0 0 0 0 0 0 0 0 0 -| ---> bulk
  0 0 0 0 0 0 0 0 0 0 --> lower border
  0 0 0 0 0 0 0 0 0 0 --> lower boundary
  ```
  When the block is inserted, the matrix becomes as
  ```
  0 0 0 0 0 0 0 0 0 0 --> upper boundary
  0 0 0 0 0 0 0 0 0 0 --> upper border
  0 0 0 0 0 0 0 0 0 0 -|
  0 0 0 0 1 1 0 0 0 0 -| ---> bulk
  0 0 0 0 1 1 0 0 0 0 --> lower border
  0 0 0 0 1 1 0 0 0 0 --> lower boundary
  ```
  As the block propagates across the bulk, we see the following steps:
  1. Moving one step up
  ```
  0 0 0 0 0 0 0 0 0 0 --> upper boundary
  0 0 0 0 0 0 0 0 0 0 --> upper border
  0 0 0 0 1 1 0 0 0 0 -|
  0 0 0 0 1 1 0 0 0 0 -| ---> bulk
  0 0 0 0 0 0 0 0 0 0 --> lower border
  0 0 0 0 0 0 0 0 0 0 --> lower boundary
  ```
  2. As the block reaches the upper border
  ```
  0 0 0 0 1 1 0 0 0 0 --> upper boundary
  0 0 0 0 1 1 0 0 0 0 --> upper border
  0 0 0 0 1 1 0 0 0 0 -|
  0 0 0 0 0 0 0 0 0 0 -| ---> bulk
  0 0 0 0 0 0 0 0 0 0 --> lower border
  0 0 0 0 0 0 0 0 0 0 --> lower boundary
  ```
  3. As the block leaves the domain and it enters the other processor
  ```
          Processor i                                 Processor i+1
  0 0 0 0 1 1 0 0 0 0 --> upper boundary    0 0 0 0 0 0 0 0 0 0 --> upper boundary
  0 0 0 0 1 1 0 0 0 0 --> upper border      0 0 0 0 0 0 0 0 0 0 --> upper border
  0 0 0 0 0 0 0 0 0 0 -|                    0 0 0 0 0 0 0 0 0 0 -|
  0 0 0 0 0 0 0 0 0 0 -| ---> bulk          0 0 0 0 0 0 0 0 0 0 -| ---> bulk
  0 0 0 0 0 0 0 0 0 0 --> lower border      0 0 0 0 1 1 0 0 0 0 --> lower border
  0 0 0 0 0 0 0 0 0 0 --> lower boundary    0 0 0 0 1 1 0 0 0 0 --> lower boundary
  ```
  The tasks involved are
   1. exchange of boundaries between processors
   2. update of the bulk
   3. update of the borders
  
  Here, we can see that the task 1 & 2 are completely independent of each other and they can be overlapped. Task 3 needs the exchange of the boundaries to be completed. Hence it cannot be overlapped.
 
 The problem essentailly teaches us how to overlap the communication & calculation and the challenges involved in doing a domain decomposition.
 
 
