# Exercise 6 - Matrix Multiplication 

## This exercise is to understand how to perform a matrix mulitplication parallely using MPI and Hybrid MPI-OpenMP

Matrix multiplication is one of the very common problems used for benchmarking. A parallel code for matrix multiplication is developed here. Two versions were developed - One with MPI and the other with MPI and OpenMP. We measure the times for communication and computation for sizes of matrix one small size and other big size. A small size is chosen such that the scaling doesn't occur very well. A bigger size is chose where an almost linear scaling can be achieved.

We also try to the BLAS libraries to understand how an efficient mulitplication kernel leads to lesser computation time and how it affects the scaling. 

The problem definition, experiment and results are discussed in the `Report.pdf`. The MPI code and its relevant experiments are found in `mpi_mat_mul/` folder, and the `hybrid_mat_mul` has the hybrid code and its results. The codes come with make files and the compliers used are:

    * ICC 19.1.3.304
    * Open MPI 3.1.4
 
 The `Plots/` folder contains the post processing of the results from the experiment.
