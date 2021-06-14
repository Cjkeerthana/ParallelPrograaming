# Exercise 7 - GPU Programming

## Matrix Multiplication Using CUDA

For GPU Computing, we use CUDA and OpenACC to do a matrix multiplication. In both version we use CuBLAS library for performing the matrix multiplication.

In the CUDA version, the gathering of columns of B from various processors for multiplication takes place in the CPU processors through MPI. While, with the OpenACC, the communication takes place among GPU's directly doing a MPI between GPUs. We measure the copying time between GPU & CPU, communication time for collecting the columns of B and the calculation time in the GPU. The plots are given in the python notebooks - for CUDA are found in Plots_CUDA.ipynb and for OpenACC are found in Plots_ACC.ipynb. 


In the previous exercises, we understood the communication & computation times involved in a matrix multiplication problem using MPI & Hybrid MPI-OpenMP. A Hybrid code takes a total time of almost 0.2s on 4 nodes (using 80 processors) for a matrix size of 2048, while on GPU (8 GPUs) it takes a total time of 0.1s. Moving to bigger matrix sizes, a hybrid code took a time of 60s on 4 nodes for a matrix size of 32768 while on GPU (8 GPUs) takes almost 10s. The scaling is not linear for a GPU even if we reach the size of 32768. 

For the openACC version, we obtain a total time 2 times as high as the hybrid version for the size of 2048 and almost 14s for the size of 30000. The openACC doesn't provide a good scaling even for bigger sizes, as the communication occurs between the GPUs rather than between MPI processes in the CPU nodes. 


A GPU does accelerate the execution of the code, while the speed up in a GPU and its capital cost come as disadvtanges. Also, a GPU requires a completely vectorized code to use its full potential. Porting the code from a serial code to GPU version of it also takes effort. Before considering porting to GPU, all these details have to be taken care of and evaulated from an investment point of view.
