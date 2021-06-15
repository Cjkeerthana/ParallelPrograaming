# Exercise 5 - IO in Parallel Computing and Hybrid MPI - OpenMP

## The exercise in IO to understand how writing files can create an effect in the performance of the code. 

While solving a Laplace equation, to save the intermediate results, we save the iterations into files while solving by writing out the matrices. This dumping of intermedaite results can cause a significant overhead. Saving all the intermediate steps is very inefficient in terms of memory requirements and the time to write each file. Hence we figure out a write interval (e.g. if we write out every 10th iteration, then the write interval is 10) for a fixed number of iterations such that the performance is not affected much and a better scaling is achieved.

Also the type of file written out determines the dumping time for each file. Text files take enormous time and occupy a lot of space, while understanding a text file is simpler as we can see the output and write headers. A binary file occupies less memory and written way faster than text files, while headers cannot be written on the file, hence information about the matrix cannot be stored in a single file to read it and use it later on. To make compromise with faster reading writing, memory requirements and information handling, we move to HDF5 files, which take a bit more time than the binary files but a header can be written into them, but have much faster speed and lesser memory than the text files.

We also explore parallel IO in MPI to understand how to use parallelize the IO process. 

## The exercise in Hybrid MPI - OpenMP is to understand how a hybrid code can perform.

It is an attempt to make a hybrid code and analyze its performance for solving the Laplace Equation using Jacobi method.
