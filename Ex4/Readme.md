# Exercise 4 - Parallelzing Jacobi Algorithm

## This exercise is to parallelize the Jacobi Algorithm to solve a 2D Laplace equation. 

We conduct a performance analysis to study how the size of the domain determines the scaling, how the communication and computation times varies as the size of the domain varies, how the blocking and non-blocking communication can have an effect on the scaling. Two different sizes of the square domain are considered 1200 and 12000. 

The serial code for the problem is present in folder `serial_code/` The mpi code is present in the `mpi_code`. Both come with a make file to do the compilation. The mpi code gets complied with two versions - `mpi_jacobi_blocking` for the version with blocking communcation and `mpi_jacobi_NB` for the non-blocking version.

The output from performing the experiment is present in the `output/` folder, it contains text files with the various times measured. The `scripts/` folder contain the scripts and the settings used to run the experiment on the Ulysses cluster at SISSA.

A post processing of the output is done in the python notebook `Plots-N1200.ipynb` for the size 1200 and `Plots-N12k.ipynb` for 12000

A detailed description of the problem and a disucssion on the results can be found in the `Report.pdf`

The codes were complied using

    * ICC 19.1.3.304
    * Open MPI 3.1.4
