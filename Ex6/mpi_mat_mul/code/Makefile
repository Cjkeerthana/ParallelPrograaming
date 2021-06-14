CC=mpicc

CCFLAGS= -O3

OPENMP_FLAG= -qopenmp

BLAS_FLAG= -D__DGEMM

INCLUDES= -I/opt/sissa/compiler/intel/2020.4/compilers_and_libraries_2020.4.304/linux/mkl/include

LIBRARIES_SERIAL= -L/opt/sissa/compiler/intel/2020.4/compilers_and_libraries_2020.4.304/linux/mkl/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

LIBRARIES_THREADED= -L/opt/sissa/compiler/intel/2020.4/compilers_and_libraries_2020.4.304/linux/mkl/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -qopenmp

all: mpi_mat_mul_naive mpi_mat_mul.o mpi_mat_mul_blas

mpi_mat_mul_naive:mpi_mat_mul.c
	${CC} ${CCFLAGS} mpi_mat_mul.c -o mpi_mat_mul_naive

mpi_mat_mul.o:mpi_mat_mul.c
	${CC} ${BLAS_FLAG} ${INCLUDES} -c mpi_mat_mul.c

mpi_mat_mul_blas:mpi_mat_mul.o
	${CC} ${CCFLAGS} ${BLAS_FLAG} ${LIBRARIES_SERIAL} mpi_mat_mul.o -o mpi_mat_mul_blas

clean:
	rm -r mpi_mat_mul_naive mpi_mat_mul.o mpi_mat_mul_blas
