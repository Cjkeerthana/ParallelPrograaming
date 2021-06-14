CC=mpicc

CCFLAGS= -O3

BLOCKING_FLAG= -DBLOCKING

all: mpi_jacobi_blocking mpi_jacobi_NB

mpi_jacobi_blocking:mpi_jacobi.c
        ${CC} ${CCFLAGS} ${BLOCKING_FLAG} mpi_jacobi.c -o mpi_jacobi_blocking

mpi_jacobi_NB:mpi_jacobi.c
        ${CC} ${CCFLAGS}  mpi_jacobi.c -o mpi_jacobi_NB

clean:
	rm -r mpi_jacobi_NB mpi_jacobi_blocking

