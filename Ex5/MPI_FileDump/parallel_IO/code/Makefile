CC=mpicc

CCFLAGS= -O3

INCLUDES= -I${HDF5_INC}

LIBRARIES= -L${HDF5_LIB} -lhdf5

mpi_jacobi_phdf5:mpi_jacobi_pIO.c
	${CC} ${CCFLAGS} ${INCLUDES} ${LIBRARIES} mpi_jacobi_pIO.c -o mpi_jacobi_phdf5

clean:
	rm -r mpi_jacobi_phdf5

