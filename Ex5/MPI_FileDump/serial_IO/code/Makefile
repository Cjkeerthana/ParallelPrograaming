CC=mpicc

CCFLAGS= -O3

LIBRARIES= -L/opt/ohpc/pub/libs/intel/hdf5/1.10.5/lib -lhdf5

INCLUDE= -I/opt/ohpc/pub/libs/intel/hdf5/1.10.5/include

DAT_FLAG= -DTEXT

HDF5_FLAG= -Dhdf5

all: mpi_jacobi_bin mpi_jacobi_dat mpi_jacobi_hdf5

mpi_jacobi_bin:mpi_jacobi_IO.c
	${CC} ${CCFLAGS} ${INCLUDE} ${LIBRARIES} mpi_jacobi_IO.c -o mpi_jacobi_bin

mpi_jacobi_dat:mpi_jacobi_IO.c
	${CC} ${CCFLAGS} ${INCLUDE} ${LIBRARIES} ${DAT_FLAG} mpi_jacobi_IO.c -o mpi_jacobi_dat

mpi_jacobi_hdf5:mpi_jacobi_IO.c
	${CC} ${CCFLAGS} ${INCLUDE} ${LIBRARIES} ${HDF5_FLAG} mpi_jacobi_IO.c -o mpi_jacobi_hdf5

clean:
	rm -r mpi_jacobi_bin mpi_jacobi_dat mpi_jacobi_hdf5
	
