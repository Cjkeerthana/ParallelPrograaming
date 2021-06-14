CC=mpicc

CCFLAGS= -O3

DAT_FLAG= -DTEXT

OPENMP_FLAG= -qopenmp

all: hybrid_jacobi_bin hybrid_jacobi_dat

hybrid_jacobi_bin:hybrid_jacobi.c
	${CC} ${CCFLAGS} hybrid_jacobi.c -o hybrid_jacobi_bin ${OPENMP_FLAG}

hybrid_jacobi_dat:hybrid_jacobi.c
	${CC} ${CCFLAGS} ${DAT_FLAG} hybrid_jacobi.c -o hybrid_jacobi_dat ${OPENMP_FLAG}

clean:
	rm -r hybrid_jacobi_bin hybrid_jacobi_dat
	
