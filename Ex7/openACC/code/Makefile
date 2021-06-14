PGCC=mpicc

CCFLAGS= -O3

PGI_FLAGS= -ta=tesla -Minfo=all -acc=noautopar

CUDA_FLAGS= -Mcuda -Mcudalib=cublas -lcublas

mat_mul_acc:mat_mul_acc.c
	${PGCC} ${CCFLAGS} ${PGI_FLAGS} ${CUDA_FLAGS} mat_mul_acc.c -o mat_mul_acc


clean:
	rm -r mat_mul_acc
