#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void swap(int **a, int **b){
	int *temp = *a;
	*a = *b;
	*b = temp;	
}

void print(int *a, int size){
	for(int j=0; j<size; j++)
	{
		printf("%d\t", a[j]);
	}
	printf("\n");
}
int main(int argc, char* argv[]){
	int *sum, *buf, *tmp;
	int rank, np,i,j, src=0, des=0;
	int sum1 = 0, buf1 = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Status status[2];
	MPI_Request request[2];
	if(np <= 1){
		fprintf(stderr, "Usage: mpi -np n %s no. of iterations \n", argv[0]);
		MPI_Finalize();
		exit(-1);
	}

	sum = (int *)calloc(np, sizeof(int));
	buf = (int *)calloc(np, sizeof(int));
	tmp = (int *)calloc(np, sizeof(int));
	for(i=0; i<np; i++){
		tmp[i] = rank;
	}

	src = (rank-1);
	des = (rank+1);

	if(rank == 0)	src = np-1;
	if(rank == np-1) des = 0;

	
	for(i=0; i<np; i++){
		MPI_Isend(tmp, np, MPI_INT, des, 100, MPI_COMM_WORLD, &request[0]);
		MPI_Irecv(buf, np, MPI_INT, src, 100, MPI_COMM_WORLD, &request[1]);
		for(j=0; j<np; j++){
                        sum[j] += tmp[j];
                }
		MPI_Waitall(2, request, status);
		swap(&tmp, &buf);
	}

	printf("the array in processor pid %d is:\n", rank);
	print(sum, np);
	MPI_Finalize();
	return 0;
}
