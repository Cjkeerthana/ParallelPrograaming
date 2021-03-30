#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


void insert_block(int **matrix, int row_start, int col_start, int row_end, int col_end){
	for(int i=row_start; i<row_end; i++){
		for(int j=col_start; j<col_end; j++){
			matrix[i][j] = 1;
		}
	}
}

void update(int **mat, int row_start, int row_end, int col){
	for(int i=row_end-1; i>=row_start; i--){
		for(int j=0; j<col; j++){
			mat[i][j] = mat[i-1][j];
		}
	}
}

void print(int **mat, int row_start, int row_end, int col){
	for(int i=row_start; i<row_end; i++){
		for(int j=0; j<col; j++){
			printf("%d\t", mat[i][j]);
		}
		printf("\n");
	}
}

int main(int argc, char* argv[]){
	int **bulk, *up_boundary, *down_boundary, **print_buf;
	int rank, np,i,j,src=0, des=0;
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

	int  n = (argc <= 1)? 2*np : atoll(argv[1]);

	int n_local = n/np;

	//local matrices are allocated as nlocal+2
	bulk = (int **)calloc(n_local+2, sizeof(int*));
        bulk[0] = (int *)calloc(n*(n_local+2), sizeof(int));
        for(i=1; i<n_local+2; i++){
                bulk[i] = bulk[0] + i*n;
        }

	//for printing from each processor
	if(rank == 0){
		print_buf = (int **)malloc((n_local)*sizeof(int*));
		print_buf[0] = (int *)malloc(n*(n_local)*sizeof(int));
		for(i=1; i<n_local; i++){
			print_buf[i] = print_buf[0] + i*n;
		}
	}

	up_boundary = bulk[n_local+1]; //upper ghost layer
	down_boundary = bulk[0];       //lower ghost layer

	//upper and lower ghost cells --> indices 0, nlocal+1
	//upper and lower borders --> indices 1 , nlocal
	//bulk --> indices 2 - nlocal-1
	
	//inserting a block 2*2 block of ones in the middle of first processor bulk
	if(rank == 0){
		insert_block(bulk, 1, n/2-1, 3, n/2+1);	
	}

	//update the boundaries so that it has the border values of bulk
	if(rank != np-1) update(bulk, n_local+1, n_local+2, n);  //up boundary
        if(rank != 0) update(bulk, 1, 2, n); //down boundary*/

	src = (rank-1);
	des = (rank+1);

	if(rank == 0)	src = np-1;
	if(rank == np-1) des = 0;

	MPI_Datatype MPI_ARRAYROW;
	MPI_Type_contiguous(n, MPI_INT, &MPI_ARRAYROW);
	MPI_Type_commit(&MPI_ARRAYROW);

        //printing for time = 0
	if(n < 16){
		if(rank != 0) MPI_Send(bulk[1],n_local,MPI_ARRAYROW,0,200,MPI_COMM_WORLD);
        	if(rank == 0){
        		printf("time t = 0\n");
        		print(bulk,1,n_local+1,n);
        		for(j=1; j<np; j++){
        			MPI_Recv(print_buf[0],n_local,MPI_ARRAYROW,j,200,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        			print(print_buf,  0, n_local, n);
        		 }
        	 }
	}

	for(i=1; i<n; i++){		
		MPI_Isend(up_boundary, 1, MPI_ARRAYROW, des, 100, MPI_COMM_WORLD, &request[0]);
                MPI_Irecv(down_boundary, 1, MPI_ARRAYROW, src, 100, MPI_COMM_WORLD, &request[1]);
                update(bulk, 2, n_local+1, n);    //bulk update (update 1 -> 2 -> 3 -> .. ->nlocal)
                MPI_Waitall(2, request, status);
		if(rank != np-1) update(bulk, n_local+1, n_local+2, n);  //upper border update (update nlocal -> nlocal+1)
                update(bulk, 1, 2, n); //down border update (update 0 -> 1)
		
	       	//printing
		if(n < 16){
                	if(rank != 0) MPI_Send(bulk[1],n_local,MPI_ARRAYROW,0,200,MPI_COMM_WORLD);
                	if(rank == 0){
				printf("time t = %d\n",i);
                        	print(bulk,1,n_local+1,n);
                        	for(j=1; j<np; j++){
                                	MPI_Recv(print_buf[0],n_local,MPI_ARRAYROW,j,200,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                	print(print_buf,0,n_local,n);
                        	}
                	}
		}

	}

	MPI_Finalize();
	return 0;
}
