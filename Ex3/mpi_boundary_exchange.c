/* Ex3 - MPI BOUNDARY EXCHANGE
 * - Parallel Programming 
 * - Prof. Ivan Girotto
 * - Exercise submitted by Keerthana C J 
 *
 * This program represents a boundary exchange problem which is predominantly
 * used for solving Finite Differences/Finite Element where the boundary data in
 * each processor has to be exchanged across all processors
 * 
 * To understand the boundary problem we construct a simple version here by
 * inserting a small block in one side of the domain and allowing it to propagate
 * to the other side of the domain in one direction.
 *
 * For this we take sqaure matrix of size n with zeros and a square block of dimension 2
 * is inserted in the middle of the matrix. The matrix is decomposed only in 1D across the
 * processors. The block has to propagate from the bottom edge to the top edge of the matrix.
 * 
 * Each process has a local matrix of size n*(nlocal + 2 boundaries) The boundaries are ghost
 * cells which are used for the purpose of exchanging the borders of the domain. 
 * The upper boundary has a copy of the upper border and the lower boundary has a copy of the
 * lower border
 * eg. 10*(5+2)
 * 0 0 0 0 0 0 0 0 0 0 --> upper boundary
 * 0 0 0 0 0 0 0 0 0 0 --> upper border
 * 0 0 0 0 0 0 0 0 0 0 -|
 * 0 0 0 0 0 0 0 0 0 0 -| ---> bulk
 * 0 0 0 0 0 0 0 0 0 0 --> lower border
 * 0 0 0 0 0 0 0 0 0 0 --> lower boundary
 *
 * When the block is inserted, the matrix becomes as
 * 0 0 0 0 0 0 0 0 0 0 --> upper boundary
 * 0 0 0 0 0 0 0 0 0 0 --> upper border
 * 0 0 0 0 0 0 0 0 0 0 -|
 * 0 0 0 0 1 1 0 0 0 0 -| ---> bulk
 * 0 0 0 0 1 1 0 0 0 0 --> lower border
 * 0 0 0 0 1 1 0 0 0 0 --> lower boundary
 *
 * As the block propagates across the bulk, we see the following steps:
 * 1. Moving one step up
 * 0 0 0 0 0 0 0 0 0 0 --> upper boundary
 * 0 0 0 0 0 0 0 0 0 0 --> upper border
 * 0 0 0 0 1 1 0 0 0 0 -|
 * 0 0 0 0 1 1 0 0 0 0 -| ---> bulk
 * 0 0 0 0 0 0 0 0 0 0 --> lower border
 * 0 0 0 0 0 0 0 0 0 0 --> lower boundary
 *
 * 2. As the block reaches the upper border
 * 0 0 0 0 1 1 0 0 0 0 --> upper boundary
 * 0 0 0 0 1 1 0 0 0 0 --> upper border
 * 0 0 0 0 1 1 0 0 0 0 -|
 * 0 0 0 0 0 0 0 0 0 0 -| ---> bulk
 * 0 0 0 0 0 0 0 0 0 0 --> lower border
 * 0 0 0 0 0 0 0 0 0 0 --> lower boundary
 *
 * 3. As the block leaves the domain and it enters the other processor
 * Processor i 							Processor i+1
 * 0 0 0 0 1 1 0 0 0 0 --> upper boundary		0 0 0 0 0 0 0 0 0 0 --> upper boundary
 * 0 0 0 0 1 1 0 0 0 0 --> upper border			0 0 0 0 0 0 0 0 0 0 --> upper border
 * 0 0 0 0 0 0 0 0 0 0 -|				0 0 0 0 0 0 0 0 0 0 -|
 * 0 0 0 0 0 0 0 0 0 0 -| ---> bulk			0 0 0 0 0 0 0 0 0 0 -| ---> bulk
 * 0 0 0 0 0 0 0 0 0 0 --> lower border			0 0 0 0 1 1 0 0 0 0 --> lower border
 * 0 0 0 0 0 0 0 0 0 0 --> lower boundary		0 0 0 0 1 1 0 0 0 0 --> lower boundary
 * 	
 * The tasks involved are
 *  1. exchange of boundaries between processors
 *  2. update of the bulk
 *  3. update of the borders
 *  
 * Here, we can see that the task 1 & 2 are completely independent of each other and they can be 
 * overlapped. Task 3 needs the exchange of the boundaries to be completed. Hence it cannot be
 * overlapped.
 *
 * The problem essentailly teaches us how to overlap the communication & calculation and the 
 * challenges involved in doing a domain decomposition.
 */ 



#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

//A function to insert the block where ever needed
void insert_block(int **matrix, int row_start, int col_start, int row_end, int col_end){
 for(int i=row_start; i<row_end; i++){
  for(int j=col_start; j<col_end; j++){
   matrix[i][j] = 1;
  }
 }
}

//A function to propagate the block 
void update(int **mat, int row_start, int row_end, int col){
 for(int i=row_end-1; i>=row_start; i--){
  for(int j=0; j<col; j++){
   mat[i][j] = mat[i-1][j];
  }
 }
}

// A function to print the matrix
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

 //if n is not given, it takes 2*no.of procs by default
 int  n = (argc <= 1)? 2*np : atoll(argv[1]);
 int n_local = n/np;

 //local matrices are allocated as nlocal+2
 bulk = (int **)calloc(n_local+2, sizeof(int*));
 bulk[0] = (int *)calloc(n*(n_local+2), sizeof(int));
 for(i=1; i<n_local+2; i++){
  bulk[i] = bulk[0] + i*n;
 }

 //for printing from each processor if the matrix size is less than 16
 if(rank == 0 && n < 16){
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
 if(rank != 0) update(bulk, 1, 2, n); 			  //down boundary

 src = (rank-1);
 des = (rank+1);
 if(rank == 0)	src = np-1;
 if(rank == np-1) des = 0;

 MPI_Datatype MPI_ARRAYROW;
 MPI_Type_contiguous(n, MPI_INT, &MPI_ARRAYROW);
 MPI_Type_commit(&MPI_ARRAYROW);

  //printing if matrix size is less than 16 for time = 0
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
  MPI_Isend(up_boundary, 1, MPI_ARRAYROW, des, 100, MPI_COMM_WORLD, &request[0]); //Exchange boundaries
  MPI_Irecv(down_boundary, 1, MPI_ARRAYROW, src, 100, MPI_COMM_WORLD, &request[1]); //Exchange boundaries
  update(bulk, 2, n_local+1, n);   			   //bulk update (update 1 -> 2 -> 3 -> .. ->nlocal)
  MPI_Waitall(2, request, status);  			   //Wait for communication to finish
  if(rank != np-1) update(bulk, n_local+1, n_local+2, n);  //upper border update (update nlocal -> nlocal+1)
  update(bulk, 1, 2, n);				   //down border update (update 0 -> 1)
		
  //printing if matrix size is less than 16 
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
