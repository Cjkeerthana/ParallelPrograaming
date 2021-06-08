/* Ex2 - MPI RING EXCHANGE
 * - Parallel Programming 
 * - Prof. Ivan Girotto
 * - Exercise submitted by Keerthana C J 
 *
 * This program performs a ring exchange betwee the processors
 * i.e. A variable in processor i is sent to processor i+1, and
 * processor i receives from processor i-1 in a cyclic manner.
 * In the next iteration, the variable recieved is sent to the processor 
 * i+1 and so on, until all the processors exchange the variable stored
 * in each of them with all other processors.
 * The processor sums the variable received everytime and stores
 * it in another variable. Finally, every processor has the same sum. 
 *
 * This an abstract implementation of MPI_Allgather 
 *    -- An MPI Collective Operation
 * 
 * Non-blocking send and recieve is used for the program to overlap
 * communication & calculation.
 * A conclusion by studying the program is creating a parallelism where
 * we overlap the various tasks involved.
 * The tasks for a processor i were
 *   1. sending the array to processor i+1
 *   2. receiving the array from processor i-1
 *   3. performing a sum once recieved from processor i-1
 *   4. swapping the buffers to prepare for the next send & recv operations
 *
 * Here tasks 1 & 2 can be completely overlapped as they are independent 
 * from each other. Task 3 can also be overlapped as it does a calculation
 * with the buffer which will be sent (hence it is only read rather written)
 * Task 4 cannot be overlapped as it involved exchanging of buffer which has 
 * to be sent and recieved. 
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// A function to swap the arrays. A double pointer is used to pass the array
// by reference
void swap(int **a, int **b){
 int *temp = *a;
 *a = *b;
 *b = temp;	
}

// A function to print the array
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

//Intialize the array with the rank of the processor
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

	
 for(i=0; i<np; i++){ //Time goes from 0 to no. of procs
  MPI_Isend(tmp, np, MPI_INT, des, 100, MPI_COMM_WORLD, &request[0]); //comm -- Task 1
  MPI_Irecv(buf, np, MPI_INT, src, 100, MPI_COMM_WORLD, &request[1]); //comm -- Task 2
  for(j=0; j<np; j++){						      //calc -- Task 3
   sum[j] += tmp[j];
  }
  MPI_Waitall(2, request, status); // A wait for all 3 to be completed
  swap(&tmp, &buf);						      //swap -- Task 4
 }

 printf("the array in processor pid %d is:\n", rank);
 print(sum, np);
 
 MPI_Finalize();
 return 0;
}
