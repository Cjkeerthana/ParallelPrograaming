/* Ex1 - IDENTITY MATRIX INITATION USING MPI 
 * - Parallel Programming 
 * - Prof. Ivan Girotto
 * - Exercise submitted by Keerthana C J
 *
 * This is a program an improvement of mpi_IdMatInit.c
 *
 * Improvements done
 *   - Round robin allocation of rests from n/np 
 *   - Data assimilation only in the disk
 *   		 
 */

#include <mpi.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <unistd.h> 

// A function to print the matrix
void print_matrix(int **matrix, int rows, int coloumns){
 for(int i = 0; i < rows; i++){
  for(int j =0; j < coloumns; j++){
   printf("%d\t", matrix[i][j]);
  }
  printf("\n");
 }
}

int main(int argc, char* argv[]) 
{ 
 unsigned long int n, n_local, offset=0, res, i, j; 
 int **local_matrix, **matrix;
 int pid, np, master =0, tag =1;
 double start_time, end_time;  

 MPI_Status status;
 MPI_Init(&argc, &argv); 
 MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
 MPI_Comm_size(MPI_COMM_WORLD, &np);


 // the number of processors is less than 1 
 if ( np < 1) {               			
  fprintf (stderr , " Usage : mpi -np n %s number_of_iterations \n", argv[0] ) ;
  MPI_Finalize() ;
  exit(-1) ;
 }

 if (argc <= 1 ) n = 10;  // if no args are given
 else n = atoll(argv[1]); // converting input argument to size
 
 
 //MPI Datatypes to send rows of the locally allocated matrix
 MPI_Datatype MPI_ARRAYROW;
 MPI_Type_contiguous(n, MPI_INT, &MPI_ARRAYROW);
 MPI_Type_commit(&MPI_ARRAYROW);


 start_time = MPI_Wtime(); // start of wall time to 
 n_local  = n/np;      // no. of rows to be allocated by each processor 
 res = n%np;
 if(pid >= res) offset = res;
 if(pid  < res) n_local++;

 //Each processor allocates a local matrix contiguously
 local_matrix = (int **)calloc(n_local, sizeof(int*));
 local_matrix[0] = (int *)calloc(n*n_local, sizeof(int));

 for (int i=1;i<n_local;i++) {
  local_matrix[i] = local_matrix[0] + i*n;
 }
			
 unsigned long int r_local_start = n_local * pid + offset, r_local_end = n_local * (pid+1) + offset; 		
		
 for(i=0,j=r_local_start; i<n_local, j<r_local_end; i++, j++) 
   local_matrix[i][j] = 1;   //substitute 1 at the specific index for the identity matrix
		
 //printf("The matrix in processor pid %d\n", pid);
 //print_matrix(local_matrix, n_local, n);

 if(pid != 0){		
   //Each processor sends the local matrix to the root
   MPI_Send(local_matrix[0],n_local,MPI_ARRAYROW,0,200,MPI_COMM_WORLD);
   end_time=MPI_Wtime();
 }
 else{	
   //Recieve local matrices from each processor and prints them out
  printf("The initalized matrix is :\n");
  print_matrix(local_matrix, n_local, n);
  for(j=1; j<np; j++){
   if(j == res) n_local--;
   MPI_Recv(local_matrix[0],n_local,MPI_ARRAYROW,j,200,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   print_matrix(local_matrix, n_local, n); //print the matrix
   }
   end_time=MPI_Wtime();
 }

 MPI_Barrier(MPI_COMM_WORLD);

 if( pid != 0) printf ( "\n # walltime on processor %i : %10.8f \n",pid, end_time - start_time ) ; //calculating the wall time on each processor
 else printf ( "\n # walltime on master processor : %10.8f \n", end_time - start_time ) ; // wall time on root

  MPI_Finalize();
  return 0;
}

	
