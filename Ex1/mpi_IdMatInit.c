/* Ex1 - IDENTITY MATRIX INITATION USING MPI 
 * - Parallel Programming 
 * - Prof. Ivan Girotto
 * - Exercise submitted by Keerthana C J
 *
 * This is a program which initalizes a square Identity Matrix             
 * of dimension 'n' provided as arguments.
 *				 
 * If dimension is not provided it creates a default size of 10.		 
 *
 * Each processor initiates a local matrix of dimension n * n/np, where np 
 * is the number of processors used.				 	
 *
 * After intitation, the matrix is sent to the root for assimilation.	 
 * The root allocates matrix an entire matrix of size n*n. It recieves the 
 * local matrices from the other processors and initates the rows that might    
 * be left out from n/np.
 *
 * Improvements suggested to this program
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
 unsigned long int n, n_local, i, j; 
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
 n_local  = n/(np-1);      // no. of rows to be allocated by each processor 
	
 if(pid != 0){
  //Each processor allocates a local matrix contiguously
  local_matrix = (int **)calloc(n_local, sizeof(int*));
  local_matrix[0] = (int *)calloc(n*n_local, sizeof(int));

  for (int i=1;i<n_local;i++) {
   local_matrix[i] = local_matrix[0] + i*n;
  }
		
	
  unsigned long int r_local_start = n_local * (pid-1), r_local_end = n_local * pid; 		
		
  for(i=0,j=r_local_start; i<n_local, j<r_local_end; i++, j++) 
   local_matrix[i][j] = 1;   //substitute 1 at the specific index for the identity matrix
		
   //printf("The matrix in processor pid %d\n", pid);
   //print_matrix(local_matrix, n_local, n);
		
   //Each processor sends to the root
   MPI_Send(local_matrix[0],n_local,MPI_ARRAYROW,0,pid,MPI_COMM_WORLD);
   end_time=MPI_Wtime();
 }

 if(pid == 0){
   //contiguous allocation of global matrix by the rooot
   matrix = (int **)calloc(n, sizeof(int*));
   matrix[0] = (int *)calloc(n*n, sizeof(int));
	
   for (i=1;i<n;i++) { 
     matrix[i] = matrix[0] + i*n; 
   }

   
   unsigned long int row_local_start = n_local*(np-1), row_local_end = n;
   //rows left over from n/np
   for (i=row_local_start,j=row_local_start;i<row_local_end,j<row_local_end;i++,j++) {
    matrix[i][j] = 1;
   }
		
   //Print the allocated matrix by the root if the size if less than 16
   if(n < 16) {  
    printf("The global matrix intialized in the master before recieving from other procs\n");
    print_matrix(matrix, n, n);
   }
		
   //Recieve local matrices from each processor
   for(i=1, j=0; i<np, j<row_local_start; i++,j+=n_local){
    MPI_Recv(matrix[j],n_local,MPI_ARRAYROW,i,i,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   }

   end_time=MPI_Wtime();

  //Print out the matrix if n < 16
  if(n < 16) { 
   printf("After assembling all the matrices in the global matrix\n");
   print_matrix(matrix, n, n); 
  }
  
 }

 MPI_Barrier(MPI_COMM_WORLD);
 if( pid != 0)
  printf ( "\n # walltime on processor %i : %10.8f \n",pid, end_time - start_time ) ; //calculating the wall time on each processor
 else
  printf ( "\n # walltime on master processor : %10.8f \n", end_time - start_time ) ; // wall time on root
  MPI_Finalize();
  return 0;
}

	
