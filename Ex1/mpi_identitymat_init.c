#include <mpi.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <unistd.h> 


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
	unsigned long int n_local, i, j; 
	int **local_matrix, **matrix;
	int pid, np, master =0, tag =1;
	double start_time, end_time;   // times
	MPI_Status status;
	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	if ( argc <= 1 || np <= 1) {                   			// if no arguments are given
		fprintf (stderr , " Usage : mpi -np n %s number_of_iterations \n", argv[0] ) ;
		MPI_Finalize() ;
		exit(-1) ;
	}

	unsigned long int n = atoll(argv[1]); // converting input argument to size

	start_time = MPI_Wtime(); // start of wall time
	n_local  = n/(np-1);  // no. of rows to be allocated by each processor 

                local_matrix = (int **)calloc(n_local, sizeof(int*));
                local_matrix[0] = (int *)calloc(n*n_local, sizeof(int));

                for (int i=1;i<n_local;i++) {
                        local_matrix[i] = local_matrix[0] + i*n;
                }
		
	if(pid != 0){
		unsigned long int r_local_start = n_local * (pid-1), r_local_end = n_local * pid; 		
		for(i=0,j=r_local_start; i<n_local, j<r_local_end; i++, j++) 
			local_matrix[i][j] = 1;				    //substitute 1 at the specific index for the identity matrix
		
		print_matrix(local_matrix, n_local, n);
		
		end_time=MPI_Wtime();
		printf ( "\n # walltime on processor %i : %10.8f \n",pid, end_time - start_time ) ; //calculating the wall time on each processor
	}
	        matrix = (int **)calloc(n, sizeof(int*));
	        matrix[0] = (int *)calloc(n*n, sizeof(int));
	
	        for (i=1;i<n;i++) {
	                matrix[i] = matrix[0] + i*n;
	        }

	if(pid == 0)
	{	
		unsigned long int row_local_start = n_local*(np-1), row_local_end = n;
		for (i=row_local_start,j=row_local_start;i<row_local_end,j<row_local_end;i++,j++) {
	        	matrix[i][j] = 1;
	        }

		print_matrix(matrix, n, n);
		end_time=MPI_Wtime();
		printf ( "\n # walltime on master processor : %10.8f \n", end_time - start_time ) ;
	}
	
	MPI_Datatype MPI_ARRAYROW;
    	MPI_Type_contiguous(n, MPI_INT, &MPI_ARRAYROW);
    	MPI_Type_commit(&MPI_ARRAYROW);

	int *displace;
    	displace = (int *)calloc(np,sizeof(int));
    	int *dim_list;
   	dim_list = (int *)calloc(np,sizeof(int));
   	j = 0;
  	for (int i=0; i<np; i++) {
        	displace[i] = j;
        	j += n_local;
        	dim_list[i] = n;
    	}	

    	// MPI_Gatherv call.
    	MPI_Barrier(MPI_COMM_WORLD);
    	MPI_Gatherv(local_matrix[0],n_local,MPI_ARRAYROW,matrix[0],&dim_list[pid],&displace[pid],MPI_ARRAYROW,0,MPI_COMM_WORLD);
	if(pid == 0) { print_matrix(matrix, n, n); }
	MPI_Finalize();
	return 0;
}

	
