#include <mpi.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <unistd.h> 

int main(int argc, char* argv[]) 
{ 
	unsigned long int elements_per_process, sum=0, partial_sum=0, left_over_sum=0; 
	int pid, np, master =0, tag =1;
	double start_time, end_time;   // times
	MPI_Status status;
	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	if ( argc <=1 || np <= 1) {                   			// if no arguments are given
		fprintf (stderr , " Usage : mpi -np n %s number_of_iterations \n", argv[0] ) ;
		MPI_Finalize() ;
		exit(-1) ;
	}

	unsigned long int n = atoll(argv[1]); // converting input argument to size

	start_time = MPI_Wtime(); // start of wall time
	elements_per_process = n/(np-1);  // estimating number of elements to be given for each processor

	if(pid !=0){				// slaves to generate partial sums
		partial_sum = 0, sum=0;
		unsigned long int i=0;
		unsigned long int index_start = elements_per_process * (pid-1), index_end = elements_per_process * pid; // calculating the interval in which the partial sum has to be calculated
		for(i=index_start+1;i<=index_end;i++)		// calculating the partial sum		
			partial_sum += i;
		end_time=MPI_Wtime();
		printf ( "\n # walltime on processor %i : %10.8f \n",pid, end_time - start_time ) ; //calculating the wall time on each processor
	}
	else{	//Master to integrate the partial sums and the left over sum
		unsigned long int i=0;
		left_over_sum=0;
		sum=0;			
		unsigned long int elements_left = n%(np-1);	//calculating the number of left over elements
                unsigned long int index_start = elements_per_process*(np-1), index_end = n;	// the left over elements are the last few elements of N numbers
                if(elements_left > 1)		// if there are more than 1 left over elements find the sum of the same
			for(i=index_start+1;i<=index_end;i++)
				left_over_sum += i;
		else if(elements_left == 1)	// if there is only one left over element add the last element 
			left_over_sum = n;
		else 				// if there are no left over elements the left over sum would be zero
			left_over_sum = 0;

		sum += left_over_sum;		// add the left over sum
		
//		printf("Left over sum in the master is %lu \n",sum);
		end_time=MPI_Wtime();
		printf ( "\n # walltime on master processor : %10.8f \n", end_time - start_time ) ;
	}

//	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&partial_sum, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD); // Doing an MPI_Reduce on the master
	if(pid == 0){
		sum += left_over_sum;			// Summing up with the left over sum
		printf("The total sum is %lu\n",sum);
	}

	MPI_Finalize();
	return 0;
}

	
