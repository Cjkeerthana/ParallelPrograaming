#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

/*** function declarations ***/

// save matrix to file
void save_gnuplot( double *M, size_t dim, size_t loc_dim, size_t glob_dimension_start, size_t start, size_t end);

// evolve Jacobi
void evolve( double * matrix, double *matrix_new, size_t row_start, size_t row_end, size_t dimension );

void print ( double *matrix, size_t loc_dimension, size_t dimension);
// return the elapsed time
double seconds( void );

/*** end function declaration ***/

int main(int argc, char* argv[]){
			
  //MPI Initialize
  int rank, np, prev, next;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Status status[2];
  MPI_Request request[2], request1[2];


  // timing variables
  double t_start, t_end, increment, increment_offset;

  // indexes for loops
  size_t i, j, it;
  
  // initialize matrix
  double *matrix, *matrix_new, *tmp_matrix;

  size_t dimension = 0, iterations = 0, row_peek = 0, col_peek = 0, loc_dimension = 0, res = 0, offset = 0;
  size_t byte_dimension = 0;

  // check on input parameters
  if(rank == 0 && argc != 5) {
    fprintf(stderr,"\nwrong number of arguments. Usage: ./a.out dim it n m\n");
    return 1;
  }

  dimension = atoi(argv[1]);
  loc_dimension = dimension/np;
  res = dimension%np;
  if(rank < res) loc_dimension++;
  iterations = atoi(argv[2]);
  row_peek = atoi(argv[3]);
  col_peek = atoi(argv[4]);

  if(rank == 0){
    printf("matrix size = %zu\n", dimension);
    printf("number of iterations = %zu\n", iterations);
    printf("element for checking = Mat[%zu,%zu]\n",row_peek, col_peek);

   if((row_peek > dimension) || (col_peek > dimension)){
     fprintf(stderr, "Cannot Peek a matrix element outside of the matrix dimension\n");
     fprintf(stderr, "Arguments n and m must be smaller than %zu\n", dimension);
     return 1;
   }
  }

  byte_dimension = sizeof(double) * ( dimension + 2 ) * ( loc_dimension + 2 );
  matrix = ( double* )malloc( byte_dimension );
  matrix_new = ( double* )malloc( byte_dimension );

  memset( matrix, 0, byte_dimension );
  memset( matrix_new, 0, byte_dimension );

  //fill initial values  
  for( i = 0; i <= loc_dimension+1; ++i )
    for( j = 1; j <= dimension; ++j )
      matrix[ ( i * ( dimension + 2 ) ) + j ] = 0.5;
	      
  // set up borders 
  increment = 100.0 / ( dimension+1 );
  if(rank >= res) offset = res;
  
  for( i=0; i < loc_dimension+1; ++i ){
    matrix[ i * ( dimension + 2 ) ] = i * increment + rank * loc_dimension * increment + offset * increment;
    matrix_new[ i * ( dimension + 2 ) ] = i * increment + rank * loc_dimension * increment + offset * increment;
  }

  if(rank == np-1){
   for( i=0; i<= dimension; ++i){
     matrix[(loc_dimension+1) * (dimension + 2) + i] = (dimension + 1 - i) * increment;
     matrix_new[(loc_dimension+1) * (dimension + 2) + i] = (dimension + 1 - i) * increment;
   }
  }
  else{
   matrix[(loc_dimension + 1) * (dimension + 2)] = increment + (rank+1) * loc_dimension * increment + offset * increment;
   matrix_new[(loc_dimension + 1) * (dimension + 2)] = increment + (rank+1) * loc_dimension * increment + offset * increment;
  }

  if(rank == 0){
   for( i=0; i<=dimension+1; ++i){
     matrix[i] = 0.0;
     matrix_new[i] = 0.0; 
   }
  }
  else{
   matrix[0] = loc_dimension * increment + (rank-1) * loc_dimension * increment + offset * increment;
   matrix_new[0] = loc_dimension * increment + (rank-1) * loc_dimension * increment + offset * increment;
  }


  prev = (rank - 1);
  next = (rank + 1);

  if(rank == 0) prev = MPI_PROC_NULL;
  if(rank == np-1) next = MPI_PROC_NULL;
  
  MPI_Datatype MPI_ARRAYROW;
  MPI_Type_contiguous(dimension+2, MPI_DOUBLE, &MPI_ARRAYROW);
  MPI_Type_commit(&MPI_ARRAYROW);

  // start algorithm
  t_start = seconds();
  for( it = 0; it < iterations; ++it ){
    MPI_Isend(matrix + (dimension+2), 1, MPI_ARRAYROW, prev, 100, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(matrix + loc_dimension*(dimension+2), 1, MPI_ARRAYROW, next, 200, MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(matrix, 1, MPI_ARRAYROW, prev, 200, MPI_COMM_WORLD, &request1[0]);
    MPI_Irecv(matrix + (loc_dimension+1)*(dimension+2), 1, MPI_ARRAYROW, next, 100, MPI_COMM_WORLD, &request1[1]);    
    evolve( matrix, matrix_new, 2, loc_dimension, dimension );	//bulk
    MPI_Waitall(2, request1, status);
    evolve( matrix, matrix_new, loc_dimension, loc_dimension+1, dimension);
    evolve( matrix, matrix_new, 1, 2, dimension);

    // swap the pointers
    tmp_matrix = matrix;
    matrix = matrix_new;
    matrix_new = tmp_matrix;

  }
  t_end = seconds();
  
  if(rank == 0) printf( "\nelapsed time = %f seconds\n", t_end - t_start );
  
  size_t glob_dimension_start = rank * loc_dimension + 1  +  offset;
  size_t glob_dimension_end = glob_dimension_start + loc_dimension - 1;
  if(row_peek >= glob_dimension_start && row_peek <= glob_dimension_end){
    size_t loc_dimension_row_peek = row_peek - (rank * loc_dimension) - offset;
    printf( "\nmatrix[%zu,%zu] = %f\n", row_peek, col_peek, matrix[ ( loc_dimension_row_peek + 1 ) * ( dimension + 2 ) + ( col_peek + 1 ) ] );
  }
 
  //printf("My rank is %d and my global dimension is %lu\n", rank, glob_dimension_start-1);
  /*
  if(rank == 0)save_gnuplot( matrix, dimension, loc_dimension, glob_dimension_start-1, 0, loc_dimension);
  if(rank == np-1) save_gnuplot( matrix, dimension, loc_dimension, glob_dimension_start, 1, loc_dimension+1);
  if(rank != 0 && rank != np-1)save_gnuplot( matrix, dimension, loc_dimension, glob_dimension_start, 1, loc_dimension);
  //if(rank == 0) print(matrix, loc_dimension, dimension);
  */

  MPI_File fp;
  MPI_Status status_file;

  MPI_File_open(MPI_COMM_WORLD, "mpi_solution.dat", MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
  MPI_Offset file_offset = sizeof(double)*(dimension+2)*loc_dimensio;
   

  free( matrix );
  free( matrix_new );

  MPI_Finalize();

  return 0;
}

void evolve( double * matrix, double *matrix_new, size_t row_start, size_t row_end, size_t dimension){
  
  size_t i , j;

  //This will be a row dominant program.
  for( i = row_start ; i < row_end; ++i )
    for( j = 1; j <= dimension; ++j )
      matrix_new[ ( i * ( dimension + 2 ) ) + j ] = ( 0.25 ) * 
	( matrix[ ( ( i - 1 ) * ( dimension + 2 ) ) + j ] + 
	  matrix[ ( i * ( dimension + 2 ) ) + ( j + 1 ) ] + 	  
	  matrix[ ( ( i + 1 ) * ( dimension + 2 ) ) + j ] + 
	  matrix[ ( i * ( dimension + 2 ) ) + ( j - 1 ) ] ); 
}


void save_gnuplot( double *M, size_t dimension, size_t loc_dimension, size_t glob_dimension_start, size_t row_start, size_t row_end){
  
  size_t i , j;
  const double h = 0.1;
  FILE *file;

  file = fopen( "solution_mpi.dat", "a+" );
  fseek(file, glob_dimension_start * (dimension + 2), SEEK_SET);

  for( i = row_start; i <= row_end; ++i )
    for( j = 0; j < dimension + 2; ++j )
      fprintf(file, "%f\t%f\t%f\n", h * j, -h * i, M[ ( i * ( dimension + 2 ) ) + j ] );

  fclose( file );

}

void print(double* matrix, size_t loc_dimension, size_t dimension){
  size_t i, j;
  for( i=0; i<loc_dimension+2; ++i){
   for( j=0; j<dimension+2; ++j){
    printf("%f\t", matrix[i*(dimension+2) + j]);
   }
   printf("\n");
  }
}

// A Simple timer for measuring the walltime
double seconds(){

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}

