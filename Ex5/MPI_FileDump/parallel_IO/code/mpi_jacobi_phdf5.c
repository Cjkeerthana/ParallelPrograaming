#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <hdf5.h>
/*
#define BLOCKING
#define TEXT
#define hdf5
*/
#define par_hdf5

/*** function declarations ***/

// save matrix to text file
void save_gnuplot_txt( FILE* file, double *M, size_t dim, size_t start, size_t end);

//save matrix to binary file
void save_gnuplot_bin( FILE* file, double *M, size_t dimension, size_t row_start, size_t row_end, size_t glob_dim);

// evolve Jacobi
void evolve( double * matrix, double *matrix_new, size_t row_start, size_t row_end, size_t dimension );

//print the matrix
void print ( double *matrix, size_t loc_dimension, size_t dimension);

//exchange borders between processors non blocking
void exchange_borders_nonblocking (double *matrix, size_t loc_dimension, size_t dimension, int prev, int next, MPI_Datatype MPI_ARRAYROW, MPI_Request **request);

//exchange borders between processors blocking
void exchange_borders_blocking (double *matrix, size_t loc_dimension, size_t dimension, int prev, int next, MPI_Datatype MPI_ARRAYROW);

// return the elapsed time
double seconds( void );

//print into text or bin files
void print_file(double *matrix_new, size_t loc_dimension, size_t dimension, int rank, int np, size_t res, size_t it, MPI_Datatype MPI_ARRAYROW);

//print into hdf5 files
void print_hdf5_time_evolution( double * matrix, // memeory buffer 
				size_t d_loc, // local dimension (boundaries excluded)
				size_t dimension, // global dimension (boundaries excluded)
				int rank, 
				int npes, 
				size_t rest, // rest of the data assigned to a given process
				size_t it );  // iteration  

//print into hdf5 files parallely
void print_hdf5_time_evolution_par( double * matrix, // memeory buffer 
				    size_t d_loc, // local dimension (boundaries excluded)
				    size_t dimension, // global dimension (boundaries excluded)
				    int rank, 
				    int npes, 
				    size_t rest, // rest of the data assigned to a given process
				    size_t offset, // offset regarding the rest distribution assigned to a given process     
				    size_t it );  // iteration

/*** end function declaration ***/

int main(int argc, char* argv[]){
			
  //MPI Initialize
  int rank, np, prev, next;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Status status[2];
  MPI_Request* request;
  request = (MPI_Request *)malloc(2*sizeof(MPI_Request));


  // timing variables
  double t_start, t_end, t_start_comm, t_end_comm, t_start_comp, t_end_comp, t_start_dump, t_end_dump, t_comm = 0, t_comp = 0, t_dump =0, increment, increment_offset;

  // indexes for loops
  size_t i, j, it;
  
  // initialize matrix
  double *matrix, *matrix_new, *tmp_matrix;

  size_t dimension = 0, iterations = 0, row_peek = 0, col_peek = 0, loc_dimension = 0, res = 0, offset = 0, dump = 0;
  size_t byte_dimension = 0;

  // check on input parameters
  if(rank == 0 && argc != 6) {
    fprintf(stderr,"\nwrong number of arguments. Usage: ./a.out dim it n m\n");
    MPI_Abort(MPI_COMM_WORLD,1);
    exit(-1);
  }

  dimension = atoi(argv[1]);
  loc_dimension = dimension/np;
  res = dimension%np;
  if(rank < res) loc_dimension++;
  iterations = atoi(argv[2]);
  row_peek = atoi(argv[3]);
  col_peek = atoi(argv[4]);
  dump = atoi(argv[5]);

  if(rank == 0){
    printf("matrix size = %zu\n", dimension);
    printf("number of iterations = %zu\n", iterations);
    printf("element for checking = Mat[%zu,%zu]\n",row_peek, col_peek);
    printf("write interval = %zu\n", dump);

   if((row_peek > dimension) || (col_peek > dimension)){
     fprintf(stderr, "Cannot Peek a matrix element outside of the matrix dimension\n");
     fprintf(stderr, "Arguments n and m must be smaller than %zu\n", dimension);
     return 1;
   }
  }


  //allocate local matrices
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
    
   //Exchange borders	  
#ifdef BLOCKING	  
    t_start_comm = seconds();
    exchange_borders_blocking(matrix, loc_dimension, dimension, prev, next, MPI_ARRAYROW);
    t_end_comm = seconds();   
#else
    t_start_comm = seconds();
    exchange_borders_nonblocking(matrix, loc_dimension, dimension, prev, next, MPI_ARRAYROW, &request);
#endif
    
    //evolve interior
    t_start_comp = seconds();
    evolve( matrix, matrix_new, 2, loc_dimension, dimension );	//bulk
    t_end_comp = seconds();

#ifndef BLOCKING    
    MPI_Waitall(2, request, status);
    t_end_comm = seconds();
#endif

    //evolve borders
    evolve( matrix, matrix_new, loc_dimension, loc_dimension+1, dimension);
    evolve( matrix, matrix_new, 1, 2, dimension);

    // swap the pointers
    tmp_matrix = matrix;
    matrix = matrix_new;
    matrix_new = tmp_matrix;

    //printf("%zu\n", it%dump);
    //dump into files
    t_start_dump = seconds();
    if((it+1)%dump == 0 || it == iterations-1) {

#ifdef hdf5
    print_hdf5_time_evolution(matrix_new, loc_dimension, dimension, rank, np, res, it+1);
#endif
 
#ifdef par_hdf5
    print_hdf5_time_evolution_par(matrix_new, loc_dimension, dimension, rank, np, res, offset, it+1); 
#else
    print_file(matrix_new, loc_dimension, dimension, rank,  np,  res, it, MPI_ARRAYROW);
#endif
    }
    t_end_dump = seconds();


    t_comm = t_comm + t_end_comm - t_start_comm;
    t_comp = t_comp + t_end_comp - t_start_comp;
    t_dump = t_dump + t_end_dump - t_start_dump;
  }
  t_end = seconds();
  

   printf( "\n total elapsed time on processor %d = %f seconds\n", rank, t_end - t_start );
   printf( "\n communication elapsed time on processor %d = %f seconds\n",rank, t_comm );
   printf("\n computation time for the evolving the bulk by processor %d = %f seconds\n",rank, t_comp);
   printf("\n computation time for dumping the files by processor %d = %f seconds\n", rank, t_dump);
   printf("\n");
  
  
  //peek into row and coloumn
  size_t glob_dimension_start = rank * loc_dimension + 1  +  offset;
  size_t glob_dimension_end = glob_dimension_start + loc_dimension - 1;
  if(row_peek >= glob_dimension_start && row_peek <= glob_dimension_end){
    size_t loc_dimension_row_peek = row_peek - (rank * loc_dimension) - offset;
    printf( "\nmatrix[%zu,%zu] = %f\n", row_peek, col_peek, matrix[ ( loc_dimension_row_peek + 1 ) * ( dimension + 2 ) + ( col_peek + 1 ) ] );
  }
  
  free( matrix );
  free( matrix_new );

  MPI_Finalize();

  return 0;
}

// A Simple timer for measuring the walltime
double seconds(){

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}

void exchange_borders_nonblocking (double *matrix, size_t loc_dimension, size_t dimension, int prev, int next, MPI_Datatype MPI_ARRAYROW, MPI_Request **request){
   
   MPI_Isend(matrix + (dimension+2), 1, MPI_ARRAYROW, prev, 100, MPI_COMM_WORLD, *request);
   MPI_Isend(matrix + loc_dimension*(dimension+2), 1, MPI_ARRAYROW, next, 200, MPI_COMM_WORLD, *request+1);
   MPI_Irecv(matrix, 1, MPI_ARRAYROW, prev, 200, MPI_COMM_WORLD, *request);
   MPI_Irecv(matrix + (loc_dimension+1)*(dimension+2), 1, MPI_ARRAYROW, next, 100, MPI_COMM_WORLD, *request+1);

}


void exchange_borders_blocking (double *matrix, size_t loc_dimension, size_t dimension, int prev, int next, MPI_Datatype MPI_ARRAYROW){
  MPI_Sendrecv(matrix + (dimension+2), 1, MPI_ARRAYROW, prev, 100, matrix + (loc_dimension+1)*(dimension+2), 1, MPI_ARRAYROW, next, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Sendrecv(matrix + loc_dimension*(dimension+2), 1, MPI_ARRAYROW, next, 200, matrix, 1, MPI_ARRAYROW, prev, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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


void save_gnuplot( FILE* file, double *M, size_t dimension, size_t row_start, size_t row_end, size_t glob_dim){
  
  size_t i , j;
  const double h = 0.1;

#ifdef TEXT
  for( i = row_start; i <= row_end; ++i )
    for( j = 0; j < dimension + 2; ++j )
      fprintf(file, "%f\t%f\t%f\n", h * j, -h * (i + glob_dim), M[ ( i * ( dimension + 2 ) ) + j ] );
#else
  for( i = row_start; i <= row_end; ++i )
   fwrite(M+i*(dimension+2), sizeof(double), (dimension+2), file);
#endif
}

//dump into a text file or a binary file
void print_file(double *matrix_new, size_t loc_dimension, size_t dimension, int rank, int np, size_t res, size_t it, MPI_Datatype MPI_ARRAYROW){
  FILE* file;
  const char fname[] = "mpi_solution_t";
  char time_step[10];
  size_t glob_dim = 0;
#ifdef TEXT
  char extension[] = ".dat";
#else 
  char extension[] = ".bin";
#endif
  if(rank != 0 && rank != np-1) MPI_Send(matrix_new+(dimension+2), loc_dimension, MPI_ARRAYROW, 0, 300, MPI_COMM_WORLD);
  if(rank == np-1) MPI_Send(matrix_new+(dimension+2), loc_dimension+1, MPI_ARRAYROW, 0, 300, MPI_COMM_WORLD);
  if(rank == 0){
    char filename[100] = "\0";
    strcat(filename, fname);
    sprintf(time_step, "%zu", it+1);
    strcat(filename, time_step);
    strcat(filename, extension);
#ifdef TEXT
    file = fopen( filename, "w" );
#else
    file = fopen(filename,"wb");
#endif
    save_gnuplot(file, matrix_new, dimension, 0, loc_dimension, glob_dim);
    for(size_t i=1; i<np; ++i){
      if(i == res)loc_dimension--;
      if(i == np-1)loc_dimension++;
      MPI_Recv(matrix_new+(dimension+2), loc_dimension, MPI_ARRAYROW, i, 300, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      glob_dim = glob_dim + loc_dimension;
      save_gnuplot(file, matrix_new, dimension, 1, loc_dimension, glob_dim);
     }
     fclose(file);
   }
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

void print_hdf5_time_evolution( double * matrix,  size_t d_loc, size_t dimension, int rank, int npes, size_t rest, size_t it )  
{
   int include_lower_border = 0;
  int offset = 0;

  if( !rank ) {

    if( rank == npes - 1 ) include_lower_border = 1;  // if rank 0 has the last row i.e. npes =1, must print it
    
    hid_t   file_id, group_id,dataset_id,dataspace_id;  /* identifiers */
    herr_t  status;
    
    char fname[100];
    sprintf(fname,"test_%d.h5", it);
    /* Open an existing HDF5 file or create if not existing. */
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    /* set the file dimensions = data domain dimension, in our case */
    hsize_t dims_data[ 2 ] = { dimension + 2, dimension + 2 };
    hid_t file_space = H5Screate_simple(2, dims_data, NULL);
    dataset_id = H5Dcreate( file_id, "/temp", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    /* set the dimensions of the local data dumped by process 0 */
    dims_data[ 0 ] = d_loc + include_lower_border + 1;
    dims_data[ 1 ] = dimension + 2;
    hid_t mem_space = H5Screate_simple(2, dims_data, NULL);

    /* set the global starting point of the 0 process slice in regards to the file dimensions (global domain dimensions) */ 
    hsize_t start_3d_mem[ 2 ];
    start_3d_mem[ 0 ] = 0;
    start_3d_mem[ 1 ] = 0;

    /* set the hyperslab (portion) of the data to write on the file */
    status = H5Sselect_hyperslab ( file_space, H5S_SELECT_SET, start_3d_mem, NULL, dims_data, NULL );
    
    /* write the portion local to process 0 */
    status = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, matrix );
    
    int i = 0;
    for( i = 1; i < npes; i++ ){

      if( i == rest ){
	d_loc -= 1;
	offset = rest;
      }

      if( i == npes - 1 ) include_lower_border = 1;

      MPI_Recv( matrix, ( dimension + 2 ) * ( d_loc + include_lower_border ), MPI_DOUBLE, i, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      /* set the dimensions of the data to print as received from process i */
      H5Sclose( mem_space );
      dims_data[ 0 ] = d_loc + include_lower_border;
      dims_data[ 1 ] = dimension + 2;
      mem_space = H5Screate_simple( 2, dims_data, NULL );

      /* set the global starting point of the ith process slice in regards to the file dimensions (global domain dimensions) */ 
      start_3d_mem[ 0 ] = d_loc * i + 1 + offset;
      start_3d_mem[ 1 ] = 0;

      status = H5Sselect_hyperslab ( file_space, H5S_SELECT_SET, start_3d_mem, NULL, dims_data, NULL );
      status = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, matrix );
    }
    // relase the objects to avoid memory leaks
      status = H5Sclose( mem_space );
      status = H5Sclose( file_space );
      status = H5Dclose( dataset_id );    
      status = H5Fclose( file_id );
     } else {
    if( rank == npes - 1 ) d_loc++;
    MPI_Send( matrix + ( dimension + 2 ), ( dimension + 2 ) * d_loc, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
    }
}


void print_hdf5_time_evolution_par( double * matrix, // memeory buffer 
				    size_t d_loc, // local dimension (boundaries excluded)
				    size_t dimension, // global dimension (boundaries excluded)
				    int rank, 
				    int npes, 
				    size_t rest, // rest of the data assigned to a given process
				    size_t offset, // offset regarding the rest distribution assigned to a given process     
				    size_t it )  // iteration  
{
  // Init Par HDF5
  H5Eset_current_stack (H5E_DEFAULT);
  hid_t plist_id = H5Pcreate (H5P_FILE_ACCESS);
  hid_t hdf5_status = H5Pset_fapl_mpio (plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  
  hid_t   file_id, group_id,dataset_id,dataspace_id;  /* identifiers */
  herr_t  status;
  
  char fname[100];
  sprintf(fname,"test_%d.h5", it);
  /* Open an existing HDF5 file or create if not existing. */
  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    
  /* set the file dimensions = data domain dimension, in our case */
  hsize_t dims_data[ 2 ] = { dimension + 2, dimension + 2 };
  hid_t file_space = H5Screate_simple(2, dims_data, NULL);
  dataset_id = H5Dcreate( file_id, "/temp", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  
  /* set the dimensions of the local data dumped by process 0 */
  dims_data[ 0 ] = !rank || rank == npes - 1 ? d_loc + 1 : d_loc;
  if( npes == 1 ) dims_data[ 0 ] += 1; 
  dims_data[ 1 ] = dimension + 2;

  hid_t mem_space = H5Screate_simple(2, dims_data, NULL);
  
  /* set the global starting point of the 0 process slice in regards to the file dimensions (global domain dimensions) */ 
  hsize_t start_3d_mem[ 2 ];
  start_3d_mem[ 0 ] = rank ? d_loc * rank + offset + 1 : 0;
  start_3d_mem[ 1 ] = 0;
  
  /* set the hyperslab (portion) of the data to write on the file */
  status = H5Sselect_hyperslab ( file_space, H5S_SELECT_SET, start_3d_mem, NULL, dims_data, NULL );
  
    /* set pbject for par IO */
  hid_t xfer_plist = H5Pcreate (H5P_DATASET_XFER);
  herr_t ret = H5Pset_dxpl_mpio (xfer_plist, H5FD_MPIO_COLLECTIVE);
  
  /* write the portion local to process 0 */
  status = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, mem_space, file_space, xfer_plist, rank ? matrix + (dimension+2) : matrix );
  
  // relase the objects to avoid memory leaks
  status = H5Sclose( mem_space );
  status = H5Pclose( xfer_plist );
  status = H5Pclose( plist_id );
  status = H5Sclose( file_space );
  status = H5Dclose( dataset_id );    
  status = H5Fclose(file_id);
}
