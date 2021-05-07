#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <omp.h>
#include <sched.h>
#include <arpa/inet.h>
#include "mkl_cblas.h"

//A function for printing a matrix
void print_matrix(double *M, size_t rows, size_t cols);

//A function for matrix multiplication C(m/n) = A(m/p)*B(p/n)
void vector_mul(double* A, double* B, double* C, size_t m, size_t n, size_t p);

int main(int argc, char* argv[]) 
{ 
 int rank, np;
 MPI_Status status;
 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &np);
 
 if(argc <= 1 || np < 1){
   fprintf(stderr, "Usage: mpi -np n %s number_of_iterations \n", argv[0]);
   MPI_Finalize();
   exit(-1);	
 }
 
 size_t n = atoll(argv[1]);

 size_t n_loc=0 , i=0, j=0, glob_i=0, glob_j=0, offset=0, res=0, temp=0, n_loc_max=0, n_loc_min=0;
 double start_time, end_time, start_comm_time, end_comm_time, start_calc_time, end_calc_time;
 double comm_time=0, calc_time=0, tot_time=0;
 double all_comm_time=0, all_calc_time=0, all_tot_time=0;
 n_loc = n/np;
 res = n%np;
 if(rank >= res)offset = res;
 if(rank < res) n_loc++;
 n_loc_max = (res > 0) ? n/np+1 : n/np;
 n_loc_min = n/np;
 double *A_loc, *B_loc, *C_loc, *B_temp;
 
 size_t byte_dim = 0;
 byte_dim = sizeof(double)*n*n_loc;
 //Every processsor allocates its local matrix
 A_loc = (double *) malloc(byte_dim);
 B_loc = (double *) malloc(byte_dim);
 C_loc = (double *) malloc(byte_dim);

 memset( A_loc, 0, byte_dim );
 memset( B_loc, 0, byte_dim );
 memset( C_loc, 0, byte_dim );

 srand(rank*124);
 
 for(i=0; i<n_loc; i++){
  for(j=0; j<n; j++){
     A_loc[i*n+j] = rand()%50;   // A random matrix for A
     //B_loc[i*n+j] = rand()%50; // A random matrix for B
  }
  glob_i = n_loc*rank + offset + i;
  B_loc[i*n + glob_i] = 1; //Initialize an identity matrix for B
 }

 MPI_Barrier(MPI_COMM_WORLD);

 //MPI Datatype for gathering every row
 MPI_Datatype MPI_ARRAYROW;
 MPI_Type_contiguous(n, MPI_DOUBLE, &MPI_ARRAYROW);
 MPI_Type_commit(&MPI_ARRAYROW);

 //Print the A matrix for checking 
 if(n < 16){
  if(rank != 0) MPI_Send(A_loc,n_loc,MPI_ARRAYROW,0,200,MPI_COMM_WORLD);
  if(rank == 0){
     size_t temp = n_loc;
     printf("The matrix A is \n");
     print_matrix(A_loc, temp, n);
     for(j=1; j<np; j++){
      if(j == res) temp--;
      MPI_Recv(C_loc,temp,MPI_ARRAYROW,j,200,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      print_matrix(C_loc, temp, n);
     }
     memset( C_loc, 0, byte_dim );
  } 
 }

 //Allocate a temporary column in each processor (which is stored contiguously)
 B_temp = (double *) malloc(sizeof(double)*n*n_loc_max); 
 memset(B_temp, 0, sizeof(double)*n*n_loc_max);
 
 //Datatype for gathering every column -- 2 types of blocks for handling rests
 MPI_Datatype column_type1, column_type2;
 MPI_Type_vector(n_loc, n_loc_max, n, MPI_DOUBLE, &column_type1);
 MPI_Type_vector(n_loc, n_loc_min, n, MPI_DOUBLE, &column_type2);
 MPI_Type_commit(&column_type1);
 MPI_Type_commit(&column_type2);

 int *displace1, *count1, *displace2, *count2;
 //for column_type1
 displace1 = (int *) calloc(np, sizeof(int)); 
 count1 = (int *) calloc(np, sizeof(int));
 //for column_type2
 displace2 = (int *) calloc(np, sizeof(int));
 count2 = (int *) calloc(np, sizeof(int));
 
 for(i=0; i<np; i++){
    if(i < res){
     count1[i] = n_loc_max*n_loc_max; 
     count2[i] = n_loc_min*n_loc_max;
     displace1[i] = i*count1[i];
     displace2[i] = i*count2[i];
    }
    else{
     count1[i] = n_loc_max*n_loc_min;
     count2[i] = n_loc_min*n_loc_min;
     displace1[i] = i*count1[i] + res*n_loc_max;
     displace2[i] = i*count2[i] + res*n_loc_min;
    }
 }
 
 start_time = MPI_Wtime();
 for(i=0; i<np; i++){

 //Gather every coloumn in temporary B coloumn of each processor
  start_comm_time = MPI_Wtime();
  if(i < res)  MPI_Allgatherv(B_loc+i*n_loc_max, 1, column_type1, B_temp, count1, displace1, MPI_DOUBLE, MPI_COMM_WORLD);
  else MPI_Allgatherv(B_loc+i*n_loc_min+res, 1, column_type2, B_temp, count2, displace2, MPI_DOUBLE, MPI_COMM_WORLD);
  end_comm_time = MPI_Wtime();


 //Each row of the matrix A is multiplied with the gathered column
  start_calc_time = MPI_Wtime();
#ifdef __DGEMM
  if(i < res) cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, n_loc, n_loc_max, n, 1.0, A_loc, n, B_temp, n_loc_max, 0.0, C_loc+(i*n_loc_max), n); 		
  else cblas_dgemm( CblasRowMajor, CblasNoTrans, CblasNoTrans, n_loc, n_loc_min, n, 1.0, A_loc, n, B_temp, n_loc_min, 0.0, C_loc+(i*n_loc_min)+res, n);
#else
  if(i < res)vector_mul(A_loc, B_temp, C_loc+(i*n_loc_max), n_loc, n_loc_max, n);
  else vector_mul(A_loc, B_temp, C_loc+(i*n_loc_min)+res, n_loc, n_loc_min, n);
#endif
  end_calc_time = MPI_Wtime();

  comm_time += end_comm_time - start_comm_time;
  calc_time += end_calc_time - start_calc_time;
 }
 end_time = MPI_Wtime();
 tot_time = end_time-start_time;

 //Print C Matrix for checking
 if(n < 16){
  if(rank != 0) MPI_Send(C_loc,n_loc,MPI_ARRAYROW,0,200,MPI_COMM_WORLD);
  if(rank == 0){
     printf("The matrix C is \n");
     print_matrix(C_loc, n_loc, n);
     for(j=1; j<np; j++){
     if(j == res) n_loc--;
     MPI_Recv(C_loc,n_loc,MPI_ARRAYROW,j,200,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     print_matrix(C_loc, n_loc, n);
     }
  }
 }

 //gather times from all processors
 MPI_Reduce(&tot_time, &all_tot_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
 MPI_Reduce(&comm_time, &all_comm_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
 MPI_Reduce(&calc_time, &all_calc_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

 if(rank == 0) {
  printf("\nTimes\n");
  printf("Total time for multiplication averaged on all processors is %f\n", all_tot_time/np);
  printf("Total time for communication averaged on all processors is %f\n", all_comm_time/np);
  printf("Total time for calculation averaged on all processors is %f\n", all_calc_time/np);
 }

 MPI_Finalize();
 return 0;
}

void print_matrix(double *M, size_t rows, size_t cols){
 size_t i=0, j=0;
 for(i=0; i<rows; i++){
  for(j=0; j<cols; j++){
   printf("%.2f\t", M[i*cols + j]);
  }
  printf("\n");
 }
 printf("\n");
}

//Here C takes a stride of z as the entire matrix of C has a stride of "n"
void vector_mul(double* A, double* B, double* C, size_t x, size_t y, size_t z){

#pragma omp parallel for proc_bind(master)
 for(size_t i=0; i<x; i++)
  for(size_t j=0; j<y; j++)
    for(size_t k=0; k<z; k++)
       C[i*z+j] += A[i*z+k]*B[k*y+j];     
}
