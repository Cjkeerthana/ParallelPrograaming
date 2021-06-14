#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <openacc.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <cublas.h>
//#include "cuda2acc.h"

void print_matrix(double *M, size_t rows, size_t cols);

int main(int argc, char* argv[]) 
{ 

 int rank, np, igpu, ngpu;
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
 
 double t_start, t_end, t_start_calc, t_end_calc, t_start_comm, t_end_comm, t_start_copy, t_end_copy, t_calc =0, t_comm = 0, t_copy =0, t_tot=0;
 double all_tot = 0, all_comm = 0, all_calc = 0, all_copy = 0;

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

 ngpu = acc_get_num_devices(acc_device_nvidia);
 igpu = rank % ngpu;
 acc_set_device_num(igpu, acc_device_nvidia);

t_start_copy = MPI_Wtime();
#pragma acc enter data create(A_loc[0:n_loc*n], B_loc[0:n_loc*n], C_loc[0:n_loc*n])
#pragma acc enter data create(B_temp[0:n*n_loc_max])
#pragma acc update device(A_loc[0:n_loc*n], B_loc[0:n_loc*n])
t_end_copy = MPI_Wtime();

t_copy += t_end_copy-t_start_copy;

t_start = MPI_Wtime();
for(i=0; i<np; i++){ 
t_start_comm = MPI_Wtime();
#pragma acc host_data use_device(B_loc, B_temp)
 {
  if(i < res)  MPI_Allgatherv(B_loc+i*n_loc_max, 1, column_type1, B_temp, count1, displace1, MPI_DOUBLE, MPI_COMM_WORLD);
  else MPI_Allgatherv(B_loc+i*n_loc_min+res, 1, column_type2, B_temp, count2, displace2, MPI_DOUBLE, MPI_COMM_WORLD);
 }
#pragma acc wait
t_end_comm = MPI_Wtime();

t_start_calc = MPI_Wtime();
#pragma acc host_data use_device(A_loc, B_temp, C_loc)
 {
  if(i < res) cublasDgemm('n', 'n', n_loc_max, n_loc, n, 1.0, B_temp, n_loc_max, A_loc, n, 0.0, C_loc+(i*n_loc_max), n);
  else cublasDgemm('n', 'n', n_loc_min, n_loc, n, 1.0, B_temp, n_loc_min, A_loc, n, 0.0, C_loc+(i*n_loc_min)+res, n);
 }
t_end_calc = MPI_Wtime();

t_calc += t_end_calc - t_start_calc;
t_comm += t_end_comm - t_start_comm;
}
t_end = MPI_Wtime();
t_tot = t_end - t_start;

t_start_copy = MPI_Wtime();
#pragma acc update host(C_loc[0:n_loc*n])
t_end_copy = MPI_Wtime();
 
t_copy += t_end_copy-t_start_copy;

 //gather times from all processors
 MPI_Reduce(&t_tot, &all_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
 MPI_Reduce(&t_comm, &all_comm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
 MPI_Reduce(&t_calc, &all_calc, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
 MPI_Reduce(&t_copy, &all_copy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

 if(rank == 0) {
  printf("\nTimes\n");
  printf("Total time for multiplication averaged on all processors is %f\n", (all_tot+all_copy)/np);
  printf("Total time for communication averaged on all processors is %f\n", all_comm/np);
  printf("Total time for calculation averaged on all processors is %f\n", all_calc/np);
  printf("Total time for copying averaged on all processors is %f\n", all_copy/np);
 }

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


