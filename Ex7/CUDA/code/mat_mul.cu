#include<stdio.h>
#include<cuda_runtime.h>
#include<stdlib.h>
#include<time.h>
#include<sys/time.h>
#include<helper_cuda.h>
#include<cublas_v2.h>
#include"mat_mul.h"

double *d_A, *d_B, *d_C;
size_t size_A, size_B, size_C;
cublasHandle_t handle;
double t_start_cpy, t_end_cpy;
cudaEvent_t t_start_calc, t_end_calc;
const double alf = 1.0;
const double bet = 0.0;
const double *alpha = &alf;
const double *beta = &bet;

void cuda_initialize(double* A, size_t m, size_t n, size_t p, int rank, double* t_cpy)
{

 int deviceCount;
 cudaGetDeviceCount(&deviceCount);
 int device_id = rank % deviceCount;
 cudaSetDevice(device_id);
 
 size_A = m*p*sizeof(double);
 size_B = p*n*sizeof(double);
 size_C = m*n*sizeof(double);

 cudaMalloc((void **)&d_A, size_A);
 cudaMalloc((void **)&d_B, size_B);
 cudaMalloc((void **)&d_C, size_C);
 
 cublasCreate(&handle);

 t_start_cpy = seconds();
 cudaMemcpy(d_A, A, size_A, cudaMemcpyHostToDevice);
 t_end_cpy = seconds();
 
 *t_cpy += t_end_cpy - t_start_cpy;

}

void cuda_vector_mul(double* B, double** C, size_t m, size_t n, size_t p, double* t_cpy, double* t_calc){

 t_start_cpy = seconds();
 cudaMemcpy(d_B, B, size_B, cudaMemcpyHostToDevice);
 t_end_cpy = seconds();

 *t_cpy += t_end_cpy - t_start_cpy;

 //cudaThreadSynchronize();
 //printf("Calculating in GPU\n");
 //cudaEventRecord(t_start_calc, 0);
 cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, m, p, alpha, d_B, n, d_A, p, beta, d_C, n);
 //cudaThreadSynchronize();
 //cudaEventRecord(t_end_calc, 0);
 
 //float calc_time = 0.0f;
 //cudaEventElapsedTime(&calc_time, t_start_calc, t_end_calc);
 //*t_calc += (double) calc_time;

 t_start_cpy = seconds();
 cudaMemcpy(*C, d_C, size_C, cudaMemcpyDeviceToHost); 
 t_end_cpy = seconds();

 *t_cpy += t_end_cpy - t_start_cpy;

}

void cuda_stop(){
 
 cublasDestroy(handle);
 
 cudaFree(d_A);
 cudaFree(d_B);
 cudaFree(d_C);
}
