#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <openacc.h>

void saxpy(int n, double a, double* x, double* restrict y){

#pragma acc parallel loop
for(int i=0; i<n; i++){
  y[i] = a*x[i] + y[i];
 }

}


int main(int argc, char** argv){

int N = 1<<20;

if(argc > 1)
 N = atoi(argv[1]);

double* x = (double*) malloc (N*sizeof(double));
double* y = (double*) malloc (N*sizeof(double));

for(int i=0; i<N; i++){
 x[i] = 2.0;
 y[i] = 1.0;
}

saxpy(N, 3.0, x, y);

return 0;

}


