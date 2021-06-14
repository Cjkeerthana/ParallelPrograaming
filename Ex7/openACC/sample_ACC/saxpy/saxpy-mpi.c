#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <openacc.h>

void saxpy(double *x, double* restrict y, int a, int len)
{
    int i;
#pragma acc declare present(x, y)
#pragma acc parallel loop //deviceptr(x, y)
    for(i = 0 ; i < len; i++)
    {
        y[i] = a * x[i] + y[i];
    }
}

int main(int argc, char** argv)
{
    int len_glb, len_loc, i, a;
    len_glb = 1024;
    a = 2;
    int irank, nrank, igpu, ngpu, ierr;
    double *x_loc, *x_glb, *y_glb, *y_ref;
    double *y_glb_gpu;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);

    len_loc = len_glb / nrank;
    ngpu = acc_get_num_devices(acc_device_nvidia);
    igpu = irank % ngpu;
    acc_set_device_num(igpu, acc_device_nvidia);

    y_glb_gpu = malloc(len_glb * sizeof(double));
    x_loc = malloc(len_loc * sizeof(double));
    double* y_loc = malloc(len_loc * sizeof(double));
    
    if (irank == 0){
    x_glb = malloc(len_glb * sizeof(double));
    y_glb = malloc(len_glb * sizeof(double));
    y_ref = malloc(len_glb * sizeof(double));
    }

    if (irank == 0){
        for(i = 0; i < len_glb; i++){
            x_glb[i] = i;
            y_glb[i] = i + len_glb;
            y_ref[i] = a * i + (i+len_glb);
        }
    }

    MPI_Scatter(x_glb, len_loc, MPI_DOUBLE, x_loc, len_loc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(y_glb, len_loc, MPI_DOUBLE, y_loc, len_loc, MPI_DOUBLE, 0, MPI_COMM_WORLD);


#pragma acc enter data create(x_loc[0:len_loc], y_loc[0:len_loc])
#pragma acc enter data create(y_glb_gpu[0:len_glb])
#pragma acc update device(x_loc[0:len_loc], y_loc[0:len_loc])
#pragma acc update device(y_glb_gpu[0:len_glb])

saxpy(x_loc, y_loc, a, len_loc);

#pragma acc host_data use_device(y_glb_gpu, y_loc)
    MPI_Gather(y_loc, len_loc, MPI_DOUBLE, y_glb_gpu, len_loc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#pragma acc wait

#pragma acc update host(y_glb_gpu[0:len_glb])   

    if(irank == 0)
    {
            printf("y_glb[0] = %f\n", y_glb[0]);
            printf("y_glb[1] = %f\n", y_glb[1]);
            printf("y_glb[1023] = %f\n", y_glb[1023]);

            
             printf("y_glb_gpu[0] = %f\n", y_glb_gpu[0]);
             printf("y_glb_gpu[1] = %f\n", y_glb_gpu[1]);
             printf("y_glb_gpu[1023] = %f\n", y_glb_gpu[1023]);
             

            printf("y_ref[0] = %f\n", y_ref[0]);
            printf("y_ref[1] = %f\n", y_ref[1]);
            printf("y_ref[1023] = %f\n", y_ref[1023]);
    }

    free(x_loc);

    if(irank == 0){
        free(x_glb);
        free(y_glb);
        free(y_ref);
    }

    MPI_Finalize();
}
