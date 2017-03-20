/**********************    mat_mul.cu    ******************************/
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "cublas_v2.h"

/*
#define M  256
#define P  128
#define N   64
#define thread_block_size 16
 */

#define M  1024
#define P  1024
#define N  1024

#define WARP (32)
#define thread_block_size 32

#define KERNEL_1

__global__ void beastDgemm( int K, double *dA, double *dB, double *dC);
__global__ void mat_mul(double *Ad, double *Bd, double *Cd);

double get_cur_time() {
    struct timeval   tv;
    struct timezone  tz;
    double cur_time;

    gettimeofday(&tv, &tz);
    cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;

    return cur_time;
}

/*
 * generate a random matrix with double value -1~1
 * init_type
 *  1: init with random number (-1~1)
 *  0: init with 0
 */
int init_matrix(double ** m, int n, int init_type){
    int i, j;
    double tmp;
    srand(time(NULL));
    // allocate space
    *m = (double*)malloc(sizeof(double)*n*n);
    if(*m ==NULL){
        perror("allocate error");
        exit(-1);
    }

    // initiate with random number
    for( i = 0; i< n; i++){
        for(j = 0; j < n ; j++){

            if(init_type == 1){
                tmp = (double)rand()/RAND_MAX*2.0-1.0;
            }
            else{
                tmp = 0;
            }
            (*m)[i*n + j] = tmp;
        }
    }
    return 0;
}

// verify results usig blas
int verify(double *C, double *C_blas, int n){
    int i, j;

    for(i = 0; i < n; i++){
        for(j = 0; j< n; j++){
            if(abs(C[i*n +j] -  C_blas[i*n+j]) > 0.000001){
                printf("element(%d, %d)%3.12lf != %3.12lf\n",i,j, C[i*n +j], C_blas[i*n +j]);
                return -1;
            }
        }
    }
    return 0;
}
void gpu_blas_mmul(double *C, double *A, double *B, int n) {
    //int lda=n,ldb=n,ldc=n;
    int m,k;
    m = n; k = n;
    const double alf = 1;
    const double bet = 0;
    const double *alpha = &alf;
    const double *beta = &bet;

    // Create a handle for CUBLAS
    cublasHandle_t handle;
    cublasCreate(&handle);

    // Do the actual multiplication
    //cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, m, k, alpha, B, n, A, k, beta, C, n);

    // Destroy the handle
    cublasDestroy(handle);
}

int main()
{
    double  *A, *B, *C, *D;
    double  *Ad, *Bd, *Cd;

#ifndef KERNEL_1
    dim3 dimGrid(M/WARP, N/WARP);
    dim3 dimBlock(WARP);
#else
    dim3   dimBlock(thread_block_size,thread_block_size);
    dim3   dimGrid(M/thread_block_size,N/thread_block_size);
#endif

    // timer
    double t1, t2;
    double t_transfer_1, t_transfer_2, t_comp, t_alloc, t_total;
    double t_comp_cublas;



    printf("start!\n");
    if(init_matrix(&B, N, 1) == 0){
        printf("\trandom matrix B is generated \n");
    }
    if(init_matrix(&A, N, 1) == 0){
        printf("\trandom matrix A is generated \n");
    }

    // cuda results
    if(init_matrix(&C, N, 0) == 0){
        printf("\t C is init with 0\n");
    }

    // naiive results
    if(init_matrix(&D, N, 0) == 0){
        printf("\t D is init with 0\n");
    }


    /* First, compute D=AB on the host CPU. */
    t1 = get_cur_time();
    /*
    for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
            D[i*N+j]=0;
            for(k=0;k<N;k++) {
                D[i*N+j] += A[i*N+k]*B[k*N+j];
            }
        }
    }
    */

    gpu_blas_mmul(D, A, B, N);
    t2 = get_cur_time();
    t_comp_cublas = t2-t1;

    
    //gpu_blas_mmul(D, A, B, N);

    // allocate device memory 
    t1 = get_cur_time();
    cudaMalloc(&Ad,(size_t)(M*P*sizeof(double)));
    cudaMalloc(&Bd,(size_t)(P*N*sizeof(double)));
    cudaMalloc(&Cd,(size_t)(M*N*sizeof(double)));
    t2 = get_cur_time();
    t_alloc = t2-t1;
    printf("\tdeivce memory allocated\n");

    // copy content to device memory
    t1 = get_cur_time();
    cudaMemcpy(Ad,A,M*P*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(Bd,B,P*N*sizeof(double),cudaMemcpyHostToDevice);
    t2 = get_cur_time();
    t_transfer_1 = t2-t1;
    printf("\tdata copied to deivce memory\n");

    t1 = get_cur_time();
    // start kernel 
#ifndef KERNEL_1
    beastDgemm<<<dimGrid,dimBlock>>>(N, Ad,Bd,Cd);
#else
    mat_mul<<<dimGrid,dimBlock>>>(Ad,Bd,Cd);
#endif
    t2 = get_cur_time();
    t_comp = t2 -t1;
    printf("\tcuda computation completed\n");

    // copy C back
    t1 = get_cur_time();
    cudaMemcpy(C,Cd,M*N*sizeof(double),cudaMemcpyDeviceToHost);
    t2 = get_cur_time();
    t_transfer_2 = t2-t1;
    printf("\tdata copied back to host memory\n");

    // free divice memory
    cudaFree(Ad);
    cudaFree(Bd);
    cudaFree(Cd);

    // verify results
    printf("***********************\n");
    printf("finished\n");
    if(verify(C, D, N) == 0){
        double n = N;
        double scalar = 2*n*n*n*(1E-9);
        double gflops = scalar/t_comp;
        t_total = t_transfer_1 + t_transfer_2 + t_comp;
        printf("\tt_alloc\tt_transfer_1\tt_transfer_2\tt_comp\tt_total\tgflops\n");
        printf("\t%f\t%f\t%f\t%f\t%f\t%f\n",t_alloc, t_transfer_1, t_transfer_2,t_comp,t_total, gflops);
        printf("\tcublas computation time %f\n", t_comp_cublas);
    }
    else{
        printf("\tnot correct\n");
    }
}

__global__ void beastDgemm( int K, double *dA, double *dB, double *dC){
    int n, k, kk;
    int mm = blockIdx.x*WARP;
    int nn = blockIdx.y*WARP;
    int idx = threadIdx.x;
    __shared__ float sB[WARP];
    float rA, rC[WARP];

#pragma unroll
    for(n=0; n < WARP; n++){
        rC[n] = 0;
    }

    // each thread computes 1 row
    for(kk = 0; kk < K; kk += WARP){
#pragma unroll
        for(k = 0; k < WARP; k++){
            rA = dA[mm+idx + (k+kk)*K];
            sB[idx] = dB[nn+idx + (k+kk)*K];

#pragma unroll
            for(n = 0 ; n < WARP; n++)
                rC[n] += rA*sB[n];
        }
    }

#pragma unroll
    for(n = 0; n< WARP;n++)
        dC[mm+idx + (n+nn)*K] = rC[n];
}

/*
 * SUMMA using CUDA
 * this function is originally from Bigred website
 */
__global__ void mat_mul(double *Ad, double *Bd, double *Cd) {
    int    m = blockIdx.x;
    int    n = blockIdx.y;
    int    i = threadIdx.x;
    int    j = threadIdx.y;
    int    k,p;
    double  c = 0.0;

    // each thread block will share the local A and B
    __shared__  double As[thread_block_size][thread_block_size];
    __shared__  double Bs[thread_block_size][thread_block_size];

    for(p=0;p<P/thread_block_size;p++) {
        // this is the broadcast step,
        // p-th colblock will broadcast in row direction;
        // p-th rowblock will broadcast in column direction
        As[i][j] = Ad[(m*thread_block_size+i)*P+(p*thread_block_size+j)];
        Bs[i][j] = Bd[(p*thread_block_size+i)*N+(n*thread_block_size+j)];

        // make sure local A and B are filled
        __syncthreads();

        // local matrix multiplication
        for(k=0; k<thread_block_size; k++) {
            c += As[i][k] * Bs[k][j];
        }

        // this is must!
        __syncthreads();
    }
    // fill the device memory
    Cd[(m*thread_block_size+i)*N+(n*thread_block_size+j)] = c;
}
/**********************************************************************/
