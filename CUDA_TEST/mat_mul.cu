/**********************    mat_mul.cu    ******************************/
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "myCuBlas.h"

/*
#define M  256
#define P  128
#define N   64
#define thread_block_size 16
 */
#ifndef N
#error "matrix size required"
#endif
#ifndef thread_block_size
#error "thread block size required"
#endif

#define M (N) 
#define P (N) 

#define NUM_EXP (1)

#define ENABLE_CUBLAS

    dim3   dimBlock;
    dim3   dimGrid;


// timer
double t1, t2;
double t_upld, t_comp, t_dnld, t_total, t_comp_cublas;



//__global__ void beastDgemm( int K, double *dA, double *dB, double *dC);
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

/*run cuda code
* type
*   0: cuda suma   
*   1: cublas
*/
void run_cuda_summa(double *C, double *A, double *B, int n){
        double  *Ad, *Bd, *Cd;
        cudaEvent_t start_comp, stop_comp;
        cudaEvent_t start_upld, stop_upld;
        cudaEvent_t start_dnld, stop_dnld;
        cudaEventCreate(&start_comp);
        cudaEventCreate(&stop_comp);
        cudaEventCreate(&start_upld);
        cudaEventCreate(&stop_upld);
        cudaEventCreate(&start_dnld);
        cudaEventCreate(&stop_dnld);


        // allocate device memory 
        cudaMalloc(&Ad,(size_t)(M*P*sizeof(double)));
        cudaMalloc(&Bd,(size_t)(P*N*sizeof(double)));
        cudaMalloc(&Cd,(size_t)(M*N*sizeof(double)));

        // copy content to device memory
        cudaEventRecord(start_upld);
        cudaMemcpy(Ad,A,M*P*sizeof(double),cudaMemcpyHostToDevice);
        cudaMemcpy(Bd,B,P*N*sizeof(double),cudaMemcpyHostToDevice);
        cudaEventRecord(stop_upld);

        // start kernel 
        cudaEventRecord(start_comp);
        mat_mul<<<dimGrid,dimBlock>>>(Ad,Bd,Cd);
        cudaEventRecord(stop_comp);

        // copy C back
        cudaEventRecord(start_dnld);
        cudaMemcpy(C,Cd,M*N*sizeof(double),cudaMemcpyDeviceToHost);
        cudaEventRecord(stop_dnld);

        cudaEventSynchronize(stop_dnld);
        float time_eclapsed;
        cudaEventElapsedTime(&time_eclapsed, start_upld, stop_upld);
        t_upld += time_eclapsed/1000;
        printf("\tuplad completed in %fs \n", time_eclapsed/1000);

        cudaEventElapsedTime(&time_eclapsed, start_comp, stop_comp);
        t_comp += time_eclapsed/1000;
        printf("\tsumma computation in %fs \n", time_eclapsed/1000);

        cudaEventElapsedTime(&time_eclapsed, start_dnld, stop_dnld);
        t_dnld += time_eclapsed/1000;
        printf("\tdownload complete in  %fs \n", time_eclapsed/1000);

        cudaEventElapsedTime(&time_eclapsed, start_upld, stop_dnld);
        t_total += time_eclapsed/1000;
        printf("\t all complete in  %fs \n", time_eclapsed/1000);


        cudaFree(Ad);
        cudaFree(Bd);
        cudaFree(Cd);
}

void run_cublas(double *C, double *A, double *B, int n){
        double  *Ad, *Bd, *Cd;
        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);

        cublasHandle_t handle; 
        init_cublas(&handle);

        // allocate device memory 
        cudaMalloc(&Ad,(size_t)(M*P*sizeof(double)));
        cudaMalloc(&Bd,(size_t)(P*N*sizeof(double)));
        cudaMalloc(&Cd,(size_t)(M*N*sizeof(double)));
        printf("\tdeivce memory allocated \n");

        // copy content to device memory
        cudaMemcpy(Ad,A,M*P*sizeof(double),cudaMemcpyHostToDevice);
        cudaMemcpy(Bd,B,P*N*sizeof(double),cudaMemcpyHostToDevice);

    
       // cublas version */
        cudaEventRecord(start);
        gpu_blas_mmul(&handle, Cd, Ad, Bd, N);
        cudaEventRecord(stop);
        // copy C back
        cudaMemcpy(C,Cd,M*N*sizeof(double),cudaMemcpyDeviceToHost);

        cudaEventSynchronize(stop);
        float time_eclapsed;
        cudaEventElapsedTime(&time_eclapsed, start, stop);
        t_comp_cublas += time_eclapsed/1000;
        printf("\tcublas computation completed in %f\n", time_eclapsed/1000);

        cudaFree(Ad);
        cudaFree(Bd);
        cudaFree(Cd);
        finalize_cublas(&handle);
        printf("\tdevice buffer freed");
}


int main()
{
    int exp;

    double  *A, *B, *C, *D;

    // configuration here

    dimBlock = dim3(thread_block_size,thread_block_size);
    dimGrid = dim3(M/thread_block_size,N/thread_block_size);

    t_comp = 0;
    t_comp_cublas = 0;
    t_total = 0;
    
    for(exp =0; exp <NUM_EXP; exp++){
        printf("start experiemnt %d!\n",  exp);
        printf("matrix size %d, thread block size %d\n",N, thread_block_size);
        // cublas

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
        
        //gpu_blas_mmul(D, A, B, N);
        run_cublas(D, A, B, N);
        run_cuda_summa(C, A, B, N);

        if(verify(C, D, N) != 0){
            exit(-1);
        }
    }


    // verify results

    printf("***********************\n");
    printf("%d experiment finished \n", NUM_EXP);
    //t_transfer_1 = t_transfer_1/NUM_EXP; 
    //t_transfer_2 = t_transfer_2/NUM_EXP;
    t_comp = t_comp/NUM_EXP;
    //t_alloc = t_alloc/NUM_EXP;
    //t_free = t_free/NUM_EXP;
    t_comp_cublas  = t_comp_cublas/NUM_EXP;

    //t_total = t_alloc+t_transfer_1 + t_transfer_2 + t_comp + t_free;

    double n = N;
    double scalar = 2*n*n*n*(1E-9);
    double gflops = scalar/t_comp;
    double gflops_cublas = scalar/t_comp_cublas;
    //printf("\tm_size thrd_blk_size t_alloc t_transfer_1 t_transfer_2 t_comp t_free t_total gflops t_comp_blas blas_gflops\n");
    //printf("\t%d %d %f %f %f %f %f %f %f %f %f\n",N, thread_block_size, t_alloc, t_transfer_1, t_transfer_2,t_comp,t_free,t_total, gflops, t_comp_cublas, gflops_cublas);

    printf("\tm_size thrd_blk_size t_comp t gflops t_comp_blas blas_gflops t_uplad t_dnld t_total\n");
    printf("\t%d %d %lf %lf %lf %lf %lf %lf %lf\n",N, thread_block_size,t_comp, gflops, t_comp_cublas, gflops_cublas, t_upld, t_dnld,t_total);
    return 0;
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

        // p-th colblock will broadcast in row direction
        As[i][j] = Ad[(m*thread_block_size+i)*P+(p*thread_block_size+j)];
        // p-th rowblock will broadcast in column direction
        Bs[i][j] = Bd[(p*thread_block_size+i)*N+(n*thread_block_size+j)];

        // make sure local A and B are filled
        __syncthreads();

        // local matrix muliplication
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
