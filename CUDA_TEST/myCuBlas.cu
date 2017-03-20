#include "myCuBlas.h"
// use cublas to calculate dgemm 
//Note,
//  this function is orginally from:
//  https://solarianprogrammer.com/2012/05/31/matrix-multiplication-cuda-cublas-curand-thrust/
void gpu_blas_mmul(cublasHandle_t *p_handle, double *C, double *A, double *B, int n) {
    //int lda=n,ldb=n,ldc=n;
    int m,k;
    m = n; k = n;
    const double alf = 1;
    const double bet = 0;
    const double *alpha = &alf;
    const double *beta = &bet;

    // Do the actual multiplication
    //cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    cublasDgemm(*p_handle, CUBLAS_OP_N, CUBLAS_OP_N, n, m, k, alpha, B, n, A, k, beta, C, n);

}

void init_cublas(cublasHandle_t * p_handle){
    cublasCreate(p_handle);
    
}
void finalize_cublas(cublasHandle_t *p_handle){
    // Destroy the handle
    cublasDestroy(*p_handle);
}
    


