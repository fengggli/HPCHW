#include "cublas_v2.h"

void gpu_blas_mmul(cublasHandle_t *handle, double *C, double *A, double *B, int n);

void init_cublas(cublasHandle_t * p_handle);

void finalize_cublas(cublasHandle_t *p_handle);

