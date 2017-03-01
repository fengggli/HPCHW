#ifndef MATRIX_OP_H
#define MATRIX_OP_H
#define DGEMM_ROWMAJOR(A,B,C,m,n,k,alpha,beta,transf_A,transf_B, lda, ldb, ldc) \
            dgemm_(transf_B, transf_A, n, m, k, alpha, B, ldb, A, lda, beta, C, ldc)
extern void dgemm_(char * transa, char * transb, int * m, int * n, int * k,
              double * alpha, double * A, int * lda,
              double * B, int * ldb, double * beta,
              double *, int * ldc);



/*
 * matrix operation functions
 */

// get the l-th col from size-n matrx A, saved that col in col;
void get_col(double *A, double*col, int n, int l );

// get the i-th row from size-n matrix A, save into vector row
void get_row(double *A, double*row, int n, int r );

// update local C
void update_local_C(double *C, double *col, double*row, int n);

#endif
