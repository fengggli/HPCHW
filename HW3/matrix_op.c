/*
 * matrix operation collection
 */
#include "matrix_op.h"

inline void get_col(double *A, double*col, int n, int l, int b ){
    int i,j;
    for(i = 0; i <n; i++){
        for(j = 0; j < n; j++){
            if(j/b == l)
                *(col++) = A[i*n+j];
        }
    }
}

inline void get_row(double *A, double*row, int n, int r, int b){
    int j;
    for(j = 0; j < n*b; j++){
        *(row++) = A[r*n*b+j];
    }
}


inline void update_local_C(double *C, double *col, double*row, int m, int n, int b){

    /*
    int i, j;
    for(i = 0; i<n;i++){
        for(j = 0; j <n ;j++){
            C[i*n+j] += col[i]*row[j];
        }
    }
    */
        int k;
        k = b;

        double  alpha =1;
        double  beta = 1;
        char transfer_A = 'N';
        char transfer_B = 'N';

        dgemm_(&transfer_B, &transfer_A, &n, &m,&k, &alpha, row, &n, col, &k, &beta, C,&n);
}



