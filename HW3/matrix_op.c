/*
 * matrix operation collection
 */
#include "matrix_op.h"

void get_col(double *A, double*col, int n, int l ){
    int i,j;
    for(i = 0; i <n; i++){
        for(j = 0; j < n; j++){
            if(j == l)
                *(col++) = A[i*n+j];
        }
    }
}

void get_row(double *A, double*row, int n, int r ){
    int j;
    for(j = 0; j < n; j++){
        *(row++) = A[r*n+j];
    }
}


void update_local_C(double *C, double *col, double*row, int n){
    int i, j;
    for(i = 0; i<n;i++){
        for(j = 0; j <n ;j++){
            C[i*n+j] += col[i]*row[j];
        }
    }
    // cal blas here
}



