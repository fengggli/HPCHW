#ifndef MATRIX_OP_H
#define MATRIX_OP_H
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
