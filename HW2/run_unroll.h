#ifndef SEQ_MATRIX_MUL_H
#define SEQ_MATRIX_MUL_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <papi.h>
#include <cblas.h>

#include "papi_timer.h"
#define NUM_EXP (3)

#define NUM_EVENTS (5)


//#define VERBOSE


/*
 * calculate matrix multiplication using order: iki with loop rurolling
 * input:
 *  matrix C, A, B, all in flat format
 *  n: matrix size
 *  l: loop times
 *
 * output:
 *  matrix C
 */



void seq_cal_ikj_unroll(double *C, double *A, double *B, int n, int b, int l);
#endif
