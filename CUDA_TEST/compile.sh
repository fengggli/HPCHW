#!/bin/bash

for N in 256 512 1024 2048 4096 8192
do
    for thread_block_size in 1 4 8 16 32
    do
	    nvcc -D N=${N} -D thread_block_size=${thread_block_size} -o bin/mat_mul_N_${N}_thrd_${thread_block_size} mat_mul.cu myCuBlas.cu -lcublas 
    done
done
