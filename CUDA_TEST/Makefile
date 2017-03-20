all:
	nvcc -o mat_mul mat_mul.cu myCuBlas.cu -lcublas
test:
	nvcc -o my_test dot.cu run_cuda.cu
