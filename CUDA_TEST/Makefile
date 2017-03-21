all:
	nvcc -D N=1024 -D thread_block_size=16 -o mat_mul mat_mul.cu myCuBlas.cu -lcublas
test:
	nvcc  -o mat_mul mat_mul.cu myCuBlas.cu -lcublas
