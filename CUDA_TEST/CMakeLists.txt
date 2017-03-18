cmake_minimum_required(VERSION 2.8)
find_package(CUDA QUIET REQUIRED)
project(cuda_example)

cuda_add_executable(
   dot.cu
   run_cuda.c
    )
