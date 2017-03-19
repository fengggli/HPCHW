cmake_minimum_required(VERSION 2.8)
project(cuda_example)
find_package(CUDA QUIET REQUIRED)

cuda_add_executable(run_cuda dot.cu)
