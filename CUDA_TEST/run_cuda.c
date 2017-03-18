#include <stdio.h>

#define N (2048*2048)
#define THREADS_PER_BLOCK 512

int main(){
    int *a, *b, *c;
    int *dev_a, *dev_b, *dev_c;
    int size = N*sizeof(int);

    cudaMalloc((void **)dev_a, size);
    cudaMalloc((void **)dev_b, size);
    cudaMalloc((void **)dev_c, sizeof(int));

    a = (int *)malloc(size);
    b = (int *)malloc(size);
    c = (int *)malloc(sizeof(int));

    // initialize a, b, c
    cudaMemcpy(dev_a, a, size, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, size, cudaMemcpyHostToDevice);

    // lanch
    dot<<<N/THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(dev_a, dev_b, dev_c);

    // copy back
    cudaMemcpy(c, dev_c, sizeof(int), cudaMemcpyDeviceToHost);

    free(a);
    free(b);
    free(c);

    cudaFree(dev_a);
    cudaFree(dev_b);
    cudeFree(dev_c);

    return 0;
}
