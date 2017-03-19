#define N (2048*2048)
#define THREADS_PER_BLOCK 512
__global__ void dot(int *a, int *b, int *c){
    __shared__ int tmp[THREADS_PER_BLOCK];
    int index = threadIdx.x + blockIdx.x*blockDim.x;
    tmp[threadIdx.x] = a[index] *b[index];
    __syncthreads();

    if(0 == threadIdx.x){
        int sum = 0;
        for(int i =0; i< THREADS_PER_BLOCK; i++){
            sum += tmp[i];
        }
        atomicAdd(c , sum);
    }
}
