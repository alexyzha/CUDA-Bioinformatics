#include <stdio.h>

__global__ void hello() {
    printf("B: %u, T: %u\n",blockIdx.x,threadIdx.x);
}

int main(int argc, char* argv[]) {

    printf("INIT\n");
    
    hello<<<7,7>>>();
    
    cudaError_t err = cudaDeviceSynchronize();
    if(err != cudaSuccess) {
        fprintf(stderr, "CUDA ERR: %s\n", cudaGetErrorString(err));
        return -1;
    }

    return 0;
}