#include <stdio.h>

__global__ hello() {
    printf("B: %u, T: %u\n",blockIdx.x,threadIdx.x);
}

int main(int argc, char argv[]) {

    hello<<<7,7>>>();
    cudaDeviceSynchronize();

    return 0;
}