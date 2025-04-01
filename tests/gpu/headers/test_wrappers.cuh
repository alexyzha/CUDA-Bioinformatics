#include "test_util.cuh"

__global__ void cu_max(int a, int b, int* out);

__global__ void cu_base2bit(char* out);

__global__ void cu_hashing(uint64_t* out, uint64_t* key, size_t LEN);