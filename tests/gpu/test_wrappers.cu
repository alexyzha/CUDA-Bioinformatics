#include "headers/test_wrappers.cuh"

__global__ void cu_max(int a, int b, int* out) {
    *out = __max(a, b);
}

__global__ void cu_base2bit(char* out) {
    int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
    if(INDEX > 255) {
        return;
    }
    out[INDEX] = __base_to_bit(static_cast<unsigned char>(INDEX));
}

__global__ void cu_hashing(uint64_t* out, uint64_t* key, size_t LEN) {
    int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
    if(INDEX >= LEN) {
        return;
    }
    out[INDEX] = __kh_hash(key[INDEX]);
}

__global__ void cu_uf_insert(cu_union_find* UF, uint32_t EDGE_LIST, size_t EDGE_COUNT) {
    int INDEX = (blockIdx.x * blockDim.x + threadIdx.x) * 2;
    if(INDEX >= EDGE_COUNT) {
        return;
    }
    __cu_uf_join(UF, EDGE_LIST[INDEX], EDGE_LIST[INDEX + 1]);
}