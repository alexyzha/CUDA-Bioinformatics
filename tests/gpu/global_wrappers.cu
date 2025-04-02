#include "headers/global_wrappers.cuh"

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

__global__ void cu_uf_join(cu_union_find* UF, uint32_t* EDGE_LIST, size_t EDGE_COUNT) {
    int INDEX = (blockIdx.x * blockDim.x + threadIdx.x) * 2;
    if(INDEX + 1 >= EDGE_COUNT) {
        return;
    }
    __cu_uf_join(UF, EDGE_LIST[INDEX], EDGE_LIST[INDEX + 1]);
}

__global__ void cu_uf_find(cu_union_find* UF, uint32_t* INDEX_LIST, uint32_t* ROOTS, size_t LEN) {
    int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
    if(INDEX >= LEN) {
        return;
    }
    ROOTS[INDEX] = __cu_uf_find(UF, INDEX_LIST[INDEX]);
}

__global__ void cu_uf_con(cu_union_find* UF, uint32_t* EDGE_LIST, char* RET, size_t EDGE_COUNT, size_t RET_LEN) {
    int EDGE_INDEX = (blockIdx.x * blockDim.x + threadIdx.x) * 2;
    int RET_INDEX = blockIdx.x * blockDim.x + threadIdx.x;
    if(EDGE_INDEX + 1 >= EDGE_COUNT || RET_INDEX >= RET_LEN) {
        return;
    }
    RET[RET_INDEX] = __cu_uf_con(
        UF,
        EDGE_LIST[EDGE_INDEX],
        EDGE_LIST[EDGE_INDEX + 1]
    );
}