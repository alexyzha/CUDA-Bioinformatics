#include "test_util.h"

__global__ void cu_max(int a, int b, int* out);

__global__ void cu_base2bit(char* out);

__global__ void cu_hashing(uint64_t* out, uint64_t* key, size_t LEN);

__global__ void cu_uf_join(cu_union_find* UF, uint32_t* EDGE_LIST, size_t EDGE_COUNT);

__global__ void cu_uf_find(cu_union_find* UF, uint32_t* INDEX_LIST, uint32_t* ROOTS, size_t LEN);

__global__ void cu_uf_con(cu_union_find* UF, uint32_t* EDGE_LIST, char* RET, size_t EDGE_COUNT, size_t RET_LEN);