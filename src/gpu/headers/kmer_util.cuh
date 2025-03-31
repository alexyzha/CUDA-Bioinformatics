#include "util_structs.cuh"

#ifndef CU_TYPE_MACROS
#define CU_TYPE_MACROS
#define ULL unsigned long long int
#endif

#ifndef CU_KMER_MACROS
#define CU_KMER_MACROS
#define K_MAX 32
#define MAP_MAX_INDICES 128
#define MAX_CLUSTER_SIZE 128
#endif

__device__ void cu_count_kmers(
    kh_pair<uint64_t>* MAP,
    char* ALL_SEQ,
    size_t* OFFSETS,
    size_t K,
    size_t LEN,
    size_t MAP_LEN
);

__device__ void cu_index_kmers(
    kh_pair<uint32_t[MAP_MAX_INDICES + 1]>* MAP,
    char* ALL_SEQ,
    size_t* OFFSETS,
    size_t K,
    size_t LEN,
    size_t MAP_LEN
);

__device__ void cu_get_kmer_overlaps(
    kh_pair<uint32_t[MAP_MAX_INDICES + 1]>* MAP,
    size_t MAP_LEN,
    uint32_t* EDGE_LIST,
    uint32_t* EDGE_COUNT,
    uint32_t MAX_EDGES
);

__device__ void cu_get_uf(
    cu_union_find* UF,
    size_t LEN,                         // Size of UF arrays, number of nodes^2
    size_t NODES,
    uint32_t* EDGE_LIST,
    uint32_t EDGE_COUNT
);

__device__ void cur_get_clusters{
    cu_union_find* UF,
    kh_pair<uint32_t[MAX_CLUSTER_SIZE + 1]> * MAP,
    size_t LEN,                         // Size of UF arrays, number of nodes^2
    size_t MAP_LEN,                     // = 2 * exp number of elements in MAP
    size_t K
};