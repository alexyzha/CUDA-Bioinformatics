#include "util.cuh"
#include "util_structs.cuh"

__global__ void cu_count_kmers(
    kh_pair<uint64_t>* MAP,
    char* ALL_SEQ,
    size_t* OFFSETS,
    size_t K,
    size_t LEN,
    size_t MAP_LEN
);

__global__ void cu_index_kmers(
    kh_pair<uint32_t[MAP_MAX_INDICES + 1]>* MAP,
    char* ALL_SEQ,
    size_t* OFFSETS,
    size_t K,
    size_t LEN,
    size_t MAP_LEN
);

__global__ void cu_get_kmer_overlaps(
    kh_pair<uint32_t[MAP_MAX_INDICES + 1]>* MAP,
    size_t MAP_LEN,
    uint32_t* EDGE_LIST,
    uint32_t* EDGE_COUNT,
    uint32_t MAX_EDGES
);

__global__ void cu_get_uf(
    cu_union_find* UF,
    size_t LEN,                         // Size of UF arrays, number of nodes^2
    size_t NODES,
    uint32_t* EDGE_LIST,
    uint32_t EDGE_COUNT
);

__global__ void cu_get_clusters(
    cu_union_find* UF,
    kh_pair<uint32_t[MAX_CLUSTER_SIZE + 1]>* MAP,
    size_t LEN,                         // Size of UF arrays, number of nodes^2
    size_t MAP_LEN,                     // = 2 * exp number of elements in MAP
    size_t K
);