#include "cu_util.h"
#include "cu_kmer_hash.h"

#ifndef CU_TYPE_MACROS
#define CU_TYPE_MACROS
#define ULL unsigned long long int
#endif

#ifndef CU_KMER_MACROS
#define CU_KMER_MACROS
#define K_MAX 32
#define MAP_MAX_INDEXES 128
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
    kh_pair<uint32_t[MAP_MAX_INDEXES + 1]>* MAP,
    char* ALL_SEQ,
    size_t* OFFSETS,
    size_t K,
    size_t LEN,
    size_t MAP_LEN
);