#include "fq_filter.cuh"
#include "kmer_util.cuh"
#include "local_alignment.cuh"
#include "util_structs.cuh"
#include "util.cuh"

__global__ gpu_filter_reads(
    char* ALL_SEQ,
    uint32_t* OFFSETS,
    size_t LEN,
    char FILTER_MODE,
    char THRESH,
    uint64_t* FILTER_MASK,
    double PROPORTION = 0.0
);

