#include "fq_filter.cuh"
#include "kmer_util.cuh"
#include "local_alignment.cuh"
#include "util_structs.cuh"
#include "util.cuh"
#include "../../src/cpu/headers/fq_read.h"

#ifndef CU_KERNEL_PARAMS
#define CU_KERNEL_PARAMS
#define MAX_THREADS 256
#define MAX_BLOCKS 0
#endif

std::vector<fq_read*> cu_filter_fq(
    std::vector<fq_read*> READS,
    char FILTER_MODE,
    char THRESH,
    size_t K,
    double PROPORTION
);