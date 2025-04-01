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
    const std::vector<fq_read*>& READS,
    char FILTER_MODE,
    char THRESH,
    size_t K,
    double PROPORTION
);

std::unordered_map<uint64_t, uint64_t> cu_count_kmers(
    const std::vector<fq_read*>& READS,
    size_t K
);

std::unordered_map<uint64_t, std::unordered_set<int>> cu_index_kmers(
    const std::vector<fq_read*>& READS,
    size_t K
);

std::vector<std::unordered_set<int>*> cu_cluster_by_kmer(
    const std::vector<fq_read*>& READS, 
    size_t K, 
    size_t THRESH
);

