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

/**
 *  Wrapper for fastq filtering kernel.
 *  @param READS `vector<fq_read*>&` fastq reads
 *  @param FILTER_MODE `char` pick from `FASTX_FILTER_MACROS`
 *  @param THRESH `char` PHRED score threshold for discard
 *  @param K `size_t` size of sliding window if chosen as filter mode
 *  @param PROPORTION `double` additional threshold for proportion if proportion-based filter mode is chosen.
 */
std::vector<fq_read*> cu_filter_fq(
    const std::vector<fq_read*>& READS,
    char FILTER_MODE,
    char THRESH,
    size_t K,
    double PROPORTION
);

/**
 *  Wrapper for kmer counting kernel.
 *  @param READS `vector<fq_read*>&` fastq reads
 *  @param K `size_t` kmer length
 *  @return `unordered_map<uint64_t, uint64_t>` map<kmer, count>.
 */
std::unordered_map<uint64_t, uint64_t> cu_count_kmers(
    const std::vector<fq_read*>& READS,
    size_t K
);

/**
 *  Wrapper for kmer indexing kernel.
 *  @param READS `vector<fq_read*>&` fastq reads
 *  @param K `size_t` kmer length
 *  @return `unordered_map<uint64_t, unordered_set<int>>` map<kmer, list of indexes>.
 */
std::unordered_map<uint64_t, std::unordered_set<int>> cu_index_kmers(
    const std::vector<fq_read*>& READS,
    size_t K
);

/**
 *  Wrapper for kmer-based clustering kernel.
 *  @param READS `vector<fq_read*>&` fastq reads
 *  @param K `size_t` kmer length
 *  @param THRESH `size_t` number of overlaps needed to cluster
 *  @return `vector<unordered_set<int>*>` list (see below)
 *  @note `list[i]` contains of sets of indexes rooted at `i`.
 */
std::vector<std::unordered_set<int>*> cu_cluster_by_kmer(
    const std::vector<fq_read*>& READS, 
    size_t K, 
    size_t THRESH
);

/**
 *  Wrapper for local alignment kernel.
 *  @param REF `string&` reference sequence
 *  @param READS `vector<fq_read*>&` reads to align to reference sequence
 *  @return `vector<cu_alignment*>`
 *  @note `cu_alignment` holds string* cigar instead of aligned_ref/read like the CPU version, `alignment`.
 */
std::vector<cu_alignment*> cu_local_align(
    const std::string& REF, 
    const std::vector<fq_read*>& READS
);