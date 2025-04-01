#include "util_structs.cuh"

/*
 *  Returns bitmask [1 = keep][0 = discard] based on `FILTER_MODE`
 *  @param ALL_SEQ `char*` all `LEN` sequences in order, flattened
 *  @param OFFSETS `uint32_t*` array of offsets, `OFFSETS[i]` = end index of seq `i`
 *  @param LEN `size_t` number of sequences
 *  @param FILTER_MODE `char` how to filter
 *  @param THRESH `char` filter threshold
 *  @param FILTER_MASK `uint64_t*` results
 *  @param PROPORTION `double` default = 0.0, used for filter modes with proportions
 *  @return `void`
 */
__global__ void cu_filter_reads(
    char* ALL_SEQ,
    uint32_t* OFFSETS,
    size_t LEN,
    char FILTER_MODE,
    char THRESH,
    uint64_t* FILTER_MASK,
    double PROPORTION = 0.0
);

/*
 *  Modifies ALL_SEQ, denotes end of trimmed seq with '\0`
 *  @param ALL_SEQ `char*` all `LEN` sequences in order, flattened
 *  @param OFFSETS `uint32_t*` array of offsets, `OFFSETS[i]` = end index of seq `i`
 *  @param LEN `size_t` number of sequences
 *  @param K `size_t` size of sliding window
 *  @param THRESH `char` filter threshold
 *  @param PROPORTION `double` default = 0.0
 *  @return `void`
 */
__global__ void cu_filter_reads_sw(
    char* ALL_SEQ,
    uint32_t* OFFSETS,
    size_t LEN,
    size_t K,
    char THRESH,
    double PROPORTION = 0.5
);