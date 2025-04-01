#include "util_structs.cuh"

/*
 *  Smith waterman local alignment kernel
 *  @param ALL_SEQ `char*` flattened char array containing all sequences including reference
 *  @param CIGAR_BUF `char*` alloc buffer for trackback string rebuild to cigar. Exp length = `LEN * MAX_CIGAR_LEN`
 *  @param CACHE `int*` flattened DP cache
 *  @param OFFSETS `size_t*` size_t array containing end indices of all sequences
 *  @param LEN `size_t` number of sequences in `ALL_SEQ` including the reference sequence
 *  @param RET `int*` Exp length: `LEN * 3`
 *  @return Fills `int[i(+0:+2)]`, where `i` = `(blockIdx.x * blockDim.x + threadIdx.x) * 3`
 */
__global__ void cu_local_alignment(
    char* ALL_SEQ,
    char* CIGAR_BUF,
    int* CACHE,
    uint32_t* OFFSETS,
    size_t LEN,
    int* RET
);

/*
 *  Smith waterman local alignment kernel that takes cu_fq_read* as an input instead of a flattened sequence `char*` array
 *  @param REF `char*` reference sequence
 *  @param READS `cu_fq_read*` all reads in cu_fq_read format
 *  @param CIGAR_BUF `char*` alloc buffer for trackback string rebuild to cigar. Exp length = `LEN * MAX_CIGAR_LEN`
 *  @param CACHE `int*` flattened DP cache
 *  @param OFFSETS `size_t*` size_t array containing end indices of all sequences
 *  @param LEN `size_t` number of sequences in `ALL_SEQ` including the reference sequence
 *  @param REF_SIZE `size_t` reference sequence length
 *  @param RET `int*` Exp length: `LEN * 3`
 *  @return Fills `int[i(+0:+2)]`, where `i` = `(blockIdx.x * blockDim.x + threadIdx.x) * 3`
 */
__global__ void cu_fq_local_alignment(
    char* REF,
    cu_fq_read* READS,
    char* CIGAR_BUF,
    int* CACHE,
    uint32_t* OFFSETS,
    size_t LEN,
    size_t REF_SIZE,
    int* RET
);