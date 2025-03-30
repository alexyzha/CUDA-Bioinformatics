#include "headers/cu_util.h"
#include "headers/cu_alignment.h"
#include "headers/cu_fq_read.h"

#ifndef CU_LOCAL_ALIGNMENT_MACROS
#define CU_LOCAL_ALIGNMENT_MACROS
#define MAX_CIGAR_LEN 340
#define MAX_REF_LEN 200
#define MAX_READ_LEN 170
#define ROWS (MAX_REF_LEN + 1)
#define COLS (MAX_READ_LEN + 1)
#define DP(i, j) cache[(i) * COLS + (j)]
#define ALIGN_MATCH 2
#define ALIGN_MISMATCH -1
#define ALIGN_GAP -2
#endif

/*
 *  @param ALL_SEQ `char*` flattened char array containing all sequences including reference
 *  @param CIGAR_BUF `char*` alloc buffer for trackback string rebuild to cigar. Exp length = `LEN * MAX_CIGAR_LEN`
 *  @param CACHE `int*` flattened DP cache
 *  @param OFFSETS `size_t*` size_t array containing end indices of all sequences
 *  @param LEN `size_t` number of sequences in `ALL_SEQ` including the reference sequence
 *  @param RET `cu_alignment*` Exp length: `LEN`
 *  @return Fills `cu_alignment[i]`, where `i` = `blockIdx.x * blockDim.x + threadIdx.x`
 */
__device__ void local_alignment(
    char* ALL_SEQ,
    char* CIGAR_BUF,
    int* CACHE,
    size_t* OFFSETS,
    size_t LEN,
    cu_alignment* RET
);

__device__ void fq_local_alignment(
    char* REF,
    cu_fq_read* READS,
    char* CIGAR_BUF,
    int* CACHE,
    size_t* OFFSETS,
    size_t LEN,
    size_t REF_SIZE,
    cu_alignment* RET
);