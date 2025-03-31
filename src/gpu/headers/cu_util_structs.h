#include "cu_util.h"

#ifndef C_CU_ALIGNMENT
#define C_CU_ALIGNMENT

/*
 *  Prettier format to return local alignment results in
 *  @param cigar `char*` pointer to cigar string
 *  @param end_ref `int` end index of reference
 *  @param end_read `int` end index of read
 *  @param score `int` max smith waterman alignment score
 */
typedef struct {
    char* cigar;
    int end_ref;
    int end_read;
    int score;
} cu_alignment;

#endif

#ifndef C_CU_FQ_READ
#define C_CU_FQ_READ

/*
 *  1:1 container for fq_read from CPU code
 *  @param id `char*` sequence id
 *  @param seq `char*` sequence
 *  @param qual `char*` quality
 *  @param metadata `char*` metadata
 *  @param size `size_t` length of `seq` and `qual`
 */
typedef struct {
    char* id;
    char* seq;
    char* qual;
    char* metadata;
    size_t size;
} cu_fq_read;

#endif

#ifndef C_CU_UF
#define C_CU_UF

/*
 *  Union find container
 *  @param p `uint32_t*` parent
 *  @param h `uint32_t*` height of root
 *  @param LEN `size_t` length of `p` and `h`
 */
typedef struct {
    uint32_t* p;
    uint32_t* h;
    size_t LEN;
} cu_union_find;

cu_union_find* cu_uf_construct(int n);

__device__ int cu_uf_find(cu_union_find* UF, int x);

__device__ void cu_uf_join(cu_union_find* UF, int x, int y);

__device__ bool cu_uf_con(cu_union_find* UF, int x, int y);

#endif

#ifndef C_KMER_HASH_TABLE
#define C_KMER_HASH_TABLE
#define EMPTY 0

/*
 *  Hash table key-value pair for kmer counting/indexing
 *  @param key `uint64_t` kmer hash
 *  @param value `T` 
 */
template<typename T>
struct kh_pair {
    uint64_t key;
    T value;
};

__device__ uint64_t kh_hash(uint64_t key);

template<typename T>
kh_pair<T>* kh_construct(int n) {
    // Malloc
    kh_pair<T>* map;
    cudaMalloc(&map, sizeof(kh_pair<T>) * n);

    // Init all as empty
    cudaMemset(map, 0, sizeof(kh_pair<T>) * n);
    return map;
}

/*
*  This is a simple linear probling template for inserting into a kh_pair hashtable.
*/
template<typename T>
__device__ void kh_insert(kh_pair<T>* map, const kh_pair<T>* pairs, size_t LEN, size_t MAP_SIZE) {
    // Thread index OOB check
    size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= LEN) {
        return;
    }

    // Get hash
    uint64_t key = pairs[LEN].key;
    T value = pairs[LEN].value;
    uint64_t index = kh_hash(key);

    // Linear probing
    for(int i = 0; i < MAP_SIZE; ++i) {
        int cur = (index + i) % (MAP_SIZE - 1);
        int prev = atomicCAS(&map[cur].key, EMPTY, key);
        if(prev == key || prev == EMPTY) {
            map[index].value = value;
            return;
        }
    }

    // Full/heavy collision
    // WAH WAH CERR CERR THROW STDS
};

#endif