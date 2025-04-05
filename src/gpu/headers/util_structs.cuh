#include "util.cuh"

/*****************************************************
 *           ALIGNMENT STRUCT DECLARATION            *
 *****************************************************/
#ifndef C_CU_ALIGNMENT
#define C_CU_ALIGNMENT

/**
 *  Prettier format to return local alignment results in.
 *  @param cigar `string*` pointer to cigar string
 *  @param end_ref `int` end index of reference
 *  @param end_read `int` end index of read
 *  @param score `int` max smith waterman alignment score
 *  @note Modified version of struct `alignment`, which has aligned_ref/read instead of cigar.
 */
typedef struct {
public:
    std::string* cigar;
    int end_ref;
    int end_read;
    int score;

} cu_alignment;

#endif

/*****************************************************
 *            FQ_READ STRUCT DECLARATION             *
 *****************************************************/
#ifndef C_CU_FQ_READ
#define C_CU_FQ_READ

/**
 *  1:1 container for fq_read from CPU code.
 *  @param id `char*` sequence id
 *  @param seq `char*` sequence
 *  @param qual `char*` quality
 *  @param metadata `char*` metadata
 *  @param size `size_t` length of `seq` and `qual`
 *  @note All strings are represented as `char*`.
 */
typedef struct {
public:
    char* id;
    char* seq;
    char* qual;
    char* metadata;
    size_t size;

} cu_fq_read;

#endif

/*****************************************************
 *           UNION FIND STRUCT DECLARATION           *
 *****************************************************/
#ifndef C_CU_UF
#define C_CU_UF

/**
 *  Union find container.
 *  @param p `uint32_t*` parent
 *  @param h `uint32_t*` height of root
 *  @param LEN `size_t` length of `p` and `h`
 *  @note This struct is meant to be passed as a device pointer.
 */
typedef struct {
public:
    uint32_t* p;
    uint32_t* h;
    size_t LEN;

} cu_union_find;

/**
 *  Union find constructor.
 *  @param n `int` number of desired nodes
 *  @return `cu_union_find*` device pointer
 *  @note **THIS RETURNS A DEVICE POINTER** 
 */
cu_union_find* cu_uf_construct(int n);

/**
 *  Union find destructor.
 *  @param d_uf `cu_union_find*` device pointer to union find being deleted
 *  @note **`d_uf` EXPECTS A DEVICE POINTER**
 */
void cu_uf_destruct(cu_union_find* d_uf);

/**
 *  Union find find root of `x` - thread safe.
 *  @param UF `cu_union_find*` pointer to union find where find will take place
 *  @param x `int` query
 *  @return `int` root of tree containing `x`
 *  @note **`UF` EXPECTS A DEVICE POINTER**
 */
__device__ int __cu_uf_find(cu_union_find* UF, int x);

/**
 *  Union find join/"union" `x` and `y` - thread safe.
 *  @param UF `cu_union_find*` pointer to union find where find will take place
 *  @param x `int` first node to be joined
 *  @param y `int` second node to be joined
 *  @return `void`
 *  @note **`UF` EXPECTS A DEVICE POINTER**
 */
__device__ void __cu_uf_join(cu_union_find* UF, int x, int y);

/**
 *  Union find check if `x` and `y` are joined - thread safe.
 *  @param UF `cu_union_find*` pointer to union find where find will take place
 *  @param x `int` first query
 *  @param y `int` second query
 *  @return `bool` Joined ? true : false
 *  @note **`UF` EXPECTS A DEVICE POINTER**
 */
__device__ bool __cu_uf_con(cu_union_find* UF, int x, int y);

#endif

/*****************************************************
 *         KMER HASHTABLE STRUCT DECLARATION         *
 *****************************************************/
#ifndef C_KMER_HASH_TABLE
#define C_KMER_HASH_TABLE
#define EMPTY -1

/**
 *  Hash table key-value pair for kmer counting/indexing.
 *  @param key `uint64_t` kmer hash
 *  @param value `T` 
 */
template<typename T>
struct kh_pair {
    uint64_t key;
    T value;
};

/**
 *  Convenience function for hashing 64-bit keys.
 *  @param key `uint64_t` 64-bit encoded kmer as key
 *  @return `uint64_t` hashed key
 *  @note **THIS IS DEVICE CODE**
 */
__device__ uint64_t __kh_hash(uint64_t key);

/**
 *  Kmer hashtable constructor.
 *  @param n `int` max number of keys
 *  @note In order to avoid collisions, the true size of the hashtable will be `2*n`.
 */
template<typename T>
kh_pair<T>* kh_construct(int n) {
    // Malloc
    kh_pair<T>* map;
    CUDA_CHECK(cudaMalloc(&map, sizeof(kh_pair<T>) * n));

    // Init all as empty
    CUDA_CHECK(cudaMemset(map, EMPTY, sizeof(kh_pair<T>) * n));
    return map;
}

/**
 *  This is a simple linear probling template for inserting into a kh_pair hashtable.
 *  @param map `kh_pair<T>*` hashtable
 *  @param pairs `kh_pair<T>*` key-value pairs
 *  @param LEN `size_t` length of `pairs`
 *  @param MAP_SIZE `size_t` length of `map`
 *  @return `void`
 *  @note **THIS FUNCTION EXPECTS DEVICE POINTERS**
 *  @note This function is not meant to be used. It is just a template for my reference.
 */
template<typename T>
__device__ void kh_insert(kh_pair<T>* map, const kh_pair<T>* pairs, size_t LEN, size_t MAP_SIZE) {
    // Thread index OOB check
    size_t INDEX = blockIdx.x * blockDim.x + threadIdx.x;
    if(INDEX >= LEN) {
        return;
    }

    // Get hash
    uint64_t key = pairs[INDEX].key;
    T value = pairs[INDEX].value;
    uint64_t index = __kh_hash(key) % MAP_SIZE;

    // Linear probing
    for(int i = 0; i < MAP_SIZE; ++i) {
        int cur = (index + i) % MAP_SIZE;
        uint64_t prev = atomicCAS(
            (ULL*)&map[cur].key, 
            (ULL)EMPTY, 
            (ULL)key
        );
        if(prev == key || prev == (ULL)EMPTY) {
            atomicAdd(
                (ULL*)&map[cur].value, 
                (ULL)1
            );
            break;
        }
    }

    // Full/heavy collision
};

#endif