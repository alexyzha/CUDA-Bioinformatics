#include "util.h"

#ifndef C_KMER_HASH_TABLE
#define C_KMER_HASH_TABLE
#define EMPTY 0

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