#include "headers/cu_kmers.h"

__device__ void cu_count_kmers(kh_pair<uint64_t>* MAP, char* ALL_SEQ, size_t* OFFSETS, size_t K, size_t LEN, size_t MAP_LEN) {
    // OOB check for block/thread & K OOB check
    int SEQ_NUM = blockIdx.x * blockDim.x + threadIdx.x;
    if(SEQ_NUM >= LEN || K > K_MAX || K <= 0) {
        return;
    }

    // Rolling hash for O(1) kmer key creation
    uint64_t hash = 0;
    uint64_t mask = (1ULL << (K * 2)) - 1;
    size_t SEQ_BEGIN = OFFSETS[SEQ_NUM - 1] + 1;
    size_t SEQ_SIZE = OFFSETS[SEQ_NUM] + 1 - SEQ_BEGIN;
    if(SEQ_SIZE < K) {
        return;
    }

    // Rolling hash main logic
    for(int i = 0; i < SEQ_SIZE; ++i) {
        hash <<= 2;

        // __base_to_bit returns 0x80 in some cases. For ACGT bases, mask out non-2-bit values
        hash |= (static_cast<char>(0b11) & __base_to_bit(ALL_SEQ[SEQ_BEGIN + i]));

        // From index SEQ_BEGIN to SEQ_BEGIN + K - 1, hash is not yet a full kmer
        if(i >= K - 1) {
            // Get hash
            uint64_t index = kh_hash(hash);

            // Linear probing
            for(int i = 0; i < MAP_LEN; ++i) {
                int cur = (index + i) % MAP_LEN;
                int prev = atomicCAS(
                    (ULL*)&MAP[cur].key, 
                    (ULL)EMPTY, 
                    (ULL)hash
                );
                if(prev == hash || prev == EMPTY) {
                    atomicAdd(
                        (ULL*)&MAP[cur].value, 
                        (ULL)1
                    );
                    break;
                }
            }
        }
    }
}

__device__ void cu_index_kmers(kh_pair<uint32_t[MAP_MAX_INDEXES + 1]>* MAP, char* ALL_SEQ, size_t* OFFSETS, size_t K, size_t LEN, size_t MAP_LEN) {
    // OOB check for block/thread & K OOB check
    int SEQ_NUM = blockIdx.x * blockDim.x + threadIdx.x;
    if(SEQ_NUM >= LEN || K > K_MAX || K <= 0) {
        return;
    }

    // Rolling hash for O(1) kmer key creation
    uint64_t hash = 0;
    uint64_t mask = (1ULL << (K * 2)) - 1;
    size_t SEQ_BEGIN = OFFSETS[SEQ_NUM - 1] + 1;
    size_t SEQ_SIZE = OFFSETS[SEQ_NUM] + 1 - SEQ_BEGIN;
    if(SEQ_SIZE < K) {
        return;
    }

    // Rolling hash main logic
    for(int i = 0; i < SEQ_SIZE; ++i) {
        hash <<= 2;

        // __base_to_bit returns 0x80 in some cases. For ACGT bases, mask out non-2-bit values
        hash |= (static_cast<char>(0b11) & __base_to_bit(ALL_SEQ[SEQ_BEGIN + i]));

        // From index SEQ_BEGIN to SEQ_BEGIN + K - 1, hash is not yet a full kmer
        if(i >= K - 1) {
            // Get hash
            uint64_t index = kh_hash(hash);

            // Linear probing
            for(int i = 0; i < MAP_LEN; ++i) {
                int cur = (index + i) % MAP_LEN;
                uint64_t prev = atomicCAS((ULL*)&MAP[cur].key, (ULL)EMPTY, (ULL)hash);

                // Insert new or found a match in map
                if(prev == EMPTY || prev == hash) {
                    // value[0] = last free index
                    uint32_t* values = &MAP[cur].value[0];
                    uint32_t count = atomicAdd(&values[0], 1);

                    // Prevent overflow/add index to list
                    if(count < MAP_MAX_INDEXES) {
                        values[count + 1] = SEQ_NUM;
                    } else {
                        atomicSub(&values[0], 1);
                    }
                    break;
                }
            }
        }
    }
}