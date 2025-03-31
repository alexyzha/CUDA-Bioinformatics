#include "headers/kmers.cuh"

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
            uint64_t index = kh_hash(hash) % MAP_LEN;

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

__device__ void cu_index_kmers(kh_pair<uint32_t[MAP_MAX_INDICES + 1]>* MAP, char* ALL_SEQ, size_t* OFFSETS, size_t K, size_t LEN, size_t MAP_LEN) {
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
            uint64_t index = kh_hash(hash) % MAP_LEN;

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
                    if(count < MAP_MAX_INDICES) {
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

__device__ void cu_get_kmer_overlaps(kh_pair<uint32_t[MAP_MAX_INDICES + 1]>* MAP, size_t MAP_LEN, uint32_t* EDGE_LIST, uint32_t* EDGE_COUNT, uint32_t MAX_EDGES) {
    // OOB check for block/thread
    int SEQ_NUM = blockIdx.x * blockDim.x + threadIdx.x;
    if(SEQ_NUM >= MAP_LEN) {
        return;
    }

    // Get array of matching indices
    uint32_t* matches = MAP[SEQ_NUM].value;
    uint32_t count = matches[0];

    // Iterate through all pairs
    for(int i = 1; i <= count; ++i) {
        for(int j = i + 1; j <= count; ++j) {
            uint32_t src = matches[i];
            uint32_t dest = matches[j];

            // Add edges to list, watch for end of array
            int index = atomicAdd(EDGE_COUNT, 2);
            if(index + 1 >= MAX_EDGES) {
                return;
            }

            // Append pair to end
            EDGE_LIST[index] = src;
            EDGE_LIST[index + 1] = dest;
        }
    }
}

__device__ void cu_get_uf(cu_union_find* UF, size_t LEN, size_t NODES, uint32_t* EDGE_LIST, uint32_t EDGE_COUNT) {
    // Block/thread OOB checks
    int INDEX = (blockIdx.x * blockDim.x + threadIdx.x) * 2;
    if(INDEX + 1 >= EDGE_COUNT) {
        return;
    }

    // UF join
    uint32_t src = EDGE_LIST[INDEX];
    uint32_t dest = EDGE_LIST[INDEX + 1];
    cu_uf_join(UF, src, dest);
}

__device__ void cur_get_clusters{cu_union_find* UF, kh_pair<uint32_t[MAX_CLUSTER_SIZE + 1]> * MAP, size_t LEN, size_t MAP_LEN, size_t K} {
    // Block/thread OOB checks
    int INDEX = blockIdx.x * blockDim.x + threadIdx.x;
    if(INDEX >= LEN) {
        return;
    }

    // Check root size
    int px = cu_uf_find(UF, INDEX);
    if(UF->h[px] < K) {
        return;
    }

    // Add to cluster
    uint64_t index = kh_hash((ULL)INDEX) % MAP_LEN;
    for(int i = 0; i < MAP_LEN; ++i) {
        int cur = (index + i) % MAP_LEN;
        uint64_t prev = atomicCAS((ULL*)&MAP[cur].key, (ULL)EMPTY, (ULL)INDEX);

        // Insert new or found a match in map
        if(prev == EMPTY || prev == INDEX) {
            // cluster[0] = last free index
            uint32_t* cluster = &MAP[cur].value[0];
            uint32_t count = atomicAdd(&cluster[0], 1);

            // Prevent overflow/add index to list
            if(count < MAX_CLUSTER_SIZE) {
                cluster[count + 1] = INDEX;
            } else {
                atomicSub(&cluster[0], 1);
            }
            break;
        }
    }
}