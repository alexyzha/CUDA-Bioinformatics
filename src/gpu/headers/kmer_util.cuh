#include "util.cuh"
#include "util_structs.cuh"

/**
 *  @brief Counts kmers.
 *  @param MAP `kh_pair<uint64_t>*` hashtable, use encoded kmer as key
 *  @param ALL_SEQ `char*` all `LEN` sequences in order, flattened
 *  @param OFFSETS `uint32_t*` array of offsets, `OFFSETS[i]` = end index of seq `i`
 *  @param K `size_t` length of kmer
 *  @param LEN `size_t` number of sequences
 *  @param MAP_LEN `size_t` size of `MAP`
 *  @note `MAP_LEN` is expected to be 4^k * 2.
 */
__global__ void cu_kmer_count(
    kh_pair<uint64_t>* MAP,
    char* ALL_SEQ,
    uint32_t* OFFSETS,
    size_t K,
    size_t LEN,
    size_t MAP_LEN
);

/**
 *  @brief Indexes kmers.
 *  @param MAP `kh_pair<uint32_t[MAP_MAX_INDICES + 1]>*` hashtable, use encoded kmer as key
 *  @param ALL_SEQ `char*` all `LEN` sequences in order, flattened
 *  @param OFFSETS `uint32_t*` array of offsets, `OFFSETS[i]` = end index of seq `i`
 *  @param K `size_t` length of kmer
 *  @param LEN `size_t` number of sequences
 *  @param MAP_LEN `size_t` size of `MAP`
 *  @note `MAP_LEN` is expected to be 4^k * 2.
 *  @note `MAP` memory usage scales very fast with `K`.
 */
__global__ void cu_kmer_index(
    kh_pair<uint32_t[MAP_MAX_INDICES + 1]>* MAP,
    char* ALL_SEQ,
    uint32_t* OFFSETS,
    size_t K,
    size_t LEN,
    size_t MAP_LEN
);

/**
 *  @brief Populates set-length edge list based on which kmers are indexed together.
 *  @param MAP `kh_pair<uint32_t[MAP_MAX_INDICES + 1]>*` hashtable, use encoded kmer as key
 *  @param MAP_LEN `size_t` size of `MAP`
 *  @param EDGE_LIST `uint32_t*` will be populated with pairs of joined edges
 *  @param EDGE_COUNT `uint32_t*` number of edges in `EDGE_LIST`
 *  @param MAX_EDGES `uint32_t` length of `EDGE_LIST`
 *  @note `MAP_LEN` is expected to be 4^k * 2.
 */
__global__ void cu_get_kmer_overlaps(
    kh_pair<uint32_t[MAP_MAX_INDICES + 1]>* MAP,
    size_t MAP_LEN,
    uint32_t* EDGE_LIST,
    uint32_t* EDGE_COUNT,
    uint32_t MAX_EDGES
);

/**
 *  @brief Populates union find based on edge list.
 *  @param UF `cu_union_find*` union find
 *  @param LEN `size_t` length of union find parent/height arrays
 *  @param EDGE_LIST `uint32_t*`
 *  @param EDGE_COUNT `uint32_t` number of pairs in edge list * 2
 *  @note `LEN` is expected to be the total number of reads
 */
__global__ void cu_get_uf(
    cu_union_find* UF,
    size_t LEN,
    uint32_t* EDGE_LIST,
    uint32_t EDGE_COUNT
);

/**
 *  @brief Clusters reads based on union find results.
 *  @param UF `cu_union_find*` contains trees of joined reads
 *  @param MAP `kh_pair<uint32_t[MAP_MAX_INDICES + 1]>*` results
 *  @param LEN `size_t` length of union find parent/height arrays
 *  @param MAP_LEN `size_t` length of `MAP`
 *  @param THRESH `size_t` minimum size of tree to be included
 *  @note `LEN` is expected to be the total number of reads
 *  @note `MAP_LEN` is expected to be `LEN` * 2
 *  @note `MAP` is structured like: <root, array of indexes in subtrees>
 */
__global__ void cu_get_clusters(
    cu_union_find* UF,
    kh_pair<uint32_t[MAP_MAX_INDICES + 1]>* MAP,
    size_t LEN,                         // Size of UF arrays = number of nodes
    size_t MAP_LEN,                     // = 2 * exp number of elements in MAP
    size_t THRESH
);