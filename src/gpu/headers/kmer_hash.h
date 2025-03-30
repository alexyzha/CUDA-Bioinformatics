#include "util.h"

#ifndef C_KMER_HASH_TABLE
#define C_KMER_HASH_TABLE

template<typename T>
struct khash_pair {
    uint64_t key;
    T value;
};

#endif