#include "headers/fq_read.h"

#ifndef FASTX_FILTER_MACROS
#define FASTX_FILTER_MACROS
#define AVERAGE_DISCARD_WHOLE 1
#define SINGLE_DISCARD_WHOLE 2
#define SLIDING_WINDOW 3
#define PROPORTION_DISCARD_WHOLE 4
#endif

/*
 *  - PERC = cutoff proportion when < 1.0
 *  - PERC = sliding window window size when >= 1.0
 */
std::vector<fa_read*> filter_fa(std::vector<fa_read*> reads, char FILTER_BY, char THRESH, double PERC = 0.0);

/*
 *  - Input can be either vectors of fx_read
 */
std::vector<double> gc_per_read(std::vector<fa_read*> reads);

/*
 *  - Input can be either vectors of fx_read
 */
double gc_global(std::vector<fa_read*> reads);

/*
 *  - K_MAX = 32
 *  - Returns map with bitmasked/"hashed" keys
 *  - May consider using pair<u64, u64> for higher K_MAX later
 *  - A = 00, C = 01, G = 10, T = 11
 *  - B/c bit inversion swaps to complementary
 */
std::unordered_map<uint64_t, uint64_t> count_kmer(std::vector<fa_read*> reads, size_t k);