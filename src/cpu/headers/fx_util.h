#include "headers/fq_read.h"
#include "headers/util_structs.h"

#ifndef FASTX_FILTER_MACROS
#define FASTX_FILTER_MACROS
#define AVERAGE_DISCARD_WHOLE 1
#define SINGLE_DISCARD_WHOLE 2
#define SLIDING_WINDOW 3
#define PROPORTION_DISCARD_WHOLE 4
#endif

#ifndef FASTX_ALIGNMENT_MACROS
#define FASTX_ALIGNMENT_MACROS
#define ALN_MATCH 2
#define ALN_MISMATCH -1
#define ALN_GAP -2
#endif

/*
 *  - Convenience function for vector of fa/fq reads
 */
void fa_to_file(const std::vector<fq_read*>& reads, std::string file_path);

/*
 *  - PERC = cutoff proportion when < 1.0
 *  - PERC = sliding window window size when >= 1.0
 */
std::vector<fq_read*> filter_fq(const std::vector<fq_read*>& reads, char FILTER_BY, char THRESH, double PERC = 0.0);

/*
 *  - Input can be either vectors of fx_read
 */
std::vector<double> gc_per_read(const std::vector<fq_read*>& reads);

/*
 *  - Input can be either vectors of fx_read
 */
double gc_global(const std::vector<fq_read*>& reads);

/*
 *  - K_MAX = 32
 *  - Returns map with bitmasked/"hashed" keys
 *  - May consider using pair<u64, u64> for higher K_MAX later
 *  - A = 00, C = 01, G = 10, T = 11
 *  - B/c bit inversion swaps to complementary
 */
std::unordered_map<uint64_t, uint64_t> count_kmer(const std::vector<fq_read*>& reads, size_t k);

/*
 *  - Returns map<hash, set<indexes of all reads containing hash>>
 */
std::unordered_map<uint64_t, std::unordered_set<int>> index_kmer(const std::vector<fq_read*>& reads, size_t k);

/*
 *  - Ref = fasta, read = fastq
 */
alignment local_align(const std::string& ref, const std::string& read);

/*
 *  - THRESH = number of overlaps needed to be placed in an island. DEFAULT = 1
 *  - READS = number of unique reads in kmer_map 
 *  - Return format: vector of size READS. index i holds ptr to set of all reads rooted at i 
 */
std::vector<std::unordered_set<int>*> cluster_by_kmer(std::unordered_map<uint64_t, std::unordered_set<int>>& kmer_map, int READS, int THRESH = 1);

