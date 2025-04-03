#include "headers/fq_read.h"
#include "headers/util_structs.h"

/**
 *  @brief Describes filter mode for fasta/fastq filtering.
 */
#ifndef FASTX_FILTER_MACROS
#define FASTX_FILTER_MACROS
#define AVERAGE_DISCARD_WHOLE 1
#define SINGLE_DISCARD_WHOLE 2
#define SLIDING_WINDOW 3
#define PROPORTION_DISCARD_WHOLE 4
#endif

/**
 *  @brief Describes alignment scoring rules.
 */
#ifndef FASTX_ALIGNMENT_MACROS
#define FASTX_ALIGNMENT_MACROS
#define ALN_MATCH 2
#define ALN_MISMATCH -1
#define ALN_GAP -2
#endif

/**
 *  @brief Describes maximum kmer length when using uint64_t as a key.
 */
#ifndef K_MAX
#define K_MAX 32
#endif

/**
 *  @brief Convenience function for streaming multiple fq_reads to file.
 *  @param reads `const vector<fq_read*>&`
 *  @param file_path `string` path to file to stream reads to
 *  @return `void`
 */
void fq_to_file(const std::vector<fq_read*>& reads, std::string file_path);

/**
 *  @brief Filters vector of fq_read.
 *  @param reads `const vector<fq_read*>&` vector with reads to filter
 *  @param FILTER_BY `char` expects value from `FASTX_FILTER_MACROS`
 *  @param THRESH `char` threshold used to filter reads by quality
 *  @param PERC `double` cutoff proportion when < 1.0
 *  @param PERC `static_cast<int>` sliding window size when >= 1.0
 *  @return `vector<fq_read*>` vector with filtered fq_reads
 *  @note `PERC` mode is decided by `FILTER_BY` value.
 *  @note Sliding window mode: `PERC` is converted to `int` and used as sliding window size (expects >= 1.0, whole number).
 *  @note Proportion discard whole mode: `PERC` is used as cutoff proportion (expects < 1.0).
 */
std::vector<fq_read*> filter_fq(const std::vector<fq_read*>& reads, char FILTER_BY, char THRESH, double PERC = 0.0);

/**
 *  @brief Calculates GC% for all fq_reads in a list.
 *  @param reads `const vector<fq_read*>&`
 *  @return `vector<double>` list of GC%
 *  @note return[i] = GC% of read `i`, or reads[i].
 */
std::vector<double> gc_per_read(const std::vector<fq_read*>& reads);

/**
 *  @brief Calculates total GC% for all nucleotides across all reads.
 *  @param reads `const vector<fq_read*>&`
 *  @return `double`
 */
double gc_global(const std::vector<fq_read*>& reads);

/**
 *  @brief Counts kmers in a list of reads. "Hashed" using one-hot encoding.
 *  @param reads `const vector<fq_reads*>&`
 *  @param k `size_t` kmer size
 *  @return `unordered_map<uint64_t, uint64_t>` = [kmer hash, occurrences of that hash]
 *  @note If using uint64_t as key, `k` cannot exceed 32.
 *  @note For hashing: any non-ACGT base returns (char)0x80 but is masked with (char)0b11 to produce (char)0b00, the bitwise representation of base `A`
 */
std::unordered_map<uint64_t, uint64_t> count_kmer(const std::vector<fq_read*>& reads, size_t k);

/**
 *  @brief Indexes reads based on kmer occurrence.
 *  @param reads `const vector<fq_read*>&` 
 *  @param k `size_t` kmer size
 *  @return `unordered_map<uint64_t, unordered_set<int>>` = [kmer hash, set<indexes of all reads containing kmer>]
 *  @note If using uint64_t as key, `k` cannot exceed 32.
 *  @note For hashing: any non-ACGT base returns (char)0x80 but is masked with (char)0b11 to produce (char)0b00, the bitwise representation of base A.
 */
std::unordered_map<uint64_t, std::unordered_set<int>> index_kmer(const std::vector<fq_read*>& reads, size_t k);

/**
 *  @brief Smith-Waterman local alignment.
 *  @param ref `string&` usually fasta read 
 *  @param read `string&` usually fastq read
 *  @return `alignment` definition in `/src/cpu/headers/util_structs.h`
 *  @note Scores based on `FASTX_ALIGNMENT_MACROS`.
 */
alignment local_align(const std::string& ref, const std::string& read);

/**
 *  @brief Groups fastq reads by kmer similarity.
 *  @param kmer_map `unordered_map<uint64_t, unordered_set<int>>&` output of `index_kmer(...)`
 *  @param READS `int` number of reads
 *  @param THRESH `int` number of similar kmers to connect 2 reads, default value = 1
 *  @return `vector<unordered_set<int>*>` index `i` holds ptr to set of all reads clustered with read `i`
 *  @note Internal implementation uses union find.
 */
std::vector<std::unordered_set<int>*> cluster_by_kmer(std::unordered_map<uint64_t, std::unordered_set<int>>& kmer_map, int READS, int THRESH = 1);

