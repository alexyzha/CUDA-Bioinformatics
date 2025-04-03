#include "sam_container.h"
#include "util_structs.h"
#include "fx_util.h"

/**
 *  @brief Streams a single sam read to a file.
 *  @param out `ofstream&` out file stream
 *  @param read `sam_read&` read to be streamed
 *  @return `void`
 */
void sam_read_to_file(std::ofstream& out, sam_read& read);

/**
 *  @brief Streams all the contents of a `sam_container` to a file.
 *  @param sam `sam_container&`
 *  @param file_path `string`
 *  @return `void`
 *  @note Format: <all headers> <all reads>
 */
void sam_to_file(sam_container& sam, std::string file_path);

/**
 *  @brief Creates a cigar string using an `alignment`.
 *  @param align `alignment&`
 *  @return `string`
 *  @note Does not compress \"1*\" -> \"*\".
 */
std::string make_cigar(alignment& align);

/**
 *  @brief Maps fastq reads to indices of a reference sequence.
 *  @param ref `string&` reference sequence
 *  @param ref_id `string` reference sequence ID
 *  @param reads `vector<fq_read*>&`
 *  @param k `size_t` kmer length
 *  @return `vector<sam_read*>`
 *  @note If using uint64_t as key, `k` cannot exceed 32.
 */
std::vector<sam_read*> map_reads_to_ref(std::string& ref, std::string ref_id, std::vector<fq_read*>& reads, size_t k);

/**
 *  @brief Given a reference sequence, streams sam reads to a vcf file.
 *  @param file_path `string`
 *  @param ref `string&` reference sequence
 *  @param ref_id `string` reference sequence ID
 *  @param reads `vector<sam_read*>&`
 *  @param CHROMO `int` chromosome
 *  @return `void`
 */
void sam_to_vcf(std::string file_path, std::string& ref, std::string ref_id, std::vector<sam_read*>& reads, int CHROMO);