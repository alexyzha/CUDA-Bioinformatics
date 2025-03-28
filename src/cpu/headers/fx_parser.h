#include "fa_read.h"
#include "fq_read.h"
#include "util.h"

/*
 *  Reads all fasta reads from valid .fasta file.
 *  @param file_path `string` path to file to read fasta reads from
 *  @return `vector<fa_read*>` vector of pointers to all fasta reads from file
 */
std::vector<fa_read*> read_fasta(std::string file_path);

/*
 *  Reads all fastq reads from valid .fastq file.
 *  @param file_path `string` path to file to read fastq reads from
 *  @return `vector<fq_read*>` vector of pointers to all fastq reads from file
 */
std::vector<fq_read*> read_fastq(std::string file_path);