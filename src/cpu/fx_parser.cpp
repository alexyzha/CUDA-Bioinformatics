#include "headers/fx_parser.h"

std::vector<fa_read*> read_fasta(std::string file_path) {
    std::ifstream file(file_path);
    std::vector<fa_read*> reads = {};
    if(file.is_open()) {
        std::string cur_line = "";
        std::string cur_id = "";
        std::string cur_seq = "";
        std::string cur_meta = "";
        while(std::getline(file, cur_line)) {
            trim(cur_line);
            if(cur_line.empty()) {                                                  // Last line in FASTA could be empty
                break;
            }
            if(cur_line[0] == '>') {                                                // Start of a new read section
                if(!cur_seq.empty()) {
                    reads.push_back(
                        new fa_read(cur_id,
                                    cur_seq.size(),
                                    cur_seq,
                                    cur_meta)
                    );
                }
                cur_seq.clear();                                                
                std::vector<std::string> header = split_by(cur_line, ' ');
                cur_id = header[0].substr(1);                                       // Trim '>' from header begin
                if(header.size() > 2) {
                    cur_meta = header[1];
                }
            } else {
                cur_seq += cur_line;                                                // No need to trim '\n' as it is done before
            }
        }
        if(!cur_seq.empty()) {
            reads.push_back(                                                        // FASTA will never end with a header
                new fa_read(cur_id,
                            cur_seq.size(),
                            cur_seq,
                            cur_meta)
            );
        }
        file.close();
    } else {
        std::cerr << "[" << file_path << "] DOES NOT EXIST OR ISN'T ABLE TO BE OPENED" << std::endl;
    }
    return reads;
}

/*
 *  FASTQ FORMAT:
 *  1. \@ header
 *  2. seq
 *  3. \+ [optional header]
 *  4. qscores
 */
std::vector<fq_read*> read_fastq(std::string file_path) {
    std::ifstream file(file_path);
    std::vector<fq_read*> reads = {};
    int line_num = 0;
    if(file.is_open()) {
        std::string cur_line = "";
        std::string cur_id = "";
        std::string cur_seq = "";
        std::string cur_qual = "";
        std::string cur_meta = "";
        while(std::getline(file, cur_line)) {
            trim(cur_line);
            switch(line_num % 4) {
                case 0: {                                                           // New read block
                    if(!cur_seq.empty()) {
                        reads.push_back(
                            new fq_read(cur_id,
                                        cur_seq.size(),
                                        cur_seq,
                                        cur_qual,
                                        cur_meta)
                        );
                    }
                    cur_seq.clear();
                    if(cur_line.empty()) {                                          // Last line of a FASTQ file could be empty
                        break;
                    }
                    std::vector<std::string> header = split_by(cur_line, ' ');
                    cur_id = header[0].substr(1);
                    if(header.size() > 2) {
                        cur_meta = header[1];
                    }
                    break;
                }
                case 1:                                                             // Sequence
                    cur_seq = cur_line;
                    break;
                case 2:                                                             // Optional header
                    break;
                case 3:                                                             // Quality scores
                    cur_qual = cur_line;
                    break;
                default:
                    throw std::invalid_argument("IMPOSSIBLE MODULO RESULT");
                    break;
            }
            ++line_num;
        }
        if(!cur_seq.empty()) {
            reads.push_back(
                new fq_read(cur_id,
                            cur_seq.size(),
                            cur_seq,
                            cur_qual,
                            cur_meta)
            );
        }
        file.close();
    } else {
        std::cerr << "[" << file_path << "] DOES NOT EXIST OR ISN'T ABLE TO BE OPENED" << std::endl;
    }
    return reads;
}