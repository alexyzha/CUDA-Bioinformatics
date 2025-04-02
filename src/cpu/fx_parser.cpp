#include "headers/fx_parser.h"

std::vector<fa_read*> read_fasta(std::string file_path) {
    // File stream handling + create output vector
    std::ifstream file(file_path);
    std::vector<fa_read*> reads = {};
    if(!file.is_open()) {
        throw std::invalid_argument("FILE DOES NOT EXIST OR ISN'T ABLE TO BE OPENED");
    }

    // File stream/temp variables
    std::string cur_line = "";
    std::string cur_id = "";
    std::string cur_seq = "";
    std::string cur_meta = "";

    // Main processing loop
    while(std::getline(file, cur_line)) {
        // Get rid of any trailing/leading whitespace in line
        trim(cur_line);

        // Last line in FASTA could be empty
        if(cur_line.empty()) {
            break;
        }

        // Check if current line is the start of a new read section
        if(cur_line[0] == '>') {
            if(!cur_seq.empty()) {
                reads.push_back(new fa_read(
                    cur_id,
                    cur_seq.size(),
                    cur_seq,
                    cur_meta
                ));
            }
            cur_seq.clear();                                                
            std::vector<std::string> header = split_by(cur_line, ' ');
            cur_id = header[0].substr(1);                                       // Trim '>' from header begin
            if(header.size() > 2) {
                cur_meta = header[1];
            }
        } else {
            cur_seq += cur_line;
        }
    }

    // FASTA will never end with a header section, therefore we need push the last read
    if(!cur_seq.empty()) {
        reads.push_back(new fa_read(
            cur_id,
            cur_seq.size(),
            cur_seq,
            cur_meta
        ));
    }
    file.close();
    return reads;
}

std::vector<fq_read*> read_fastq(std::string file_path) {
    
    /*
     *  REMINDER [FASTQ FORMAT]
     *  1. \@ header
     *  2. seq
     *  3. \+ [optional header]
     *  4. qscores
     */

    // File stream handling + create output 
    std::ifstream file(file_path);
    std::vector<fq_read*> reads = {};
    if(!file.is_open()) {
        throw std::invalid_argument("FILE DOES NOT EXIST OR ISN'T ABLE TO BE OPENED");
    }
    
    // File stream/temp variables
    int line_num = 0;
    std::string cur_line = "";
    std::string cur_id = "";
    std::string cur_seq = "";
    std::string cur_qual = "";
    std::string cur_meta = "";

    // Main processing loop
    while(std::getline(file, cur_line)) {
        trim(cur_line);
        switch(line_num % 4) {
            // New read block case/header case
            case 0: {
                // First block being read case
                if(!cur_seq.empty()) {
                    reads.push_back(new fq_read(
                        cur_id,
                        cur_seq.size(),
                        cur_seq,
                        cur_qual,
                        cur_meta
                    ));
                }
                cur_seq.clear();
                if(cur_line.empty()) {                                          // Last line of a FASTQ file could be empty
                    break;
                }

                // Parse header, get rip of header marker
                std::vector<std::string> header = split_by(cur_line, ' ');
                cur_id = header[0].substr(1);
                if(header.size() > 2) {
                    cur_meta = header[1];
                }
                break;
            }

            // Sequence line
            case 1:
                cur_seq = cur_line;
                break;
            
            // Optional additional header line
            case 2:
                break;
            
            // Quality line
            case 3:
                cur_qual = cur_line;
                break;

            // (line_num % 4) should not return anything other than 0, 1, 2, or 3
            default:
                throw std::invalid_argument("IMPOSSIBLE MODULO RESULT");
                break;
        }
        ++line_num;
    }

    // Since fastq files don't end with a header, need to push last read manually
    if(!cur_seq.empty()) {
        reads.push_back(new fq_read(
            cur_id,
            cur_seq.size(),
            cur_seq,
            cur_qual,
            cur_meta
        ));
    }
    file.close();
    return reads;
}