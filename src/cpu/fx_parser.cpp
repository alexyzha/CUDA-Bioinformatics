#include "fx_parser.h"

std::vector<std::string> split_by(const std::string& s, char c) {
    std::vector<std::string> substrings = {};
    std::string temp = "";
    for(auto& ch : s) {
        if(ch == c) {
            if(!temp.empty()) {
                substrings.push_back(temp);
            }
            temp.clear();
        } else {
            temp.push_back(ch);
        }
    }
    if(!temp.empty()) {
        substrings.push_back(temp);
    }
    return substrings;
}

std::vector<fa_read*> read_fasta(std::string file_path) {
    std::ifstream file(file_path);
    std::vector<fa_read*> reads = {};
    if(file.is_open()) {
        std::string cur_line = "";
        std::string cur_id = "";
        std::string cur_seq = "";
        std::string cur_meta = "";
        int cur_len = 0;
        while(std::getline(file, cur_line)) {
            if(cur_line[0] == '>') {
                if(cur_len) {
                    reads.push_back(new fa_read(cur_id, cur_len, cur_seq, cur_meta));
                }
                cur_seq.clear();
                std::vector<std::string> header = split_by(cur_line, ' ');
                cur_id = header[0].substr(1);                                   // Trim '>' from header begin
                if(header.size() > 2) {
                    cur_meta = header[1];
                }
                int temp = 0;                                                   // Extract number from last part of header
                for(auto& c : header.back()) {                                  // format: 'length=XXX'
                    if(c >= '0' && c <= '9') {
                        temp *= 10;
                        temp += (c - '0');
                    }
                }
                cur_len = temp;
            } else {
                cur_seq += cur_line;
                cur_seq.pop_back();                                             // Trim '\n'
            }
        }
        reads.push_back(new fa_read(cur_id, cur_len, cur_seq, cur_meta));       // FASTA will never end with a header
        file.close();
    } else {
        std::cerr << "[" << file_path << "] DOES NOT EXIST OR ISN'T ABLE TO BE OPENED" << std::endl;
    }
    return reads;
}

/*
 *  1. header
 *  2. seq
 *  3. + [optional]
 *  4. qscores
 */
std::vector<fq_read*> read_fastq(std::string file_path) {
    std::ifstream file(file_path);
    std::vector<fq_read*> reads = {};
    if(file.is_open()) {
        /*

            erm
            thinkin about exact implementation rn

        */
    } else {
        std::cerr << "[" << file_path << "] DOES NOT EXIST OR ISN'T ABLE TO BE OPENED" << std::endl;
    }
    return reads;
}