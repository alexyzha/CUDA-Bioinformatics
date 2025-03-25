#include "headers/xam_parser.h"

std::vector<sam_read*> read_sam(std::string file_path) {
    std::ifstream file(file_path);
    std::vector<sam_read*> reads = {};
    if(!file.is_open()) {
        throw std::invalid_argument("FILE DOES NOT EXIST OR ISN'T ABLE TO BE OPENED");
    }
    std::string cur_line = "";
    while(std::getline(file, cur_line)) {
        trim(cur_line);
        if(cur_line.empty() || cur_line[0] == '@') {
            continue;
        }
        std::vector<std::string> read_line = split_by(cur_line, '\n');
        reads.push_back(new sam_read{
            [&](){                                                              // Lambda to clearly separate this section
                return read_line.size() >= 12 
                    ? std::vector<std::string>(read_line.begin() + 12, read_line.end())
                    : std::vector<std::string>();
            }(), 
            read_line[0],
            read_line[2],
            read_line[5],
            read_line[6],
            read_line[9],
            read_line[10],
            static_cast<uint16_t>(std::stoi(read_line[1])),
            static_cast<size_t>(std::stoi(read_line[3])),
            static_cast<size_t>(std::stoi(read_line[7])),
            static_cast<size_t>(std::stoi(read_line[8])),
            read_line[4][0]
        });
    }
    return reads;
}

std::unordered_map<std::string, std::vector<std::string>> read_sam_headers(std::string file_path) {
    std::ifstream file(file_path);
    std::unordered_map<std::string, std::vector<std::string>> headers;
    if(!file.is_open()) {
        throw std::invalid_argument("FILE DOES NOT EXIST OR ISN'T ABLE TO BE OPENED");
    }
    std::string cur_line = "";
    while(std::getline(file, cur_line)) {
        trim(cur_line);
        if(cur_line.empty() || cur_line[0] != '@') {
            break;
        }
        std::string marker = cur_line.substr(1, 2);
        cur_line = cur_line.substr(3);
        trim(cur_line);
        headers[marker].push_back(cur_line);
    }
    return headers;
}