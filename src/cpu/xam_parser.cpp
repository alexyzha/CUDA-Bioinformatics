#include "headers/xam_parser.h"

std::vector<sam_read*> read_sam(std::string file_path) {
    // File handling
    std::ifstream file(file_path);
    std::vector<sam_read*> reads = {};
    if(!file.is_open()) {
        throw std::invalid_argument("FILE DOES NOT EXIST OR ISN'T ABLE TO BE OPENED");
    }

    // File handling variables
    std::string cur_line = "";

    // Main processing loop
    while(std::getline(file, cur_line)) {
        trim(cur_line);
        
        // Skip headers and empty lines
        if(cur_line.empty() || cur_line[0] == '@') {
            continue;
        }

        // Sam reads are tab-separated
        std::vector<std::string> read_line = split_by(cur_line, '\t');
        reads.push_back(new sam_read{
            [&](){                    
                // Instantiate with an empty vector if there are no tags in the line
                return read_line.size() >= 12 
                    ? std::vector<std::string>(read_line.begin() + 11, read_line.end())
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
    file.close();
    return reads;
}

std::unordered_map<std::string, std::vector<std::string>> read_sam_headers(std::string file_path) {
    // File handling
    std::ifstream file(file_path);
    std::unordered_map<std::string, std::vector<std::string>> headers;
    if(!file.is_open()) {
        throw std::invalid_argument("FILE DOES NOT EXIST OR ISN'T ABLE TO BE OPENED");
    }

    // File handling variables
    std::string cur_line = "";
    
    // Main processing loop
    while(std::getline(file, cur_line)) {
        trim(cur_line);

        // Skip non header lines. After header section, break
        if(cur_line.empty() || cur_line[0] != '@') {
            break;
        }

        // Header marker is always 3 chars: @**
        std::string marker = cur_line.substr(1, 2);
        
        // Header contents begin at header[3]
        cur_line = cur_line.substr(3);
        trim(cur_line);
        headers[marker].push_back(cur_line);
    }
    file.close();
    return headers;
}