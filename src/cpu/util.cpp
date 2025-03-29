#include "headers/util.h"

std::vector<std::string> split_by(const std::string& s, char c) {
    std::vector<std::string> substrings = {};
    std::string temp = "";

    // Main processing loop
    for(auto& ch : s) {
        if(ch == c) {
            
            // Only push if there is a substring to push. That is: don't push empty strings
            if(!temp.empty()) {
                substrings.push_back(temp);
            }
            temp.clear();
        } else {
            temp.push_back(ch);
        }
    }

    // Sometimes last substring isn't pushed 
    if(!temp.empty()) {
        substrings.push_back(temp);
    }
    return substrings;
}

void ltrim(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char c) {
        return !std::isspace(c);
    }));
}

void rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char c) {
        return !std::isspace(c);
    }).base(), s.end());
}

void trim(std::string& s) {
    rtrim(s);
    ltrim(s);
}

int sum(const std::string& s, char offset) {
    int ret = 0;
    for(auto& c : s) {
        ret += (c - offset);
    }
    return ret;
}

char base_to_bit(char base) {
    switch(base) {
        case 'A':
            return 0b00;
        case 'C':
            return 0b01;
        case 'G':
            return 0b10;
        case 'T':
            return 0b11;
        default:
            // std::cerr << "INPUT INCLUDES A NON-ACGT BASE" << std::endl;
            return 0x80;
    }
}