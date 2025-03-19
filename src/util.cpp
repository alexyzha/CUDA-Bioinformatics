#include "util.h"

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

inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char c) {
        return !std::isspace(c);
    }));
}

inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char c) {
        return !std::isspace(c);
    }).base(), s.end());
}

inline void trim(std::string &s) {
    rtrim(s);
    ltrim(s);
}