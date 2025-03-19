#pragma once
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <locale>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

std::vector<std::string> split_by(const std::string& s, char c);

inline void ltrim(std::string &s);

inline void rtrim(std::string &s);

inline void trim(std::string &s);