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

#ifndef ANSI_ESC_COMMON
#define ANSI_ESC_COMMON
#define BLACK "\x1B[30m"
#define LBLACK "\x1B[90m"
#define RED "\x1B[31m"
#define LRED "\x1B[91m"
#define GREEN "\x1B[32m"
#define LGREEN "\x1B[92m"
#define YELLOW "\x1B[33m"
#define LYELLOW "\x1B[93m"
#define BLUE "\x1B[34m"
#define LBLUE "\x1B[94m"
#define MAGENTA "\x1B[35m"
#define LMAGENTA "\x1B[95m"
#define CYAN "\x1B[36m"
#define LCYAN "\x1B[96m"
#define WHITE "\x1B[37m"
#define LWHITE "\x1B[97m"
#define RESET "\x1B[0m"
#endif

std::vector<std::string> split_by(const std::string& s, char c);

void ltrim(std::string& s);

void rtrim(std::string& s);

void trim(std::string& s);