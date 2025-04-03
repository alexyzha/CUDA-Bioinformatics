#pragma once
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <locale>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

/**
 *  @brief Convenience macros for printing GTEST error statements.
 */
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

/**
 *  @brief Convenience macros for PHRED q-score related tasks.
 */
#ifndef PHRED_COMMON
#define PHRED_COMMON
#define PHRED_BEGIN '!'
#define PHRED_BEGIN_INT 33
#endif

/**
 *  @brief Splits a string into substrings by `c`.
 *  @param s `const string&` string to be split
 *  @param c `char` char to split `s` by
 *  @return `vector<string>` list of substrings after splitting `s` by `c`
 */
std::vector<std::string> split_by(const std::string& s, char c);

/**
 *  @brief Trims a string of whitespace from the left.
 *  @param s `string&`
 *  @return `void`
 */
void ltrim(std::string& s);

/**
 *  @brief Trims a string of whitespace from the right.
 *  @param s `string&`
 *  @return `void`
 */
void rtrim(std::string& s);

/**
 *  @brief Trims a string of whitespace from both sides.
 *  @param s `string&`
 *  @return `void`
 */
void trim(std::string& s);

/**
 *  @brief Treats a string as a vector of numbers and returns its sum.
 *  @param s `string&` string to be summed
 *  @param offset `char` a set value to subtract from every instance of char addition, default = 33 = `PHRED_BEGIN`
 *  @return `int`
 */
int sum(const std::string& s, char offset = PHRED_BEGIN);

/**
 *  @brief Encodes nucleotides as bits.
 *  @param base `char` ACGT... etc.
 *  @return `char`
 *  @note A -> 0b00, C -> 0b01, G -> 0b10, T -> 0b11, all others -> 0x80.
 */
char base_to_bit(char base);

/**
 *  @brief Debugging helper function. Prints the bitwise representation of ints, chars, etc.
 *  @param num `*`
 *  @return `string` chainable in a GTEST error output
 */
template<typename C>
std::string meow(C num) {
    std::string ret = "";
    for(int bits = sizeof(num) * 8; bits--; num>>=1) {
        ret.push_back('0' + !!(num & 1));
    }
    std::reverse(ret.begin(), ret.end());
    return ret;
}