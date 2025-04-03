#include "util.h"

#ifndef C_ALIGNMENT
#define C_ALIGNMENT

/**
 *  @brief Holds the results from a Smith-Waterman local alignment.
 *  @note The CUDA version `cu_alignment` holds a `char*` to a CIGAR string instead of both aligned ref and read.
 */
struct alignment {
public:
    /**
     *  @brief `alignment` default constructor.
     *  @note Default constructor sets `int` to 0, `string*` to nullptr.
     */
    alignment();

    /**
     *  @brief `alignment` constructor.
     *  @param s `int` score
     *  @param erf `int` end of reference sequence
     *  @param erd `int` end of read sequence
     *  @param arf `string*` aligned reference sequence
     *  @param ard `string*` aligned read sequence
     */
    alignment(int s, int erf, int erd, std::string* arf, std::string* ard);

public:
    int score;
    int end_ref;
    int end_read;
    std::string* aligned_ref;
    std::string* aligned_read;

};

#endif

#ifndef C_UF
#define C_UF

/**
 *  @brief `union_find` only has a default constructor.
 *  @note This implementation of union find uses unordered maps instead of arrays for parent and height.
 *  @note The CUDA implementation of union find uses different find/join functions and arrays instead of maps.
 */
class union_find {
public:
    /**
     *  @brief Finds the root of `x` while performing path compression.
     *  @param x `int`
     *  @return `int`
     *  @note Different implementation than the CUDA version's find.
     */
    int find(int x);

    /**
     *  @brief Joins `x` and `y`. Performs path compression for both variables.
     *  @param x `int`
     *  @param y `int`
     *  @return `void`
     *  @note Different implementation than the CUDA version's join.
     */
    void join(int x, int y);

    /**
     *  @brief Returns true if `x` and `y` are connected. Performs path compression for both variables.
     *  @param x `int`
     *  @param y `int`
     *  @return `bool`
     */
    bool con(int x, int y);

private:
    std::unordered_map<int, int> p;
    std::unordered_map<int, int> h;
    
};

#endif