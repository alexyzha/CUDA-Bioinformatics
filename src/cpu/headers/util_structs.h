#include "util.h"

#ifndef C_ALIGNMENT
#define C_ALIGNMENT

struct alignment {
public:
    /*
     *  `alignment` default constructor.
     */
    alignment();

    /*
     *  `alignment` constructor.
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

/*
 *  `union_find` only has a default constructor.
 *  @param ._. My poor baby...
 *  @return I'm sorry but the umap was the only way... D,:
 *  
 */
class union_find {
public:
    /*
     *  Finds the root of `x` while performing path compression.
     *  @param x `int`
     *  @return `int`
     */
    int find(int x);

    /*
     *  Joins `x` and `y`. Performs path compression for both variables.
     *  @param x `int`
     *  @param y `int`
     *  @return `void`
     */
    void join(int x, int y);

    /*
     *  Returns true if `x` and `y` are connected. Performs path compression for both variables.
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