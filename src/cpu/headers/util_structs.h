#include "util.h"

#ifndef C_ALIGNMENT
#define C_ALIGNMENT

struct alignment {
public:
    int score;
    int end_ref;
    int end_read;
    std::string aligned_ref;
    std::string aligned_read;

};

#endif

#ifndef C_UF
#define C_UF

/*
 *  My poor baby...
 *  I'm sorry but the umap was the only way... D,:
 */
class union_find {
public:
    int find(int x);
    void join(int x, int y);
    bool con(int x, int y);

private:
    std::unordered_map<int, int> p;
    std::unordered_map<int, int> h;
    
};

#endif