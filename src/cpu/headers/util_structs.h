#include "util.h"

struct alignment {
public:
    int score;
    int end_ref;
    int end_read;
    std::string aligned_ref;
    std::string aligned_read;

};