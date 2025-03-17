#include "fa_read.h"
#include <cstdint>

class fq_read : public fa_read {
public:
    fq_read::fq_read(std::string id, size_t len, std::string seq_, std::string qual_, std::string meta = "");
    const std::string& get_quality() const;
    uint16_t operator[](size_t index) const;                    // Returns 16-bit bitmask | key: {LOW 0th}[read][qscore]{HIGH 15th}
private:
    std::string qual;
};