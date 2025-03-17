#include "fq_read.h"

fq_read::fq_read(std::string id, size_t len, std::string seq_, std::string qual_, std::string meta) : fa_read(id, len, seq_, meta), qual(qual_) {}

const std::string& fq_read::get_quality() const {
    return qual;
}

uint16_t fq_read::operator[](size_t index) const {
    if(index < 0 || index >= length) {
        throw std::out_of_range("INDEX OUT OF RANGE IN FQ_READ");
    }
    return (static_cast<uint16_t>(qual[index]) << 8) | static_cast<uint16_t>(seq[index]);
}