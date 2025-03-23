#include "headers/fa_read.h"

fa_read::fa_read(std::string id, size_t len, std::string seq_, std::string meta) : read_id(id), length(len), seq(seq_), metadata(meta) {};

void fa_read::to_file(std::ofstream& out) {
    out << ">" << read_id << " ";
    if(get_metadata() != "") {
        out << metadata << " ";
    }
    out << "length=" << length << "\n";
    out << seq << "\n";
}

const size_t fa_read::size() const {
    return length;
}

const std::string& fa_read::get_id() const {
    return read_id;
}

const std::string& fa_read::get_seq() const {
    return seq;
}

const std::string& fa_read::get_metadata() const {
    return metadata;
}

char fa_read::operator[](size_t index) const {
    if(index >= length) {
        throw std::out_of_range("INDEX OUT OF RANGE IN FA_READ");
    }
    return seq[index];
}