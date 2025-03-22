#include "sam_container.h"

const std::unordered_map<std::string, std::vector<std::string>>& sam_container::get_headers() const {
    return headers;
}

const std::vector<sam_read*>& sam_container::get_reads() const {
    return reads;
}

void sam_container::set_headers(std::string file_path) {
    headers = read_sam_headers(file_path);
}

void sam_container::set_reads(std::string file_path) {
    reads = read_sam(file_path);
}