#include "headers/xam_util.h"

void sam_read_to_file(std::ofstream& out, sam_read& read) {
    auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = read;
    out << qname << "\t" 
        << flags << "\t" 
        << rname << "\t" 
        << pos << "\t" 
        << mapq << "\t"
        << cigar << "\t" 
        << rnext << "\t"
        << posnext << "\t"
        << tlen << "\t" 
        << seq << "\t" 
        << qual << "\t";
    for(auto& tag : tags) {
        out << tag << "\t";
    }
    out << "\n";
}

void sam_to_file(sam_container& sam, std::string file_path) {
    std::ofstream file(file_path);
    if(!file.is_open()) {
        throw std::invalid_argument("FILE DOESN'T EXIST OR CANNOT OPEN FILE");
    }
    const auto& headers = sam.get_headers();
    for(auto& [tag, lines] : headers) {
        for(auto& line : lines) {
            file << "@" << tag << "\t" << line << "\n";
        }
    }
    const auto& reads = sam.get_reads();
    for(auto& read : reads) {
        sam_read_to_file(file, *read);
    }
    file.close();
}