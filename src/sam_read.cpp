#include "sam_read.h"



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