#include "util.h"

#ifndef C_SAM_READ
#define C_SAM_READ

struct sam_read {
public:
    std::vector<std::string> tags;
    std::string qname;
    std::string rname;
    std::string cigar;
    std::string rnext;
    std::string seq;
    std::string qual;
    uint16_t flags;
    size_t pos;
    size_t posnext;
    size_t tlen;
    char mapq;

};

#endif

/*
 *  1. Query name (QNAME): read id
 *  2. FLAG: bitmask for flags, req uint16_t
 *  3. Reference name (RNAME): which contig in the reference genome the read is aligned to (should be found in headers in one of the @SQ lines)
 *  4. Position: leftmost mapping pos of the first matching base within the read
 *  5. Mapping quality (MAPQ): phred qscore (255 = placeholder)
 *  6. CIGAR str
 *  7. Reference name for mate (RNEXT): "=" if identical to RNAME, describes paired-end mate of the read
 *  8. Position of mate (PNEXT): = 4 (pos) or follows same rules
 *  9. Template length (TLEN): length of template seq to which the read maps
 *  10. Sequence (SEQ): actual read sequence
 *  11. Quality string (QUAL): quality of read
 *  12. Predefined tags: additional information for alignment or read
 *  
 *  Tab separated
 *  Headers begin with @**
 *  Simple aggregate used, structured binding syntax:
 *  - auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = read;
 *  Decl/construct:
 *  - sam_read read = {vec<str>, str*6, uin16, sizet*3, char}
 */