#include "util.h"

#ifndef C_SAM_READ
#define C_SAM_READ

/**
 *  @brief Contains all data from a singluar SAM read.
 *  @note Built specifically with structured bindings in mind: no functions.
 *  @note See declaration for more detailed description of data layout.
 */
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
 *  DESCRIPTION OF SAM READ DATA LAYOUT IN FILE (1-12)
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
 *  Simple aggregate used, structured binding syntax (13, 14):
 *  13. auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = read;
 *  Decl/construct:
 *  14. sam_read read = {vec<str>, str, str, str, str, str, str, uint16_t, sizet, size_t, size_t, char};
 */