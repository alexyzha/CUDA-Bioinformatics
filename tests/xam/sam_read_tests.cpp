#include <gtest/gtest.h>
#include "../../src/sam_read.h"

TEST(SAM_READ, CONSTRUCTOR_NORMAL) {
    sam_read read;
    EXPECT_NO_THROW([&]() {
        read = sam_read({}, "A", "B", "C", "D", "E", "F", 7, 8, 9, 10, 'K');
    }()) << RED << "NORMAL CONSTRUCTOR FAILED" << std::endl;
    EXPECT_TRUE(read.tags.empty()) << RED << "TAGS SOMEHOW NOT EMPTY" << std::endl;
    EXPECT_EQ(read.qname, "A") << RED << "QNAME MISMATCH" << std::endl;
    EXPECT_EQ(read.rname, "B") << RED << "RNAME MISMATCH" << std::endl;
    EXPECT_EQ(read.cigar, "C") << RED << "CIGAR MISMATCH" << std::endl;
    EXPECT_EQ(read.rnext, "D") << RED << "RNEXT MISMATCH" << std::endl;
    EXPECT_EQ(read.seq, "E") << RED << "SEQ MISMATCH" << std::endl;
    EXPECT_EQ(read.qual, "F") << RED << "QUALITY MISMATCH" << std::endl;
    EXPECT_EQ(read.flags, 7) << RED << "FLAG MISMATCH" << std::endl;
    EXPECT_EQ(read.pos, 8) << RED << "POS MISMATCH" << std::endl;
    EXPECT_EQ(read.posnext, 9) << RED << "POSNEXT MISMATCH" << std::endl;
    EXPECT_EQ(read.tlen, 10) << RED << "TLEN MISMATCH" << std::endl;
    EXPECT_EQ(read.mapq, 'K') << RED << "MAPQ MISTMATCH" << std::endl;
}

TEST(SAM_READ, CONSTRUCTOR_AGG_EV) {
    sam_read read;
    EXPECT_NO_THROW([&](){
        read = {{}, "A", "B", "C", "D", "E", "F", 7, 8, 9, 10, 'K'};
    }()) << RED << "AGG CONSTRUCTOR FAILED" << std::endl;
    EXPECT_TRUE(read.tags.empty()) << RED << "TAGS SOMEHOW NOT EMPTY" << std::endl;
    EXPECT_EQ(read.qname, "A") << RED << "QNAME MISMATCH" << std::endl;
    EXPECT_EQ(read.rname, "B") << RED << "RNAME MISMATCH" << std::endl;
    EXPECT_EQ(read.cigar, "C") << RED << "CIGAR MISMATCH" << std::endl;
    EXPECT_EQ(read.rnext, "D") << RED << "RNEXT MISMATCH" << std::endl;
    EXPECT_EQ(read.seq, "E") << RED << "SEQ MISMATCH" << std::endl;
    EXPECT_EQ(read.qual, "F") << RED << "QUALITY MISMATCH" << std::endl;
    EXPECT_EQ(read.flags, 7) << RED << "FLAG MISMATCH" << std::endl;
    EXPECT_EQ(read.pos, 8) << RED << "POS MISMATCH" << std::endl;
    EXPECT_EQ(read.posnext, 9) << RED << "POSNEXT MISMATCH" << std::endl;
    EXPECT_EQ(read.tlen, 10) << RED << "TLEN MISMATCH" << std::endl;
    EXPECT_EQ(read.mapq, 'K') << RED << "MAPQ MISTMATCH" << std::endl;
}

TEST(SAM_READ, CONSTRUCTOR_AGG_FV) {
    sam_read read;
    EXPECT_NO_THROW([&](){
        read = {{"L", "M", "N"}, "A", "B", "C", "D", "E", "F", 7, 8, 9, 10, 'K'};
    }()) << RED << "NESTED AGG CONSTRUCTOR FAILED" << std::endl;
    EXPECT_EQ(read.tags.size(), 3) << RED << "TAGS HAS INCORRECT NUMBER OF ITEMS" << std::endl;
    EXPECT_EQ(read.tags, {"L", "M", "N"}) << RED << "TAGS ITEM MISMATCH" << std::endl;
    EXPECT_EQ(read.qname, "A") << RED << "QNAME MISMATCH" << std::endl;
    EXPECT_EQ(read.rname, "B") << RED << "RNAME MISMATCH" << std::endl;
    EXPECT_EQ(read.cigar, "C") << RED << "CIGAR MISMATCH" << std::endl;
    EXPECT_EQ(read.rnext, "D") << RED << "RNEXT MISMATCH" << std::endl;
    EXPECT_EQ(read.seq, "E") << RED << "SEQ MISMATCH" << std::endl;
    EXPECT_EQ(read.qual, "F") << RED << "QUALITY MISMATCH" << std::endl;
    EXPECT_EQ(read.flags, 7) << RED << "FLAG MISMATCH" << std::endl;
    EXPECT_EQ(read.pos, 8) << RED << "POS MISMATCH" << std::endl;
    EXPECT_EQ(read.posnext, 9) << RED << "POSNEXT MISMATCH" << std::endl;
    EXPECT_EQ(read.tlen, 10) << RED << "TLEN MISMATCH" << std::endl;
    EXPECT_EQ(read.mapq, 'K') << RED << "MAPQ MISTMATCH" << std::endl;
}

TEST(SAM_READ, SB_ACCESS) {
    sam_read read;
    EXPECT_NO_THROW([&](){
        read = {{"L", "M", "N"}, "A", "B", "C", "D", "E", "F", 7, 8, 9, 10, 'K'};
    }()) << RED << "NESTED AGG CONSTRUCTOR FAILED" << std::endl;
    auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = read;
    EXPECT_EQ(tags.size(), 3) << RED << "SB TAGS HAS INCORRECT NUMBER OF ITEMS" << std::endl;
    EXPECT_EQ(tags, {"L", "M", "N"}) << RED << "SB TAGS ITEM MISMATCH" << std::endl;
    EXPECT_EQ(qname, "A") << RED << "SB QNAME MISMATCH" << std::endl;
    EXPECT_EQ(rname, "B") << RED << "SB RNAME MISMATCH" << std::endl;
    EXPECT_EQ(cigar, "C") << RED << "SB CIGAR MISMATCH" << std::endl;
    EXPECT_EQ(rnext, "D") << RED << "SB RNEXT MISMATCH" << std::endl;
    EXPECT_EQ(seq, "E") << RED << "SB SEQ MISMATCH" << std::endl;
    EXPECT_EQ(qual, "F") << RED << "SB QUALITY MISMATCH" << std::endl;
    EXPECT_EQ(flags, 7) << RED << "SB FLAG MISMATCH" << std::endl;
    EXPECT_EQ(pos, 8) << RED << "SB POS MISMATCH" << std::endl;
    EXPECT_EQ(posnext, 9) << RED << "SB POSNEXT MISMATCH" << std::endl;
    EXPECT_EQ(tlen, 10) << RED << "SB TLEN MISMATCH" << std::endl;
    EXPECT_EQ(mapq, 'K') << RED << "SB MAPQ MISTMATCH" << std::endl;
}