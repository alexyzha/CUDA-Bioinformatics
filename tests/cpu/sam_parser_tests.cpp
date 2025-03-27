#include <gtest/gtest.h>
#include "../../src/cpu/headers/xam_parser.h"

TEST(SAM_PARSER, EMPTY) {
    std::string file = "infiles/empty.txt";
    std::vector<sam_read*> reads = read_sam(file);
    std::unordered_map<std::string, std::vector<std::string>> headers = read_sam_headers(file);
    EXPECT_TRUE(reads.empty()) << RED << "READ FROM EMPTY FILE" << RESET << std::endl;
    EXPECT_TRUE(headers.empty()) << RED << "HEADER FROM EPMTY FILE" << RESET << std::endl;
}

TEST(SAM_PARSER, READS_ONLY_NO_TAGS) {
    std::string file = "infiles/sam_parser_read_no_tags_in.txt";
    std::vector<sam_read*> reads = read_sam(file);
    std::unordered_map<std::string, std::vector<std::string>> headers = read_sam_headers(file);
    EXPECT_TRUE(headers.empty()) << RED << "HEADER FROM FILE WITH NO HEADER" << std::endl;
    EXPECT_EQ(reads.size(), 1) << RED << "GOT != 1 READ FROM FILE WITH SINGULAR READ" << std::endl;
    auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = (*reads[0]);
    EXPECT_EQ(tags, std::vector<std::string>{}) << RED << "TAGS CAPTURED WHEN EXP NO TAGS" << RESET << std::endl;
    EXPECT_EQ(qname, "ID1") << RED << "WRONG QNAME" << RESET << std::endl;
    EXPECT_EQ(flags, 0) << RED << "WRONG FLAGS" << RESET << std::endl;
    EXPECT_EQ(rname, "REF1") << RED << "WRONG REFNAME" << RESET << std::endl;
    EXPECT_EQ(pos, 1) << RED << "WRONG POS" << RESET << std::endl;
    EXPECT_EQ(mapq, 'Q') << RED << "WRONG MAPQ" << RESET << std::endl;
    EXPECT_EQ(cigar, "CIGAR") << RED << "WRONG CIGAR STRING" << RESET << std::endl;
    EXPECT_EQ(rnext, "=") << RED << "WRONG RNEXT" << RESET << std::endl;
    EXPECT_EQ(posnext, 2) << RED << "WRONG POSNEXT" << RESET << std::endl;
    EXPECT_EQ(tlen, 4) << RED << "WRONG TLEN" << RESET << std::endl;
    EXPECT_EQ(seq, "ACTG") << RED << "WRONG TLEN" << RESET << std::endl;
    EXPECT_EQ(qual, "????") << RED << "WRONG QUALITY STRING" << RESET << std::endl;
}

TEST(SAM_PARSER, READS_ONLY_YES_TAGS) {
    std::string file = "infiles/sam_parser_read_yes_tags_in.txt";
    std::vector<sam_read*> reads = read_sam(file);
    std::unordered_map<std::string, std::vector<std::string>> headers = read_sam_headers(file);
    EXPECT_TRUE(headers.empty()) << RED << "HEADER FROM FILE WITH NO HEADER" << std::endl;
    EXPECT_EQ(reads.size(), 1) << RED << "GOT != 1 READ FROM FILE WITH SINGULAR READ" << std::endl;
    auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = (*reads[0]);
    EXPECT_EQ(tags, std::vector<std::string>({"TAG1", "TAG2", "TAG3"})) << RED << "WRONG/NO TAGS CAPTURED WHEN EXP \'TAG[1:3]\'" << RESET << std::endl;
    EXPECT_EQ(qname, "ID1") << RED << "WRONG QNAME" << RESET << std::endl;
    EXPECT_EQ(flags, 0) << RED << "WRONG FLAGS" << RESET << std::endl;
    EXPECT_EQ(rname, "REF1") << RED << "WRONG REFNAME" << RESET << std::endl;
    EXPECT_EQ(pos, 1) << RED << "WRONG POS" << RESET << std::endl;
    EXPECT_EQ(mapq, 'Q') << RED << "WRONG MAPQ" << RESET << std::endl;
    EXPECT_EQ(cigar, "CIGAR") << RED << "WRONG CIGAR STRING" << RESET << std::endl;
    EXPECT_EQ(rnext, "=") << RED << "WRONG RNEXT" << RESET << std::endl;
    EXPECT_EQ(posnext, 2) << RED << "WRONG POSNEXT" << RESET << std::endl;
    EXPECT_EQ(tlen, 4) << RED << "WRONG TLEN" << RESET << std::endl;
    EXPECT_EQ(seq, "ACTG") << RED << "WRONG TLEN" << RESET << std::endl;
    EXPECT_EQ(qual, "????") << RED << "WRONG QUALITY STRING" << RESET << std::endl;
}

TEST(SAM_PARSER, MULTI_READS) {
    // EXP
    std::vector<std::vector<std::string>> exp_tags = {
        {"OOPY"}, {"TAG1", "TAG2", "TAG3"}, {"BRUH", "BRUH_AGAIN"}
    };
    std::vector<std::string> exp_qname = {"ID1", "ID2", "ID3"};
    std::vector<uint16_t> exp_flags = {0, 6, 19};
    std::vector<std::string> exp_rname = {"REF1", "R2", "D2"};
    std::vector<size_t> exp_pos = {1, 109, 4200};
    std::vector<char> exp_mapq = {'Q', 'P', 'R'};
    std::vector<std::string> exp_cigar = {"CIGAR", "STRINGY", "YIPPEE"};
    std::vector<std::string> exp_rnext = {"=", "*", "BRUH"};
    std::vector<size_t> exp_posnext = {2400, 6900, 420};
    std::vector<size_t> exp_tlen = {12034, 102398, 120};
    std::vector<std::string> exp_seq = {"LKJAHG", "UPASDGBH", "GACTCGA"};
    std::vector<std::string> exp_qual = {"????&1234", "!09123", "&@*1010"};
    // Test
    std::string file = "infiles/sam_parser_multi_read_in.txt";
    std::vector<sam_read*> reads = read_sam(file);
    std::unordered_map<std::string, std::vector<std::string>> headers = read_sam_headers(file);
    EXPECT_TRUE(headers.empty()) << RED << "HEADER FROM FILE WITH NO HEADER" << std::endl;
    EXPECT_EQ(reads.size(), 3) << RED << "GOT != 3 READ FROM FILE WITH 3 READS" << std::endl;
    for(int i = 0; i < 3; ++i) {
        auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = (*reads[i]);
        EXPECT_EQ(tags, exp_tags[i]) << RED << "WRONG/NO TAGS CAPTURED FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(qname, exp_qname[i]) << RED << "WRONG QNAME FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(flags, exp_flags[i]) << RED << "WRONG FLAGS FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(rname, exp_rname[i]) << RED << "WRONG REFNAME FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(pos, exp_pos[i]) << RED << "WRONG POS FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(mapq, exp_mapq[i]) << RED << "WRONG MAPQ FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(cigar, exp_cigar[i]) << RED << "WRONG CIGAR STRING FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(rnext, exp_rnext[i]) << RED << "WRONG RNEXT FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(posnext, exp_posnext[i]) << RED << "WRONG POSNEXT FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(tlen, exp_tlen[i]) << RED << "WRONG TLEN FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(seq, exp_seq[i]) << RED << "WRONG TLEN FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(qual, exp_qual[i]) << RED << "WRONG QUALITY STRING FOR READ [" << i << "]" << RESET << std::endl;
    }
}

TEST(SAM_PARSER, HEADER_ONLY_NO_REPEAT) {
    // EXP
    std::unordered_map<std::string, std::vector<std::string>> exp_headers {
        {"AA", {"SOMETHING_IN_HEADER"}},
        {"AB", {"SOMETHING ELSE IN HEADER"}},
        {"AC", {"ANOTHER SOMETHING IN HEADER"}}
    };
    // Test
    std::string file = "infiles/sam_parser_only_header_in.txt";
    std::vector<sam_read*> reads = read_sam(file);
    std::unordered_map<std::string, std::vector<std::string>> headers = read_sam_headers(file);
    EXPECT_TRUE(reads.empty()) << RED << "READ FROM FILE WITH ONLY HEADERS/NO READS" << RESET << std::endl;
    EXPECT_EQ(headers.size(), 3) << RED << "GOT != 3 HEADER CATEGORIES FROM FILE WITH 3 CATS" << RESET << std::endl;
    for(auto& [key, vec] : headers) {
        EXPECT_NE(exp_headers.find(key), exp_headers.end()) << RED << "UNEXPECTED HEADER DELIM PARSED" << RESET << std::endl;
        EXPECT_EQ(vec.size(), exp_headers[key].size()) << RED << "KEY [" << key << "] HAS != 1 HEADER ASSOC, EXP 1" << RESET << std::endl; 
        EXPECT_EQ(vec[0], exp_headers[key][0]) << RED << "HEADER CONTENT MISMATCH AT KEY [" << key << "]" << RESET << std::endl;
    }
}

TEST(SAM_PARSER, HEADER_ONLY_YES_REPEAT) {
    // EXP
    std::unordered_map<std::string, std::vector<std::string>> exp_headers {
        {"AA", {"SOMETHING_IN_HEADER", "SPOOKY SAME HEADER"}},
        {"AB", {"BUT THIS IS DIFFERENT D:", "ANOTHER AB!"}},
        {"AC", {"LONELY AC D:"}}
    };
    // Test
    std::string file = "infiles/sam_parser_dupe_header_in.txt";
    std::vector<sam_read*> reads = read_sam(file);
    std::unordered_map<std::string, std::vector<std::string>> headers = read_sam_headers(file);
    EXPECT_TRUE(reads.empty()) << RED << "READ FROM FILE WITH ONLY HEADERS/NO READS" << RESET << std::endl;
    EXPECT_EQ(headers.size(), 3) << RED << "GOT != 3 HEADER CATEGORIES FROM FILE WITH 3 CATS" << RESET << std::endl;
    for(auto& [key, vec] : headers) {
        EXPECT_NE(exp_headers.find(key), exp_headers.end()) << RED << "UNEXPECTED HEADER DELIM PARSED" << RESET << std::endl;
        EXPECT_EQ(vec.size(), exp_headers[key].size()) << RED << "KEY [" << key << "] HAS MISMATCH HEADER ASSOC" << RESET << std::endl; 
        for(int i = 0; i < vec.size(); ++i) {
            EXPECT_EQ(vec[i], exp_headers[key][i]) << RED << "HEADER CONTENT MISMATCH AT KEY [" << key << "][" << i << "]" << RESET << std::endl;
        }
    }
}

TEST(SAM_PARSER, BOTH_READS_AND_HEADER) {
    // EXP
    std::unordered_map<std::string, std::vector<std::string>> exp_headers {
        {"CA", {"HEADER1", "HEADER3"}},
        {"TZ", {"HEADER2", "HEADER4"}},
        {"FG", {"HEADER5"}}
    };
    std::vector<std::vector<std::string>> exp_tags = {
        {"TAG1"}, {"TAG2", "TAG2A"}, {"TAG3", "TAG3A", "TAG3B"}
    };
    std::vector<std::string> exp_qname = {"ID1", "ID2", "ID3"};
    std::vector<uint16_t> exp_flags = {1, 2, 3};
    std::vector<std::string> exp_rname = {"REF1", "REF2", "REF3"};
    std::vector<size_t> exp_pos = {10, 20, 30};
    std::vector<char> exp_mapq = {'A', 'B', 'C'};
    std::vector<std::string> exp_cigar = {"AIGAR", "BIGAR", "CIGAR"};
    std::vector<std::string> exp_rnext = {"=", "*", "!"};
    std::vector<size_t> exp_posnext = {100, 200, 300};
    std::vector<size_t> exp_tlen = {1000, 2000, 3000};
    std::vector<std::string> exp_seq = {"AAA", "BBB", "CCC"};
    std::vector<std::string> exp_qual = {"???", "?*!", "!!!"};
    // Test header
    std::string file = "infiles/sam_parser_all_allowed_in.txt";
    std::unordered_map<std::string, std::vector<std::string>> headers = read_sam_headers(file);
    EXPECT_EQ(headers.size(), 3) << RED << "GOT != 3 HEADER CATEGORIES FROM FILE WITH 3 CATS" << RESET << std::endl;
    for(auto& [key, vec] : headers) {
        EXPECT_NE(exp_headers.find(key), exp_headers.end()) << RED << "UNEXPECTED HEADER DELIM PARSED" << RESET << std::endl;
        EXPECT_EQ(vec.size(), exp_headers[key].size()) << RED << "KEY [" << key << "] HAS MISMATCH HEADER ASSOC" << RESET << std::endl; 
        for(int i = 0; i < vec.size(); ++i) {
            EXPECT_EQ(vec[i], exp_headers[key][i]) << RED << "HEADER CONTENT MISMATCH AT KEY [" << key << "][" << i << "]" << RESET << std::endl;
        }
    }
    // Test read
    std::vector<sam_read*> reads = read_sam(file);
    EXPECT_EQ(reads.size(), 3) << RED << "GOT != 3 READ FROM FILE WITH 3 READS" << std::endl;
    for(int i = 0; i < 3; ++i) {
        auto [tags, qname, rname, cigar, rnext, seq, qual, flags, pos, posnext, tlen, mapq] = (*reads[i]);
        EXPECT_EQ(tags, exp_tags[i]) << RED << "WRONG/NO TAGS CAPTURED FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(qname, exp_qname[i]) << RED << "WRONG QNAME FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(flags, exp_flags[i]) << RED << "WRONG FLAGS FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(rname, exp_rname[i]) << RED << "WRONG REFNAME FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(pos, exp_pos[i]) << RED << "WRONG POS FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(mapq, exp_mapq[i]) << RED << "WRONG MAPQ FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(cigar, exp_cigar[i]) << RED << "WRONG CIGAR STRING FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(rnext, exp_rnext[i]) << RED << "WRONG RNEXT FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(posnext, exp_posnext[i]) << RED << "WRONG POSNEXT FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(tlen, exp_tlen[i]) << RED << "WRONG TLEN FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(seq, exp_seq[i]) << RED << "WRONG TLEN FOR READ [" << i << "]" << RESET << std::endl;
        EXPECT_EQ(qual, exp_qual[i]) << RED << "WRONG QUALITY STRING FOR READ [" << i << "]" << RESET << std::endl;
    }
}
