#include <gtest/gtest.h>
#include "../../src/cpu/headers/sam_container.h"

TEST(SAM_CONTAINER, CONSTRUCTOR_EMPTY) {
    // Class being tested
    sam_container cont("infiles/empty.txt");
    const auto& headers = cont.get_headers();
    const auto& reads = cont.get_reads();

    // Validate class contents
    EXPECT_TRUE(headers.empty()) << RED << "HEADER FROM EMPTY FILE" << RESET << std::endl;
    EXPECT_TRUE(reads.empty()) << RED << "READ FROM EMPTY FILE" << RESET << std::endl;
}

TEST(SAM_CONTAINER, CONSTRUCTOR_READ_ONLY) {
    // Create testing variables
    sam_container cont("infiles/sam_parser_read_yes_tags_in.txt");
    const auto& headers = cont.get_headers();
    const auto& reads = cont.get_reads();
    
    // Prelim validations
    EXPECT_TRUE(headers.empty()) << RED << "HEADER FROM NO HEADER FILE" << RESET << std::endl;
    EXPECT_EQ(reads.size(), 1) << RED << "EXPECTED 1 READ" << RESET << std::endl;
    
    // Read validation
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

TEST(SAM_CONTAINER, CONSTRUCTOR_MULTI_READ) {
    // Create testing variables
    sam_container cont("infiles/sam_parser_multi_read_in.txt");
    const auto& headers = cont.get_headers();
    const auto& reads = cont.get_reads();
    
    // Prelim validations
    EXPECT_TRUE(headers.empty()) << RED << "HEADER FROM NO HEADER FILE" << RESET << std::endl;
    EXPECT_EQ(reads.size(), 3) << RED << "EXPECTED 3 READS" << RESET << std::endl;
    
    // Read expected
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
    
    // Read validation
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

TEST(SAM_CONTAINER, HEADER_ONLY_NO_DUPES) {
    // Create testing variables
    sam_container cont("infiles/sam_parser_only_header_in.txt");
    const auto& headers = cont.get_headers();
    const auto& reads = cont.get_reads();
    
    // Prelim validations
    EXPECT_TRUE(reads.empty()) << RED << "GOT READS FROM FILE WITH NO READS" << RESET << std::endl;
    EXPECT_EQ(headers.size(), 3) << RED << "EXP 3 HEADERS" << RESET << std::endl;
    
    // Header expected
    std::unordered_map<std::string, std::vector<std::string>> exp_headers {
        {"AA", {"SOMETHING_IN_HEADER"}},
        {"AB", {"SOMETHING ELSE IN HEADER"}},
        {"AC", {"ANOTHER SOMETHING IN HEADER"}}
    };
    
    // Header validation
    for(auto& [key, vec] : headers) {
        EXPECT_NE(exp_headers.find(key), exp_headers.end()) << RED << "UNEXPECTED HEADER DELIM PARSED" << RESET << std::endl;
        EXPECT_EQ(vec.size(), exp_headers[key].size()) << RED << "KEY [" << key << "] HAS != 1 HEADER ASSOC, EXP 1" << RESET << std::endl; 
        EXPECT_EQ(vec[0], exp_headers[key][0]) << RED << "HEADER CONTENT MISMATCH AT KEY [" << key << "]" << RESET << std::endl;
    }    
}

TEST(SAM_CONTAINER, HEADER_ONLY_WITH_DUPES) {
    // Create testing variables
    sam_container cont("infiles/sam_parser_dupe_header_in.txt");
    const auto& headers = cont.get_headers();
    const auto& reads = cont.get_reads();
    
    // Prelim validations
    EXPECT_TRUE(reads.empty()) << RED << "GOT READS FROM FILE WITH NO READS" << RESET << std::endl;
    EXPECT_EQ(headers.size(), 3) << RED << "EXP 3 UNIQUE HEADERS" << RESET << std::endl;
    
    // Header expected
    std::unordered_map<std::string, std::vector<std::string>> exp_headers {
        {"AA", {"SOMETHING_IN_HEADER", "SPOOKY SAME HEADER"}},
        {"AB", {"BUT THIS IS DIFFERENT D:", "ANOTHER AB!"}},
        {"AC", {"LONELY AC D:"}}
    };

    // Header validation
    for(auto& [key, vec] : headers) {
        EXPECT_NE(exp_headers.find(key), exp_headers.end()) << RED << "UNEXPECTED HEADER DELIM PARSED" << RESET << std::endl;
        EXPECT_EQ(vec.size(), exp_headers[key].size()) << RED << "KEY [" << key << "] HAS MISMATCH HEADER ASSOC" << RESET << std::endl; 
        for(int i = 0; i < vec.size(); ++i) {
            EXPECT_EQ(vec[i], exp_headers[key][i]) << RED << "HEADER CONTENT MISMATCH AT KEY [" << key << "][" << i << "]" << RESET << std::endl;
        }
    }
}

TEST(SAM_CONTAINER, BOTH_READS_AND_HEADER) {
    // Create testing variables
    sam_container cont("infiles/sam_parser_all_allowed_in.txt");
    const auto& headers = cont.get_headers();
    const auto& reads = cont.get_reads();
    
    // Prelim validations
    EXPECT_EQ(reads.size(), 3) << RED << "EXPECTED 3 READS" << RESET << std::endl;
    EXPECT_EQ(headers.size(), 3) << RED << "EXP 3 UNIQUE HEADERS" << RESET << std::endl;
    
    // Headers expected
    std::unordered_map<std::string, std::vector<std::string>> exp_headers {
        {"CA", {"HEADER1", "HEADER3"}},
        {"TZ", {"HEADER2", "HEADER4"}},
        {"FG", {"HEADER5"}}
    };

    // Reads expected
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
    
    // Header validation
    for(auto& [key, vec] : headers) {
        EXPECT_NE(exp_headers.find(key), exp_headers.end()) << RED << "UNEXPECTED HEADER DELIM PARSED" << RESET << std::endl;
        EXPECT_EQ(vec.size(), exp_headers[key].size()) << RED << "KEY [" << key << "] HAS MISMATCH HEADER ASSOC" << RESET << std::endl; 
        for(int i = 0; i < vec.size(); ++i) {
            EXPECT_EQ(vec[i], exp_headers[key][i]) << RED << "HEADER CONTENT MISMATCH AT KEY [" << key << "][" << i << "]" << RESET << std::endl;
        }
    }
    // Read validation
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

TEST(SAM_CONTAINER, SORTING) {
    // Create testing variables
    sam_container cont("infiles/sam_parser_multi_read_in.txt");
    
    // Function being tested
    cont.sort([&](sam_read* a, sam_read* b) {
        return a->posnext > b->posnext;
    });
    const auto& headers = cont.get_headers();
    const auto& reads = cont.get_reads();
    
    // Prelim validations
    EXPECT_TRUE(headers.empty()) << RED << "GOT HEADERS FROM FILE WITH NO HEADERS" << RESET << std::endl;
    EXPECT_EQ(reads.size(), 3) << RED << "GOT != 3 READS FROM A FILE WITH 3 READS" << RESET << std::endl;
    
    // Expected ID orders
    std::vector<std::string> exp_id = {"ID2", "ID1", "ID3"};
    
    // Validate function results
    for(int i = 0; i < 3; ++i) {
        EXPECT_EQ(reads[i]->qname, exp_id[i]) << RED << "READS NOT IN SORTED ORDER" << RESET << std::endl;
    }
}