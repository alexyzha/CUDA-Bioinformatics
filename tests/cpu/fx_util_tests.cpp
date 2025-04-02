#include <gtest/gtest.h>
#include "../../src/cpu/headers/fx_util.h"

TEST(FX_UTIL, FILTER_FQ_AVERAGE_DISCARD_ALL) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::vector<fq_read*> ret;
    reads.push_back(new fq_read("ID1", 5, "ACGTA", "?????", ""));
    reads.push_back(new fq_read("ID2", 5, "GATCA", "???!!", ""));
    
    // Function being tested
    ret = filter_fq(reads, AVERAGE_DISCARD_WHOLE, '?');

    // Test against manually entered expected values
    EXPECT_EQ(ret.size(), 1) << RED << "EXP ONLY 1 READ" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_id(), "ID1") << RED << "WRONG ID" << RESET << std::endl;
    EXPECT_EQ(ret[0]->size(), 5) << RED << "WRONG SIZE" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_seq(), "ACGTA") << RED << "WRONG SEQ" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_quality(), "?????") << RED << "WRONG QUALITY" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_metadata(), "") << RED << "METADATA" << RESET << std::endl;
}

TEST(FX_UTIL, FILTER_FQ_SINGLE_DISCARD_ALL) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::vector<fq_read*> ret;
    reads.push_back(new fq_read("ID1", 5, "ACGTA", "?????", ""));
    reads.push_back(new fq_read("ID2", 5, "GATCA", "???!?", ""));
    reads.push_back(new fq_read("ID3", 5, "CATGA", "??!??", ""));
    reads.push_back(new fq_read("ID4", 5, "ATCGA", "?!???", ""));
    reads.push_back(new fq_read("ID5", 5, "GCTAG", "!???!", ""));
    
    // Function being tested
    ret = filter_fq(reads, SINGLE_DISCARD_WHOLE, '?');

    // Test against manually entered expected values
    EXPECT_EQ(ret.size(), 1) << RED << "EXP ONLY 1 READ" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_id(), "ID1") << RED << "WRONG ID" << RESET << std::endl;
    EXPECT_EQ(ret[0]->size(), 5) << RED << "WRONG SIZE" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_seq(), "ACGTA") << RED << "WRONG SEQ" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_quality(), "?????") << RED << "WRONG QUALITY" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_metadata(), "") << RED << "METADATA" << RESET << std::endl;
}

TEST(FX_UTIL, FILTER_FQ_SLIDING_WINDOW) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::vector<fq_read*> ret;
    reads.push_back(new fq_read("ID1", 5, "GTAGC", "??!!?", ""));
    
    // Function being tested
    ret = filter_fq(reads, SLIDING_WINDOW, '?', 2.0);
    
    // Test against manually entered expected values
    EXPECT_EQ(ret.size(), 1) << RED << "EXP ONLY 1 READ" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_id(), "ID1") << RED << "WRONG ID" << RESET << std::endl;
    EXPECT_EQ(ret[0]->size(), 2) << RED << "WRONG SIZE" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_seq(), "GT") << RED << "WRONG SEQ" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_quality(), "??") << RED << "WRONG QUALITY" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_metadata(), "") << RED << "METADATA" << RESET << std::endl;
}

TEST(FX_UTIL, FILTER_FQ_PROPORTION_DISCARD_WHOLE) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::vector<fq_read*> ret;
    reads.push_back(new fq_read("ID1", 5, "AGTCA", "??!!!", ""));
    reads.push_back(new fq_read("ID2", 5, "GTAGC", "??!??", ""));
    
    // Function being tested
    ret = filter_fq(reads, PROPORTION_DISCARD_WHOLE, '?', 0.5);

    // Test against manually entered expected values
    EXPECT_EQ(ret.size(), 1) << RED << "EXP ONLY 1 READ" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_id(), "ID2") << RED << "WRONG ID" << RESET << std::endl;
    EXPECT_EQ(ret[0]->size(), 5) << RED << "WRONG SIZE" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_seq(), "GTAGC") << RED << "WRONG SEQ" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_quality(), "??!??") << RED << "WRONG QUALITY" << RESET << std::endl;
    EXPECT_EQ(ret[0]->get_metadata(), "") << RED << "METADATA" << RESET << std::endl;
}

TEST(FX_UTIL, GC_PER_READ) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::vector<double> ret;
    std::vector<double> exp = {
        static_cast<double>(2) / 4,
        static_cast<double>(3) / 5,
        static_cast<double>(4) / 6,
        static_cast<double>(5) / 7,
        static_cast<double>(5) / 8
    };
    reads.push_back(new fq_read("ID1", 4, "AGTC", "??!!", ""));
    reads.push_back(new fq_read("ID2", 5, "GTAGC", "??!??", ""));
    reads.push_back(new fq_read("ID3", 6, "GGAGCT", "??!?!?", ""));
    reads.push_back(new fq_read("ID4", 7, "GGAGCTG", "??!?!??", ""));
    reads.push_back(new fq_read("ID5", 8, "GGAGCTGA", "??!?!??!", ""));
    
    // Function being tested
    ret = gc_per_read(reads);
    
    // Validate function return values 
    for(int i = 0; i < 5; ++i) {
        EXPECT_EQ(ret[i], exp[i]) << RED << "EXP[" << i << "] != RET[" << i << "]" << RESET << std::endl;
    }
}

TEST(FX_UTIL, GC_GLOBAL) {
    // Create testing variables
    std::vector<fq_read*> reads;
    double ret;
    double exp = static_cast<double>(2 + 3 + 4 + 5 + 5) / (4 + 5 + 6 + 7 + 8);
    reads.push_back(new fq_read("ID1", 4, "AGTC", "??!!", ""));
    reads.push_back(new fq_read("ID2", 5, "GTAGC", "??!??", ""));
    reads.push_back(new fq_read("ID3", 6, "GGAGCT", "??!?!?", ""));
    reads.push_back(new fq_read("ID4", 7, "GGAGCTG", "??!?!??", ""));
    reads.push_back(new fq_read("ID5", 8, "GGAGCTGA", "??!?!??!", ""));
    
    // Function being tested
    ret = gc_global(reads);
    
    // Validate function return value
    EXPECT_EQ(ret, exp) << RED << "GLOBAL GC COUNT OFF" << RESET << std::endl;
}

TEST(FX_UTIL, COUNT_KMER_UNDER_LENGTH) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::unordered_map<uint64_t, uint64_t> res;
    for(int i = 0; i < 5; ++i) {
        reads.push_back(new fq_read(
            "ID" + std::to_string(i),
            4,
            [&]() {
                std::string ret = "";
                for(int j = 0; j < 4; ++j) {
                    ret.push_back(static_cast<char>(rand() % 0xff));
                }
                return ret;
            }(),
            "????"
        ));
    }

    // Function being tested
    res = count_kmer(reads, 5);

    // Validate function return value
    EXPECT_TRUE(res.empty()) << RED << "READ A K-MER UNDER LENGTH 5" << RESET << std::endl;
}

TEST(FX_UTIL, COUNT_KMER_K_OOB) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::unordered_map<uint64_t, uint64_t> res;
    
    // Expects function [count_kmer()] to throw with k out of bound
    EXPECT_THROW([&](){
        res = count_kmer(reads, 40);
    }(), std::invalid_argument) << RED << "EXP THROW WITH K = 40 > 32" << RESET << std::endl;
    EXPECT_THROW([&](){
        res = count_kmer(reads, -1);
    }(), std::invalid_argument) << RED << "EXP THROW WITH K = -1 <= 0" << RESET << std::endl;
    EXPECT_THROW([&](){
        res = count_kmer(reads, 0);
    }(), std::invalid_argument) << RED << "EXP THROW WITH K = 0 <= 0" << RESET << std::endl;
}

TEST(FX_UTIL, COUNT_KMER_NORMAL_1) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::vector<std::string> seqs = {
        "ACG", "GAT", "CAT", "CAG", "ACG"
    };
    std::unordered_map<uint64_t, uint64_t> exp = {
        {0b000110, 2}, {0b100011, 1},
        {0b010011, 1}, {0b010010, 1}
    };
    std::unordered_map<uint64_t, uint64_t> res;
    for(int i = 0; i < 5; ++i) {
        reads.push_back(new fq_read(
            "ID" + std::to_string(i),
            3,
            seqs[i],
            "???"    
        ));
    }

    // Function being tested
    res = count_kmer(reads, 3);

    // Validate function output
    EXPECT_EQ(res.size(), 4) << RED << "MORE THAN 4 UNIQUE 3-MERS INDEXED" << RESET << std::endl;
    for(auto& [key, val] : exp) {
        EXPECT_NE(res.find(key), res.end()) << RED << "RES DIDN'T INDEX: " << meow(key) << RESET << std::endl;
        EXPECT_EQ(res[key], val) << RED << "RES INDEXED WRONG COUNT OF: " << meow(key) << RESET << std::endl;
    }
}

TEST(FX_UTIL, COUNT_KMER_NORMAL_2) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::vector<std::string> seqs = {
        "ACG", "GAT", "CAT", "CAG", "ACG"
    };
    std::unordered_map<uint64_t, uint64_t> exp = {
        {0b0001, 2}, {0b0110, 2}, {0b1000, 1},
        {0b0011, 2}, {0b0100, 2}, {0b0010, 1}
    };
    std::unordered_map<uint64_t, uint64_t> res;
    for(int i = 0; i < 5; ++i) {
        reads.push_back(new fq_read(
            "ID" + std::to_string(i),
            3,
            seqs[i],
            "???"
        ));
    }

    // Function being tested
    res = count_kmer(reads, 2);

    // Validate function output
    EXPECT_EQ(res.size(), 6) << RED << "NOT 6 UNIQUE 2-MERS INDEXED" << RESET << std::endl;
    for(auto& [key, val] : exp) {
        EXPECT_NE(res.find(key), res.end()) << RED << "RES DIDN'T INDEX: " << meow(key) << RESET << std::endl;
        EXPECT_EQ(res[key], val) << RED << "RES INDEXED WRONG COUNT OF: " << meow(key) << RESET << std::endl;
    }
}

TEST(FX_UTIL, INDEX_KMER_UNDER_LENGTH) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::unordered_map<uint64_t, std::unordered_set<int>> res;
    for(int i = 0; i < 5; ++i) {
        reads.push_back(new fq_read(
            "ID" + std::to_string(i),
            4,
            [&]() {
                std::string ret = "";
                for(int j = 0; j < 4; ++j) {
                    ret.push_back(static_cast<char>(rand() % 0xff));
                }
                return ret;
            }(),
            "????"
        ));
    }

    // Function being tested
    res = index_kmer(reads, 5);

    // Validate function output
    EXPECT_TRUE(res.empty()) << RED << "READ A K-MER UNDER LENGTH 5" << RESET << std::endl;
}

TEST(FX_UTIL, INDEX_KMER_K_OOB) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::unordered_map<uint64_t, std::unordered_set<int>> res;
    
    // Expects function [index_kmer()] to throw with k out of bound
    EXPECT_THROW([&](){
        res = index_kmer(reads, 69);
    }(), std::invalid_argument) << RED << "EXP THROW WITH K = 69 > 32" << RESET << std::endl;
    EXPECT_THROW([&](){
        res = index_kmer(reads, -1);
    }(), std::invalid_argument) << RED << "EXP THROW WITH K = -1 <= 0" << RESET << std::endl;
    EXPECT_THROW([&](){
        res = index_kmer(reads, 0);
    }(), std::invalid_argument) << RED << "EXP THROW WITH K = 0 <= 0" << RESET << std::endl;
}

TEST(FX_UTIL, INDEX_KMER_NORMAL_1) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::vector<std::string> seqs = {
        "ACG", "GAT", "CAT", "CAG", "ACG"
    };
    std::unordered_map<uint64_t, std::unordered_set<int>> exp = {
        {0b000110, {0, 4}}, {0b100011, {1}},
        {0b010011, {2}}, {0b010010, {3}}
    };
    std::unordered_map<uint64_t, std::unordered_set<int>> res;
    for(int i = 0; i < 5; ++i) {
        reads.push_back(new fq_read(
            "ID" + std::to_string(i),
            3,
            seqs[i],
            "???"
        ));
    }
    
    // Function being tested
    res = index_kmer(reads, 3);

    // Validate function output
    EXPECT_EQ(res.size(), 4) << RED << "MORE THAN 4 UNIQUE 3-MERS INDEXED" << RESET << std::endl;
    for(auto& [key, val] : exp) {
        EXPECT_NE(res.find(key), res.end()) << RED << "RES DIDN'T INDEX: " << meow(key) << RESET << std::endl;
        EXPECT_EQ(val.size(), res[key].size()) << RED << "RES INDEXED A DIFFERENT COUNT OF KEY: " << meow(key) << RESET << std::endl;
        for(auto& index : val) {
            EXPECT_TRUE(res[key].count(index)) << RED << "RES AT KEY: " << meow(key) << " DIDN'T GET INDEX: " << index << RESET << std::endl;
        }
    }
}

TEST(FX_UTIL, INDEX_KMER_NORMAL_2) {
    // Create testing variables
    std::vector<fq_read*> reads;
    std::vector<std::string> seqs = {
        "ACG", "GAT", "CAT", "CAG", "ACG"
    };
    std::unordered_map<uint64_t, std::unordered_set<int>> exp = {
        {0b0001, {0, 4}}, {0b0110, {0, 4}}, {0b1000, {1}},
        {0b0011, {1, 2}}, {0b0100, {2, 3}}, {0b0010, {3}}
    };
    std::unordered_map<uint64_t, std::unordered_set<int>> res;
    for(int i = 0; i < 5; ++i) {
        reads.push_back(new fq_read(
            "ID" + std::to_string(i),
            3,
            seqs[i],
            "???"
        ));
    }

    // Function being tested
    res = index_kmer(reads, 2);

    // Validate function output
    EXPECT_EQ(res.size(), 6) << RED << "MORE THAN 6 UNIQUE 2-MERS INDEXED" << RESET << std::endl;
    for(auto& [key, val] : exp) {
        EXPECT_NE(res.find(key), res.end()) << RED << "RES DIDN'T INDEX: " << meow(key) << RESET << std::endl;
        EXPECT_EQ(val.size(), res[key].size()) << RED << "RES INDEXED A DIFFERENT COUNT OF KEY: " << meow(key) << RESET << std::endl;
        for(auto& index : val) {
            EXPECT_TRUE(res[key].count(index)) << RED << "RES AT KEY: " << meow(key) << " DIDN'T GET INDEX: " << index << RESET << std::endl;
        }
    }
}

TEST(FX_UTIL, LOCAL_ALIGN_MATCH) {
    // Create testing variables
    std::string ref = "GATCATC";
    std::string read = "GATC";

    // Function being tested
    alignment align = local_align(ref, read);

    // Validate function output
    EXPECT_EQ(align.score, 8) << RED << "CHECK ALIGNMENT SCORE CALC" << RESET << std::endl;
    EXPECT_EQ(align.end_ref, 3) << RED << "END REFERENCE INDEX INCORRECT" << RESET << std::endl;
    EXPECT_EQ(align.end_read, 3) << RED << "END READ INDEX INCORRECT" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_ref), read) << RED << "PROCESSED REFERENCE SEQ INCORRECT" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_read), read) << RED << "PROCESSED READ SEQ INCORRECT" << RESET << std::endl;
}

TEST(FX_UTIL, LOCAL_ALIGN_MATCH_MIDDLE) {
    // Create testing variables
    std::string ref = "GAAAAGTC";
    std::string read = "AAAA";

    // Function being tested
    alignment align = local_align(ref, read);

    // Validate function output
    EXPECT_EQ(align.score, 8) << RED << "CHECK ALIGNMENT SCORE CALC" << RESET << std::endl;
    EXPECT_EQ(align.end_ref, 4) << RED << "END REFERENCE INDEX INCORRECT" << RESET << std::endl;
    EXPECT_EQ(align.end_read, 3) << RED << "END READ INDEX INCORRECT" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_ref), read) << RED << "PROCESSED REFERENCE SEQ INCORRECT" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_read), read) << RED << "PROCESSED READ SEQ INCORRECT" << RESET << std::endl;
}

TEST(FX_UTIL, LOCAL_ALIGN_MISMATCH) {
    // Create testing variables
    std::string ref = "GAAAAGTC";
    std::string read = "GAACA";

    // Function being tested
    alignment align = local_align(ref, read);

    // Validate function output
    EXPECT_EQ(align.score, 7) << RED << "CHECK ALIGNMENT SCORE CALC" << RESET << std::endl;
    EXPECT_EQ(align.end_ref, 4) << RED << "END REFERENCE INDEX INCORRECT" << RESET << std::endl;
    EXPECT_EQ(align.end_read, 4) << RED << "END READ INDEX INCORRECT" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_ref), "GAAAA") << RED << "PROCESSED REFERENCE SEQ INCORRECT" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_read), read) << RED << "PROCESSED READ SEQ INCORRECT" << RESET << std::endl;
}

TEST(FX_UTIL, LOCAL_ALIGN_GAP_IN_READ) {
    // Create testing variables
    std::string ref = "GCTAGCT";
    std::string read = "GCTGCT";

    // Function being tested
    alignment align = local_align(ref, read);

    // Validate function output
    EXPECT_EQ(align.score, 10) << RED << "CHECK ALIGNMENT SCORE CALC" << RESET << std::endl;
    EXPECT_EQ(align.end_ref, 6) << RED << "END REFERENCE INDEX INCORRECT" << RESET << std::endl;
    EXPECT_EQ(align.end_read, 5) << RED << "END READ INDEX INCORRECT" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_ref), ref) << RED << "PROCESSED REFERENCE SEQ INCORRECT" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_read), "GCT-GCT") << RED << "PROCESSED READ SEQ INCORRECT" << RESET << std::endl;
}

TEST(FX_UTIL, LOCAL_ALIGN_GAP_IN_REF) {
    // Create testing variables
    std::string ref = "ATAGCT";
    std::string read = "ATATGCT";

    // Function being tested
    alignment align = local_align(ref, read);

    // Validate function output
    EXPECT_EQ(align.score, 10) << RED << "CHECK ALIGNMENT SCORE CALC" << RESET << std::endl;
    EXPECT_EQ(align.end_ref, 5) << RED << "END REFERENCE INDEX INCORRECT" << RESET << std::endl;
    EXPECT_EQ(align.end_read, 6) << RED << "END READ INDEX INCORRECT" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_ref), "ATA-GCT") << RED << "PROCESSED REFERENCE SEQ INCORRECT" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_read), read) << RED << "PROCESSED READ SEQ INCORRECT" << RESET << std::endl;
}

TEST(FX_UTIL, LOCAL_ALIGN_MULTI_ALIGN_GAP) {
    // Create testing variables
    std::string ref = "ATAGCT";
    std::string read = "ATAAGCT";

    // Function being tested
    alignment align = local_align(ref, read);

    // Validate function output
    EXPECT_EQ(align.score, 10) << RED << "CHECK ALIGNMENT SCORE CALC" << RESET << std::endl;
    EXPECT_EQ(align.end_ref, 5) << RED << "END REFERENCE INDEX INCORRECT" << RESET << std::endl;
    EXPECT_EQ(align.end_read, 6) << RED << "END READ INDEX INCORRECT" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_ref), "AT-AGCT") << RED << "PROCESSED REFERENCE SEQ INCORRECT" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_read), read) << RED << "PROCESSED READ SEQ INCORRECT" << RESET << std::endl;
}

TEST(FX_UTIL, LOCAL_ALIGN_EMPTY) {
    // Create testing variables
    std::string ref = "";
    std::string read = "";

    // Function being tested
    alignment align = local_align(ref, read);

    // Validate function output
    EXPECT_EQ(align.score, 0) << RED << "EMPTY ALIGN SCORE EXP 0" << RESET << std::endl;
    EXPECT_EQ(align.end_ref, -1) << RED << "EMPTY ALIGN END REF EXP 0" << RESET << std::endl;
    EXPECT_EQ(align.end_read, -1) << RED << "EMPTY ALIGN END READ EXP 0" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_ref), "") << RED << "EMPTY ALIGN REF EXP EMPTY STRING" << RESET << std::endl;
    EXPECT_EQ((*align.aligned_read), "") << RED << "EMPTY ALIGN READ EXP EMPTY STRING" << RESET << std::endl;
}

TEST(FX_UTIL, CLUSTER_BY_KMER) {
    // This is an integration test technically
    // Create testing variables
    std::vector<fq_read*> reads;
    std::vector<std::string> seqs = {
        "ACG", "GAT", "CAT", "CAG", "ACG"
    };
    std::unordered_map<uint64_t, std::unordered_set<int>> kmer_map;
    for(int i = 0; i < 5; ++i) {
        reads.push_back(new fq_read(
            "ID" + std::to_string(i),
            3,
            seqs[i],
            "???"
        ));
    }
    kmer_map = index_kmer(reads, 2);

    // Function being tested
    std::vector<std::unordered_set<int>*> res = cluster_by_kmer(kmer_map, reads.size(), 2);

    // Validate function output
    EXPECT_NE([&]() {
        return res[0] ? res[0] : res[4];
    }(), nullptr) << RED << "EXP SET IN [0] OR [4], NULLPTR EVERYWHERE ELSE" << RESET << std::endl;
    for(int i = 0; i <= 3; ++i) {
        EXPECT_EQ(res[i], nullptr) << RED << "EXP NULLPTR AT RES[" << i << "]" << RESET << std::endl;
    }
    if(res[0]) {
        EXPECT_EQ(res[0]->size(), 1) << RED << "MORE THAN 1 CHILD INDEXED AT RES[0]" << RESET << std::endl;
        EXPECT_TRUE(res[0]->count(4)) << RED << "4 NOT INDEXED AS CHILD AT RES[0]" << RESET << std::endl;
    } else if(res[4]) {
        EXPECT_EQ(res[4]->size(), 1) << RED << "MORE THAN 1 CHILD INDEXED AT RES[4]" << RESET << std::endl;
        EXPECT_TRUE(res[4]->count(0)) << RED << "0 NOT INDEXED AS CHILD AT RES[4]" << RESET << std::endl;
    }
}