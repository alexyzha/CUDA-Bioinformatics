#include <gtest/gtest.h>
#include "../../src/cpu/headers/fx_util.h"

TEST(FX_UTIL, FILTER_FQ_AVERAGE_DISCARD_ALL) {
    std::vector<fq_read*> reads;
    std::vector<fq_read*> ret;
    reads.push_back(new fq_read("ID1", 5, "ACGTA", "?????", ""));
    reads.push_back(new fq_read("ID2", 5, "GATCA", "???!!", ""));
    ret = filter_fq(reads, AVERAGE_DISCARD_WHOLE, '?');
    EXPECT_EQ(ret.size(), 1) << RED << "EXP ONLY 1 READ" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_id(), "ID1") << RED << "WRONG ID" << RESET << std::endl;
    EXPECT_EQ(reads[0]->size(), 5) << RED << "WRONG SIZE" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_seq(), "ACGTA") << RED << "WRONG SEQ" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_quality(), "?????") << RED << "WRONG QUALITY" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_metadata(), "") << RED << "METADATA" << RESET << std::endl;
}

TEST(FX_UTIL, FILTER_FQ_SINGLE_DISCARD_ALL) {
    std::vector<fq_read*> reads;
    std::vector<fq_read*> ret;
    reads.push_back(new fq_read("ID1", 5, "ACGTA", "?????", ""));
    reads.push_back(new fq_read("ID2", 5, "GATCA", "???!?", ""));
    reads.push_back(new fq_read("ID3", 5, "CATGA", "??!??", ""));
    reads.push_back(new fq_read("ID4", 5, "ATCGA", "?!???", ""));
    reads.push_back(new fq_read("ID5", 5, "GCTAG", "!???!", ""));
    ret = filter_fq(reads, SINGLE_DISCARD_WHOLE, '?');
    EXPECT_EQ(ret.size(), 1) << RED << "EXP ONLY 1 READ" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_id(), "ID1") << RED << "WRONG ID" << RESET << std::endl;
    EXPECT_EQ(reads[0]->size(), 5) << RED << "WRONG SIZE" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_seq(), "ACGTA") << RED << "WRONG SEQ" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_quality(), "?????") << RED << "WRONG QUALITY" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_metadata(), "") << RED << "METADATA" << RESET << std::endl;
}

TEST(FX_UTIL, FILTER_FQ_SLIDING_WINDOW) {
    std::vector<fq_read*> reads;
    std::vector<fq_read*> ret;
    reads.push_back(new fq_read("ID1", 5, "GTAGC", "??!!?", ""));
    ret = filter_fq(reads, SLIDING_WINDOW, '?', 2.0);
    EXPECT_EQ(ret.size(), 1) << RED << "EXP ONLY 1 READ" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_id(), "ID1") << RED << "WRONG ID" << RESET << std::endl;
    EXPECT_EQ(reads[0]->size(), 2) << RED << "WRONG SIZE" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_seq(), "GT") << RED << "WRONG SEQ" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_quality(), "??") << RED << "WRONG QUALITY" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_metadata(), "") << RED << "METADATA" << RESET << std::endl;
}

TEST(FX_UTIL, FILTER_FQ_PROPORTION_DISCARD_WHOLE) {
    std::vector<fq_read*> reads;
    std::vector<fq_read*> ret;
    reads.push_back(new fq_read("ID1", 5, "AGTCA", "??!!!", ""));
    reads.push_back(new fq_read("ID2", 5, "GTAGC", "??!??", ""));
    ret = filter_fq(reads, PROPORTION_DISCARD_WHOLE, '?', 0.5);
    EXPECT_EQ(ret.size(), 1) << RED << "EXP ONLY 1 READ" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_id(), "ID2") << RED << "WRONG ID" << RESET << std::endl;
    EXPECT_EQ(reads[0]->size(), 5) << RED << "WRONG SIZE" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_seq(), "GTAGC") << RED << "WRONG SEQ" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_quality(), "??!??") << RED << "WRONG QUALITY" << RESET << std::endl;
    EXPECT_EQ(reads[0]->get_metadata(), "") << RED << "METADATA" << RESET << std::endl;
}

TEST(FX_UTIL, GC_PER_READ) {
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
    ret = gc_per_read(reads);
    for(int i = 0; i < 5; ++i) {
        EXPECT_EQ(ret[i], exp[i]) << RED << "EXP[" << i << "] != RET[" << i << "]" << RESET << std::endl;
    }
}

TEST(FX_UTIL, GC_GLOBAL) {
    std::vector<fq_read*> reads;
    double ret;
    double exp = static_cast<double>(2 + 3 + 4 + 5 + 5) / (4 + 5 + 6 + 7 + 8);
    reads.push_back(new fq_read("ID1", 4, "AGTC", "??!!", ""));
    reads.push_back(new fq_read("ID2", 5, "GTAGC", "??!??", ""));
    reads.push_back(new fq_read("ID3", 6, "GGAGCT", "??!?!?", ""));
    reads.push_back(new fq_read("ID4", 7, "GGAGCTG", "??!?!??", ""));
    reads.push_back(new fq_read("ID5", 8, "GGAGCTGA", "??!?!??!", ""));
    ret = gc_global(reads);
    EXPECT_EQ(ret, exp) << RED << "GLOBAL GC COUNT OFF" << RESET << std::endl;
}

TEST(FX_UTIL, COUNT_KMER_UNDER_LENGTH) {
    std::vector<fq_read*> reads;
    std::unordered_map<uint64_t, uint64_t> res;
    for(int i = 0; i < 5; ++i) {
        reads.push_back(
            new fq_read("ID" + std::to_string(i),
            4,
            [&]() {
                std::string ret = "";
                for(int j = 0; j < 4; ++j) {
                    ret.push_back(static_cast<char>(rand() % 0xff));
                }
                return ret;
            }(),
            "????",
            "")
        );
    }
    res = count_kmer(reads, 5);
    EXPECT_TRUE(res.empty()) << RED << "READ A K-MER UNDER LENGTH 5" << RESET << std::endl;
}

TEST(FX_UTIL, COUNT_KMER_K_OOB) {
    std::vector<fq_read*> reads;
    std::unordered_map<uint64_t, uint64_t> res;
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
        reads.push_back(
            new fq_read("ID" + std::to_string(i),
            3,
            seqs[i],
            "???",
            "")
        );
    }
    res = count_kmer(reads, 3);
    EXPECT_EQ(res.size(), 4) << RED << "MORE THAN 4 UNIQUE 3-MERS INDEXED" << RESET << std::endl;
    for(auto& [key, val] : exp) {
        EXPECT_NE(res.find(key), res.end()) << RED << "RES DIDN'T INDEX: " << key << RESET << std::endl;
        EXPECT_EQ(res[key], val) << RED << "RES INDEXED WRONG COUNT OF: " << key << RESET << std::endl;
    }
}

TEST(FX_UTIL, COUNT_KMER_NORMAL_2) {
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
        reads.push_back(
            new fq_read("ID" + std::to_string(i),
            3,
            seqs[i],
            "???",
            "")
        );
    }
    res = count_kmer(reads, 2);
    EXPECT_EQ(res.size(), 6) << RED << "NOT 6 UNIQUE 2-MERS INDEXED" << RESET << std::endl;
    for(auto& [key, val] : exp) {
        EXPECT_NE(res.find(key), res.end()) << RED << "RES DIDN'T INDEX: " << key << RESET << std::endl;
        EXPECT_EQ(res[key], val) << RED << "RES INDEXED WRONG COUNT OF: " << key << RESET << std::endl;
    }
}

TEST(FX_UTIL, INDEX_KMER_UNDER_LENGTH) {
    std::vector<fq_read*> reads;
    std::unordered_map<uint64_t, std::unordered_set<int>> res;
    for(int i = 0; i < 5; ++i) {
        reads.push_back(
            new fq_read("ID" + std::to_string(i),
            4,
            [&]() {
                std::string ret = "";
                for(int j = 0; j < 4; ++j) {
                    ret.push_back(static_cast<char>(rand() % 0xff));
                }
                return ret;
            }(),
            "????",
            "")
        );
    }
    res = index_kmer(reads, 5);
    EXPECT_TRUE(res.empty()) << RED << "READ A K-MER UNDER LENGTH 5" << RESET << std::endl;
}

TEST(FX_UTIL, INDEX_KMER_K_OOB) {
    std::vector<fq_read*> reads;
    std::unordered_map<uint64_t, std::unordered_set<int>> res;
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
        reads.push_back(
            new fq_read("ID" + std::to_string(i),
            3,
            seqs[i],
            "???",
            "")
        );
    }
    res = index_kmer(reads, 3);
    EXPECT_EQ(res.size(), 4) << RED << "MORE THAN 4 UNIQUE 3-MERS INDEXED" << RESET << std::endl;
    for(auto& [key, val] : exp) {
        EXPECT_NE(res.find(key), res.end()) << RED << "RES DIDN'T INDEX: " << key << RESET << std::endl;
        EXPECT_EQ(val.size(), res[key].size()) << RED << "RES INDEXED A DIFFERENT COUNT OF KEY: " << key << RESET << std::endl;
        for(auto& index : val) {
            EXPECT_TRUE(res[key].count(index)) << RED << "RES AT KEY: " << key << " DIDN'T GET INDEX: " << index << RESET << std::endl;
        }
    }
}

TEST(FX_UTIL, INDEX_KMER_NORMAL_2) {
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
        reads.push_back(
            new fq_read("ID" + std::to_string(i),
            3,
            seqs[i],
            "???",
            "")
        );
    }
    res = index_kmer(reads, 2);
    EXPECT_EQ(res.size(), 6) << RED << "MORE THAN 6 UNIQUE 2-MERS INDEXED" << RESET << std::endl;
    for(auto& [key, val] : exp) {
        EXPECT_NE(res.find(key), res.end()) << RED << "RES DIDN'T INDEX: " << key << RESET << std::endl;
        EXPECT_EQ(val.size(), res[key].size()) << RED << "RES INDEXED A DIFFERENT COUNT OF KEY: " << key << RESET << std::endl;
        for(auto& index : val) {
            EXPECT_TRUE(res[key].count(index)) << RED << "RES AT KEY: " << key << " DIDN'T GET INDEX: " << index << RESET << std::endl;
        }
    }
}

