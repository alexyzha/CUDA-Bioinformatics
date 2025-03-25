#include <gtest/gtest.h>
#include "../../src/cpu/headers/util_structs.h"

TEST(UTIL_S, ALIGNMENT) {
    std::string ref = "ABC";
    std::string read = "123";
    alignment test_align = {1, 2, 3, ref, read};
    EXPECT_EQ(test_align.score, 1) << RED << "WRONG ALIGNMENT SCORE" << RESET << std::endl;
    EXPECT_EQ(test_align.end_ref, 2) << RED << "WRONG END REF" << RESET << std::endl;
    EXPECT_EQ(test_align.end_read, 3) << RED << "WRONG END READ" << RESET << std::endl;
    EXPECT_EQ(test_align.aligned_ref, ref) << RED << "WRONG REF SEQ" << RESET << std::endl;
    EXPECT_EQ(test_align.aligned_read, read) << RED << "WRONG READ SEQ" << RESET << std::endl;
}

TEST(UTIL_S, UNION_FIND) {
    std::vector<std::pair<int, int>> pairs = {
        {0, 1}, {2, 5}, {8, 11}, {10, 14},
        {3, 4}, {6, 7}, {9, 15}, {12, 13},
        {1, 2}, {5, 8}, {11, 10}, {4, 6}
    };
    std::vector<std::pair<int, int>> exp_con = {
        {1, 2}, {5, 1}, {8, 1}, {11, 2},
        {14, 1}, {14, 2}, {14, 8}, {3, 7}
    };
    std::vector<std::pair<int, int>> exp_sep = {
        {9, 1}, {9, 2}, {15, 3}, {15, 4},
        {12, 5}, {12, 6}, {13, 7}, {13, 8}
    };
    union_find uf;
    for(auto& [x, y] : pairs) {
        uf.join(x, y);
    }
    for(auto& [x, y] : exp_con) {
        EXPECT_TRUE(uf.con(x, y)) << RED << "EXP CON, GOT NOT CON" << RESET << std::endl;
    }
    for(auto& [x, y] : exp_sep) {
        EXPECT_FALSE(uf.con(x, y)) << RED << "EXP NOT CON, GOT CON" << RESET << std::endl;
    }
}