#include <gtest/gtest.h>
#include "../../src/cpu/headers/util.h"

TEST(UTIL, SPLIT_BY_EMPTY) {
    std::string test_str = "";
    std::vector<std::string> ret;
    EXPECT_NO_THROW([&](){
        ret = split_by(test_str, ' ');
    }()) << RED << "BAD ASSIGNMENT, EXP NO THROW" << RESET << std::endl;
    EXPECT_TRUE(ret.empty()) << RED << "EXP EMPTY VECTOR" << std::endl;
}

TEST(UTIL, SPLIT_BY_NORMAL) {
    std::string test_str = "ABC DEF GHI JKL";
    std::vector<std::string> exp = {"ABC", "DEF", "GHI", "JKL"};
    std::vector<std::string> ret;
    EXPECT_NO_THROW([&](){
        ret = split_by(test_str, ' ');
    }()) << RED << "BAD ASSIGNMENT, EXP NO THROW" << RESET << std::endl;
    EXPECT_EQ(exp, ret) << RED << "SPLIT != EXP RESULT" << RESET << std::endl;
}

TEST(UTIL, SPLIT_BY_DNE) {
    std::string test_str = "MNO PQR STU";
    std::vector<std::string> ret;
    EXPECT_NO_THROW([&](){
        ret = split_by(test_str, '\t');
    }()) << RED << "BAD ASSIGNMENT, EXP NO THROW" << RESET << std::endl;
    EXPECT_EQ({test_str}, ret) << RED << "SPLIT != EXP RESULT" << RESET << std::endl;
}

TEST(UTIL, TRIM) {
    std::string test_str = " BHKGYLKJ<BJMHGJ \t\n";
    std::string exp = "BHKGYLKJ<BJMHGJ";
    trim(test_str);
    EXPECT_EQ(test_str, exp) << RED << "TRIMMED WRONG" << RESET << std::endl;
}

TEST(UTIL, BASE2BIT) {
    EXPECT_EQ(base_to_bit('A'), 0x00) << RED << "BASE->BIT CONVERTS A WRONG" << RESET << std::endl;
    EXPECT_EQ(base_to_bit('C'), 0x01) << RED << "BASE->BIT CONVERTS C WRONG" << RESET << std::endl;
    EXPECT_EQ(base_to_bit('G'), 0x10) << RED << "BASE->BIT CONVERTS G WRONG" << RESET << std::endl;
    EXPECT_EQ(base_to_bit('T'), 0x11) << RED << "BASE->BIT CONVERTS T WRONG" << RESET << std::endl;
    EXPECT_EQ(base_to_bit('U'), 0x80) << RED << "BASE->BIT CONVERTS *[!=ACGT] WRONG" << RESET << std::endl;
}