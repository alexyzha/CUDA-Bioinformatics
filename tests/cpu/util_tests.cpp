#include <gtest/gtest.h>
#include "../../src/cpu/headers/util.h"

TEST(UTIL, SPLIT_BY_EMPTY) {
    // Create testing variables
    std::string test_str = "";
    std::vector<std::string> ret;

    // Validate results
    EXPECT_NO_THROW([&](){
        ret = split_by(test_str, ' ');
    }()) << RED << "BAD ASSIGNMENT, EXP NO THROW" << RESET << std::endl;
    EXPECT_TRUE(ret.empty()) << RED << "EXP EMPTY VECTOR" << std::endl;
}

TEST(UTIL, SPLIT_BY_NORMAL) {
    // Create testing variables
    std::string test_str = "ABC DEF GHI JKL";
    std::vector<std::string> exp = {"ABC", "DEF", "GHI", "JKL"};
    std::vector<std::string> ret;
    
    // Validate results
    EXPECT_NO_THROW([&](){
        ret = split_by(test_str, ' ');
    }()) << RED << "BAD ASSIGNMENT, EXP NO THROW" << RESET << std::endl;
    EXPECT_EQ(exp, ret) << RED << "SPLIT != EXP RESULT" << RESET << std::endl;
}

TEST(UTIL, SPLIT_BY_DNE) {
    // Create testing variables
    std::string test_str = "MNO PQR STU";
    std::vector<std::string> ret;
    
    // Validate results
    EXPECT_NO_THROW([&](){
        ret = split_by(test_str, '\t');
    }()) << RED << "BAD ASSIGNMENT, EXP NO THROW" << RESET << std::endl;
    EXPECT_EQ(std::vector<std::string>{test_str}, ret) << RED << "SPLIT != EXP RESULT" << RESET << std::endl;
}

TEST(UTIL, TRIM) {
    // Create testing variables
    std::string test_str = " BHKGYLKJ<BJMHGJ \t\n";
    std::string exp = "BHKGYLKJ<BJMHGJ";

    // Function to be tested
    trim(test_str);

    // Validate function results
    EXPECT_EQ(test_str, exp) << RED << "TRIMMED WRONG" << RESET << std::endl;
}

TEST(UTIL, BASE2BIT) {
    // Validate function results
    EXPECT_EQ(base_to_bit('A'), 0b00) << RED << "BASE->BIT CONVERTS A WRONG" << RESET << std::endl;
    EXPECT_EQ(base_to_bit('C'), 0b01) << RED << "BASE->BIT CONVERTS C WRONG" << RESET << std::endl;
    EXPECT_EQ(base_to_bit('G'), 0b10) << RED << "BASE->BIT CONVERTS G WRONG" << RESET << std::endl;
    EXPECT_EQ(base_to_bit('T'), 0b11) << RED << "BASE->BIT CONVERTS T WRONG" << RESET << std::endl;

    // Convert both masks to unsigned b/c only bits matter
    EXPECT_EQ(static_cast<unsigned char>(base_to_bit('U')), static_cast<unsigned char>(0x80)) << RED << "BASE->BIT CONVERTS *[!=ACGT] WRONG" << RESET << std::endl;
}

TEST(UTIL, MEOW) {
    // Create testing variables
    std::string exp = "";
    int in = 0;
    for(int i = 0; i < 32; ++i) {
        int lsb = rand() % 2;
        in <<= 1;
        in |= lsb;
        exp.push_back('0' + lsb);
    }

    // Validate results
    EXPECT_EQ(exp, meow(in));
}