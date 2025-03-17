#include <gtest/gtest.h>
#include "../src/fa_read.h"

TEST(FA_READ, CONSTRUCTOR_1) {
    fa_read test_read("ABCDE", 10, "ACGTACGTAC");
    EXPECT_EQ("ABCDE", test_read.get_id()) << "ID MISMATCH" << std::endl;
    EXPECT_EQ(10, test_read.size()) << "READ SIZE MISMATCH" << std::endl;
    EXPECT_EQ("ACGTACGTAC", test_read.get_seq()) << "READ SEQ MISMATCH" << std::endl;
    EXPECT_EQ("", test_read.get_metadata()) << "METADATA MISMATCH" << std::endl;
}

TEST(FA_READ, CONSTRUCTOR_2) {
    fa_read test_read("ALHGKHIJKL", 21, "AOALJKBNQJPIHUOBJNKLM", "SUPER_META_DATA_AHH");
    EXPECT_EQ("ALHGKHIJKL", test_read.get_id()) << "ID MISMATCH" << std::endl;
    EXPECT_EQ(21, test_read.size()) << "READ SIZE MISMATCH" << std::endl;
    EXPECT_EQ("AOALJKBNQJPIHUOBJNKLM", test_read.get_seq()) << "READ SEQ MISMATCH" << std::endl;
    EXPECT_EQ("SUPER_META_DATA_AHH", test_read.get_metadata()) << "METADATA MISMATCH" << std::endl;
}

TEST(FA_READ, OPERATOR_IN_BOUND) {
    std::string seq = "CFHGVJHBIUOJIKLNMBGHFGDTAGCVHGJHKJLKMKNJGUYFGCVHBJ";
    fa_read test_read("TEST", 50, seq);
    for(int i = 0; i < 50; ++i) {
        EXPECT_NO_THROW(test_read[i]) << "VALID INDEX ACCESS THROWS ERROR" << std::endl;
        EXPECT_EQ(test_read[i], seq[i]) << "READ[" << i << "] DOESN'T MATCH TEMPLATE[" << i << "]" << std::endl;
    }
}

TEST(FA_READ, OPERATOR_OOB) {
    fa_read test_read("TEST", 17, "MNBHUYTGFGHJKLMJH");
    EXPECT_THROW(test_read[-1], std::out_of_range) << "EXPECTED OOB EXCEPTION USING -1" << std::endl;
    EXPECT_THROW(test_read[17], std::out_of_range) << "EXPECTED OOB EXCEPTION USING 17" << std::endl;
}

/* TO DO */
TEST(FA_READ, FILE_OUT) {
    
}

TEST(FA_READ, RANDOM_STRESS) {
    std::string seq = "";
    int len = rand() % 100000;
    for(int i = 0; i < len; ++i) {
        seq.push_back(static_cast<char>(rand() % 0xff));
    }
    fa_read test_read("STRESS TEST", len, seq);
    for(int i = 0; i < len; ++i) {
        EXPECT_NO_THROW(test_read[i]) << "VALID INDEX ACCESS THROWS ERROR" << std::endl;
        EXPECT_EQ(test_read[i], seq[i]) << "READ[" << i << "] DOESN'T MATCH TEMPLATE[" << i << "]" << std::endl;
    }
    EXPECT_THROW(test_read[-1], std::out_of_range) << "EXPECTED OOB EXCEPTION USING -1" << std::endl;
    EXPECT_THROW(test_read[len], std::out_of_range) << "EXPECTED OOB EXCEPTION USING 17" << std::endl;
}