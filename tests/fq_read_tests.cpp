#include <gtest/gtest.h>
#include "../src/fq_read.h"

TEST(FQ_READ, CONSTRUCTOR_1) {
    fq_read test_read("TESTY CAL", 35, "KLJKHJGHFTUYHJKUIYHUGHFGHJHKJHGHGCF", "QQQLKAJBHHUJIKSLAJNBHKJGUIHJNLASJKH");
    EXPECT_EQ("TESTY CAL", test_read.get_id()) << "ID MISMATCH" << std::endl;
    EXPECT_EQ(35, test_read.size()) << "READ SIZE MISMATCH" << std::endl;
    EXPECT_EQ("KLJKHJGHFTUYHJKUIYHUGHFGHJHKJHGHGCF", test_read.get_seq()) << "READ SEQ MISMATCH" << std::endl;
    EXPECT_EQ("QQQLKAJBHHUJIKSLAJNBHKJGUIHJNLASJKH", test_read.get_quality()) << "READ QUALITY MISMATCH" << std::endl;
    EXPECT_EQ("", test_read.get_metadata()) << "METADATA MISMATCH" << std::endl;
}

TEST(FQ_READ, CONSTRUCTOR_2) {
    fq_read test_read("T", 18, "LJKHJGUYFDTGCHGJHL", "ASJHQUAL QUAL LLLL", "META DATA...");
    EXPECT_EQ("T", test_read.get_id()) << "ID MISMATCH" << std::endl;
    EXPECT_EQ(18, test_read.size()) << "READ SIZE MISMATCH" << std::endl;
    EXPECT_EQ("LJKHJGUYFDTGCHGJHL", test_read.get_seq()) << "READ SEQ MISMATCH" << std::endl;
    EXPECT_EQ("ASJHQUAL QUAL LLLL", test_read.get_quality()) << "READ QUALITY MISMATCH" << std::endl;
    EXPECT_EQ("META DATA...", test_read.get_metadata()) << "METADATA MISMATCH" << std::endl;
}

TEST(FQ_READ, OPERATOR_IN_BOUND) {
    std::string seq = "JKIYUTYFTDTRYUGHJKJUIHGUFHGTYTUYHIJHKJ";
    std::string qual = "MKHJGHFYTGHJKLKJKHJGHFGCVGJHKJLKMNBVBN";
    fq_read test_read("TEST", 38, seq, qual);
    for(int i = 0; i < 38; ++i) {
        EXPECT_NO_THROW(test_read[i]) << "VALID INDEX ACCESS THROWS ERROR" << std::endl;
        uint16_t exp_val = qual[i];
        exp_val <<= 8;
        exp_val |= seq[i];
        EXPECT_EQ(test_read[i], exp_val) << "READ[" << i << "] DOESN'T MATCH TEMPLATE[" << i << "]" << std::endl;
    }
}

TEST(FQ_READ, OPERATOR_OOB) {
    fq_read test_read("TEST", 20, "MUYTRFDSDERTGYHJKKJH", "KJHGFDEDSDFGHYUIOKJM");
    EXPECT_THROW(test_read[-1], std::out_of_range) << "EXPECTED OOB EXCEPTION USING -1" << std::endl;
    EXPECT_THROW(test_read[20], std::out_of_range) << "EXPECTED OOB EXCEPTION USING 20" << std::endl;
}

/* TO DO */
TEST(FQ_READ, RANDOM_STRESS) {

}