#include <gtest/gtest.h>
#include "../../src/cpu/fq_read.h"

TEST(FQ_READ, CONSTRUCTOR_1) {
    // Create read (default metadata = "")
    fq_read test_read("TESTY CAL", 35, "KLJKHJGHFTUYHJKUIYHUGHFGHJHKJHGHGCF", "QQQLKAJBHHUJIKSLAJNBHKJGUIHJNLASJKH");
    
    // Validate read constructor (no metadata)
    EXPECT_EQ("TESTY CAL", test_read.get_id()) << RED << "ID MISMATCH" << RESET << std::endl;
    EXPECT_EQ(35, test_read.size()) << RED << "READ SIZE MISMATCH" << RESET << std::endl;
    EXPECT_EQ("KLJKHJGHFTUYHJKUIYHUGHFGHJHKJHGHGCF", test_read.get_seq()) << RED << "READ SEQ MISMATCH" << RESET << std::endl;
    EXPECT_EQ("QQQLKAJBHHUJIKSLAJNBHKJGUIHJNLASJKH", test_read.get_quality()) << RED << "READ QUALITY MISMATCH" << RESET << std::endl;
    EXPECT_EQ("", test_read.get_metadata()) << RED << "METADATA MISMATCH" << RESET << std::endl;
}

TEST(FQ_READ, CONSTRUCTOR_2) {
    // Create read
    fq_read test_read("T", 18, "LJKHJGUYFDTGCHGJHL", "ASJHQUAL QUAL LLLL", "META DATA...");
    
    // Validate read constructor
    EXPECT_EQ("T", test_read.get_id()) << RED << "ID MISMATCH" << RESET << std::endl;
    EXPECT_EQ(18, test_read.size()) << RED << "READ SIZE MISMATCH" << RESET << std::endl;
    EXPECT_EQ("LJKHJGUYFDTGCHGJHL", test_read.get_seq()) << RED << "READ SEQ MISMATCH" << RESET << std::endl;
    EXPECT_EQ("ASJHQUAL QUAL LLLL", test_read.get_quality()) << RED << "READ QUALITY MISMATCH" << RESET << std::endl;
    EXPECT_EQ("META DATA...", test_read.get_metadata()) << RED << "METADATA MISMATCH" << RESET << std::endl;
}

TEST(FQ_READ, OPERATOR_IN_BOUND) {
    // Create read
    std::string seq = "JKIYUTYFTDTRYUGHJKJUIHGUFHGTYTUYHIJHKJ";
    std::string qual = "MKHJGHFYTGHJKLKJKHJGHFGCVGJHKJLKMNBVBN";
    fq_read test_read("TEST", 38, seq, qual);
    
    // Test [] operator with 16 bit bitmask that contains both [seq/qual] at index i
    for(int i = 0; i < 38; ++i) {
        EXPECT_NO_THROW(test_read[i]) << RED << "VALID INDEX ACCESS THROWS ERROR" << RESET << std::endl;
        uint16_t exp_val = qual[i];
        exp_val <<= 8;
        exp_val |= seq[i];
        EXPECT_EQ(test_read[i], exp_val) << RED << "READ[" << i << "] DOESN'T MATCH TEMPLATE[" << i << "]" << RESET << std::endl;
    }
}

TEST(FQ_READ, OPERATOR_OOB) {
    // Create read
    fq_read test_read("TEST", 20, "MUYTRFDSDERTGYHJKKJH", "KJHGFDEDSDFGHYUIOKJM");
    
    // Test out of bound index access
    EXPECT_THROW(test_read[-1], std::out_of_range) << RED << "EXPECTED OOB EXCEPTION USING -1" << RESET << std::endl;
    EXPECT_THROW(test_read[20], std::out_of_range) << RED << "EXPECTED OOB EXCEPTION USING 20" << RESET << std::endl;
}

/*
 *  to_file takes ostream as input, file_io burden put on programmer
 */
TEST(FQ_READ, FILE_OUT) {
    // File actual content vs expected content comparison
    std::vector<std::string> file_actual = {};
    std::vector<std::string> file_exp {
        "@SRR32254469.1 VH01851:45:AAG3H73M5:1:1101:21602:1000 length=26",
        "GGTCCGTCCAGGCTGCCTGCAATGAT",
        "+SRR32254469.1 VH01851:45:AAG3H73M5:1:1101:21602:1000 length=26",
        "??????????????????????????"
    };
    std::string temp = "";

    // Create read
    fq_read test_read("SRR32254469.1", 26, "GGTCCGTCCAGGCTGCCTGCAATGAT", "??????????????????????????", "VH01851:45:AAG3H73M5:1:1101:21602:1000");
    
    // fstream handle + streaming read to file
    std::ofstream outfile("outfiles/fq_read_out.txt");
    EXPECT_TRUE(outfile.is_open()) << RED << "OUTFILE NOT CONFIGURED PROPERLY" << RESET << std::endl;
    test_read.to_file(outfile);
    outfile.close();

    // fstream handle + checking read output to file
    std::ifstream infile("outfiles/fq_read_out.txt");
    EXPECT_TRUE(infile.is_open()) << RED << "UNABLE TO OPEN INFILE" << RESET << std::endl;
    while(std::getline(infile, temp)) {
        file_actual.push_back(temp);
    }
    EXPECT_EQ(file_actual.size(), file_exp.size()) << RED << "DIFFERENT NUMBER OF OUTPUT LINES THAN EXPECTED" << RESET << std::endl;
    for(int i = 0; i < file_exp.size(); ++i) {
        EXPECT_EQ(file_actual[i], file_exp[i]) << RED << "FILE OUTPUT & EXPECTED OUTPUT MISMATCH AT LINE [" << i << "]" << RESET << std::endl;
    }
    infile.close();
}

TEST(FQ_READ, RANDOM_STRESS) {
    // Generating sequence + quality scores
    std::string seq = "";
    std::string qual = "";
    int len = rand() % 100000;
    for(int i = 0; i < len; ++i) {
        seq.push_back(static_cast<char>(rand() % 0xff));
        qual.push_back(static_cast<char>(rand() % 0xff));
    }

    // Creating read with generated sequence + quality scores
    fq_read test_read("STRESS TEST", len, seq, qual);

    // Testing [] operator against seq/qual bitmask as the template
    for(int i = 0; i < len; ++i) {
        EXPECT_NO_THROW(test_read[i]) << RED << "VALID INDEX ACCESS THROWS ERROR" << RESET << std::endl;
        uint16_t exp_val = qual[i];
        exp_val <<= 8;
        exp_val |= seq[i];
        EXPECT_EQ(test_read[i], exp_val) << RED << "READ[" << i << "] DOESN'T MATCH TEMPLATE[" << i << "]" << RESET << std::endl;
    }
    EXPECT_THROW(test_read[-1], std::out_of_range) << RED << "EXPECTED OOB EXCEPTION USING -1" << RESET << std::endl;
    EXPECT_THROW(test_read[len], std::out_of_range) << RED << "EXPECTED OOB EXCEPTION USING 17" << RESET << std::endl;
}