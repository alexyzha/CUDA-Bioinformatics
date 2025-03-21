#include <gtest/gtest.h>
#include "../../src/fa_read.h"

TEST(FA_READ, CONSTRUCTOR_1) {
    // Create read (default metadata = "")
    fa_read test_read("ABCDE", 10, "ACGTACGTAC");

    // Validate read constructor (no metadata)
    EXPECT_EQ("ABCDE", test_read.get_id()) << RED << "ID MISMATCH" << RESET << std::endl;
    EXPECT_EQ(10, test_read.size()) << RED << "READ SIZE MISMATCH" << RESET << std::endl;
    EXPECT_EQ("ACGTACGTAC", test_read.get_seq()) << RED << "READ SEQ MISMATCH" << RESET << std::endl;
    EXPECT_EQ("", test_read.get_metadata()) << RED << "METADATA MISMATCH" << RESET << std::endl;
}

TEST(FA_READ, CONSTRUCTOR_2) {
    // Create read
    fa_read test_read("ALHGKHIJKL", 21, "AOALJKBNQJPIHUOBJNKLM", "SUPER_META_DATA_AHH");

    // Validate read constructor
    EXPECT_EQ("ALHGKHIJKL", test_read.get_id()) << RED << "ID MISMATCH" << RESET << std::endl;
    EXPECT_EQ(21, test_read.size()) << RED << "READ SIZE MISMATCH" << RESET << std::endl;
    EXPECT_EQ("AOALJKBNQJPIHUOBJNKLM", test_read.get_seq()) << RED << "READ SEQ MISMATCH" << RESET << std::endl;
    EXPECT_EQ("SUPER_META_DATA_AHH", test_read.get_metadata()) << RED << "METADATA MISMATCH" << RESET << std::endl;
}

TEST(FA_READ, OPERATOR_IN_BOUND) {
    // Create read
    std::string seq = "CFHGVJHBIUOJIKLNMBGHFGDTAGCVHGJHKJLKMKNJGUYFGCVHBJ";

    // Validate [] operator with seq as the template
    fa_read test_read("TEST", 50, seq);
    for(int i = 0; i < 50; ++i) {
        EXPECT_NO_THROW(test_read[i]) << RED << "VALID INDEX ACCESS THROWS ERROR" << RESET << std::endl;
        EXPECT_EQ(test_read[i], seq[i]) << RED << "READ[" << i << "] DOESN'T MATCH TEMPLATE[" << i << "]" << RESET << std::endl;
    }
}

TEST(FA_READ, OPERATOR_OOB) {
    // Create read
    fa_read test_read("TEST", 17, "MNBHUYTGFGHJKLMJH");

    // Expect throw from out of bound accesses
    EXPECT_THROW(test_read[-1], std::out_of_range) << RED << "EXPECTED OOB EXCEPTION USING -1" << RESET << std::endl;
    EXPECT_THROW(test_read[17], std::out_of_range) << RED << "EXPECTED OOB EXCEPTION USING 17" << RESET << std::endl;
}

/*
 *  to_file takes ostream as input, file_io burden put on programmer
 */
TEST(FA_READ, FILE_OUT) {
    // File actual content vs expected content comparison
    std::vector<std::string> file_actual = {};
    std::vector<std::string> file_exp {
        ">SRR32254469.1 VH01851:45:AAG3H73M5:1:1101:21602:1000 length=26",
        "GGTCCGTCCAGGCTGCCTGCAATGAT"
    };
    std::string temp = "";

    // Create read
    fa_read test_read("SRR32254469.1", 26, "GGTCCGTCCAGGCTGCCTGCAATGAT", "VH01851:45:AAG3H73M5:1:1101:21602:1000");
    
    // fstream handle + streaming read to file
    std::ofstream outfile("outfiles/fa_read_out.txt");
    EXPECT_TRUE(outfile.is_open()) << RED << "OUTFILE NOT CONFIGURED PROPERLY" << RESET << std::endl;
    test_read.to_file(outfile);
    outfile.close();
    
    // fstream handle + checking read output to file
    std::ifstream infile("outfiles/fa_read_out.txt");
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

TEST(FA_READ, RANDOM_STRESS) {
    // Generating sequence
    std::string seq = "";
    int len = rand() % 100000;
    for(int i = 0; i < len; ++i) {
        seq.push_back(static_cast<char>(rand() % 0xff));
    }

    // Creating read with generated sequence
    fa_read test_read("STRESS TEST", len, seq);
    
    // Testing [] operator against seq as the template
    for(int i = 0; i < len; ++i) {
        EXPECT_NO_THROW(test_read[i]) << RED << "VALID INDEX ACCESS THROWS ERROR" << RESET << std::endl;
        EXPECT_EQ(test_read[i], seq[i]) << RED << "READ[" << i << "] DOESN'T MATCH TEMPLATE[" << i << "]" << RESET << std::endl;
    }
    EXPECT_THROW(test_read[-1], std::out_of_range) << RED << "EXPECTED OOB EXCEPTION USING -1" << RESET << std::endl;
    EXPECT_THROW(test_read[len], std::out_of_range) << RED << "EXPECTED OOB EXCEPTION USING 17" << RESET << std::endl;
}