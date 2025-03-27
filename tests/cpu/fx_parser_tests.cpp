#include <gtest/gtest.h>
#include "../../src/cpu/headers/fx_parser.h"

TEST(FA_PARSER_CPU, SINGLE) {
    // Expected fa_read variables
    std::string exp_id = "SRR32254469.1";
    std::string exp_meta = "VH01851:45:AAG3H73M5:1:1101:21602:1000";
    std::string exp_seq = "GGTCCGTCCAGGCTGCCTGCAATGATGTCCCCTTTCTACAA";
    int exp_len = 41;
    std::vector<fa_read*> reads;

    // Validate file_io
    EXPECT_NO_THROW([&]() {
        reads = read_fasta("infiles/fa_parser_single_in.txt");
    }()) << RED << "INFILE CONFIG ERROR" << RESET << std::endl;
    EXPECT_EQ(reads.size(), 1) << RED << "EXPECTED 1 READ, GOT " << reads.size() << RESET << std::endl;

    // Validate file read contents
    EXPECT_EQ(exp_id, reads[0]->get_id()) << RED << "ID PARSING ERROR" << RESET << std::endl;
    EXPECT_EQ(exp_meta, reads[0]->get_metadata()) << RED << "METADATA PARSING ERROR" << RESET << std::endl;
    EXPECT_EQ(exp_seq, reads[0]->get_seq()) << RED << "SEQ PARSING ERROR" << RESET << std::endl;
    EXPECT_EQ(exp_len, reads[0]->size()) << RED << "LENGTH CALCULATION ERROR" << RESET << std::endl;
}

TEST(FA_PARSER_CPU, MULTIPLE) {
    // Expected fa_read variables
    std::vector<std::string> exp_id {"SRR32254469.1", "SRR32254469.1", "SRR32254469.2"};
    std::vector<std::string> exp_meta {
        "VH01851:45:AAG3H73M5:1:1101:21602:1000",
        "VH01851:45:AAG3H73M5:1:1101:21602:1000",
        "VH01851:45:AAG3H73M5:1:1101:31581:1019"
    };
    std::vector<std::string> exp_seq {
        "GGTCCGTCCAGGCTGCCTGCAATGATGTCCCCTTTCTACAAGGGACACTGTGAGCAGGAATGATGTCTGGGAGCAGAATTCTCGCAGTCTTTCCACTGAATTCTCTGTCCCTGTCCTTTGCCTGGCTCTCTGTTGCAGATGGAAGCTG",
        "GTCACAATTGTCTCCTTGATGCCATGAGCCCCTGCAGTGTCACAATGGCTCCTTGCTTCCATGAGGAACACAAAAGCTTTGCATCAACTCTGAGCAACAGCTGCTGAACACACCTGGGAGCATCTGCTTGCATCACACGAGGGCTTCAAG",
        "GNATCAGGCCCCACTGTGTCACAATGGCTCCTTGGCTTCATGTGCATCACCTGTGTCACCATGGTCTCCTTGATTCCATGAGGTTCCACTGTGTCACAATGGCTCCTTGGCTTCATGTGCATCACCTGTGTCACCATGGCCTCCTTGATT"
    };
    std::vector<int> exp_len {148, 150, 150};
    std::vector<fa_read*> reads;

    // Validate file_io
    EXPECT_NO_THROW([&](){
        reads = read_fasta("infiles/fa_parser_multi_in.txt");
    }()) << RED << "INFILE CONFIG ERROR" << RESET << std::endl;
    EXPECT_EQ(reads.size(), 3) << RED << "EXPECTED 3 READS, GOT " << reads.size() << RESET << std::endl;

    // Validate file read contents
    for(int i = 0; i < 3; ++i) {
        EXPECT_EQ(reads[i]->get_id(), exp_id[i]) << RED << "ID PARSING ERROR" << RESET << std::endl;
        EXPECT_EQ(reads[i]->get_metadata(), exp_meta[i]) << RED << "METADATA PARSING ERROR" << RESET << std::endl;
        EXPECT_EQ(reads[i]->get_seq(), exp_seq[i]) << RED << "SEQ PARSING ERROR" << RESET << std::endl;
        EXPECT_EQ(reads[i]->size(), exp_len[i]) << RED << "LENGTH CALCULATION ERROR" << RESET << std::endl;
    }
}

/*
 *  Header sometimes only contains:
 *  [ID]
 *  As opposed to:
 *  [ID] [METADATA] [LENGTH=X] or [ID] [LENGTH=X]
 */
TEST(FA_PARSER_CPU, ID_ONLY) {
    // Expected fa_read variables
    std::string exp_id = "SRR32254469.2";
    std::string exp_seq = "GTGACACAGCAAGACCTCGTGGAAGCAAGGAGCCTTTGTGATACTCCAGATGCGCATGGAAGCCAGGCTGCACTGAGACACAGCGGGGCCTTGTTGGAACCAAAGACACCTCTGCAAAAAGCAGGGTAGCTGACAGTACCAAGGAGTCCAT";
    int exp_len = 151;
    std::vector<fa_read*> reads;

    // Validate file_io
    EXPECT_NO_THROW([&]() {
        reads = read_fasta("infiles/fa_parser_only_id_in.txt");
    }()) << RED << "INFILE CONFIG ERROR" << RESET << std::endl;
    EXPECT_EQ(reads.size(), 1) << RED << "EXPECTED 1 READ, GOT " << reads.size() << RESET << std::endl;

    // Validate file read contents
    EXPECT_EQ(exp_id, reads[0]->get_id()) << RED << "ID PARSING ERROR" << RESET << std::endl;
    EXPECT_EQ("", reads[0]->get_metadata()) << RED << "METADATA PARSING ERROR, EXPECTED NONE" << RESET << std::endl;
    EXPECT_EQ(exp_seq, reads[0]->get_seq()) << RED << "SEQ PARSING ERROR" << RESET << std::endl;
    EXPECT_EQ(exp_len, reads[0]->size()) << RED << "LENGTH CALCULATION ERROR" << RESET << std::endl;
}

TEST(FA_PARSER_CPU, EMPTY) {
    // Expect reads to be empty
    std::vector<fa_read*> reads;
    
    // Validate file_io
    EXPECT_NO_THROW([&](){
        reads = read_fasta("infiles/empty.txt");
    }()) << RED << "INFILE CONFIG ERROR" << RESET << std::endl;

    // Validate output
    EXPECT_TRUE(reads.empty()) << RED << "EXPECTED READS TO BE EMPTY" << RESET << std::endl;
}

TEST(FQ_PARSER_CPU, SINGLE) {
    // Expected fq_read variables
    std::string exp_id = "SRR32254469.1";
    std::string exp_meta = "VH01851:45:AAG3H73M5:1:1101:21602:1000";
    std::string exp_seq = "GTCACAATTGTCTCCTTGATGCCATGAGCCCCTGCAGTGTCACAATGGCTCCTTGCTTCCATGAGGAACACAAAAGCTTTGCATCAACTCTGAGCAACAGCTGCTGAACACACCTGGGAGCATCTGCTTGCATCACACGAGGGCTTCAAG";
    std::string exp_qual = "??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????";
    int exp_len = 150;
    std::vector<fq_read*> reads;

    // Validate file_io
    EXPECT_NO_THROW([&]() {
        reads = read_fastq("infiles/fq_parser_single_in.txt");
    }()) << RED << "INFILE CONFIG ERROR" << RESET << std::endl;
    EXPECT_EQ(reads.size(), 1) << RED << "EXPECTED 1 READ, GOT " << reads.size() << RESET << std::endl;

    // Validate file read contents
    EXPECT_EQ(exp_id, reads[0]->get_id()) << RED << "ID PARSING ERROR" << RESET << std::endl;
    EXPECT_EQ(exp_meta, reads[0]->get_metadata()) << RED << "METADATA PARSING ERROR" << RESET << std::endl;
    EXPECT_EQ(exp_seq, reads[0]->get_seq()) << RED << "SEQ PARSING ERROR" << RESET << std::endl;
    EXPECT_EQ(exp_len, reads[0]->size()) << RED << "LENGTH CALCULATION ERROR" << RESET << std::endl;
    EXPECT_EQ(exp_qual, reads[0]->get_quality()) << "QUALITY PARSING ERROR" << RESET << std::endl;
}

TEST(FQ_PARSER_CPU, MULTIPLE) {
    // Expected fq_read variables
    std::vector<std::string> exp_id {"SRR32254469.2", "SRR32254469.2", "SRR32254469.3"};
    std::vector<std::string> exp_meta {
        "VH01851:45:AAG3H73M5:1:1101:31581:1019",
        "VH01851:45:AAG3H73M5:1:1101:31581:1019",
        "VH01851:45:AAG3H73M5:1:1101:26298:1038"
    };
    std::vector<std::string> exp_seq {
        "GNATCAGGCCCCACTGTGTCACAATGGCTCCTTGGCTTCATGTGCATCACCTGTGTCACCATGGTCTCCTTGATTCCATGAGGTTCCACTGTGTCACAATGGCTCCTTGGCTTCATGTGCATCACCTGTGTCACCATGGCCTCCTTGATT",
        "GTGACACAGCAAGACCTCGTGGAAGCAAGGAGCCTTTGTGATACTCCAGATGCGCATGGAAGCCAGGCTGCACTGAGACACAGCGGGGCCTTGTTGGAACCAAAGACACCTCTGCAAAAAGCAGGGTAGCTGACAGTACCAAGGAGTCCAT",
        "GTGCAAGAGATGATATTGTAAAAAAGATCCACTGGGTTTATTTGGCAGGGATGGATGGTTCAGGACTATGCCCAGAGCACTTTTCCATTTATCAGGCAAAGAGACCTCCCACTGATGGGTATAGTTTTGGGGACTGGTTGCCTAATAGAAA"
    };
    std::vector<std::string> exp_qual {
        "??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????",
        "???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????",
        "???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????"
    };
    std::vector<int> exp_len {150, 151, 151};
    std::vector<fq_read*> reads;

    // Validate file_io
    EXPECT_NO_THROW([&](){
        reads = read_fastq("infiles/fq_parser_multi_in.txt");
    }()) << RED << "INFILE CONFIG ERROR" << RESET << std::endl;
    EXPECT_EQ(reads.size(), 3) << RED << "EXPECTED 3 READS, GOT " << reads.size() << RESET << std::endl;

    // Validate file read contents
    for(int i = 0; i < 3; ++i) {
        EXPECT_EQ(reads[i]->get_id(), exp_id[i]) << RED << "ID PARSING ERROR" << RESET << std::endl;
        EXPECT_EQ(reads[i]->get_metadata(), exp_meta[i]) << RED << "METADATA PARSING ERROR" << RESET << std::endl;
        EXPECT_EQ(reads[i]->get_seq(), exp_seq[i]) << RED << "SEQ PARSING ERROR" << RESET << std::endl;
        EXPECT_EQ(reads[i]->size(), exp_len[i]) << RED << "LENGTH CALCULATION ERROR" << RESET << std::endl;
        EXPECT_EQ(reads[i]->get_quality(), exp_qual[i]) << RED << "QUALITY PARSING ERROR" << RESET << std::endl;
    }
}

/*
 *  Header sometimes only contains:
 *  [ID]
 *  As opposed to:
 *  [ID] [METADATA] [LENGTH=X] or [ID] [LENGTH=X]
 */
TEST(FQ_PARSER_CPU, ID_ONLY) {
    // Expected fq_read variables
    std::string exp_id = "SRR32254469.1";
    std::string exp_seq = "GGTCCGTCCAGGCTGCCTGCAATGATGTCCCCTTTCTACAAGGGACACTGTGAGCAGGAATGATGTCTGGGAGCAGAATTCTCGCAGTCTTTCCACTGAATTCTCTGTCCCTGTCCTTTGCCTGGCTCTCTGTTGCAGATGGAAGCTG";
    std::string exp_qual = "????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????";
    int exp_len = 148;
    std::vector<fq_read*> reads;

    // Validate file_io
    EXPECT_NO_THROW([&]() {
        reads = read_fastq("infiles/fq_parser_only_id_in.txt");
    }()) << RED << "INFILE CONFIG ERROR" << RESET << std::endl;
    EXPECT_EQ(reads.size(), 1) << RED << "EXPECTED 1 READ, GOT " << reads.size() << RESET << std::endl;

    // Validate file read contents
    EXPECT_EQ(exp_id, reads[0]->get_id()) << RED << "ID PARSING ERROR" << RESET << std::endl;
    EXPECT_EQ("", reads[0]->get_metadata()) << RED << "METADATA PARSING ERROR, EXPECTED NONE" << RESET << std::endl;
    EXPECT_EQ(exp_seq, reads[0]->get_seq()) << RED << "SEQ PARSING ERROR" << RESET << std::endl;
    EXPECT_EQ(exp_len, reads[0]->size()) << RED << "LENGTH CALCULATION ERROR" << RESET << std::endl;
    EXPECT_EQ(exp_qual, reads[0]->get_quality()) << RED << "QUALITY PARSING ERROR" << RESET << std::endl;
}

TEST(FQ_PARSER_CPU, EMPTY) {
    // Expect reads to be empty
    std::vector<fq_read*> reads;
    
    // Validate file_io
    EXPECT_NO_THROW([&](){
        reads = read_fastq("infiles/empty.txt");
    }()) << RED << "INFILE CONFIG ERROR" << RESET << std::endl;

    // Validate output
    EXPECT_TRUE(reads.empty()) << RED << "EXPECTED READS TO BE EMPTY" << RESET << std::endl;
}