#include <gtest/gtest.h>
#include "../../src/cpu/headers/xam_parser.h"

TEST(SAM_PARSER, EMPTY) {
    std::string file = "infiles/empty.txt";
    std::vector<sam_read*> reads = read_sam(file);
    std::unordered_map<std::string, std::vector<std::string>> headers = read_sam_headers(file);
    EXPECT_TRUE(reads.empty()) << RED << "READ FROM EMPTY FILE" << RESET << std::endl;
    EXPECT_TRUE(headers.empty()) << RED << "HEADER FROM EPMTY FILE" << RESET << std::endl;
}

TEST(SAM_PARSER, READS_ONLY_NO_TAGS) {

}

TEST(SAM_PARSER, READS_ONLY_YES_TAGS) {

}

TEST(SAM_PARSER, MULTI_READS) {

}

TEST(SAM_PARSER, HEADER_ONLY_NO_REPEAT) {

}

TEST(SAM_PARSER, HEADER_ONLY_YES_REPEAT) {

}

TEST(SAM_PARSER, BOTH_READS_AND_HEADER) {

}
