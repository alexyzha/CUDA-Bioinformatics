#include <gtest/gtest.h>
#include "../../src/cpu/headers/xam_util.h"
#include "../../src/cpu/headers/xam_parser.h"

TEST(XAM_UTIL, SAM_TO_FILE_SINGLE) {
    std::ofstream file("outfiles/sam_to_file_out.txt");
    EXPECT_TRUE(file.is_open()) << RED << "UNABLE TO OPEN SAM2F_OUT" << RESET << std::endl;

}