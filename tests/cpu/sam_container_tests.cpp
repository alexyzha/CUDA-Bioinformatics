#include <gtest/gtest.h>
#include "../../src/cpu/headers/sam_container.h"

TEST(SAM_CONTAINER, CONSTRUCTOR_EMPTY) {
    sam_container cont("infiles/empty.txt");
    const auto& headers = cont.get_headers();
    const auto& reads = cont.get_reads();
    EXPECT_TRUE(headers.empty());
    EXPECT_TRUE(reads.empty());
}

TEST(SAM_CONTAINER, CONSTRUCTOR_READ_ONLY) {

}

TEST(SAM_CONTAINER, CONSTRUCTOR_HEADER_ONLY) {

}