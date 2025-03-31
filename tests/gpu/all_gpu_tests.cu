#include "headers/test_prototypes.cuh"

int main(int argc, char* argv[]) {
    std::vector<TEST_RESULT*> results;

    // Run all tests
    CU_UTIL_TESTS(results);


    // Sort results
    std::sort(results.begin(), results.end(), [&](TEST_RESULT* a, TEST_RESULT* b) {
        if(a->TEST_SUITE == b->TEST_SUITE) {
            return a->TEST_NAME < b->TEST_NAME;
        }
        return a->TEST_SUITE < b->TEST_SUITE;
    });

    // Output to XML


    return 0;
}