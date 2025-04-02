#include "headers/test_prototypes.cuh"

int main(int argc, char* argv[]) {
    std::vector<TEST_RESULT*> results;

    // Run all tests
    CU_UTIL_TESTS(results);
    CU_UTIL_STRUCTS_TESTS(results);
    CU_WRAPPER_TESTS(results);

    // Sort results
    std::sort(results.begin(), results.end(), [&](TEST_RESULT* a, TEST_RESULT* b) {
        if(a->TEST_SUITE == b->TEST_SUITE) {
            return a->TEST_NAME < b->TEST_NAME;
        }
        return a->TEST_SUITE < b->TEST_SUITE;
    });

    // Output header
    std::cout << GREEN << "[==========] Running " << results.size() << " tests." << RESET << std::endl;
    std::string prev = "";
    int passed = 0;
    int failed = 0;
    for(auto& result : results) {
        if(result->TEST_SUITE != prev) {
            std::cout << GREEN << "[----------] " << "Test suite: " << result->TEST_SUITE << RESET << std::endl;
            prev = result->TEST_SUITE;
        }
        std::string status = result->TEST_PASSED ? "       OK" : "***FAILED";
        std::cout << (result->TEST_PASSED ? GREEN : RED)
                  << "[" << status << " ] "
                  << result->TEST_SUITE << "." << result->TEST_NAME
                  << " (" << result->TIME_TAKEN << " ms)"
                  << RESET << std::endl;
        if(result->TEST_PASSED) {
            ++passed;
        } else {
            ++failed;
        }
    }

    // Output summary
    std::cout << GREEN << "[==========] " << results.size() << " tests ran." << RESET << std::endl;
    std::cout << GREEN << "[  PASSED  ] " << passed << RESET << std::endl;
    if(failed) {
        std::cout << RED << "[  FAILED  ] " << failed << RESET << std::endl;
    }

    // Cleanup
    for(auto& result : results) {
        delete result;
    }
    return failed ? -1 : 0;
}