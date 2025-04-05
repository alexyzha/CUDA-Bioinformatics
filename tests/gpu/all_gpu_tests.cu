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

    // Output to XML
    if(argc > 1 && !failed) {
        // File handling
        std::ofstream file(argv[1]);
        if(!file.is_open()) {
            std::cerr << RED << "Invalid file path [" << argv[1] 
                      << "].\nRun ./all_gpu_tests [file_path]." 
                      << RESET << std::endl;
        }
        
        // Groups
        int total = results.size();
        double total_time = 0.0;
        std::unordered_map<std::string, std::vector<TEST_RESULT*>> suites;
        for(auto& result : results) {
            suites[result->TEST_SUITE].push_back(result);
            total_time += (double)result->TIME_TAKEN / 1000.0;
        }
        std::string time = []() {
            std::time_t now = std::time(nullptr);
            std::tm tm_now;
            char buf[32];
            localtime_r(&now, &tm_now);
            std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm_now);
            std::stringstream ss;
            ss << buf << "." << std::setw(3) << std::setfill('0') << (now % 1000);
            return ss.str();
        }();

        // Stream to file
        file << std::fixed << std::setprecision(3);
        file << "<testsuites tests=\"" << total << "\" failures=\"" << 0 
             << "\" disabled=\"0\" errors=\"0\" time=\"" << total_time
             << "\" timestamp=\"" << time << "\" name=\"AllTests\">\n";
        for(auto& [suite, tests] : suites) {
            double suite_time = 0.0;
            for(auto& test : tests) {
                suite_time += (double)test->TIME_TAKEN / 1000.0;
            }
            file << "\t<testsuite name=\"" << suite << "\" tests=\"" << tests.size()
                 << "\" failures=\"" << 0 << "\" disabled=\"0\" skipped=\"0\" errors=\"0\" time = \""
                 << suite_time << "\" timestamp=\"" << time << "\">\n";

            // Stream all tests to file
            for(auto& test : tests) {
                file << "\t\t<testcase name=\"" << test->TEST_NAME
                    << "\" file=\"" << suite << "_tests.cpp\" line=\"0\" status=\"run\" result=\""
                    << "completed" << "\" time=\"" << (double)test->TIME_TAKEN / 1000.0
                    << "\" timestamp=\"" << time
                    << "\" classname=\"" << suite << "\"/>\n"; 
            }
            file << "\t</testsuite>\n";
        }
        file << "</testsuites>\n";
    }

    // Cleanup
    for(auto& result : results) {
        delete result;
    }
    return failed ? -1 : 0;
}