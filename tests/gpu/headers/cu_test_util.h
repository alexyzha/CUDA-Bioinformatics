#pragma once
#include <iostream>
#include <string>
#include <chrono>
#include <stdexcept>

#ifndef ANSI_ESC_COMMON
#define ANSI_ESC_COMMON
#define BLACK "\x1B[30m"
#define LBLACK "\x1B[90m"
#define RED "\x1B[31m"
#define LRED "\x1B[91m"
#define GREEN "\x1B[32m"
#define LGREEN "\x1B[92m"
#define YELLOW "\x1B[33m"
#define LYELLOW "\x1B[93m"
#define BLUE "\x1B[34m"
#define LBLUE "\x1B[94m"
#define MAGENTA "\x1B[35m"
#define LMAGENTA "\x1B[95m"
#define CYAN "\x1B[36m"
#define LCYAN "\x1B[96m"
#define WHITE "\x1B[37m"
#define LWHITE "\x1B[97m"
#define RESET "\x1B[0m"
#endif

#ifndef C_CU_TEST_RESULT
#define C_CU_TEST_RESULT

struct TEST_RESULT {
    std::string TEST_SUITE;
    std::string TEST_NAME;
    bool TEST_PASSED;
    int TIME_TAKEN;             // ms
};

#endif

template<typename T>
TEST_RESULT* TEST(std::string TEST_SUITE, std::string TEST_NAME, T TEST_FXN) {
    // Take start time
    auto START = std::chrono::high_resolution_clock::now();

    try {
        // Expects throw on test fail
        TEST_FXN();
        
        // Take end time, return
        auto END = std::chrono::high_resolution_clock::now();
        int TIME_TAKEN = duration_cast<milliseconds>(END - START).count();
        return new TEST_RESULT{TEST_SUITE, TEST_NAME, true, TIME_TAKEN};
    } catch(const std::exception& e) {
        // Print error message
        std::cerr << RED 
                  << TEST_SUITE << '.' << TEST_NAME 
                  << "\n> TEST FAILED\n> " << e.what() 
                  << RESET << std::endl;

        // Take end time, return
        auto END = std::chrono::high_resolution_clock::now();
        int TIME_TAKEN = duration_cast<milliseconds>(END - START).count();
        return new TEST_RESULT{TEST_SUITE, TEST_NAME, false, TIME_TAKEN};
    } catch(...) {
        // Print error message
        std::cerr << RED 
                  << TEST_SUITE << '.' << TEST_NAME 
                  << "\n> TEST FAILED\n> UNKNOWN EXCEPTION"
                  << RESET << std::endl;

        // Take end time, return
        auto END = std::chrono::high_resolution_clock::now();
        int TIME_TAKEN = duration_cast<milliseconds>(END - START).count();
        return new TEST_RESULT{TEST_SUITE, TEST_NAME, false, TIME_TAKEN};
    }
}