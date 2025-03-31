#include "headers/test_prototypes.cuh"

void CU_UTIL_TESTS(std::vector<TEST_RESULT*>& RESULTS) {

    // Simple __max test
    RESULTS.push_back(TEST("CU_UTIL", "MAX_A<B", [](){
        
        // __global

    }));

    // Simple __max test
    RESULTS.push_back(TEST("CU_UTIL", "MAX_A>B", [](){
        
        // __global

    }));

    // Test output for __base_to_bit given ACGT
    RESULTS.push_back(TEST("CU_UTIL", "BASE_TO_BIT_ACGT", [](){
        
        // __global

    }));
    
    // Test output for __base_to_bit given all chars aside from ACGT
    RESULTS.push_back(TEST("CU_UTIL", "BASE_TO_BIT_NON_ACGT", [](){
        for(char i = 0; char <= static_cast<char>(0xff); ++i) {
            if(i == 'A' || i == 'C' || i == 'G' || i == 'T') {
                continue;
            }
            
            // __global
            
        }
    }));

}