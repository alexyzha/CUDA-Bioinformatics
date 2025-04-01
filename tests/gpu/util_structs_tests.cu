#include "headers/test_prototypes.cuh"

void CU_UTIL_STRUCTS_TESTS(std::vector<TEST_RESULT*>& RESULTS) {

    // Kmer hashtable constuctor test basic
    RESULTS.push_back(new TEST("CU_UTIL_STRUCTS", "KMER_HASHTABLE_CONSTRUCTOR_BASIC", [](){
        // Create test variables
        kh_pair<uint64_t>* d_map = kh_construct<uint64_t>(100);
        kh_pair<uint64_t> h_map[100];
    
        // Copy device -> host
        CUDA_CHECK(cudaMemcpy(h_map, d_map, sizeof(h_map), cudaMemcpyDeviceToHost));
        for(int i = 0; i < 100; ++i) {
            EXPECT_EQ((uint64_t)0, h_map[i].key);
            EXPECT_EQ((uint64_t)0, h_map[i].value);
        }

        // Clean up
        CUDA_CHECK(cudaFree(d_map));
    }));
    
    // Kmer hashtable constuctor test array
    RESULTS.push_back(new TEST("CU_UTIL_STRUCTS", "KMER_HASHTABLE_CONSTRUCTOR_ARRAY", [](){
        // Create test variables
        kh_pair<uint64_t[100]>* d_map = kh_construct<uint64_t[100]>(10);
        kh_pair<uint64_t[100]> h_map[10];
    
        // Copy device -> host
        CUDA_CHECK(cudaMemcpy(h_map, d_map, sizeof(h_map), cudaMemcpyDeviceToHost));
        for(int i = 0; i < 10; ++i) {
            EXPECT_EQ((uint64_t)0, h_map[i].key);
            for(int j = 0; j < 100; ++j) {
                EXPECT_EQ((uint64_t)0, h_map[i].value[j]);
            }
        }

        // Clean up
        CUDA_CHECK(cudaFree(d_map));
    }));
    
    // Hashing test
    RESULTS.push_back(new TEST("CU_UTIL_STRUCTS", "KMER_HASHING", [](){
        
        // __global__ kernel done

    }));

    /*
     *  Need __global__ wrappers for union find tests
     */

}