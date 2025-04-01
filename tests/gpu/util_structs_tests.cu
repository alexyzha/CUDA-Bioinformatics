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
        // Create test variables
        uint64_t* d_key;
        uint64_t* d_value;
        std::vector<uint64_t> h_key(1024, rand());
        auto validate_hash = [](int key, int hash) -> bool {
            key = (~key) + (key << 21);
            key = key ^ (key >> 24);
            key = (key + (key << 3)) + (key << 8);
            key = key ^ (key >> 14);
            key = (key + (key << 2)) + (key << 4);
            key = key ^ (key >> 28);
            key = key + (key << 31);
            return (key == hash);
        };

        // Allocate mem, set/copy mem
        CUDA_CHECK(cudaMalloc(&d_key, sizeof(uint64_t) * 1024));
        CUDA_CHECK(cudaMalloc(&d_value, sizeof(uint64_t) * 1024));
        CUDA_CHECK(cudaMemset(d_value, 0, sizeof(uint64_t) * 1024));
        CUDA_CHECK(cudaMemcpy(d_key, h_key.data(), sizeof(uint64_t) * 1024, cudaMemcpyHostToDevice));

        // Run kernel
        cu_hashing <<<32, 32>>> (d_value, d_key, 1024);
        CUDA_CHECK(cudaDeviceSynchronize());

        // Copy mem back & validate results
        std::vector<uint64_t> h_value(1024);
        CUDA_CHECK(cudaMemcpy(h_value, d_value, sizeof(uint64_t) * 1024, cudaMemcpyDeviceToHost));
        for(int i = 0; i < 1024; ++i) {
            EXPECT_TRUE(validate_hash(h_key[i], h_value[i]));
        }

        // Clean up
        CUDA_CHECK(cudaFree(d_value));
        CUDA_CHECK(cudaFree(d_key));
    }));

    /*
     *  Need __global__ wrappers for union find tests
     */

}