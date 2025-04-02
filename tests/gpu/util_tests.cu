#include "headers/test_prototypes.cuh"

void CU_UTIL_TESTS(std::vector<TEST_RESULT*>& RESULTS) {
    srand(0x56709AAE);

    // Simple __max test
    RESULTS.push_back(TEST("CU_UTIL", "MAX_A<B", [](){
        // Device vars
        int* d_ret;
        CUDA_CHECK(cudaMalloc(&d_ret, sizeof(int)));
        
        // Run kernel
        cu_max <<<1, 1>>> (1, 2, d_ret);
        CUDA_CHECK(cudaDeviceSynchronize());
        
        // Check results
        int h_ret;
        CUDA_CHECK(cudaMemcpy(&h_ret, d_ret, sizeof(int), cudaMemcpyDeviceToHost));
        EXPECT_EQ(h_ret, 2);
        
        // Clean up
        CUDA_CHECK(cudaFree(d_ret));
    }));

    // Simple __max test
    RESULTS.push_back(TEST("CU_UTIL", "MAX_A>B", [](){
        // Device vars
        int* d_ret;
        CUDA_CHECK(cudaMalloc(&d_ret, sizeof(int)));
        
        // Run kernel
        cu_max <<<1, 1>>> (420, 69, d_ret);
        CUDA_CHECK(cudaDeviceSynchronize());
        
        // Check results
        int h_ret;
        CUDA_CHECK(cudaMemcpy(&h_ret, d_ret, sizeof(int), cudaMemcpyDeviceToHost));
        EXPECT_EQ(h_ret, 420);
        
        // Clean up
        CUDA_CHECK(cudaFree(d_ret));
    }));

    // Test output for __base_to_bit given all chars
    RESULTS.push_back(TEST("CU_UTIL", "BASE_TO_BIT_ALL", [](){
        // Device vars
        char* d_ret;
        CUDA_CHECK(cudaMalloc(&d_ret, sizeof(char) * 255));
        CUDA_CHECK(cudaMemset(d_ret, 0, sizeof(char) * 255));
        
        // Run kernel
        cu_base2bit <<<16, 16>>> (d_ret);
        CUDA_CHECK(cudaDeviceSynchronize());

        // Check results
        std::vector<char> h_ret(255);
        CUDA_CHECK(cudaMemcpy(h_ret.data(), d_ret, sizeof(char) * 255, cudaMemcpyDeviceToHost));
        for(unsigned char i = 0; i < 255; ++i) {
            switch(i) {
                case 'A':
                    EXPECT_EQ(h_ret[i], static_cast<char>(0b00));
                    break;
                case 'C':
                    EXPECT_EQ(h_ret[i], static_cast<char>(0b01));
                    break;
                case 'G':
                    EXPECT_EQ(h_ret[i], static_cast<char>(0b10));
                    break;
                case 'T':
                    EXPECT_EQ(h_ret[i], static_cast<char>(0b11));
                    break;
                default:
                    EXPECT_EQ(static_cast<unsigned char>(h_ret[i]), static_cast<unsigned char>(0x80));
            }
        }
        
        // Clean up
        CUDA_CHECK(cudaFree(d_ret));
    }));

}