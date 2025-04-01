#include "headers/wrappers.cuh"

std::vector<fq_read*> cu_filter_fq(std::vector<fq_read*> READS, char FILTER_MODE, char THRESH, size_t K, double PROPORTION) {
    // Size variables
    int READ_LEN = READS.size();
    int FILTER_LENGTH = (READ_LEN / 64) + 1;
    
    // Host variables
    std::string temp = "";
    std::vector<uint32_t> h_offsets(READ_LEN);
    for(int i = 0; i < READ_LEN; ++i) {
        temp += READS[i]->get_quality();
        h_offsets[i] = temp.size() - 1;
    }
    int SEQ_LEN = temp.size() + 1;                  // temp + '\0'
    const char* h_allseq = temp.c_str();
    
    // Device variables
    char* d_allseq;
    uint32_t* d_offsets;    
    uint64_t* d_filter;

    // Allocate mem for device variables
    CUDA_CHECK(cudaMalloc(&d_allseq, sizeof(char) * SEQ_LEN));
    CUDA_CHECK(cudaMalloc(&d_offsets, sizeof(uint32_t) * READ_LEN));
    CUDA_CHECK(cudaMalloc(&d_filter, sizeof(uint64_t) * FILTER_LENGTH));
    CUDA_CHECK(cudaMemset(d_filter, 0, sizeof(uint64_t) * FILTER_LENGTH));

    // Copy mem host -> device
    CUDA_CHECK(cudaMemcpy(d_allseq, h_allseq, sizeof(char) * SEQ_LEN, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_offsets, h_offsets.data(), sizeof(uint32_t) * READ_LEN, cudaMemcpyHostToDevice));

    // Run kernels
    int THREADS = MAX_THREADS;
    int BLOCKS = (READ_LEN / THREADS) + 1;
    if(FILTER_MODE == SLIDING_WINDOW) {
        cu_filter_reads_sw <<<BLOCKS, THREADS>>> (
            d_allseq, 
            d_offsets, 
            READ_LEN, 
            K, 
            THRESH, 
            PROPORTION
        );
    } else {
        cu_filter_reads <<<BLOCKS, THREADS>>> (
            d_allseq, 
            d_offsets, 
            READ_LEN, 
            FILTER_MODE, 
            THRESH, 
            d_filter, 
            PROPORTION
        );
    }
    
    // Sync
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy mem back
    std::vector<fq_read*> filtered_reads;
    if(FILTER_MODE != SLIDING_WINDOW) {
        // Copy mem back
        std::vector<uint64_t> ret_filter(FILTER_LENGTH, 0);
        CUDA_CHECK(cudaMemcpy(ret_filter.data(), d_filter, sizeof(uint64_t) * FILTER_LENGTH, cudaMemcpyDeviceToHost));

        // Filter and return
        for(int i = 0; i < READ_LEN; ++i) {
            int f_index = i / 64;
            int f_bit = i % 64;
            if(ret_filter[f_index] & (1ULL << f_bit)) {
                filtered_reads.push_back(READS[i]);
            }
        }
    } else {
        // Copy mem back
        std::vector<char> temp_allseq(SEQ_LEN);
        CUDA_CHECK(cudaMemcpy(temp_allseq.data(), d_allseq, sizeof(char) * SEQ_LEN, cudaMemcpyDeviceToHost));
        const char* ret_allseq = temp_allseq.data();

        // Filter, trim, and return
        int seq_begin = 0;
        for(int i = 0; i < READ_LEN; ++i) {
            int og_size = READS[i]->size();
            
            // Get seq substring from seq buffer
            std::string cur_seq(ret_allseq + seq_begin);
            if(cur_seq.size() >= og_size) {
                filtered_reads.push_back(READS[i]);
            } else if(cur_seq.size()) {
                filtered_reads.push_back(new fq_read(
                    READS[i]->get_id(),
                    cur_seq.size(),
                    READS[i]->get_seq().substr(0, cur_seq.size()),
                    cur_seq,
                    READS[i]->get_metadata()
                ));
            }

            // Move to next seq begin
            seq_begin = h_offsets[i] + 1;
        }
    }
    
    // Cleanup and return
    CUDA_CHECK(cudaFree(d_allseq));
    CUDA_CHECK(cudaFree(d_offsets));
    CUDA_CHECK(cudaFree(d_filter));
    return filtered_reads;
}