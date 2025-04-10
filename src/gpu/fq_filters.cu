#include "headers/fq_filter.cuh"

__global__ void cu_filter_reads(char* ALL_SEQ, uint32_t* OFFSETS, size_t LEN, char FILTER_MODE, char THRESH, uint64_t* FILTER_MASK, double PROPORTION) {
    // Block/thread OOB checks
    int SEQ_NUM = blockIdx.x * blockDim.x + threadIdx.x;
    if(SEQ_NUM >= LEN) {
        return;
    }

    // Seq variables
    uint32_t SEQ_BEGIN = (SEQ_NUM ? OFFSETS[SEQ_NUM - 1] + 1 : 0);
    uint32_t SEQ_SIZE = OFFSETS[SEQ_NUM] - SEQ_BEGIN + 1;

    // Check seq
    switch(FILTER_MODE) {
        case AVERAGE_DISCARD_WHOLE: {
            int sum = 0;
            for(int i = 0; i < SEQ_SIZE; ++i) {
                sum += ALL_SEQ[SEQ_BEGIN + i];
            }

            // Do not discard if average >= THRESH
            if(sum >= (SEQ_SIZE * THRESH)) {
                atomicOr((ULL*)&FILTER_MASK[SEQ_NUM / 64], (1ULL << (SEQ_NUM % 64)));
            }
            break;
        }

        case SINGLE_DISCARD_WHOLE: {
            for(int i = 0; i < SEQ_SIZE; ++i) {
                if(ALL_SEQ[SEQ_BEGIN + i] < THRESH) {
                    // Discard
                    return;
                }
            }
            
            // Do not discard
            atomicOr((ULL*)&FILTER_MASK[SEQ_NUM / 64], (1ULL << (SEQ_NUM % 64)));
            break;
        }

        // (% of quality < THRESH) > PERC = discard
        case PROPORTION_DISCARD_WHOLE: {
            // Count under thresh
            uint64_t count = 0;
            for(int i = 0; i < SEQ_SIZE; ++i) {
                if(ALL_SEQ[SEQ_BEGIN + i] < THRESH) {
                    ++count;
                }
            }

            // Discard if proportion under thresh > PROPORTION
            if((static_cast<double>(count) / SEQ_SIZE) <= PROPORTION) {
                atomicOr((ULL*)&FILTER_MASK[SEQ_NUM / 64], (1ULL << (SEQ_NUM % 64)));
            }
            break;
        }
    }
}

__global__ void cu_filter_reads_sw(char* ALL_SEQ, uint32_t* OFFSETS, size_t LEN, size_t K, char THRESH, double PROPORTION) {
    // Block/thread OOB checks
    int SEQ_NUM = blockIdx.x * blockDim.x + threadIdx.x;
    if(SEQ_NUM >= LEN) {
        return;
    }

    // Seq variables
    uint32_t SEQ_BEGIN = (SEQ_NUM ? OFFSETS[SEQ_NUM - 1] + 1 : 0);
    uint32_t SEQ_SIZE = OFFSETS[SEQ_NUM] - SEQ_BEGIN + 1;
    if(K > SEQ_SIZE) {
        return;
    }

    // Check and modify seq
    uint64_t sum = 0;
    int trim_size = SEQ_SIZE;
    
    // Set up first window
    for(int i = 0; i < K; ++i) {
        sum += ALL_SEQ[SEQ_BEGIN + i];
    }

    // Validate first window, iff valid, search for largest trim size
    if(sum < (THRESH * K)) {
        trim_size = 0;
    } else {
        for(int i = K; i < SEQ_SIZE; ++i) {
            sum -= ALL_SEQ[SEQ_BEGIN + i - K];
            sum += ALL_SEQ[SEQ_BEGIN + i];
            if(sum < (THRESH * K)) {
                trim_size = i;
                break;
            }
        }
    }

    // Trim read
    if(trim_size < SEQ_SIZE) {
        ALL_SEQ[SEQ_BEGIN + trim_size] = '\0';
    }
}