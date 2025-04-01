#include "headers/local_alignment.cuh"

__global__ void cu_local_alignment(char* ALL_SEQ, char* CIGAR_BUF, int* CACHE, uint32_t* OFFSETS, size_t LEN, int* RET) {
    // OOB check for block/thread
    int SEQ_NUM = blockIdx.x * blockDim.x + threadIdx.x;
    if(SEQ_NUM >= LEN || !SEQ_NUM) {
        return;
    }

    // Extract ref and read offsets for convenience
    size_t REF_BEGIN = 0;
    size_t REF_SIZE = OFFSETS[0] + 1;
    size_t READ_BEGIN = OFFSETS[SEQ_NUM - 1] + 1;
    size_t READ_SIZE = OFFSETS[SEQ_NUM] + 1 - READ_BEGIN;
    size_t CIGAR_BEGIN = SEQ_NUM * MAX_CIGAR_LEN;

    // Find preallocated space from CACHE
    int* cache = &CACHE[(MAX_REF_LEN + 1) * (MAX_READ_LEN + 1) * SEQ_NUM];
    int MAX_SCORE = 0;
    int MAX_I = 0;
    int MAX_J = 0;

    // Run DP
    for(int i = 1; i <= REF_SIZE; ++i) {
        for(int j = 1; j <= READ_SIZE; ++j) {
            int dgnl = (
                DP(i - 1, j - 1)
                + (ALL_SEQ[REF_BEGIN + i - 1] == ALL_SEQ[READ_BEGIN + j - 1] 
                    ? ALIGN_MATCH 
                    : ALIGN_MISMATCH)
            );
            int up = DP(i - 1, j) + ALIGN_GAP;
            int left = DP(i, j - 1) + ALIGN_GAP;
            DP(i, j) = __max(__max(0, dgnl), __max(up, left));
            if(DP(i, j) > MAX_SCORE) {
                MAX_SCORE = DP(i, j);
                MAX_I = i;
                MAX_J = j;
            }
        }
    }

    // Traceback cigar string from DP results
    char* cigar = &CIGAR_BUF[CIGAR_BEGIN];
    char prev = '\0';
    int cigar_iter = 0;
    int count = 0;
    int i = MAX_I;
    int j = MAX_J;

    // Main traceback loop
    while(i > 0 && j > 0 && DP(i, j) > 0) {
        int score = DP(i, j);
        int match = (
            ALL_SEQ[REF_BEGIN + i - 1] == ALL_SEQ[READ_BEGIN + j - 1]
                ? ALIGN_MATCH
                : ALIGN_MISMATCH
        );
        char cur;
        
        // Main logic
        if(score == DP(i - 1, j - 1) + match) {
            cur = (match == ALIGN_MATCH ? 'M' : 'X');
            --i;
            --j;
        } else if(score == DP(i - 1, j) + ALIGN_GAP) {
            cur = 'D';
            --i;
        } else if(score == DP(i, j - 1) + ALIGN_GAP) {
            cur = 'I';
            --j;
        } else {
            break;
        }

        // Add to cigar string
        if(cur == prev) {
            ++count;
        } else {
            if(count) {
                cigar_iter += sprintf(
                    cigar + cigar_iter, 
                    "%d%c",
                    count,
                    prev 
                );
            }
            count = 1;
            prev = cur;
        }
    }

    // Check for leftovers
    if(count) {
        cigar_iter += sprintf(
            cigar + cigar_iter, 
            "%d%c",
            count,
            prev 
        );
    }
    cigar[cigar_iter] = '\0';

    // Cigar string is built in reverse order
    for(int k = 0; k < cigar_iter / 2; ++k) {
        char temp = cigar[k];
        cigar[k] = cigar[cigar_iter - k - 1];
        cigar[cigar_iter - k - 1] = temp;
    }

    // "Return" in RET
    int RET_BEGIN = SEQ_NUM * 3;
    RET[RET_BEGIN] = MAX_I - 1;
    RET[RET_BEGIN + 1] = MAX_J - 1;
    RET[RET_BEGIN + 2] = MAX_SCORE;
}

__global__ void cu_fq_local_alignment(char* REF, cu_fq_read* READS, char* CIGAR_BUF, int* CACHE, uint32_t* OFFSETS, size_t LEN, size_t REF_SIZE, int* RET) {
    // OOB check for block/thread
    int SEQ_NUM = blockIdx.x * blockDim.x + threadIdx.x;
    if(SEQ_NUM >= LEN || !SEQ_NUM) {
        return;
    }

    // Extract ref and read offsets for convenience
    size_t READ_SIZE = READS[SEQ_NUM].size;
    size_t CIGAR_BEGIN = SEQ_NUM * MAX_CIGAR_LEN;

    // Find preallocated space from CACHE
    int* cache = &CACHE[(MAX_REF_LEN + 1) * (MAX_READ_LEN + 1) * SEQ_NUM];
    int MAX_SCORE = 0;
    int MAX_I = 0;
    int MAX_J = 0;

    // Run DP
    for(int i = 1; i <= REF_SIZE; ++i) {
        for(int j = 1; j <= READ_SIZE; ++j) {
            int dgnl = (
                DP(i - 1, j - 1)
                + (REF[i - 1] == READS[SEQ_NUM].seq[j - 1] 
                    ? ALIGN_MATCH 
                    : ALIGN_MISMATCH)
            );
            int up = DP(i - 1, j) + ALIGN_GAP;
            int left = DP(i, j - 1) + ALIGN_GAP;
            DP(i, j) = __max(__max(0, dgnl), __max(up, left));
            if(DP(i, j) > MAX_SCORE) {
                MAX_SCORE = DP(i, j);
                MAX_I = i;
                MAX_J = j;
            }
        }
    }

    // Traceback cigar string from DP results
    char* cigar = &CIGAR_BUF[CIGAR_BEGIN];
    char prev = '\0';
    int cigar_iter = 0;
    int count = 0;
    int i = MAX_I;
    int j = MAX_J;

    // Main traceback loop
    while(i > 0 && j > 0 && DP(i, j) > 0) {
        int score = DP(i, j);
        int match = (
            REF[i - 1] == READS[SEQ_NUM].seq[j - 1] 
                ? ALIGN_MATCH
                : ALIGN_MISMATCH
        );
        char cur;
        
        // Main logic
        if(score == DP(i - 1, j - 1) + match) {
            cur = (match == ALIGN_MATCH ? 'M' : 'X');
            --i;
            --j;
        } else if(score == DP(i - 1, j) + ALIGN_GAP) {
            cur = 'D';
            --i;
        } else if(score == DP(i, j - 1) + ALIGN_GAP) {
            cur = 'I';
            --j;
        } else {
            break;
        }

        // Add to cigar string
        if(cur == prev) {
            ++count;
        } else {
            if(count) {
                cigar_iter += sprintf(
                    cigar + cigar_iter, 
                    "%d%c",
                    count,
                    prev 
                );
            }
            count = 1;
            prev = cur;
        }
    }

    // Check for leftovers
    if(count) {
        cigar_iter += sprintf(
            cigar + cigar_iter, 
            "%d%c",
            count,
            prev 
        );
    }
    cigar[cigar_iter] = '\0';

    // Cigar string is built in reverse order
    for(int k = 0; k < cigar_iter / 2; ++k) {
        char temp = cigar[k];
        cigar[k] = cigar[cigar_iter - k - 1];
        cigar[cigar_iter - k - 1] = temp;
    }

    // "Return" in RET
    int RET_BEGIN = SEQ_NUM * 3;
    RET[RET_BEGIN] = MAX_I - 1;
    RET[RET_BEGIN + 1] = MAX_J - 1;
    RET[RET_BEGIN + 2] = MAX_SCORE;
}