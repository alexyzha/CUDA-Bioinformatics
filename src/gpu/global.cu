#include "headers/global.cuh"

__global__ gpu_filter_reads(char* ALL_SEQ, uint32_t* OFFSETS, size_t LEN, char FILTER_MODE, char THRESH, uint64_t* FILTER_MASK, double PROPORTION = 0.0) {
    cu_filter_reads(ALL_SEQ, OFFSETS, LEN, FILTER_MODE, THRESH, FILTER_MASK, PROPORTION);
}