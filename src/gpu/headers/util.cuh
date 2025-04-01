#pragma once
#pragma once
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <locale>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <stdio.h>
#include <inttypes.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#ifndef CU_TYPE_MACROS
#define CU_TYPE_MACROS
#define ULL unsigned long long int
#endif

#ifndef FASTX_FILTER_MACROS
#define FASTX_FILTER_MACROS
#define AVERAGE_DISCARD_WHOLE 1
#define SINGLE_DISCARD_WHOLE 2
#define SLIDING_WINDOW 3
#define PROPORTION_DISCARD_WHOLE 4
#endif

#ifndef CU_KMER_MACROS
#define CU_KMER_MACROS
#define K_MAX 32
#define MAP_MAX_INDICES 128
#define MAX_CLUSTER_SIZE 128
#endif

#ifndef CU_LOCAL_ALIGNMENT_MACROS
#define CU_LOCAL_ALIGNMENT_MACROS
#define MAX_CIGAR_LEN 340
#define MAX_REF_LEN 200
#define MAX_READ_LEN 170
#define ROWS (MAX_REF_LEN + 1)
#define COLS (MAX_READ_LEN + 1)
#define DP(i, j) cache[(i) * COLS + (j)]
#define ALIGN_MATCH 2
#define ALIGN_MISMATCH -1
#define ALIGN_GAP -2
#endif

#ifndef CU_CALL_CHECK
#define CU_CALL_CHECK
#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if(err != cudaSuccess) { \
            fprintf(stderr, "CUDA Error: %s at %s:%d\n", \
                cudaGetErrorString(err), __FILE__, __LINE__); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)
#endif

__device__ int __max(int a, int b);

__device__ char __base_to_bit(char base);