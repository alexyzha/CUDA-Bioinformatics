#include "headers/cu_util.h"

__device__ int __max(int a, int b) {
    return a > b ? a : b;
}