#include "headers/util.cuh"

__device__ int __max(int a, int b) {
    return a > b ? a : b;
}

__device__ char __base_to_bit(char base) {
    switch(base) {
        case 'A':
            return 0b00;
        case 'C':
            return 0b01;
        case 'G':
            return 0b10;
        case 'T':
            return 0b11;
        default:
            return 0x80;
    }
    return -1;
}

__device__ int __reverse_to_cigar(char* str, int count, char op) {
    if(!count || !op) {
        return 0;
    }

    // Write op then digits in reverse
    int i = 1;
    str[0] = op;
    while(count) {
        str[i] = '0' + (count % 10);
        count /= 10;
        ++i;
    }
    return i;
}