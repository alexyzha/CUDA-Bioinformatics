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
    int i = 1;
    int temp = count;

    // Special 0 case
    if(temp == 0) {
        str[1] = '0';
        str[2] = op;
        return 3;
    }

    // Write op then digits in reverse
    str[0] = op;
    while (temp > 0) {
        str[i + 1] = '0' + (temp % 10);
        temp /= 10;
        ++i;
    }
    return i;
}