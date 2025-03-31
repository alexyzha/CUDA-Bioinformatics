#pragma once
#include <stdio.h>
#include <inttypes.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#ifndef C_CU_UF
#define C_CU_UF

typedef struct {
    uint32_t* p;
    uint32_t* h;
    size_t LEN;
} cu_union_find;

cu_union_find* uf_construct(int n);

__device__ int cu_uf_find(cu_union_find* UF, int x);

__device__ void cu_uf_join(cu_union_find* UF, int x, int y);

__device__ bool cu_uf_con(cu_union_find* UF, int x, int y);

#endif

__device__ int __max(int a, int b);

__device__ char __base_to_bit(char base);