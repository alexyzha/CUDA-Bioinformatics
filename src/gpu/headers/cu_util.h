#pragma once
#include <stdio.h>
#include <inttypes.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

__device__ int __max(int a, int b);

__device__ char __base_to_bit(char base);