cmake_minimum_required(VERSION 3.18)
project(GPU_TESTS LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_CUDA_ARCHITECTURES 86)

include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/cpu
    ${CMAKE_SOURCE_DIR}/src/cpu/headers
    ${CMAKE_SOURCE_DIR}/src/gpu
    ${CMAKE_SOURCE_DIR}/src/gpu/headers
    ${CMAKE_SOURCE_DIR}/tests/gpu/headers
)

set(CPU_SOURCES
    src/cpu/fa_read.cpp
    src/cpu/fq_read.cpp
    src/cpu/fx_parser.cpp
    src/cpu/fx_util.cpp
    src/cpu/sam_container.cpp
    src/cpu/util_structs.cpp
    src/cpu/util.cpp
    src/cpu/xam_parser.cpp
    src/cpu/xam_util.cpp
    tests/gpu/test_util.cpp
)

set(GPU_SOURCES
    src/gpu/fq_filters.cu
    src/gpu/kmer_util.cu
    src/gpu/local_alignment.cu
    src/gpu/util_structs.cu
    src/gpu/util.cu
    src/gpu/wrappers.cu
)

set(TEST_SOURCES
    tests/gpu/all_gpu_tests.cu
    tests/gpu/global_wrappers.cu
    tests/gpu/util_structs_tests.cu
    tests/gpu/util_tests.cu
    tests/gpu/wrapper_tests.cu
)

set(ALL_SOURCES
    ${CPU_SOURCES}
    ${GPU_SOURCES}
    ${TEST_SOURCES}
)

add_executable(all_gpu_tests ${ALL_SOURCES})
